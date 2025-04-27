//===----------------------------------------------------------------------===//
//
// Created by jianglimin on 23-8-23.
//
//===----------------------------------------------------------------------===//
//
// Assumption: All functions are forced to be inlined
//
//
//===----------------------------------------------------------------------===//

#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstrTypes.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Intrinsics.h"
#include "llvm/IR/Value.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/CompareNat.h"
#include "llvm/Support/Debug.h"
#include "llvm/Transforms/Utils/VenusIntrinSplit.h"
#include <map>
#include <stack>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <queue>
using namespace llvm;

namespace opts {
extern llvm::cl::opt<unsigned> VenusNrRows;
extern llvm::cl::opt<unsigned> VenusNrLanes;
extern llvm::cl::opt<unsigned> VenusCRBaseAddr;
extern llvm::cl::opt<std::string> VenusTaskName;

cl::opt<unsigned> VenusInputStructAddr(
    "venus-istruct-baseaddr",
    cl::desc("Specify the input address of the scalar struct of the Venus architecture."),
    cl::value_desc("value"),
    cl::init(0x8000a000)
);

unsigned VenusScatterID = 0;
} // namespace opts

struct RowInfo {
  unsigned VElements;       // Total vector elements in ELEN
  unsigned ElementsPerRow;
  unsigned NumRowsCeil;
  RowInfo() = default;
  RowInfo(unsigned velem, unsigned e, unsigned r)
      : VElements(velem), ElementsPerRow(e), NumRowsCeil(r) {}
};

RowInfo getIntrinsicRowInfo(Type *T) {
  unsigned const ElementTypeInBits = T->getScalarSizeInBits();
  unsigned const VElements = T->getPrimitiveSizeInBits() >>
                             llvm::Log2_32(ElementTypeInBits);
  // Get Venus architecture-specific information
  unsigned const NrLanes = opts::VenusNrLanes;
  unsigned ElementsPerRowInVenusArch;
  switch (ElementTypeInBits) {
  case 8:
    ElementsPerRowInVenusArch = 4/* NrBanks */ * NrLanes * 2/* SPM */;
    break;
  case 16:
    ElementsPerRowInVenusArch = 4/* NrBanks */ * NrLanes;
    break;
  default:
    llvm_unreachable("Unsupported Vector Element Type!");
  }
  unsigned NumRowsCeil = VElements >>
                         llvm::Log2_32(ElementsPerRowInVenusArch);
  if (VElements & (ElementsPerRowInVenusArch - 1))
    NumRowsCeil++;

  RowInfo rowInfo(VElements, ElementsPerRowInVenusArch, NumRowsCeil);
  return rowInfo;
}

bool inline intrinsicFilterNormal(unsigned Opcode) {
  const bool isLoad = (Opcode == Instruction::Load);
  const bool isVenus = ((Opcode >= Intrinsic::venus_add_ivv) &&
                        (Opcode <= Intrinsic::venus_xor_ivx));
  const bool isSpecial = (Opcode == Intrinsic::venus_mnot_misc) ||
                         (Opcode == Intrinsic::venus_barrier) ||
                         (Opcode == Intrinsic::venus_delimit) ||
                         (Opcode == Intrinsic::venus_sink) ||
                         (Opcode == Intrinsic::venus_bind) ||
                         (Opcode == Intrinsic::venus_scalar_bind) ||
                         (Opcode == Intrinsic::venus_setcsr);
   return (isLoad | isVenus) & (!isSpecial);
}


bool isNeedSplit(Value* Op) {
  RowInfo rowInfo = getIntrinsicRowInfo(Op->getType());
  const bool isOneRow = (rowInfo.NumRowsCeil == 1);
  const bool isFullRow = (Op->getType()->getPrimitiveSizeInBits() ==
                          opts::VenusNrLanes * 4 /* banks per lane */
                              * 2 /* bytes per bank */ * 8 /* bits per byte */);
  return (isOneRow & isFullRow);
}

bool isNeedSink(Value *Op) {
    Type *OpType = Op->getType();

    bool NeedSinkNormal = OpType->isVectorTy() && Op->use_empty();

    bool NeedSinkPHI = false;

    for (User *U : Op->users()) {
        if (Instruction *In = dyn_cast<Instruction>(U)) {
            if (In->getOpcode() == Instruction::PHI) {
                Value *PHIRet = dyn_cast<Value>(In);
                NeedSinkPHI = PHIRet->getType()->isVectorTy() && PHIRet->use_empty();
            }
        }
    }
    // Return whether sinking is needed, true if either condition is met
    return (NeedSinkNormal | NeedSinkPHI);
}

// Locate PHI nodes, and check if any of their operands contain undefined values.
// If such operands are found, replace them with the corresponding actual operands in VMap.

void VenusIntrinSplitPass::FindAndReplacePHIUndef(Value *Opd) {
  Type *OpdType = Opd->getType();
  StructType *OpdStructType = dyn_cast<StructType>(OpdType);
  if (OpdStructType && OpdStructType->getNumElements() == 2) {
      OpdType = OpdStructType->getElementType(0);
  }
   // find whether the used phi node has an undef value
  for (User *U : Opd->users()) {
    if (Instruction *In = dyn_cast<Instruction>(U)) {
      if (In->getOpcode() == Instruction::PHI) {
        Value *PHIRet = dyn_cast<Value>(In);
        // Skip irrelevant phi nodes
        if (VMap[PHIRet].empty()) continue;

        RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);
        for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
          Instruction *SplitIn = VMap[PHIRet].at(R);
          bool SameOperand = false;
          for (unsigned Op = 0; Op < SplitIn->getNumOperands(); Op++) {
            // Do not glue the operands that are already existed
            // SplitIn->dump() ;
            if (SplitIn->getOperand(Op) == VMap[Opd].at(R)) {
              SameOperand = true;
            }
          }
          for (unsigned Op = 0;
               !SameOperand && Op < SplitIn->getNumOperands(); Op++) {
            if (isa<UndefValue>(SplitIn->getOperand(Op))) {
              LLVM_DEBUG(dbgs() << "Found undef values in " << *SplitIn << "\n");
              // Glue the loop
              SplitIn->setOperand(Op, VMap[Opd].at(R));
              break; // only glue once
            }
          }
        }

      }
    }
  }
}

void VenusIntrinSplitPass::VenusSplitNormal(Instruction *Inst,
                                           Intrinsic::ID intrinsicID) {
  // Check whether the variables are accessed for the first time
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Value *Op1 = Inst->getOperand(0);
  Type *Op1Type = Op1->getType();
  Value *Op2 = Inst->getOperand(1);
  Type *Op2Type = Op2->getType();
  Value *OpVEW = Inst->getOperand(2);
  Type *OpVEWType = OpVEW->getType();
  Value *OpVMR = Inst->getOperand(3); // vmask_read
  Value *OpAVL = Inst->getOperand(4); // AVL
  Type *OpAVLType = OpAVL->getType();
  Value *OpRows = Inst->getOperand(5);
  // Only consider three source operands first
  const bool is3SrcOperands = (Inst->getNumOperands() == 8);
  Value *Op3;
  Value *Op4;
  if (intrinsicID == Intrinsic::venus_cmxmul_ivv) {
    Op3 = Inst->getOperand(2);
    Op4 = Inst->getOperand(3);
    OpVEW = Inst->getOperand(4);
    OpVEWType = OpVEW->getType();
    OpVMR = Inst->getOperand(5);
    OpAVL = Inst->getOperand(6);
    OpAVLType = OpAVL->getType();
    OpRows = Inst->getOperand(7);
    if (VMap.find(Op3) == VMap.end())
      VMap.insert(std::pair<Value *, std::vector<Instruction *>>(Op3, NULL));
    if (VMap.find(Op4) == VMap.end())
      VMap.insert(std::pair<Value *, std::vector<Instruction *>>(Op4, NULL));
  }
  else if (is3SrcOperands) {
    Op3 = Inst->getOperand(2);
    OpVEW = Inst->getOperand(3);
    OpVEWType = OpVEW->getType();
    OpVMR = Inst->getOperand(4);
    OpAVL = Inst->getOperand(5);
    OpAVLType = OpAVL->getType();
    OpRows = Inst->getOperand(6);
    if (VMap.find(Op3) == VMap.end())
      VMap.insert(std::pair<Value *, std::vector<Instruction *>>(Op3, NULL));
  }
  
  RowInfo const rowInfo = getIntrinsicRowInfo(Op1Type);
  // No need to split if there is only one row
  if (isNeedSplit(Op1)) return;

  LLVM_DEBUG(dbgs() << "Intrinsic Found: "
                    << Intrinsic::getBaseName(intrinsicID) << "\n");
  LLVM_DEBUG(dbgs() << "Operand has " << rowInfo.VElements << " elements\n");

  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {},
                          {Builder.getInt32(rowInfo.NumRowsCeil)});
  for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));
    Value *SubElement2 = Op2Type->isVectorTy()
                             ? UndefValue::get(VectorType::get(Op2Type->getScalarType(),
                                                               rowInfo.ElementsPerRow, false))
                             : Op2;
    Value *SubElement3;
    Value *SubElement4;
    if (intrinsicID == Intrinsic::venus_cmxmul_ivv){
      SubElement3 = UndefValue::get(VectorType::get(
          Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));    
      SubElement4 = UndefValue::get(VectorType::get(
          Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));
      }
    else if (is3SrcOperands)
      SubElement3 = UndefValue::get(VectorType::get(
          Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan
    if (VMap[Op1].empty() && ((intrinsicID != Intrinsic::venus_brdcst_ivx) ||
        (!isa<UndefValue>(Op1) && intrinsicID == Intrinsic::venus_brdcst_ivx))) {
      return;
    }
    if (VMap[Op2].empty() && Op2Type->isVectorTy()) {
      return;
    }

     if (intrinsicID == Intrinsic::venus_cmxmul_ivv && (VMap[Op3].empty() || VMap[Op4].empty())) {
       return;
     }

    if (is3SrcOperands && VMap[Op3].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert((!VMap[Op1].empty() || (intrinsicID == Intrinsic::venus_brdcst_ivx))
           && (!VMap[Op2].empty() || !Op2Type->isVectorTy())
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty()) SubElement1 = VMap[Op1].at(R);
    if (!VMap[Op2].empty()) SubElement2 = VMap[Op2].at(R);
    llvm::SmallVector<llvm::Value *, 5>
        NewInstOperands = {SubElement1, SubElement2,
                           OpVEW, OpVMR, OpAVL, OpRows};
    if (intrinsicID == Intrinsic::venus_cmxmul_ivv) {
      assert(!VMap[Op3].empty()
             && "Undef source operand when splitting instructions!");
      assert(!VMap[Op4].empty()
             && "Undef source operand when splitting instructions!");
      SubElement3 = VMap[Op3].at(R);
      SubElement4 = VMap[Op4].at(R);
      NewInstOperands = {SubElement1, SubElement2, SubElement3, SubElement4,
                         OpVEW, OpVMR, OpAVL, OpRows};
    }
    else if (is3SrcOperands) {
      assert(!VMap[Op3].empty()
             && "Undef source operand when splitting instructions!");
      SubElement3 = VMap[Op3].at(R);
      NewInstOperands = {SubElement1, SubElement2, SubElement3,
                         OpVEW, OpVMR, OpAVL, OpRows};
    }
    if (intrinsicID == Intrinsic::venus_cmxmul_ivv) {
      CallInst *NewInst = Builder.CreateIntrinsic(
      intrinsicID,
      {SubElement1->getType(), OpVEWType},
        NewInstOperands);
      VMap[Opd].push_back(NewInst);
    }    
    else{
      CallInst *NewInst = Builder.CreateIntrinsic(
      intrinsicID,
      {SubElement1->getType(), OpVEWType, OpAVLType},
        NewInstOperands);
      VMap[Opd].push_back(NewInst);
      LLVM_DEBUG(dbgs() << "opd: \n" << *Opd << "\n");
      LLVM_DEBUG({
          dbgs() << "VMap: \n";
          for (auto *Inst : VMap[Opd]) {
              dbgs() << *Inst << "\n";
          }
      });
    }
  } // row-wise
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});
  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitVaddr(Instruction *Inst,
                                             Intrinsic::ID intrinsicID) {
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Value *Op1 = Inst->getOperand(0);
  Type *Op1Type = Op1->getType();

  RowInfo const rowInfo = getIntrinsicRowInfo(Op1Type);
  // No need to split if there is only one row
  if (isNeedSplit(Op1)) return;

  // Add venus_sink to avoid incorrect dest register allocation
  for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty()
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty())
      SubElement1 = VMap[Op1].at(R);
    Builder.CreateIntrinsic(Intrinsic::venus_sink,
                            {SubElement1->getType()}, {SubElement1});
  }

  // We only emit one venus load intrinsic, and we will have only one
  // scalar register to be allocated.
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(1)});
  for (unsigned R = 0; R < 1; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty()
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty())
      SubElement1 = VMap[Op1].at(R);
    CallInst *NewInst = Builder.CreateIntrinsic(
        intrinsicID,
        {SubElement1->getType()},
        {SubElement1});
    VMap[Opd].push_back(NewInst);
  } // row-wise

   for (unsigned R = 1; R < rowInfo.NumRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty()
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty())
      SubElement1 = VMap[Op1].at(R);
    Builder.CreateIntrinsic(Intrinsic::venus_sink,
                            {SubElement1->getType()}, {SubElement1});
  }

  // The return value is a scalar value. The first operand of
  // the following instruction  will turn out to be `undef`
  // because the "unsplit" intrinsic has been removed.
  // The following codes use `user iterator` to find every instruction
  // which uses the return value of the "unsplit" intrinsic,
  // then replace them with the first return value of "split" intrinsic.
  for (User *U : Opd->users())
    if (Instruction *I = dyn_cast<Instruction>(U)) {
      LLVM_DEBUG(dbgs() << "Get Users: " << *I << "\t");
      unsigned Op;
      for (Op = 0; Op < I->getNumOperands(); Op++){
        if (I->getOperand(Op) == Opd) break;
      }
      LLVM_DEBUG(dbgs() << "@ Operand No: " << Op << "\n");
      I->setOperand(Op, VMap[Opd].at(0));
    }
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});
  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitShuffle(Instruction *Inst,
                                           Intrinsic::ID intrinsicID) {
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  Value *Op1 = Inst->getOperand(0);
  Type *Op1Type = Op1->getType();
  Value *Op2 = Inst->getOperand(1);
  Type *Op2Type = Op2->getType();
  Value *Op3 = Inst->getOperand(2);
  Type *Op3Type = Op3->getType();
  Value *OpVEW = Inst->getOperand(3);
  Type *OpVEWType = OpVEW->getType();
  Value *OpVMR = Inst->getOperand(4); // vmask_read
  Value *OpAVL = Inst->getOperand(5); // AVL
  Type *OpAVLType = OpAVL->getType();
  Value *OpRows = Inst->getOperand(6);

  RowInfo const rowInfod = getIntrinsicRowInfo(OpdType);
  RowInfo const rowInfo1 = getIntrinsicRowInfo(Op1Type);

  // No need to split if there is only one row
  if (isNeedSplit(Op1)) return;

  // Add venus_sink to avoid incorrect dest register allocation
  for (unsigned R = 0; R < rowInfo1.NumRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo1.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty()
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty())
      SubElement1 = VMap[Op1].at(R);
    Builder.CreateIntrinsic(Intrinsic::venus_sink,
                            {SubElement1->getType()}, {SubElement1});
  }

  if (Op2Type->isVectorTy()) {
    // Add venus_sink to avoid incorrect dest register allocation
    RowInfo const rowInfo2 = getIntrinsicRowInfo(Op2Type);
    for (unsigned R = 0; R < rowInfo2.NumRowsCeil; R++) {
      Value *SubElement1 = UndefValue::get(VectorType::get(
          Op1Type->getScalarType(), rowInfo2.ElementsPerRow, false));

      // If the variables appear for the first time,
      // we need to split at the next scan.
      if (VMap[Op2].empty()) {
        return;
      }

      // Create a new instruction
      // If this variable are related to the previous instructions,
      // use them instead.
      assert(!VMap[Op2].empty() &&
             "Undef source operand when splitting instructions!");
      if (!VMap[Op2].empty())
        SubElement1 = VMap[Op2].at(R);
      Builder.CreateIntrinsic(Intrinsic::venus_sink, {SubElement1->getType()},
                              {SubElement1});
    }
  }

  // Add venus_sink to avoid incorrect dest register allocation
  RowInfo const rowInfo3 = getIntrinsicRowInfo(Op3Type);
  for (unsigned R = 0; R < rowInfo1.NumRowsCeil; R++) {
    Value *SubElement3 = UndefValue::get(VectorType::get(
        Op3Type->getScalarType(), rowInfo1.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op3].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op3].empty()
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op3].empty())
      SubElement3 = VMap[Op3].at(R);
    Builder.CreateIntrinsic(Intrinsic::venus_sink,
                            {SubElement3->getType()}, {SubElement3});
  }

  Value *OpRet = UndefValue::get(llvm::VectorType::get(
      Builder.getIntNTy(OpdType->getScalarSizeInBits()),
      rowInfod.ElementsPerRow, false));
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {},
                          {Builder.getInt32(rowInfod.NumRowsCeil)});
  const unsigned MinRowsCeil =
      std::min(rowInfod.NumRowsCeil, rowInfo1.NumRowsCeil);
  unsigned RemainRowsCeil = 0;
  std::vector<CallInst *> RemainInputs;
  if (rowInfod.NumRowsCeil > rowInfo1.NumRowsCeil) {
    // Scatter needs extra vector claims
    RemainRowsCeil = rowInfod.NumRowsCeil - rowInfo1.NumRowsCeil;
  }

  for (unsigned R = 0; R < MinRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfod.ElementsPerRow, false));
    Value *SubElement2 = Op2Type->isVectorTy() ?
                                               UndefValue::get(VectorType::get(
                                                   Op2Type->getScalarType(),
                                                   rowInfod.ElementsPerRow,
                                                   false)) : Op2;
    Value *SubElement3 = UndefValue::get(VectorType::get(
        Op3Type->getScalarType(), rowInfod.ElementsPerRow, false));

    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }
    if (VMap[Op2].empty() && Op2Type->isVectorTy()) {
      return;
    }
    if (VMap[Op3].empty()) {
      return;
    }

    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty() && (!VMap[Op2].empty() || !Op2Type->isVectorTy())
           && "Undef source operand when splitting instructions!");
    if (!VMap[Op1].empty())
      SubElement1 = VMap[Op1].at(R);
    if (!VMap[Op2].empty())
      SubElement2 = VMap[Op2].at(R);
    if (!VMap[Op3].empty())
      SubElement3 = VMap[Op3].at(R);
    CallInst *NewInst = Builder.CreateIntrinsic(
        intrinsicID,
        {OpRet->getType(),
         SubElement1->getType(), SubElement2->getType(),
         OpVEWType, OpAVLType},
        {SubElement1, SubElement2, SubElement3, OpVEW, OpVMR, OpAVL, OpRows});
    VMap[Opd].push_back(NewInst);
  } // row-wise

  // Emit extra intrinsics if the original intrinsic is a scatter one.
  for (unsigned R = 0; R < RemainRowsCeil; R++) {
    ConstantInt *CI = dyn_cast<ConstantInt>(OpVMR);
    unsigned const SGIndicator = CI->getValue().getZExtValue();
    if (SGIndicator)
      report_fatal_error("Gather mode does not allow longer Opd!");
    Value *newSubID = Builder.getInt32(opts::VenusScatterID++);
    CallInst *NewInst = Builder.CreateIntrinsic(
        Intrinsic::venus_scatter,
        {OpRet->getType(), OpVEWType, OpAVLType},
        {OpVEW, OpVMR, newSubID, OpAVL, OpRows});
    VMap[Opd].push_back(NewInst);
  }
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});
  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitShuffleTest(Instruction *Inst,
                                           Intrinsic::ID intrinsicID) {
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  Value *Op1 = Inst->getOperand(0);
  Type *Op1Type = Op1->getType();
  Value *Op2 = Inst->getOperand(1);
  Type *Op2Type = Op2->getType();
  Value *Op3 = Inst->getOperand(2);
  Type *Op3Type = Op3->getType();
  Value *OpVEW = Inst->getOperand(3);
  Type *OpVEWType = OpVEW->getType();
  Value *OpVMR = Inst->getOperand(4); // vmask_read
  Value *OpAVL = Inst->getOperand(5); // AVL
  Type *OpAVLType = OpAVL->getType();
  Value *OpRows = Inst->getOperand(6);
                                            
  RowInfo const rowInfod = getIntrinsicRowInfo(OpdType);
  RowInfo const rowInfo1 = getIntrinsicRowInfo(Op1Type);
  RowInfo const rowInfo2 = getIntrinsicRowInfo(Op2Type);
  RowInfo const rowInfo3 = getIntrinsicRowInfo(Op3Type);
  // No need to split if there is only one row
  if (isNeedSplit(Op1)) return;
  if (!isNeedSplit(Opd)){
    Builder.CreateIntrinsic(Intrinsic::venus_delimit, {},
      {Builder.getInt32(rowInfod.NumRowsCeil)});
  }                                       
  Value *OpRet = UndefValue::get(llvm::VectorType::get(
      Builder.getIntNTy(OpdType->getScalarSizeInBits()),
      rowInfod.ElementsPerRow, false));

  const unsigned MinRowsCeil =
      std::min(rowInfo1.NumRowsCeil, rowInfo2.NumRowsCeil);
    
  for (unsigned R = 0; R < rowInfod.NumRowsCeil; R++) {
    Value *SubElement1 = UndefValue::get(VectorType::get(
        Op1Type->getScalarType(), rowInfo1.ElementsPerRow, false));
    Value *SubElement2 = Op2Type->isVectorTy() ?
                                               UndefValue::get(VectorType::get(
                                                   Op2Type->getScalarType(),
                                                   rowInfo2.ElementsPerRow,
                                                   false)) : Op2;
    Value *SubElement3 = UndefValue::get(VectorType::get(
        Op3Type->getScalarType(), rowInfo3.ElementsPerRow, false));
                                            
    // If the variables appear for the first time,
    // we need to split at the next scan.
    if (VMap[Op1].empty()) {
      return;
    }
    if (VMap[Op2].empty() && Op2Type->isVectorTy()) {
      return;
    }
    if (VMap[Op3].empty()) {
      return;
    }
                                            
    // Create a new instruction
    // If this variable are related to the previous instructions,
    // use them instead.
    assert(!VMap[Op1].empty() && (!VMap[Op2].empty() || !Op2Type->isVectorTy())
           && "Undef source operand when splitting instructions!");
    if(R < MinRowsCeil){
      if (!VMap[Op1].empty())
        SubElement1 = VMap[Op1].at(R);
      if (!VMap[Op2].empty())
        SubElement2 = VMap[Op2].at(R);
      if (!VMap[Op3].empty())
        SubElement3 = VMap[Op3].at(R);
      CallInst *NewInst = Builder.CreateIntrinsic(
          intrinsicID,
          {OpRet->getType(),
          SubElement1->getType(), SubElement2->getType(),
          OpVEWType, OpAVLType},
          {SubElement1, SubElement2, SubElement3, OpVEW, OpVMR, OpAVL, OpRows});
      VMap[Opd].push_back(NewInst);
    } // row-wise
    else{
      if((R < rowInfo1.NumRowsCeil) && (R < rowInfo2.NumRowsCeil)){
        if (!VMap[Op2].empty())
          SubElement2 = VMap[Op2].at(R);
        if (!VMap[Op3].empty())
          SubElement3 = VMap[Op3].at(R);
        CallInst *NewInst = Builder.CreateIntrinsic(
          intrinsicID,
          {OpRet->getType(),
          SubElement1->getType(), SubElement2->getType(),
          OpVEWType, OpAVLType},
          {SubElement1, SubElement2, 
          SubElement3, OpVEW, OpVMR, OpAVL, OpRows});
        VMap[Opd].push_back(NewInst);
      }
      else if((R >= rowInfo1.NumRowsCeil) && (R < rowInfo2.NumRowsCeil)){
        if (!VMap[Op1].empty())
          SubElement1 = VMap[Op1].at(R);
        if (!VMap[Op3].empty())
          SubElement3 = VMap[Op3].at(R);
        CallInst *NewInst = Builder.CreateIntrinsic(
          intrinsicID,
          {OpRet->getType(),
          SubElement1->getType(), SubElement2->getType(),
          OpVEWType, OpAVLType},
          {UndefValue::get(SubElement1->getType()), 
           SubElement2,
          SubElement3, OpVEW, OpVMR, OpAVL, OpRows});
        VMap[Opd].push_back(NewInst);
      }
      else if((R >= rowInfo1.NumRowsCeil) && (R >= rowInfo2.NumRowsCeil)){
        if (!VMap[Op1].empty())
          SubElement1 = VMap[Op1].at(R);
        if (!VMap[Op2].empty())
          SubElement2 = VMap[Op2].at(R);
        if (!VMap[Op3].empty())
          SubElement3 = VMap[Op3].at(R);
        CallInst *NewInst = Builder.CreateIntrinsic(
          intrinsicID,
          {OpRet->getType(),
          SubElement1->getType(), SubElement2->getType(),
          OpVEWType, OpAVLType},
          {UndefValue::get(SubElement1->getType()), 
            UndefValue::get(SubElement2->getType()), 
          SubElement3, OpVEW, OpVMR, OpAVL, OpRows});
        VMap[Opd].push_back(NewInst);
      }
  }
  }
                                            
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});

  for (unsigned R = 0; R < rowInfo1.NumRowsCeil; R++) {
          Value *SubElement1 = UndefValue::get(VectorType::get(
              Op1Type->getScalarType(), rowInfo1.ElementsPerRow, false));
                                                
          // If the variables appear for the first time,
          // we need to split at the next scan.
          if (VMap[Op1].empty()) {
            return;
          }
                                                
          // Create a new instruction
          // If this variable are related to the previous instructions,
          // use them instead.
          assert(!VMap[Op1].empty() &&
                "Undef source operand when splitting instructions!");
          if (!VMap[Op1].empty())
            SubElement1 = VMap[Op1].at(R);
          Builder.CreateIntrinsic(Intrinsic::venus_sink, {SubElement1->getType()},
                                  {SubElement1});                                                
        }
      
        for (unsigned R = 0; R < rowInfo2.NumRowsCeil; R++) {
          Value *SubElement2 = UndefValue::get(VectorType::get(
              Op2Type->getScalarType(), rowInfo2.ElementsPerRow, false));
                                                
          // If the variables appear for the first time,
          // we need to split at the next scan.
          if (VMap[Op2].empty()) {
            return;
          }
                                                
          // Create a new instruction
          // If this variable are related to the previous instructions,
          // use them instead.
          assert(!VMap[Op2].empty() &&
                "Undef source operand when splitting instructions!");
          if (!VMap[Op2].empty())
            SubElement2 = VMap[Op2].at(R);
          Builder.CreateIntrinsic(Intrinsic::venus_sink, {SubElement2->getType()},
                                  {SubElement2});
                                                
        }

  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitClaim(Instruction *Inst,
                                           Intrinsic::ID intrinsicID) {
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  Value *OpVEW = Inst->getOperand(0);
  Type *OpVEWType = OpVEW->getType();
  Value *OpVMR = Inst->getOperand(1); // vmask_read
  Value *OpSubID = Inst->getOperand(2);
  ConstantInt *CI = dyn_cast<ConstantInt>(OpSubID);
  unsigned const SubID = CI->getValue().getZExtValue();
  Value *OpAVL = Inst->getOperand(3); // AVL
  Type *OpAVLType = OpAVL->getType();
  Value *OpRows = Inst->getOperand(4);

  RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);
  // No need to split if there is only one row
  if (isNeedSplit(Opd)) return;

  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {},
                          {Builder.getInt32(rowInfo.NumRowsCeil)});
  Value *OpRet = UndefValue::get(llvm::VectorType::get(
      Builder.getIntNTy((intrinsicID == Intrinsic::venus_range_misc) ? 16 :
         OpdType->getScalarSizeInBits()),
      rowInfo.ElementsPerRow, false));
  for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
    Value *newSubID = Builder.getInt32(SubID + R);
//    LLVM_DEBUG(dbgs() << "Venus range type subID: " << SubID + R << "\n");
    CallInst *NewInst = Builder.CreateIntrinsic(
        intrinsicID,
        {OpRet->getType(), OpVEWType, OpAVLType},
        {OpVEW, OpVMR, newSubID, OpAVL, OpRows});
    VMap[Opd].push_back(NewInst);
  } // row-wise
  Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});
  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitIntrinsic (Instruction *Inst,
                                               Intrinsic::ID intrinsicID,
                                               bool NeedSink) {
  IRBuilder<> Builder(Inst);
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  StructType *OpdStructType = dyn_cast<StructType>(OpdType);
  if (OpdStructType && OpdStructType->getNumElements() == 2) {
      OpdType = OpdStructType->getElementType(0);
  }
  // return here if this intrinsic has already been splatted
  if (std::find(VListIntrinsics.begin(), VListIntrinsics.end(), Inst)
      != VListIntrinsics.end())
    return;

  if ((intrinsicID == Intrinsic::venus_range_misc) ||
      (intrinsicID == Intrinsic::venus_claim_misc)) {
    VenusSplitClaim(Inst, intrinsicID);
  } else if ((intrinsicID == Intrinsic::venus_shuffle_ivv) ||
             (intrinsicID == Intrinsic::venus_shuffle_ivx)){
    VenusSplitShuffle(Inst, intrinsicID);
  } else if (intrinsicID == Intrinsic::venus_shuffle_test){
    VenusSplitShuffleTest(Inst, intrinsicID);
  } else if ((intrinsicID == Intrinsic::venus_vaddr_misc)) {
    VenusSplitVaddr(Inst, intrinsicID);
  } else if (intrinsicFilterNormal(intrinsicID)) {
    VenusSplitNormal(Inst, intrinsicID);
  } 

  // Add venus_sink to avoid incorrect dest register allocation
  if (NeedSink) {
    RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);
    if (rowInfo.NumRowsCeil > 1) {
      for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
        if (VMap[Opd].empty()) continue;
        Value *OpSink = VMap[Opd].at(R);
        Builder.CreateIntrinsic(Intrinsic::venus_sink, {OpSink->getType()},
                                {OpSink});
      }
    }
  }
}

void VenusIntrinSplitPass::VenusSplitExtractValue(Instruction *Inst) {
  Value *Op1 = Inst->getOperand(0);
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  IRBuilder<> Builder(Inst);
  RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);

  // return here if this intrinsic has already been splatted
  if (std::find(VListIntrinsics.begin(), VListIntrinsics.end(), Inst)
      != VListIntrinsics.end())
    return;

  if (ExtractValueInst *EVInst = dyn_cast<ExtractValueInst>(Inst)) {
    Value *AggregateOperand = EVInst->getAggregateOperand();
    ArrayRef<unsigned> Indices = EVInst->getIndices();
    
    for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
      if (VMap[Op1].empty()) {
          return;
      }
      Value *SubElement = UndefValue::get(VectorType::get(
         OpdType->getScalarType(), rowInfo.ElementsPerRow, false
         ));
      assert(!VMap[Op1].empty()
             && "Undef source operand when splitting instructions!");
      SubElement = VMap[Op1].at(R);
      Value *ExtractValue_split = Builder.CreateExtractValue(SubElement,Indices[0]);
      VMap[Opd].push_back(dyn_cast<ExtractValueInst>(ExtractValue_split));
    }
  }     
  VListIntrinsics.push_back(Inst);
}

void VenusIntrinSplitPass::VenusSplitPHI(Instruction *Inst,
                                         bool SkipUndef = true) {
  Value *Opd = Inst;
  Type *OpdType = Opd->getType();
  IRBuilder<> Builder(Inst);

  // return here if this intrinsic has already been splatted
  if (std::find(VListIntrinsics.begin(), VListIntrinsics.end(), Inst)
      != VListIntrinsics.end())
    return;

  PHINode *CurrPN = dyn_cast<PHINode>(Inst);
  unsigned NumValuePN = CurrPN->getNumIncomingValues();
  llvm::SmallVector<BasicBlock *, 4> IncomingBB;
  LLVM_DEBUG(dbgs() << "Vector PHI Node @ " << Inst->getParent()->getName() <<
             " has " << NumValuePN << " incoming values\n");
  if (SkipUndef && Opd->use_empty()) {
    LLVM_DEBUG(dbgs() << " Dead Vector PHI Node\n");
    VListIntrinsics.push_back(Inst);
    return;
  }
  for (unsigned I = 0; I < NumValuePN; I++)
    IncomingBB.push_back(CurrPN->getIncomingBlock(I));

  RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);
  // No need to split if there is only one row
  if (isNeedSplit(Opd)) return;

  for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++) {
    for (unsigned Op = 0; Op < NumValuePN; Op++) {
      if (VMap[Inst->getOperand(Op)].empty()) {
        if (SkipUndef) {
          LLVM_DEBUG(dbgs() << "There are undef values in PHI Nodes: @\t"
                            << Inst->getParent()->getName() << "\tThis"
                               " should be handled in the next sub-pass.\n");
          return;
        }
      }
    }

    PHINode *PHISplit = Builder.CreatePHI(
        VectorType::get(OpdType->getScalarType(),
                        rowInfo.ElementsPerRow, false)
            , NumValuePN);
    for (unsigned Op = 0; Op < NumValuePN; Op++) {
      Value *SubElement = UndefValue::get(VectorType::get(
          OpdType->getScalarType(), rowInfo.ElementsPerRow, false
          ));
      // Create a new instruction
      // If this variable are related to the previous instructions,
      // use them instead.
//      assert(!VMap[Inst->getOperand(Op)].empty()
//             && "Undef source operand when splitting instructions!");
      if (!VMap[Inst->getOperand(Op)].empty())
        SubElement = VMap[Inst->getOperand(Op)].at(R);
      PHISplit->addIncoming(SubElement, IncomingBB[Op]);
    }
    VMap[Opd].push_back(PHISplit);
  } // row-wise
  VListIntrinsics.push_back(Inst);
}

PreservedAnalyses VenusIntrinSplitPass::run(Function &F,
                                            FunctionAnalysisManager &AM) {
  unsigned VirtualInputCounter = 0;
  unsigned StructInputAddress = opts::VenusInputStructAddr;
  // unsigned ClaimInputCounter = 0;
  LLVM_DEBUG(dbgs() << "VenusNrRows: " << opts::VenusNrRows << "\n");
  LLVM_DEBUG(dbgs() << "VenusNrLanes: " << opts::VenusNrLanes << "\n");
  // The task input in a DAG should only be venus vectors and struct
  // These input are represented in opaque pointer
  std::map<Value *, std::pair<bool, unsigned>> InputMap;
  std::ofstream debug_input_file(opts::VenusTaskName + "_input_scalar_num.txt");
  int real_input_scalar_NumOperands = 0;
  LLVM_DEBUG(dbgs() << "VenusTaskName: " << opts::VenusTaskName << "\n");

  for (Function::arg_iterator AI = F.arg_begin(), E = F.arg_end();
       AI != E; AI++) {
    Argument *Arg = &*AI;

    if (!Arg->getType()->isPointerTy())
      report_fatal_error("Inputs should only be venus vectors and structs!");
    // if (Arg->getType()->getNonOpaquePointerElementType()->isStructTy()){
    real_input_scalar_NumOperands ++;
    LLVM_DEBUG(dbgs() << "Function: " << opts::VenusTaskName << "\n");
    LLVM_DEBUG(dbgs() << "real_input_scalar_NumOperands: " << real_input_scalar_NumOperands << "\n");
    

    if (InputMap.find(Arg) == InputMap.end())
      InputMap[Arg] = std::make_pair(false, NULL);
  }
  
  if (debug_input_file.is_open()) {
      debug_input_file << real_input_scalar_NumOperands << "\n";
      debug_input_file.close();
    } else {
      report_fatal_error("Cannot write input_scalar_Num value to a file!");
    }
    
  // Split normal intrinsics
  LLVM_DEBUG(dbgs() << "********** SPLIT NORMAL INTRINSICS **********\n");
//  LLVM_DEBUG(dbgs() << "Function body:\n" << F << "\n");
  for (Function::iterator BB = F.begin(); BB != F.end(); BB++) {
//    LLVM_DEBUG(dbgs() << "Basic Block:\n" << *BB << "\n");
    for (BasicBlock::iterator I = BB->begin(); I != BB->end(); I++) {
      Instruction *Inst = dyn_cast<Instruction>(I);

      
      LLVM_DEBUG(dbgs() << "Instruction: " << *I << "\n");

      if (Inst->getOpcode() == Instruction::Ret) continue;

      LLVMContext &Context = Inst->getContext();
      Value *Op1 = Inst->getOperand(0);
      Value *Opd = Inst;
      Type *OpdType = Opd->getType();
      
      // Cmxmul instruction has two return values that form a structure.
      // To determine if the return value is a vector
      // Need to deconstruct the structure.
      StructType *OpdStructType = dyn_cast<StructType>(OpdType);
      if (OpdStructType && OpdStructType->getNumElements() == 2) {
          OpdType = OpdStructType->getElementType(0);
      }
      bool NeedSink = isNeedSink(Opd);
      IRBuilder<> Builder(Inst);

      if (NeedSink)
        LLVM_DEBUG(dbgs() << "This instruction need sinking\n");

      // If the number of uses of the operand is 1 (aka. function parameter)
      // in the load instruction, this instruction is called a
      // "bind" instruction in Venus.
      // Only vector ret type needs bind instruction.
      if (Inst->getOpcode() == Instruction::Load &&
          Op1->hasOneUse() && OpdType->isVectorTy()) {
        // Inst->dump();
        // If the dest operand has no uses, drop it otherwise the regalloc would
        // run wild.
        if (Opd->use_empty()) {
          VListIntrinsics.push_back(Inst);
          continue;
        }

        RowInfo const rowInfo = getIntrinsicRowInfo(OpdType);

        if (VMap.find(Op1) == VMap.end())
          VMap.insert(std::pair<Value *, std::vector<Instruction *>>(Op1, NULL));
        if (VMap.find(Opd) == VMap.end())
          VMap.insert(std::pair<Value *, std::vector<Instruction *>>(Opd, NULL));

        Builder.CreateIntrinsic(Intrinsic::venus_delimit, {},
                                {Builder.getInt32(rowInfo.NumRowsCeil)});
        for (unsigned R = 0; R < rowInfo.NumRowsCeil; R++){
          Value *SubElement1 = UndefValue::get(VectorType::get(
              OpdType->getScalarType(), rowInfo.ElementsPerRow, false));
          // Virtual input ensures that the consecutive bind instructions will
          // not be optimized
          ConstantInt *VirtualInput = ConstantInt::get(Type::getInt32Ty(Context),
                                                       VirtualInputCounter);
          VirtualInputCounter++;

          CallInst *NewInst = Builder.CreateIntrinsic(Intrinsic::venus_bind,
                                                      {SubElement1->getType(), VirtualInput->getType()},
                                                      {Op1, VirtualInput});
          VMap[Opd].push_back(NewInst);
          // If there is only one row, just replace all the uses
          if (isNeedSplit(Opd)) {
            Opd->replaceAllUsesWith(NewInst);
          }
        }
        Builder.CreateIntrinsic(Intrinsic::venus_delimit, {}, {Builder.getInt32(0)});
        VListIntrinsics.push_back(Inst);
      }

      // Everything except venus vector should be packed into structs
      if (Inst->getOpcode() == Instruction::GetElementPtr) {
        GetElementPtrInst *GEP = dyn_cast<GetElementPtrInst>(Inst);
        const DataLayout &DL = F.getParent()->getDataLayout();
        const unsigned DataSize = DL.getTypeAllocSize(GEP->getSourceElementType());
        const bool isGEPStructTy = GEP->getSourceElementType()->isStructTy();
        // This pointer should be found in the task argument

        // LLVM_DEBUG(dbgs() << "Function: " << F.getName() << "\n");
        // LLVM_DEBUG(dbgs() << "DataSize: " << DataSize << "\n");

        for (Function::arg_iterator AI = F.arg_begin(), E = F.arg_end();
              AI != E; AI++) {
           Argument *Arg = &*AI;

          // Fetch argument name (fallback for unnamed variables)
        // std::string VarName = Arg->hasName() ? Arg->getName().str() : "<unnamed>";

        // Attempt to extract debug information
        // if (auto *Dbg = llvm::FindDbgDeclare(Arg)) {
        //     if (auto *Var = Dbg->getVariable()) {
        //         VarName = Var->getName().str();
        //         LLVM_DEBUG(dbgs() << "Debug Info - Variable: " << VarName
        //                           << ", File: " << Var->getFile()->getFilename()
        //                           << ", Line: " << Var->getLine() << "\n");
        //     }
        // }

        // LLVM_DEBUG(dbgs() << "Argument: " << VarName << ", Type: ";
        //            Arg->getType()->print(dbgs());
        //            dbgs() << "\n");
                   
          if (Inst->getOperand(0) == Arg) {
            if (isGEPStructTy) {
              // Insert prolog
              Builder.CreateIntrinsic(Intrinsic::venus_scalar_bind,
                                      {}, {Inst->getOperand(0)});
//              llvm::Value *add = Builder.CreateAdd(
//                  Builder.CreatePtrToInt(Inst->getOperand(0), Builder.getIntNTy(32)),
//                  Builder.getInt32(InputMap[Arg].first ?
//                                                       InputMap[Arg].second :
//                                                       StructInputAddress));
              llvm::Value *add = Builder.getInt32(InputMap[Arg].first ?
                        InputMap[Arg].second : StructInputAddress);
              Inst->setOperand(0, Builder.CreateIntToPtr(add, Builder.getPtrTy()));
              if (!InputMap[Arg].first) {
                InputMap[Arg].first = true;
                InputMap[Arg].second = StructInputAddress;
                // TODO: 64 should be parameterized. It is the bus width
                StructInputAddress += alignTo(DataSize, 64);
              }
            } else {
              report_fatal_error("Scalar part must be packed into a struct");
            }
          }
        }
      }

      // Check whether the instruction is an intrinsic
      if (IntrinsicInst* intrinsicInst = dyn_cast<IntrinsicInst>(I)) {
        Intrinsic::ID intrinsicID = intrinsicInst->getIntrinsicID();
        VenusSplitIntrinsic(Inst, intrinsicID, NeedSink);
      } // check intrinsic
      else if (Inst->getOpcode() == Instruction::ExtractValue ){
        VenusSplitExtractValue(Inst);
      } // check intrinsic
    } // for Instruction
  } // for BasicBlock

  // Split PHI Node
  LLVM_DEBUG(dbgs() << "********** SPLIT PHI NODES **********\n");
  for (Function::iterator BB = F.begin(); BB != F.end(); BB++) {
    for (BasicBlock::iterator I = BB->begin(); I != BB->end(); I++) {
      Instruction *Inst = dyn_cast<Instruction>(I);
      Value *Opd = Inst;
      Type *OpdType = Opd->getType();
      IRBuilder<> Builder(Inst);

      // phi nodes should also be checked
      if (Inst->getOpcode() == Instruction::PHI
          && OpdType->isVectorTy()) {
        VenusSplitPHI(Inst);
      } // Check PHI
    } // for Instruction
  } // for BasicBlock

  // Scan the remaining unsplatted intrinsics
  LLVM_DEBUG(dbgs() << "********** SPLIT THE REMAINING NODES **********\n");
  for (Function::iterator BB = F.begin(); BB != F.end(); BB++) {
    for (BasicBlock::iterator I = BB->begin(); I != BB->end(); I++) {
      Instruction *Inst = dyn_cast<Instruction>(I);
      Value *Opd = Inst;
      bool NeedSink = isNeedSink(Opd);
      if (Inst->getOpcode() == Instruction::Ret) continue;

      // Check whether the instruction is an intrinsic
      if (IntrinsicInst* intrinsicInst = dyn_cast<IntrinsicInst>(I)) {
        Intrinsic::ID intrinsicID = intrinsicInst->getIntrinsicID();
        VenusSplitIntrinsic(Inst, intrinsicID, NeedSink);
      } // check intrinsic
      else if (Inst->getOpcode() == Instruction::ExtractValue ){
        VenusSplitExtractValue(Inst);
      } // check intrinsic
    } // for Instruction
  } // for BasicBlock

  LLVM_DEBUG(dbgs() << "********** BREAK PHI LOOPS **********\n");
  for (Function::iterator BB = F.begin(); BB != F.end(); BB++) {
    for (BasicBlock::iterator I = BB->begin(); I != BB->end(); I++) {
      Instruction *Inst = dyn_cast<Instruction>(I);
      Value *Opd = Inst;
      Type *OpdType = Opd->getType();
      IRBuilder<> Builder(Inst);

      // This time we should break the loop, to insert an undef value
      if (Inst->getOpcode() == Instruction::PHI
          && OpdType->isVectorTy()) {
        VenusSplitPHI(Inst, false);
      } // Check PHI
    } // for Instruction
  } // for BasicBlock


  LLVM_DEBUG(dbgs() << "********** CLOSE PHI LOOPS **********\n");
  for (Function::iterator BB = F.begin(); BB != F.end(); BB++) {
    for (BasicBlock::iterator I = BB->begin(); I != BB->end(); I++) {
      Instruction *Inst = dyn_cast<Instruction>(I);
      Value *Opd = Inst;
      Type *OpdType = Opd->getType();
      StructType *OpdStructType = dyn_cast<StructType>(OpdType);
      if (OpdStructType && OpdStructType->getNumElements() == 2) {
          OpdType = OpdStructType->getElementType(0);
      }
      bool NeedSink = isNeedSink(Opd);
      if (Inst->getOpcode() == Instruction::Ret) continue;
      // Check whether the instruction is an intrinsic
      if (IntrinsicInst* intrinsicInst = dyn_cast<IntrinsicInst>(I)) {
        Intrinsic::ID intrinsicID = intrinsicInst->getIntrinsicID();
        VenusSplitIntrinsic(Inst, intrinsicID, NeedSink);
        FindAndReplacePHIUndef(Opd);
      }
      else if (Inst->getOpcode() == Instruction::ExtractValue ){
        VenusSplitExtractValue(Inst);
        FindAndReplacePHIUndef(Opd);
      } // check intrinsic

      // the dest value of a phi node can also be the source of another phi node
      // By this time, all the phi nodes should be splatted
      if (Inst->getOpcode() == Instruction::PHI
          && OpdType->isVectorTy()) {
        if (isNeedSplit(Opd))
          continue;
        FindAndReplacePHIUndef(Opd);
      }
    } // for Instruction
  } // for BasicBlock

  // Remove old intrinsics bottom-up
  LLVM_DEBUG(dbgs() << "********** REMOVE SPLATTED NODES **********\n");
  while (!VListIntrinsics.empty()) {
    Instruction *I = VListIntrinsics.back();
    LLVM_DEBUG(dbgs() << "Erasing: \n"; I->dump());
    Value *IRet = dyn_cast<Value>(I);
    IRet->replaceAllUsesWith(UndefValue::get(I->getType()));
    I->eraseFromParent();
    VListIntrinsics.pop_back();
  }

  return PreservedAnalyses::all();
}
