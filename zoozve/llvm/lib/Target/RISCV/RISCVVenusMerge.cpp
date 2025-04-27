//
// Created by jianglimin on 23-8-31.
//

#include "RISCV.h"
#include "RISCVSubtarget.h"
#include "RISCVInstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/CompareNat.h"
#include <stack>


using namespace llvm;

// The vector register to be merged is incremental
#define REG_CONSEC_INCR 1
#define GET_REG(id) \
  ((REG_CONSEC_INCR) ? (CurrentInfo.Reg.at(id) + 1 - NumConsecutive) : \
                     (CurrentInfo.Reg.at(id)))

#define DEBUG_TYPE "venus"
#define RISCV_VENUS_MERGE_NAME \
  "Venus WiPU Extension Post Register Allocation Instruction Merging Pass"

namespace opts {
extern llvm::cl::opt<unsigned> VenusNrRows;
extern llvm::cl::opt<unsigned> VenusNrLanes;
} // namespace opts

namespace {

class RISCVVenusMerge : public MachineFunctionPass {
public:
  static char ID;
  static unsigned VenusIntrinsicFirst;
  static unsigned VenusIntrinsicLast;
  static unsigned VenusRegisterFirst;
  static unsigned VenusRegisterLast;

  struct VenusIntrinsicInfo{
    unsigned Opcode;            // Name of the intrinsic
    unsigned NumOp;             // Num of operands of the intrinsic
    unsigned NumDataOp;         // Num of data operands (exclude imm operands)
    unsigned Vew;               // Vector element width
    unsigned Vmask_read;        // Whether this intrinsic needs to read mask reg
    unsigned SubID;             // Sub-Intrinsic ID
    unsigned AVL;               // Application vector length
    unsigned Rows;              // How many consecutive rows in middle-end
    bool isVXInst;              // Whether this instruction has a scalar operand
    std::vector<unsigned> Reg;  // Register ID of the operands
    VenusIntrinsicInfo() = default;
  };

  RISCVVenusMerge() : MachineFunctionPass(ID) {
    initializeRISCVVenusMergePass(*PassRegistry::getPassRegistry());
  }

  bool runOnMachineFunction(MachineFunction &MF) override;

  static bool isVenusRegister(unsigned Reg);

  static bool isNormalVenusMI(unsigned Opcode);

  static bool isMaskWriteVenusMI(unsigned Opcode);

  VenusIntrinsicInfo getVenusIntrinsicInfo(MachineInstr* MI);

  VenusIntrinsicInfo getVenusBroadcastInfo(MachineInstr* MI);

  VenusIntrinsicInfo getVenusRangeInfo(MachineInstr* MI);

  VenusIntrinsicInfo getVenusCmxmulInfo(MachineInstr* MI);

  VenusIntrinsicInfo getVenusScatterInfo(MachineInstr* MI);

  StringRef getPassName() const override {
    return RISCV_VENUS_MERGE_NAME;
  }
};

char RISCVVenusMerge::ID = 0;
unsigned RISCVVenusMerge::VenusIntrinsicFirst = RISCV::venus_add_ivv;
unsigned RISCVVenusMerge::VenusIntrinsicLast = RISCV::venus_xor_ivx;
unsigned RISCVVenusMerge::VenusRegisterFirst = RISCV::VNS0;
unsigned RISCVVenusMerge::VenusRegisterLast = RISCV::VNS2047;

inline bool isThreeSourceOperands(unsigned OpNum) {
  return (OpNum == 8);
}

RISCVVenusMerge::VenusIntrinsicInfo
RISCVVenusMerge::getVenusIntrinsicInfo(MachineInstr* MI) {
  VenusIntrinsicInfo IntrinsicInfo;
  // return an empty info if current machine instruction is not an intrinsic
  if ((MI->getOpcode() == RISCV::venus_brdcst_ivx)) return getVenusBroadcastInfo(MI);
  if ((MI->getOpcode() == RISCV::venus_range_misc)) return getVenusRangeInfo(MI);
  if ((MI->getOpcode() == RISCV::venus_scatter)) return getVenusScatterInfo(MI);
  if ((MI->getOpcode() == RISCV::venus_cmxmul_ivv)) return getVenusCmxmulInfo(MI);   
  if (!isNormalVenusMI(MI->getOpcode())) return IntrinsicInfo;

  const unsigned is3Ops = isThreeSourceOperands(MI->getNumOperands());
  const unsigned isMW = isMaskWriteVenusMI(MI->getOpcode());
  const unsigned VewOperandNo = 3 - isMW + is3Ops;
  const unsigned VmaskReadOperandNo = 4 - isMW + is3Ops;
  const unsigned AVLOperandNo = 5 - isMW + is3Ops;
  const unsigned RowsOperandNo = 6 - isMW + is3Ops;
  IntrinsicInfo.Opcode = MI->getOpcode();
  IntrinsicInfo.NumOp = MI->getNumOperands();
  IntrinsicInfo.NumDataOp = IntrinsicInfo.NumOp - 4;
  IntrinsicInfo.Vew = MI->getOperand(VewOperandNo).getImm();
  IntrinsicInfo.Vmask_read = MI->getOperand(VmaskReadOperandNo).getImm();
  IntrinsicInfo.AVL = MI->getOperand(AVLOperandNo).getReg();
  IntrinsicInfo.Rows = MI->getOperand(RowsOperandNo).getImm();
  for (unsigned I = 0; I < IntrinsicInfo.NumDataOp; I++) {
    IntrinsicInfo.Reg.push_back(MI->getOperand(I).getReg());
    LLVM_DEBUG(dbgs() << "Opcode: " << MI->getOpcode()
                      << ", Operand: " << MI->getOperand(I).getReg() << "\n");
  }
  IntrinsicInfo.isVXInst = !isVenusRegister(MI->getOperand(2 - isMW).getReg());

  return IntrinsicInfo;
}

RISCVVenusMerge::VenusIntrinsicInfo
RISCVVenusMerge::getVenusBroadcastInfo(MachineInstr *MI) {
  VenusIntrinsicInfo IntrinsicInfo;

  IntrinsicInfo.Opcode = MI->getOpcode();
  IntrinsicInfo.NumOp = MI->getNumOperands();
  IntrinsicInfo.NumDataOp = 2;
  IntrinsicInfo.Vew = MI->getOperand(3).getImm();
  IntrinsicInfo.Vmask_read = MI->getOperand(4).getImm();
  IntrinsicInfo.AVL = MI->getOperand(5).getReg();
  IntrinsicInfo.Rows = MI->getOperand(6).getImm();
  IntrinsicInfo.Reg.push_back(MI->getOperand(0).getReg());
  IntrinsicInfo.Reg.push_back(MI->getOperand(2).getReg());
  IntrinsicInfo.isVXInst = true;
  LLVM_DEBUG(dbgs() << "Opcode: " << MI->getOpcode()
                    << " is a venus_brdcst.\n");

  return IntrinsicInfo;
}

RISCVVenusMerge::VenusIntrinsicInfo
RISCVVenusMerge::getVenusRangeInfo(MachineInstr *MI) {
  VenusIntrinsicInfo IntrinsicInfo;

  IntrinsicInfo.Opcode = MI->getOpcode();
  IntrinsicInfo.NumOp = MI->getNumOperands();
  IntrinsicInfo.NumDataOp = 1;
  IntrinsicInfo.Vew = MI->getOperand(1).getImm();
  IntrinsicInfo.Vmask_read = MI->getOperand(2).getImm();
  IntrinsicInfo.SubID = MI->getOperand(3).getImm();
  IntrinsicInfo.AVL = MI->getOperand(4).getReg();
  IntrinsicInfo.Rows = MI->getOperand(5).getImm();
  IntrinsicInfo.Reg.push_back(MI->getOperand(0).getReg());
  IntrinsicInfo.isVXInst = false;
  LLVM_DEBUG(dbgs() << "Opcode: " << MI->getOpcode()
                    << " is a venus_range.\n");

  return IntrinsicInfo;
}

RISCVVenusMerge::VenusIntrinsicInfo
RISCVVenusMerge::getVenusCmxmulInfo(MachineInstr *MI) {
  VenusIntrinsicInfo IntrinsicInfo;

  IntrinsicInfo.Opcode = MI->getOpcode();
  IntrinsicInfo.NumOp = MI->getNumOperands();
  IntrinsicInfo.NumDataOp = 6;

  for (unsigned I = 0; I < IntrinsicInfo.NumDataOp; I++) {
    IntrinsicInfo.Reg.push_back(MI->getOperand(I).getReg());
    LLVM_DEBUG(dbgs() << "CmxmulOpcode: " << MI->getOpcode()
                      << ", Operand"<< I <<": " << MI->getOperand(I).getReg() << "\n");
  }
  
  IntrinsicInfo.Vew = MI->getOperand(6).getImm();
  IntrinsicInfo.Vmask_read = MI->getOperand(7).getImm();
  IntrinsicInfo.AVL = MI->getOperand(8).getReg();    
  IntrinsicInfo.Rows = MI->getOperand(9).getImm();
  
  IntrinsicInfo.isVXInst = false;
  LLVM_DEBUG(dbgs() << "Opcode: " << MI->getOpcode()
                    << " is a venus_cmxmul.\n");

  return IntrinsicInfo;
}

RISCVVenusMerge::VenusIntrinsicInfo
RISCVVenusMerge::getVenusScatterInfo(llvm::MachineInstr *MI) {
  VenusIntrinsicInfo IntrinsicInfo;

  IntrinsicInfo.Opcode = MI->getOpcode();
  IntrinsicInfo.Rows = MI->getOperand(2).getImm();
  IntrinsicInfo.isVXInst = false;
  LLVM_DEBUG(dbgs() << "Venus Scatter.\n");

  return IntrinsicInfo;
}


inline bool RISCVVenusMerge::isVenusRegister(unsigned Reg) {
  return (Reg >= VenusRegisterFirst && Reg <= VenusRegisterLast);
}

inline bool RISCVVenusMerge::isNormalVenusMI(unsigned Opcode) {
  const bool inRange =  Opcode >= VenusIntrinsicFirst &&
                        Opcode <= VenusIntrinsicLast;
  const bool isScatter = Opcode == RISCV::venus_scatter;
  const bool isSpecial = Opcode == RISCV::venus_bind ||
                         Opcode == RISCV::venus_scalar_bind ||
                         Opcode == RISCV::venus_mnot_misc ||
                         Opcode == RISCV::venus_claim_misc ||
                         Opcode == RISCV::venus_sink ||
                         Opcode == RISCV::venus_delimit ||
                         Opcode == RISCV::venus_barrier ||
                         Opcode == RISCV::venus_setcsr;
  return ((inRange | isScatter) & (!isSpecial));
}

inline bool RISCVVenusMerge::isMaskWriteVenusMI(unsigned Opcode) {
#define VENUS_BUILTIN(ID, TYPE, ATTRS, FEATURE, INTRINSIC, MASK_WRITE) \
  if (MASK_WRITE) {                                                    \
    if (Opcode >= RISCV::INTRINSIC + (1 << 1) &&                       \
        Opcode <  RISCV::INTRINSIC + (1 << 2)) {                       \
      return true;                                                     \
    }                                                                  \
  }
#include "../../clang/include/clang/Basic/BuiltinsRISCVVenus.def"
  return false;
}

/// Check whether the current and the upcoming intrinsic can be combined.
/// The opcode and the number of the operands should be identical.
/// The register ID of the consecutive intrinsics should be continuous if it is
/// a vector register, otherwise the ID should remain the same.
/// There is an exception that when it comes a scatter shuffle intrinsic
bool isConsecutiveIntrinsic(RISCVVenusMerge::VenusIntrinsicInfo current,
                            RISCVVenusMerge::VenusIntrinsicInfo next) {
  if ((current.Opcode == RISCV::venus_shuffle_ivv ||
       current.Opcode == RISCV::venus_shuffle_ivx ||
       current.Opcode == RISCV::venus_scatter)
      && next.Opcode == RISCV::venus_scatter) {
    return true;
  }

  if (current.Opcode != next.Opcode) return false;
  if (current.NumOp != next.NumOp) return false;
  if (current.Opcode == RISCV::venus_shuffle_test) {
    for (unsigned I = 0; I < current.NumDataOp; I++) { 
      const bool isVenusReg = RISCVVenusMerge::isVenusRegister(current.Reg.at(I));
      if (!isVenusReg && next.Reg.at(I) != current.Reg.at(I)) {
          return false;
      }

      if (isVenusReg) {
          if (next.Reg.at(I) != RISCV::VNS0 && (next.Reg.at(I) - current.Reg.at(I) != 1)){
              return false;
          }
      }
  }
  return true;
}
  // TODO: start value of I is 0 or 1?
  for (unsigned I = 0; I < current.NumDataOp; I++) {
    const bool isVenusReg = RISCVVenusMerge::isVenusRegister(current.Reg.at(I));
    if (!isVenusReg && (next.Reg.at(I) != current.Reg.at(I)) )
      return false;
#if REG_CONSEC_INCR == 1
    if (isVenusReg && (next.Reg.at(I) - current.Reg.at(I) != 1))
#else
    if (isVenusReg && (current.Reg.at(I) - next.Reg.at(I) != 1))
#endif
      return false;
  }



  return true;
}

bool RISCVVenusMerge::runOnMachineFunction(MachineFunction &MF) {
  bool Modified = false;
  // Skip if the vector extension is not enabled.
  const RISCVSubtarget &ST = MF.getSubtarget<RISCVSubtarget>();
  if (!ST.hasExtVenus())
    return Modified;

  const TargetInstrInfo &TII = *MF.getSubtarget().getInstrInfo();
  std::stack<MachineInstr *> VStackIntrinsics;

  for (MachineFunction::iterator MBBI = MF.begin();
       MBBI != MF.end(); MBBI++) {
//    LLVM_DEBUG(dbgs() << "Contents of MachineBasicBlock:\n");
//    LLVM_DEBUG(dbgs() << *MBBI << "\n");
//    const BasicBlock *BB = MBBI->getBasicBlock();
    MachineBasicBlock &MBB = *MBBI;
//    LLVM_DEBUG(dbgs() << "Contents of BasicBlock corresponding to the MBB:\n");
//    LLVM_DEBUG(if(BB != NULL) dbgs() << *BB << "\n";
//               else dbgs() << "is NULL\n");
    for (MachineBasicBlock::iterator MII = MBBI->begin();
         MII != MBBI->end(); MII++) {
      LLVM_DEBUG(dbgs() << "Contents of MachineInstruction:\n");
      LLVM_DEBUG(dbgs() << *MII << "\n");
      DebugLoc const DL = MII->getDebugLoc();

      // Specify the range of the machine instructions to be merged
      if (isNormalVenusMI(MII->getOpcode())) {
        Modified = true;
        bool IsConsecutive;
        unsigned NumConsecutive = 0;
        VenusIntrinsicInfo CurrentInfo, NextInfo;
        MachineBasicBlock::iterator It = MII;
        do {
          NumConsecutive++;
          LLVM_DEBUG(dbgs() << "\nNumConsecutive: " << NumConsecutive << "\n");
          MachineInstr &CurrentMachInst = *It;
          VStackIntrinsics.push(&CurrentMachInst);
          LLVM_DEBUG(dbgs() << "CurrentInfo: \n");
          CurrentInfo = getVenusIntrinsicInfo(&CurrentMachInst);
          It++;
          if (It == MBBI->end()) break;

          // The redundant LUI & ADDI must be sandwiched by consecutive intrinsics
          std::vector<MachineInstr*> Land;
          Land.clear();

          while (true) {
            switch (It->getOpcode()) {
              case RISCV::ADDI: 
              case RISCV::LUI: 
              case RISCV::LW:
              case RISCV::COPY:
                Land.push_back(&*It); It++; break;
              default: goto end_loop;
            }
          }
end_loop:
          if (It == MBBI->end()) break;

          MachineInstr &NextMachInst = *It;
          LLVM_DEBUG(dbgs() << "NextInfo: \n");
          NextInfo = getVenusIntrinsicInfo(&NextMachInst);
          IsConsecutive = isConsecutiveIntrinsic(CurrentInfo, NextInfo);
          if (IsConsecutive) {
             if (!Land.empty()) {
              for (auto *Instr : Land) {
                LLVM_DEBUG(dbgs() << "Delete Sandwiched Insn: "; Instr->dump());
                VStackIntrinsics.push(Instr);
              }
            }
          }
        } while (IsConsecutive);
        LLVM_DEBUG(dbgs() << "Number of Consecutive Intrinsics: " <<
                   NumConsecutive << "\n");

        // Replace by a merged intrinsic
        MachineInstr *NewMI;
        if ((MII->getOpcode() == RISCV::venus_cmxmul_ivv)) {
          assert(CurrentInfo.Rows == NumConsecutive && "Wrong inst merge!");
          unsigned int vd2_value = GET_REG(2) - VenusRegisterFirst; // 40 = the num of gpr
          unsigned int vd1_value = GET_REG(3) - VenusRegisterFirst;
          unsigned int vs2_value = GET_REG(4) - VenusRegisterFirst;
          unsigned int vs1_value = GET_REG(5) - VenusRegisterFirst;
          unsigned int high_vd2_3bits = (vd2_value >> 7) & 0b111;
          unsigned int high_vd1_2bits = (vd1_value >> 8) & 0b11;
          unsigned int high_vs2_2bits = (vs2_value >> 8) & 0b11;
          unsigned int high_vs1_2bits = (vs1_value >> 8) & 0b11;
          
          unsigned int combined_bits = (high_vd2_3bits << 6) | (high_vd1_2bits << 4) | (high_vs2_2bits << 2) | high_vs1_2bits;
          NewMI = BuildMI(MBB, MII, DL,
                          TII.get(RISCV::venus_setcsrimm))
                      .addImm(combined_bits)
                      .addImm(1);
          Register DestRegs[] = {GET_REG(0), GET_REG(1)};
          NewMI = BuildMI(MBB, MII, DL,
                          TII.get(CurrentInfo.Opcode),DestRegs)
                      .addReg(GET_REG(2))
                      .addReg(GET_REG(3))
                      .addReg(GET_REG(4))
                      .addReg(GET_REG(5))
                      .addImm(CurrentInfo.Vew)
                      .addImm(CurrentInfo.Vmask_read)
                      .addReg(CurrentInfo.AVL)
                      .addImm(CurrentInfo.Rows);

          goto VENUS_MERGE_BUMP_ITER;
        }
        if ((MII->getOpcode() == RISCV::venus_brdcst_ivx)) {
          assert(CurrentInfo.Rows == NumConsecutive && "Wrong inst merge!");
          NewMI = BuildMI(MBB, MII, DL,
                          TII.get(CurrentInfo.Opcode),
                          GET_REG(0))
                      .addReg(GET_REG(0))
                      .addReg(CurrentInfo.Reg.at(1))
                      .addImm(CurrentInfo.Vew)
                      .addImm(CurrentInfo.Vmask_read)
                      .addReg(CurrentInfo.AVL)
                      .addImm(CurrentInfo.Rows);
          goto VENUS_MERGE_BUMP_ITER;
        }
        if ((MII->getOpcode() == RISCV::venus_range_misc)) {

          assert(CurrentInfo.Rows == NumConsecutive && "Wrong inst merge!");
          NewMI = BuildMI(MBB, MII, DL,
                          TII.get(CurrentInfo.Opcode),
                          GET_REG(0))
                      .addImm(CurrentInfo.Vew)
                      .addImm(CurrentInfo.Vmask_read)
                      .addImm(CurrentInfo.SubID)
                      .addReg(CurrentInfo.AVL)
                      .addImm(CurrentInfo.Rows);
          goto VENUS_MERGE_BUMP_ITER;
        }
        {
        const unsigned is3Ops = isThreeSourceOperands(MII->getNumOperands());
        const unsigned isMW = isMaskWriteVenusMI(MII->getOpcode());
        const unsigned DestOperandNo = 0;
        const unsigned SourceOperand1No = 1 - isMW;
        const unsigned SourceOperand2No = 2 - isMW;
        const bool isVenusReg =
            RISCVVenusMerge::isVenusRegister(CurrentInfo.Reg.at(SourceOperand2No));
        unsigned RegOp2;
        if (isVenusReg) RegOp2 = GET_REG(SourceOperand2No);
        else            RegOp2 = CurrentInfo.Reg.at(SourceOperand2No);
        LLVM_DEBUG(dbgs() << "CurrentInfo.Rows: " << CurrentInfo.Rows << "\n");
        LLVM_DEBUG(dbgs() << "NumConsecutive: " << NumConsecutive << "\n");

        // MVV/MVX instructions do not have a destination register
        if (isMW) {
          NewMI = BuildMI(MBB, MII, DL,
            TII.get(CurrentInfo.Opcode))
            .addReg(GET_REG(SourceOperand1No))
            .addReg(RegOp2)
            .addImm(CurrentInfo.Vew)
            .addImm(CurrentInfo.Vmask_read)
            .addReg(CurrentInfo.AVL)
            .addImm(CurrentInfo.Rows);
        } else {
          if (is3Ops) {
            assert(CurrentInfo.Rows == NumConsecutive && "Wrong inst merge!");
            if ((CurrentInfo.Opcode != RISCV::venus_shuffle_ivv) &&
                  (CurrentInfo.Opcode != RISCV::venus_shuffle_ivx)&&
                  (CurrentInfo.Opcode != RISCV::venus_shuffle_test)) {
                unsigned int vd2_value = GET_REG(3) - VenusRegisterFirst; // 40 = the num of gpr
                unsigned int vd1_value = GET_REG(DestOperandNo) - VenusRegisterFirst;
                unsigned int vs2_value = (isVenusReg) ? (RegOp2 - VenusRegisterFirst) : RegOp2;
                unsigned int vs1_value = GET_REG(SourceOperand1No) - VenusRegisterFirst;
                unsigned int high_vd2_3bits = (vd2_value >> 7) & 0b111;
                unsigned int high_vd1_2bits = (vd1_value >> 8) & 0b11;
                unsigned int high_vs2_2bits = (vs2_value >> 8) & 0b11;
                unsigned int high_vs1_2bits = (vs1_value >> 8) & 0b11;

                unsigned int combined_bits = (high_vd2_3bits << 6) | (high_vd1_2bits << 4) | (high_vs2_2bits << 2) | high_vs1_2bits;
                NewMI = BuildMI(MBB, MII, DL,
                                TII.get(RISCV::venus_setcsrimm))
                            .addImm(combined_bits)
                            .addImm(1);
            }
                 
            NewMI = BuildMI(MBB, MII, DL,
                            TII.get(CurrentInfo.Opcode),
                            GET_REG(DestOperandNo))
                        .addReg(GET_REG(SourceOperand1No))
                        .addReg(RegOp2)
                        .addReg(GET_REG(3))
                        .addImm(CurrentInfo.Vew)
                        .addImm(CurrentInfo.Vmask_read)
                        .addReg(CurrentInfo.AVL)
                        .addImm(CurrentInfo.Rows);
            LLVM_DEBUG(NewMI->dump());

            // Check whether shuffle need clobbering
            if (((CurrentInfo.Opcode == RISCV::venus_shuffle_ivx) ||
                 (CurrentInfo.Opcode == RISCV::venus_shuffle_ivv)) &&
                ((GET_REG(DestOperandNo) == GET_REG(SourceOperand1No)) ||
                 (GET_REG(DestOperandNo) == RegOp2))) {
              NewMI = BuildMI(MBB, MII, DL,
                              TII.get(RISCV::venus_shuffle_clobber_move),
                              GET_REG(DestOperandNo))
                          .addImm(CurrentInfo.Vew)
                          .addImm(CurrentInfo.Vmask_read)
                          .addImm(0)
                          .addReg(CurrentInfo.AVL)
                          .addImm(CurrentInfo.Rows);
            }
          } else {
            NewMI = BuildMI(MBB, MII, DL,
              TII.get(CurrentInfo.Opcode),
                            GET_REG(DestOperandNo))
              .addReg(GET_REG(SourceOperand1No))
              .addReg(RegOp2)
              .addImm(CurrentInfo.Vew)
              .addImm(CurrentInfo.Vmask_read)
              .addReg(CurrentInfo.AVL)
              .addImm(CurrentInfo.Rows);
          }
        }
        }
VENUS_MERGE_BUMP_ITER:
        LLVM_DEBUG(NewMI->dump());

        // Bump iterator

        MII = It;
      } // if [VenusIntrinsicFirst, VenusIntrinsicLast)
    } // for MachineBasicBlock::iterator
  } // for MachineFunction::iterator

  // Remove old intrinsics bottom-up
  while (!VStackIntrinsics.empty()) {
    MachineInstr *MI = VStackIntrinsics.top();
    MI->eraseFromParent();
    VStackIntrinsics.pop();
  }

  return Modified;
}

} // end of anonymous namespace

INITIALIZE_PASS(RISCVVenusMerge, "riscv-venus-merge",
                RISCV_VENUS_MERGE_NAME, false, false)

namespace llvm {

FunctionPass *createRISCVVenusMergePass() {
  return new RISCVVenusMerge();
}

} // end of namespace llvm