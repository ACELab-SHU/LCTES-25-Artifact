//
// Created by jianglimin on 23-8-23.
//

#ifndef LLVM_TRANSFORMS_UTILS_VENUSINTRINSPLIT_H
#define LLVM_TRANSFORMS_UTILS_VENUSINTRINSPLIT_H

#include "llvm/IR/PassManager.h"
#include <map>
#include <stack>
#include <utility>
#include <vector>
#include <queue>

#define DEBUG_TYPE "venus"

namespace llvm {

class VenusIntrinSplitPass : public PassInfoMixin<VenusIntrinSplitPass> {
public:
  // large vectors mapping into small subvectors
  std::map<Value *, std::vector<Instruction *>> VMap;
  // original intrinsics
  std::list<Instruction *> VListIntrinsics;

  void VenusSplitClaim(Instruction *Inst, Intrinsic::ID intrinsicID);

  void VenusSplitShuffle(Instruction *Inst, Intrinsic::ID intrinsicID);

  void VenusSplitShuffleTest(Instruction *Inst, Intrinsic::ID intrinsicID);

  void VenusSplitVaddr(Instruction *Inst, Intrinsic::ID intrinsicID);

  void VenusSplitNormal(Instruction *Inst, Intrinsic::ID intrinsicID);

  void VenusSplitPHI(Instruction *Inst, bool SkipUndef);

  void VenusSplitExtractValue(Instruction *Inst);
  
  void FindAndReplacePHIUndef(Value *Opd);

  void VenusSplitIntrinsic(Instruction *Inst,
                      Intrinsic::ID intrinsicID,
                      bool NeedSink);
  PreservedAnalyses run(Function &F, FunctionAnalysisManager &AM);
};

} // namespace llvm

#endif // LLVM_TRANSFORMS_UTILS_VENUSINTRINSPLIT_H
