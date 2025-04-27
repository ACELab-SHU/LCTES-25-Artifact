//
// Created by xusiyi on 24-1-11.
//

#include "RISCV.h"
#include "RISCVSubtarget.h"
#include "RISCVInstrInfo.h"
#include "llvm/CodeGen/LiveIntervals.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include <queue>
#include <stack>
using namespace llvm;

#define DEBUG_TYPE "vns-bind-remove"
#define RISCV_VENUS_BIND_REMOVE \
  "Venus WiPU Extension Post Register Allocation Instruction Remove Assembly Instruction Vns_Bind Pass"

namespace {


class RISCVVenusBindRemove : public MachineFunctionPass {

  const TargetInstrInfo *TII;
  MachineRegisterInfo *MRI;

public:
  static char ID;

  RISCVVenusBindRemove() : MachineFunctionPass(ID) {
    initializeRISCVVenusBindRemovePass(*PassRegistry::getPassRegistry());
  }

  bool runOnMachineFunction(MachineFunction &MF) override;

  StringRef getPassName() const override { return RISCV_VENUS_BIND_REMOVE; }

private:

  void VnsBindRemove(MachineBasicBlock &MBB);
  
};

} // end anonymous namespace

char RISCVVenusBindRemove::ID = 0;

INITIALIZE_PASS(RISCVVenusBindRemove, "riscv-venus-bind-remove", RISCV_VENUS_BIND_REMOVE,
                false, false)

bool Removed;

void RISCVVenusBindRemove::VnsBindRemove(MachineBasicBlock &MBB) {
  
  for (auto MII = MBB.begin(), E = MBB.end(); MII != E;) {
    MachineInstr &MI = *MII++;
    DebugLoc DL = MI.getDebugLoc();

    // if (MI.getOpcode() == RISCV::venus_scatter_misc) {  
    //   Removed = true;
    //   MI.eraseFromParent();
    // } 
  }
}

bool RISCVVenusBindRemove::runOnMachineFunction(MachineFunction &MF) {
  // Skip if the vector extension is not enabled.
  const RISCVSubtarget &ST = MF.getSubtarget<RISCVSubtarget>();
  if (!ST.hasExtVenus())
    return false;

  TII = ST.getInstrInfo();
  MRI = &MF.getRegInfo();

  Removed= false;

  for (MachineBasicBlock &MBB : MF){
    VnsBindRemove(MBB);
  }

  return Removed;
}

/// Returns an instance of the Insert EA pass.
FunctionPass *llvm::createRISCVVenusBindRemovePass() {
  return new RISCVVenusBindRemove();
}
