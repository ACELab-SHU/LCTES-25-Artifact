//
// Created by xusiyi on 24-1-5.
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

#define DEBUG_TYPE "venus"
#define RISCV_VENUS_INST_MODIFY_LOAD_STORE \
  "Venus WiPU Extension Post Register Allocation Instruction Modify Load and Store Pass"

namespace opts {
cl::opt<unsigned> VenusTileBaseAddr(
    "venus-tile-baseaddr",
    cl::desc("Specify the base address of a tile of the Venus architecture."),
    cl::value_desc("value"),
    cl::init(0x80000000)
);

cl::opt<unsigned> VenusVRFBaseAddr(
    "venus-vrf-baseaddr",
    cl::desc("Specify the base address of the vector register files of the Venus architecture."),
    cl::value_desc("value"),
    cl::init(0x80100000)
);

cl::opt<unsigned> VenusCRBaseAddr(
    "venus-cr-baseaddr",
    cl::desc("Specify the base address of the control registers of the Venus architecture."),
    cl::value_desc("value"),
    cl::init(0x801ff000)
);

cl::opt<std::string> VenusTaskName(
    "venus-task-name",
    cl::desc("Specify a name of the task. This is for generating a hash."),
    cl::value_desc("taskname"),
    cl::init("")
);
} // namespace opts

namespace {


class RISCVVenusInstModifyLoadStore : public MachineFunctionPass {
public:
  static char ID;

  RISCVVenusInstModifyLoadStore() : MachineFunctionPass(ID) {
    initializeRISCVVenusInstModifyLoadStorePass(*PassRegistry::getPassRegistry());
  }

  bool runOnMachineFunction(MachineFunction &MF) override;

  StringRef getPassName() const override {
    return RISCV_VENUS_INST_MODIFY_LOAD_STORE;
  }
};

char RISCVVenusInstModifyLoadStore::ID = 0;

bool RISCVVenusInstModifyLoadStore::runOnMachineFunction(MachineFunction &MF) {
  bool Modified = false;
  // Skip if the vector extension is not enabled.
  const RISCVSubtarget &ST = MF.getSubtarget<RISCVSubtarget>();
  if (!ST.hasExtVenus())
    return Modified;

  const TargetInstrInfo &TII = *MF.getSubtarget().getInstrInfo();
  std::stack<MachineInstr *> VStackIntrinsics;
  LLVM_DEBUG(dbgs() << "VIMLS: VRF base address: "
                    << format_hex(opts::VenusVRFBaseAddr, 10) << "\n");

  for (MachineFunction::iterator MBBI = MF.begin();
       MBBI != MF.end(); MBBI++) {
    MachineBasicBlock &MBB = *MBBI;
    for (MachineBasicBlock::iterator MII = MBBI->begin();
         MII != MBBI->end(); MII++) {
      DebugLoc const DL = MII->getDebugLoc();

      if (MII->getOpcode() == RISCV::venus_vaddr_misc) {
        Modified = true;
        MachineInstr &MI = *MII;
        VStackIntrinsics.push(&MI);
        LLVM_DEBUG(dbgs() << "VIMLS: Venus Vaddr Detected.\n");

        const unsigned PosReg = MII->getOperand(0).getReg();
        const unsigned VecReg = MII->getOperand(1).getReg();
        const unsigned VRFaddr = opts::VenusVRFBaseAddr +
                                 (VecReg - RISCV::VNS0) * opts::VenusNrLanes * 8 /*bytes*/;
        LLVM_DEBUG(dbgs() << "VIMLS: VRF address: "
                          << format_hex(VRFaddr, 10) << "\n");

        // Check whether the lower 12 bit is larger or equal than 0x800,
        // since the immediate number of addi is a signed number
        if (VRFaddr & 0x800) {
          // load upper immediate
          BuildMI(MBB, MII, DL, TII.get(RISCV::LUI), PosReg)
              .addImm((VRFaddr >> 12) + 1);
          // load immediate
          if (VRFaddr & 0xFFF)
            BuildMI(MBB, MII, DL, TII.get(RISCV::ADDI), PosReg)
                .addReg(PosReg).addImm((signed)((VRFaddr & 0xFFF) - 0x1000));
        } else {
          // load upper immediate
          BuildMI(MBB, MII, DL, TII.get(RISCV::LUI), PosReg)
              .addImm(VRFaddr >> 12);
          // load immediate
          if (VRFaddr & 0xFFF)
            BuildMI(MBB, MII, DL, TII.get(RISCV::ADDI), PosReg)
                .addReg(PosReg)
                .addImm(VRFaddr & 0xFFF);
        }

      } else if (MII->getOpcode() == RISCV::venus_claim_misc ||
                 MII->getOpcode() == RISCV::venus_sink ||
                 MII->getOpcode() == RISCV::venus_delimit ||
                 MII->getOpcode() == RISCV::venus_pseudo_ivv ||
                 MII->getOpcode() == RISCV::venus_pseudo_ivx) {
        Modified = true;
        MachineInstr &MI = *MII;
        VStackIntrinsics.push(&MI);
        LLVM_DEBUG(dbgs() << "VIMLS: Venus Miscs Detected: " <<
                   MII->getOpcode() << "\n");
      }
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

} // end anonymous namespace

INITIALIZE_PASS(RISCVVenusInstModifyLoadStore,
                "riscv-venus-inst-modify-load-store",
                RISCV_VENUS_INST_MODIFY_LOAD_STORE,
                false, false)

namespace llvm {

FunctionPass *createRISCVVenusInstModifyLoadStorePass() {
  return new RISCVVenusInstModifyLoadStore();
}

} // end of namespace llvm