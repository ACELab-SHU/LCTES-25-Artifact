//
// Created by jianglimin on 24-1-23.
//

#include "RISCV.h"
#include "RISCVSubtarget.h"
#include "RISCVInstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/Support/CommandLine.h"
#include <map>
#include <queue>
#include <vector>

using namespace llvm;

#define DEBUG_TYPE "venus"
#define RISCV_VENUS_REG_SHIFT_NAME \
  "Venus WiPU Extension Post Register Allocation Register Shifting Pass"

#define NUM_VENUS_INTRINSIC_OPERANDS 10

namespace opts {
extern llvm::cl::opt<unsigned> VenusNrRows;
extern llvm::cl::opt<unsigned> VenusNrLanes;
} // namespace opts

namespace {

class RISCVVenusRegShift : public MachineFunctionPass {
public:
  static char ID;

  RISCVVenusRegShift() : MachineFunctionPass(ID) {
    initializeRISCVVenusRegShiftPass(*PassRegistry::getPassRegistry());
  }

  bool runOnMachineFunction(MachineFunction &MF) override;

  StringRef getPassName() const override {
    return RISCV_VENUS_REG_SHIFT_NAME;
  }
};

char RISCVVenusRegShift::ID = 0;

bool RISCVVenusRegShift::runOnMachineFunction(llvm::MachineFunction &MF) {
  assert(0 && "Obsolete and to be tested if need re-using");
  bool Modified = false;
  // Skip if the vector extension is not enabled.
  const RISCVSubtarget &ST = MF.getSubtarget<RISCVSubtarget>();
  if (!ST.hasExtVenus()) return Modified;

  // A map recording the unordered and ordered vector registers
  std::map<std::vector<Register>, std::vector<Register>> RegRegMap;
  std::map<std::vector<Register>, std::vector<Register>>::iterator RegRegMapIter;
  Register BaseRegister = RISCV::VNS0;

  // Temporarily stores venus operands of the intrinsics
  std::vector<Register> VenusRegVec[NUM_VENUS_INTRINSIC_OPERANDS];

  // A flag indicates that whether the current intrinsic is between the
  // delimiter intrinsics
  bool IsBetweenDelimit;
  unsigned NumConsecutive;

  //
  unsigned VenusRegisterCount = 0;

//  const TargetInstrInfo &TII = *MF.getSubtarget().getInstrInfo();

  for (MachineFunction::iterator MBBI = MF.begin();
       MBBI != MF.end(); MBBI++) {
//    MachineBasicBlock &MBB = *MBBI;
    for (MachineBasicBlock::iterator MII = MBBI->begin();
         MII != MBBI->end(); MII++) {
      DebugLoc const DL = MII->getDebugLoc();
      assert(MII->getNumOperands() <= NUM_VENUS_INTRINSIC_OPERANDS &&
             "Change the macro the avoid the assertion.");

      if (MII->getOpcode() == RISCV::venus_delimit) {
        IsBetweenDelimit = !(MII->getOperand(0).getImm());
        if (IsBetweenDelimit) {
          NumConsecutive = 0;
          for (unsigned I = 0; I < NUM_VENUS_INTRINSIC_OPERANDS; I++)
            VenusRegVec[I].clear();
        }
        else {
          unsigned const CurrentVenusNumInsts = VenusRegVec[0].size();
          MachineBasicBlock::iterator MBBIter = MII;

          for (unsigned count = CurrentVenusNumInsts; count != 0; count--) {
            // Bump pointer
            MBBIter--;
            MachineInstr &CurrentMI = *MBBIter;
            for (unsigned Ops = 0; Ops < CurrentMI.getNumOperands(); Ops++) {
              if (CurrentMI.getOperand(Ops).getType() == MachineOperand::MO_Register) {
                if (CurrentMI.getOperand(Ops).getReg() >= RISCV::VNS0 &&
                    CurrentMI.getOperand(Ops).getReg() <  RISCV::VNS0 + opts::VenusNrRows) {
                  RegRegMapIter = RegRegMap.find(VenusRegVec[Ops]);
                  if (RegRegMapIter == RegRegMap.end()) {
                    std::vector<Register> NewVenusOps;
                    for (unsigned I = 0; I < CurrentVenusNumInsts; I++) {
                      NewVenusOps.push_back(BaseRegister + VenusRegisterCount);
                      VenusRegisterCount++;
                    }
                    RegRegMap.insert(
                        std::pair<std::vector<Register>, std::vector<Register>>(
                            VenusRegVec[Ops], NewVenusOps));
                    RegRegMapIter = RegRegMap.find(VenusRegVec[Ops]);
                  }
                  std::vector<Register> const NewVenusOps =
                      RegRegMapIter->second;
                  CurrentMI.getOperand(Ops).setReg(NewVenusOps.at(count - 1));
                }
              } // if is a register
            } // for operands in intrinsics
          } // for intrinsic count
        } // !IsBetweenDelimit
        continue;
      }

      if (IsBetweenDelimit) {
        for (unsigned Ops = 0; Ops < MII->getNumOperands(); Ops++) {
          if (MII->getOperand(Ops).getType() == MachineOperand::MO_Register) {
            if (MII->getOperand(Ops).getReg() >= RISCV::VNS0 &&
                MII->getOperand(Ops).getReg() <  RISCV::VNS0 + opts::VenusNrRows) {
              VenusRegVec[Ops].push_back(MII->getOperand(Ops).getReg());
            }
          }
        } // for Ops
        NumConsecutive++;
      } // IsBetweenDelimit


    }
  }

  return Modified;
}

} // end of anonymous namespace

INITIALIZE_PASS(RISCVVenusRegShift, "riscv-venus-reg-shift",
                RISCV_VENUS_REG_SHIFT_NAME, false, false)

namespace llvm {

FunctionPass *createRISCVVenusRegShiftPass() {
  return new RISCVVenusRegShift();
}

} // end of namespace llvm