//
// Created by xusiyi on 24-3-5.
//

#include "RISCV.h"
#include "RISCVSubtarget.h"
#include "RISCVInstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/Support/CommandLine.h"
#include <fstream> 
#include <sstream> 

using namespace llvm;

#define DEBUG_TYPE "venus_alu"
#define RISCV_VENUS_MERGE_NAME \
  "Venus Alu Analysis Pass"

namespace {

class RISCVVenusAlu : public MachineFunctionPass {
public:
  static char ID;
  RISCVVenusAlu() : MachineFunctionPass(ID) {
    initializeRISCVVenusAluPass(*PassRegistry::getPassRegistry());
  }

  bool runOnMachineFunction(MachineFunction &MF) override;

};

char RISCVVenusAlu::ID = 0;

bool RISCVVenusAlu::runOnMachineFunction(MachineFunction &MF) {
  std::string functionName = MF.getName().str(); // Get the function name

  std::map<std::string, int> resources {
      {"BitALU", 0},
      {"CAU", 0},
      {"SerDiv", 0}
      //,{"MAX", 256}
  };
  std::vector<std::string> bit_alu_keywords = {
      "venus_and", "venus_or", "venus_xor", "venus_brdcst", "venus_sll", "venus_srl", "venus_sra", 
      "venus_seq", "venus_sne", "venus_sltu", "venus_slt", "venus_sleu", "venus_sle", "venus_sgtu", "venus_sgt", "venus_vmnot"
  };
  std::vector<std::string> cau_keywords = {
      "venus_add", "venus_sadd", "venus_saddu", "venus_rsub", "venus_sub", "venus_ssub", "venus_ssubu", "venus_mul", "venus_mulh", "venus_mulhu", "venus_mulhsu", "venus_muladd", "venus_mulsub", "venus_addmul", "venus_submul", "venus_cmxmul"
  };
  std::vector<std::string> ser_div_keywords = {
      "venus_div", "venus_rem", "venus_divu", "venus_remu"
  };
  int BitALU = 0;
  int CAU = 0;
  int SerDiv = 0;
  const MCInstrInfo *MCII = MF.getSubtarget().getInstrInfo();

  for (MachineBasicBlock &MBB : MF) {
    for (MachineInstr &MI : MBB) {
      unsigned Opcode = MI.getOpcode();
      const MCInstrDesc &MCID = MI.getDesc();
      LLVM_DEBUG(dbgs() << "Name: " << MCII->getName(Opcode)<< "\n");
      
      for (const std::string &keyword : bit_alu_keywords) {
        if (StringRef(MCII->getName(Opcode)).startswith(keyword)) {
          BitALU++;
          break;
        }
      }
      for (const std::string &keyword : cau_keywords) {
        if (StringRef(MCII->getName(Opcode)).startswith(keyword)) {
          CAU++;
          break;
        }
      }
      for (const std::string &keyword : ser_div_keywords) {
        if (StringRef(MCII->getName(Opcode)).startswith(keyword)) {
          SerDiv++;
          break;
        }
      }
    }
  }
  resources["BitALU"] = BitALU;
  resources["CAU"] = CAU;
  resources["SerDiv"] = SerDiv;

  std::ostringstream output;
  output << "{";
  for (auto it = resources.begin(); it != resources.end(); ++it) {
      if (it != resources.begin())
          output << ", ";
      output << "\"" << it->first << "\": " << it->second;
  }
  output << "}\n";

  std::ofstream outputFile(functionName + "_resource.json");
  outputFile << output.str();
  outputFile.close();
  
  return false;
  }
} // end of anonymous namespace

INITIALIZE_PASS(RISCVVenusAlu, "riscv-venus-alu",
                "Venus Alu Analysis Pass", false, false)

namespace llvm {

FunctionPass *createRISCVVenusAluPass() {
  return new RISCVVenusAlu();
}

} // end of namespace llvm
