//===-- RISCV.h - Top-level interface for RISCV -----------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file contains the entry points for global functions defined in the LLVM
// RISC-V back-end.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_LIB_TARGET_RISCV_RISCV_H
#define LLVM_LIB_TARGET_RISCV_RISCV_H

#include "MCTargetDesc/RISCVBaseInfo.h"
#include "llvm/Target/TargetMachine.h"

namespace llvm {
class RISCVRegisterBankInfo;
class RISCVSubtarget;
class RISCVTargetMachine;
class AsmPrinter;
class FunctionPass;
class InstructionSelector;
class MCInst;
class MCOperand;
class MachineInstr;
class MachineOperand;
class PassRegistry;

FunctionPass *createRISCVCodeGenPreparePass();
void initializeRISCVCodeGenPreparePass(PassRegistry &);

bool lowerRISCVMachineInstrToMCInst(const MachineInstr *MI, MCInst &OutMI,
                                    AsmPrinter &AP);
bool lowerRISCVMachineOperandToMCOperand(const MachineOperand &MO,
                                         MCOperand &MCOp, const AsmPrinter &AP);

FunctionPass *createRISCVISelDag(RISCVTargetMachine &TM,
                                 CodeGenOpt::Level OptLevel);

FunctionPass *createRISCVMakeCompressibleOptPass();
void initializeRISCVMakeCompressibleOptPass(PassRegistry &);

FunctionPass *createRISCVGatherScatterLoweringPass();
void initializeRISCVGatherScatterLoweringPass(PassRegistry &);

FunctionPass *createRISCVSExtWRemovalPass();
void initializeRISCVSExtWRemovalPass(PassRegistry &);

FunctionPass *createRISCVMergeBaseOffsetOptPass();
void initializeRISCVMergeBaseOffsetOptPass(PassRegistry &);

FunctionPass *createRISCVExpandPseudoPass();
void initializeRISCVExpandPseudoPass(PassRegistry &);

FunctionPass *createRISCVExpandAtomicPseudoPass();
void initializeRISCVExpandAtomicPseudoPass(PassRegistry &);

FunctionPass *createRISCVInsertVSETVLIPass();
void initializeRISCVInsertVSETVLIPass(PassRegistry &);

FunctionPass *createRISCVRedundantCopyEliminationPass();
void initializeRISCVRedundantCopyEliminationPass(PassRegistry &);

FunctionPass *createRISCVVenusMergePass();
void initializeRISCVVenusMergePass(PassRegistry &);

FunctionPass *createRISCVVenusAluPass();
void initializeRISCVVenusAluPass(PassRegistry &);

FunctionPass *createRISCVVenusInstModifyLoadStorePass();
void initializeRISCVVenusInstModifyLoadStorePass(PassRegistry &);

FunctionPass *createRISCVVenusBindRemovePass();
void initializeRISCVVenusBindRemovePass(PassRegistry &);

FunctionPass *createRISCVVenusRegShiftPass();
void initializeRISCVVenusRegShiftPass(PassRegistry &);

InstructionSelector *createRISCVInstructionSelector(const RISCVTargetMachine &,
                                                    RISCVSubtarget &,
                                                    RISCVRegisterBankInfo &);
}

#endif
