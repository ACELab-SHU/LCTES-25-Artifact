//===-- RISCVDisassembler.cpp - Disassembler for RISCV --------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file implements the RISCVDisassembler class.
//
//===----------------------------------------------------------------------===//

#include "MCTargetDesc/RISCVBaseInfo.h"
#include "MCTargetDesc/RISCVMCTargetDesc.h"
#include "TargetInfo/RISCVTargetInfo.h"
#include "llvm/MC/MCContext.h"
#include "llvm/MC/MCDecoderOps.h"
#include "llvm/MC/MCDisassembler/MCDisassembler.h"
#include "llvm/MC/MCInst.h"
#include "llvm/MC/MCInstrInfo.h"
#include "llvm/MC/MCRegisterInfo.h"
#include "llvm/MC/MCSubtargetInfo.h"
#include "llvm/MC/TargetRegistry.h"
#include "llvm/Support/Endian.h"

using namespace llvm;

#define DEBUG_TYPE "riscv-disassembler"

typedef MCDisassembler::DecodeStatus DecodeStatus;

namespace {
class RISCVDisassembler : public MCDisassembler {
  std::unique_ptr<MCInstrInfo const> const MCII;

public:
  RISCVDisassembler(const MCSubtargetInfo &STI, MCContext &Ctx,
                    MCInstrInfo const *MCII)
      : MCDisassembler(STI, Ctx), MCII(MCII) {}

  DecodeStatus getInstruction(MCInst &Instr, uint64_t &Size,
                              ArrayRef<uint8_t> Bytes, uint64_t Address,
                              raw_ostream &CStream) const override;
};
} // end anonymous namespace

static MCDisassembler *createRISCVDisassembler(const Target &T,
                                               const MCSubtargetInfo &STI,
                                               MCContext &Ctx) {
  return new RISCVDisassembler(STI, Ctx, T.createMCInstrInfo());
}

extern "C" LLVM_EXTERNAL_VISIBILITY void LLVMInitializeRISCVDisassembler() {
  // Register the disassembler for each target.
  TargetRegistry::RegisterMCDisassembler(getTheRISCV32Target(),
                                         createRISCVDisassembler);
  TargetRegistry::RegisterMCDisassembler(getTheRISCV64Target(),
                                         createRISCVDisassembler);
}

static DecodeStatus DecodeGPRRegisterClass(MCInst &Inst, uint64_t RegNo,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  const FeatureBitset &FeatureBits =
      Decoder->getSubtargetInfo().getFeatureBits();
  bool IsRV32E = FeatureBits[RISCV::FeatureRV32E];

  if (RegNo >= 32 || (IsRV32E && RegNo >= 16))
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::X0 + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeVenusRegisterClass(MCInst &Inst, uint64_t RegNo,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {

  MCRegister Reg = RISCV::VNS0 + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeFPR16RegisterClass(MCInst &Inst, uint64_t RegNo,
                                             uint64_t Address,
                                             const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::F0_H + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeFPR32RegisterClass(MCInst &Inst, uint64_t RegNo,
                                             uint64_t Address,
                                             const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::F0_F + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeFPR32CRegisterClass(MCInst &Inst, uint64_t RegNo,
                                              uint64_t Address,
                                              const MCDisassembler *Decoder) {
  if (RegNo >= 8) {
    return MCDisassembler::Fail;
  }
  MCRegister Reg = RISCV::F8_F + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeFPR64RegisterClass(MCInst &Inst, uint64_t RegNo,
                                             uint64_t Address,
                                             const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::F0_D + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeFPR64CRegisterClass(MCInst &Inst, uint64_t RegNo,
                                              uint64_t Address,
                                              const MCDisassembler *Decoder) {
  if (RegNo >= 8) {
    return MCDisassembler::Fail;
  }
  MCRegister Reg = RISCV::F8_D + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeGPRNoX0RegisterClass(MCInst &Inst, uint64_t RegNo,
                                               uint64_t Address,
                                               const MCDisassembler *Decoder) {
  if (RegNo == 0) {
    return MCDisassembler::Fail;
  }

  return DecodeGPRRegisterClass(Inst, RegNo, Address, Decoder);
}

static DecodeStatus
DecodeGPRNoX0X2RegisterClass(MCInst &Inst, uint64_t RegNo, uint64_t Address,
                             const MCDisassembler *Decoder) {
  if (RegNo == 2) {
    return MCDisassembler::Fail;
  }

  return DecodeGPRNoX0RegisterClass(Inst, RegNo, Address, Decoder);
}

static DecodeStatus DecodeGPRCRegisterClass(MCInst &Inst, uint64_t RegNo,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder) {
  if (RegNo >= 8)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::X8 + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeGPRPF64RegisterClass(MCInst &Inst, uint64_t RegNo,
                                               uint64_t Address,
                                               const MCDisassembler *Decoder) {
  if (RegNo >= 32 || RegNo & 1)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::X0 + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeVRRegisterClass(MCInst &Inst, uint64_t RegNo,
                                          uint64_t Address,
                                          const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  MCRegister Reg = RISCV::V0 + RegNo;
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeVRM2RegisterClass(MCInst &Inst, uint64_t RegNo,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  if (RegNo % 2)
    return MCDisassembler::Fail;

  const RISCVDisassembler *Dis =
      static_cast<const RISCVDisassembler *>(Decoder);
  const MCRegisterInfo *RI = Dis->getContext().getRegisterInfo();
  MCRegister Reg =
      RI->getMatchingSuperReg(RISCV::V0 + RegNo, RISCV::sub_vrm1_0,
                              &RISCVMCRegisterClasses[RISCV::VRM2RegClassID]);

  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeVRM4RegisterClass(MCInst &Inst, uint64_t RegNo,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  if (RegNo % 4)
    return MCDisassembler::Fail;

  const RISCVDisassembler *Dis =
      static_cast<const RISCVDisassembler *>(Decoder);
  const MCRegisterInfo *RI = Dis->getContext().getRegisterInfo();
  MCRegister Reg =
      RI->getMatchingSuperReg(RISCV::V0 + RegNo, RISCV::sub_vrm1_0,
                              &RISCVMCRegisterClasses[RISCV::VRM4RegClassID]);

  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus DecodeVRM8RegisterClass(MCInst &Inst, uint64_t RegNo,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder) {
  if (RegNo >= 32)
    return MCDisassembler::Fail;

  if (RegNo % 8)
    return MCDisassembler::Fail;

  const RISCVDisassembler *Dis =
      static_cast<const RISCVDisassembler *>(Decoder);
  const MCRegisterInfo *RI = Dis->getContext().getRegisterInfo();
  MCRegister Reg =
      RI->getMatchingSuperReg(RISCV::V0 + RegNo, RISCV::sub_vrm1_0,
                              &RISCVMCRegisterClasses[RISCV::VRM8RegClassID]);

  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

static DecodeStatus decodeVMaskReg(MCInst &Inst, uint64_t RegNo,
                                   uint64_t Address,
                                   const MCDisassembler *Decoder) {
  MCRegister Reg = RISCV::NoRegister;
  switch (RegNo) {
  default:
    return MCDisassembler::Fail;
  case 0:
    Reg = RISCV::V0;
    break;
  case 1:
    break;
  }
  Inst.addOperand(MCOperand::createReg(Reg));
  return MCDisassembler::Success;
}

// Add implied SP operand for instructions *SP compressed instructions. The SP
// operand isn't explicitly encoded in the instruction.
static void addImplySP(MCInst &Inst, int64_t Address,
                       const MCDisassembler *Decoder) {
  if (Inst.getOpcode() == RISCV::C_LWSP || Inst.getOpcode() == RISCV::C_SWSP ||
      Inst.getOpcode() == RISCV::C_LDSP || Inst.getOpcode() == RISCV::C_SDSP ||
      Inst.getOpcode() == RISCV::C_FLWSP ||
      Inst.getOpcode() == RISCV::C_FSWSP ||
      Inst.getOpcode() == RISCV::C_FLDSP ||
      Inst.getOpcode() == RISCV::C_FSDSP ||
      Inst.getOpcode() == RISCV::C_ADDI4SPN) {
    DecodeGPRRegisterClass(Inst, 2, Address, Decoder);
  }
  if (Inst.getOpcode() == RISCV::C_ADDI16SP) {
    DecodeGPRRegisterClass(Inst, 2, Address, Decoder);
    DecodeGPRRegisterClass(Inst, 2, Address, Decoder);
  }
}

template <unsigned N>
static DecodeStatus decodeUImmOperand(MCInst &Inst, uint64_t Imm,
                                      int64_t Address,
                                      const MCDisassembler *Decoder) {
  assert(isUInt<N>(Imm) && "Invalid immediate");
  addImplySP(Inst, Address, Decoder);
  Inst.addOperand(MCOperand::createImm(Imm));
  return MCDisassembler::Success;
}

template <unsigned N>
static DecodeStatus decodeUImmNonZeroOperand(MCInst &Inst, uint64_t Imm,
                                             int64_t Address,
                                             const MCDisassembler *Decoder) {
  if (Imm == 0)
    return MCDisassembler::Fail;
  return decodeUImmOperand<N>(Inst, Imm, Address, Decoder);
}

template <unsigned N>
static DecodeStatus decodeSImmOperand(MCInst &Inst, uint64_t Imm,
                                      int64_t Address,
                                      const MCDisassembler *Decoder) {
  assert(isUInt<N>(Imm) && "Invalid immediate");
  addImplySP(Inst, Address, Decoder);
  // Sign-extend the number in the bottom N bits of Imm
  Inst.addOperand(MCOperand::createImm(SignExtend64<N>(Imm)));
  return MCDisassembler::Success;
}

template <unsigned N>
static DecodeStatus decodeSImmNonZeroOperand(MCInst &Inst, uint64_t Imm,
                                             int64_t Address,
                                             const MCDisassembler *Decoder) {
  if (Imm == 0)
    return MCDisassembler::Fail;
  return decodeSImmOperand<N>(Inst, Imm, Address, Decoder);
}

template <unsigned N>
static DecodeStatus decodeSImmOperandAndLsl1(MCInst &Inst, uint64_t Imm,
                                             int64_t Address,
                                             const MCDisassembler *Decoder) {
  assert(isUInt<N>(Imm) && "Invalid immediate");
  // Sign-extend the number in the bottom N bits of Imm after accounting for
  // the fact that the N bit immediate is stored in N-1 bits (the LSB is
  // always zero)
  Inst.addOperand(MCOperand::createImm(SignExtend64<N>(Imm << 1)));
  return MCDisassembler::Success;
}

static DecodeStatus decodeCLUIImmOperand(MCInst &Inst, uint64_t Imm,
                                         int64_t Address,
                                         const MCDisassembler *Decoder) {
  assert(isUInt<6>(Imm) && "Invalid immediate");
  if (Imm > 31) {
    Imm = (SignExtend64<6>(Imm) & 0xfffff);
  }
  Inst.addOperand(MCOperand::createImm(Imm));
  return MCDisassembler::Success;
}

static DecodeStatus decodeFRMArg(MCInst &Inst, uint64_t Imm, int64_t Address,
                                 const MCDisassembler *Decoder) {
  assert(isUInt<3>(Imm) && "Invalid immediate");
  if (!llvm::RISCVFPRndMode::isValidRoundingMode(Imm))
    return MCDisassembler::Fail;

  Inst.addOperand(MCOperand::createImm(Imm));
  return MCDisassembler::Success;
}

static DecodeStatus decodeRVCInstrSImm(MCInst &Inst, unsigned Insn,
                                       uint64_t Address,
                                       const MCDisassembler *Decoder);

static DecodeStatus decodeRVCInstrRdSImm(MCInst &Inst, unsigned Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder);

static DecodeStatus decodeRVCInstrRdRs1UImm(MCInst &Inst, unsigned Insn,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder);

static DecodeStatus decodeRVCInstrRdRs2(MCInst &Inst, unsigned Insn,
                                        uint64_t Address,
                                        const MCDisassembler *Decoder);

static DecodeStatus decodeRVCInstrRdRs1Rs2(MCInst &Inst, unsigned Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder);
                                           
static DecodeStatus decodeVenusCmxmul(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;                                       

static DecodeStatus decodeVenusBarrier(MCInst &Inst, unsigned Insn,
                                           uint32_t Address,
                                           const MCDisassembler *Decoder);

static DecodeStatus decodeVenusMnot(MCInst &Inst, unsigned Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder);

static DecodeStatus decodeVenusSetCSR(MCInst &Inst, unsigned Insn,
                                           uint32_t Address,
                                           const MCDisassembler *Decoder);

static DecodeStatus decodeVenusSetCSRIMM(MCInst &Inst, unsigned Insn,
                                           uint32_t Address,
                                           const MCDisassembler *Decoder);

static DecodeStatus decodeVenusi32ori64(MCInst &Inst, unsigned Insn,
                                            uint32_t Address,
                                            const MCDisassembler *Decoder);

static DecodeStatus decodeVenusNormalIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;
                                          
static DecodeStatus decodeVenusNormalIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;

static DecodeStatus decodeVenusShuffleIVV(MCInst &Inst, uint64_t Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder) ;

static DecodeStatus decodeVenusShuffleIVX(MCInst &Inst, uint64_t Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder) ;

static DecodeStatus decodeVenusNormalMVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;
                                          
static DecodeStatus decodeVenusNormalMVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;

static DecodeStatus decodeVenusNormalSrc3OpsIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;
                                          
static DecodeStatus decodeVenusNormalSrc3OpsIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;

static DecodeStatus decodeVenusRange(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;                        

static DecodeStatus decodeVenusShuffleClobberMove(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) ;

#include "RISCVGenDisassemblerTables.inc"

static DecodeStatus decodeRVCInstrSImm(MCInst &Inst, unsigned Insn,
                                       uint64_t Address,
                                       const MCDisassembler *Decoder) {
  uint64_t SImm6 =
      fieldFromInstruction(Insn, 12, 1) << 5 | fieldFromInstruction(Insn, 2, 5);
  DecodeStatus Result = decodeSImmOperand<6>(Inst, SImm6, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeRVCInstrRdSImm(MCInst &Inst, unsigned Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder) {
  DecodeGPRRegisterClass(Inst, 0, Address, Decoder);
  uint64_t SImm6 =
      fieldFromInstruction(Insn, 12, 1) << 5 | fieldFromInstruction(Insn, 2, 5);
  DecodeStatus Result = decodeSImmOperand<6>(Inst, SImm6, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeRVCInstrRdRs1UImm(MCInst &Inst, unsigned Insn,
                                            uint64_t Address,
                                            const MCDisassembler *Decoder) {
  DecodeGPRRegisterClass(Inst, 0, Address, Decoder);
  Inst.addOperand(Inst.getOperand(0));
  uint64_t UImm6 =
      fieldFromInstruction(Insn, 12, 1) << 5 | fieldFromInstruction(Insn, 2, 5);
  DecodeStatus Result = decodeUImmOperand<6>(Inst, UImm6, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeRVCInstrRdRs2(MCInst &Inst, unsigned Insn,
                                        uint64_t Address,
                                        const MCDisassembler *Decoder) {
  unsigned Rd = fieldFromInstruction(Insn, 7, 5);
  unsigned Rs2 = fieldFromInstruction(Insn, 2, 5);
  DecodeGPRRegisterClass(Inst, Rd, Address, Decoder);
  DecodeGPRRegisterClass(Inst, Rs2, Address, Decoder);
  return MCDisassembler::Success;
}

static DecodeStatus decodeRVCInstrRdRs1Rs2(MCInst &Inst, unsigned Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned Rd = fieldFromInstruction(Insn, 7, 5);
  unsigned Rs2 = fieldFromInstruction(Insn, 2, 5);
  DecodeGPRRegisterClass(Inst, Rd, Address, Decoder);
  Inst.addOperand(Inst.getOperand(0));
  DecodeGPRRegisterClass(Inst, Rs2, Address, Decoder);
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusMnot(MCInst &Inst, unsigned Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned Uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  unsigned Uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  decodeUImmOperand<1>(Inst, Uimm1_vew, Address, Decoder);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, Uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusSetCSR(MCInst &Inst, unsigned Insn,
                                      uint32_t Address,
                                      const MCDisassembler *Decoder) {
  unsigned addr = fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, addr, Address, Decoder);
  unsigned uimm_value = fieldFromInstruction(Insn, 15, 12);
  DecodeStatus Result = decodeUImmOperand<12>(Inst, uimm_value, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusSetCSRIMM(MCInst &Inst, unsigned Insn,
                                      uint32_t Address,
                                      const MCDisassembler *Decoder) {
  unsigned csrimm = fieldFromInstruction(Insn, 15, 12);
  decodeUImmOperand<12>(Inst, csrimm, Address, Decoder);
  unsigned csraddr = fieldFromInstruction(Insn, 7, 5);
  DecodeStatus Result = decodeUImmOperand<5>(Inst, csraddr, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusi32ori64(MCInst &Inst, unsigned Insn,
                                      uint32_t Address,
                                      const MCDisassembler *Decoder) {
  unsigned i32i64imm = fieldFromInstruction(Insn, 7, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, i32i64imm, Address, Decoder);
  (void)Result;
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusBarrier(MCInst &Inst, unsigned Insn,
                                           uint32_t Address,
                                           const MCDisassembler *Decoder) {
  return MCDisassembler::Success;
}

#define VENUS_LITTLE

#ifndef VENUS_LITTLE
static DecodeStatus decodeVenusNormalIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 49, 1) << 10 | fieldFromInstruction(Insn, 20, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 48, 1) << 10 | fieldFromInstruction(Insn, 10, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 47, 1) << 10 | fieldFromInstruction(Insn, 0, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 50, 7) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 49, 1) << 10 | fieldFromInstruction(Insn, 20, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 48, 1) << 10 | fieldFromInstruction(Insn, 10, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 47, 1) << 10 | fieldFromInstruction(Insn, 0, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 50, 7) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalMVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src1 = fieldFromInstruction(Insn, 48, 1) << 10 | fieldFromInstruction(Insn, 10, 10) ;                                          
  unsigned src2 = fieldFromInstruction(Insn, 47, 1) << 10 | fieldFromInstruction(Insn, 0, 10) ;
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 49, 8) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalMVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src1 = fieldFromInstruction(Insn, 48, 1) << 10 | fieldFromInstruction(Insn, 10, 10) ;                                          
  unsigned src2 = fieldFromInstruction(Insn, 47, 1) << 10 | fieldFromInstruction(Insn, 0, 10) ;
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 49, 8) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusCmxmul(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src3 = fieldFromInstruction(Insn, 24, 8) ;                                          
  unsigned dest = fieldFromInstruction(Insn, 16, 8) ;
  unsigned src2 = fieldFromInstruction(Insn, 8, 8) ;
  unsigned src1 = fieldFromInstruction(Insn, 0, 8) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src3, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 47, 10) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalSrc3OpsIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src3 = fieldFromInstruction(Insn, 24, 8) ;                                          
  unsigned dest = fieldFromInstruction(Insn, 16, 8) ;
  unsigned src2 = fieldFromInstruction(Insn, 8, 8) ;
  unsigned src1 = fieldFromInstruction(Insn, 0, 8) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src3, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 47, 10) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalSrc3OpsIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src3 = fieldFromInstruction(Insn, 24, 8) ;                                          
  unsigned dest = fieldFromInstruction(Insn, 16, 8) ;
  unsigned src2 = fieldFromInstruction(Insn, 8, 8) ;
  unsigned src1 = fieldFromInstruction(Insn, 0, 8) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src3, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 47, 10) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusShuffleClobberMove(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {                                         
  unsigned dest = fieldFromInstruction(Insn, 20, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);     
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  unsigned gpr_avl = fieldFromInstruction(Insn, 47, 10) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  Inst.addOperand(Op);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusRange(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {                                         
  unsigned dest = fieldFromInstruction(Insn, 20, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);     
  unsigned uimm1_vew = fieldFromInstruction(Insn, 58, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 57, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  unsigned gpr_avl = fieldFromInstruction(Insn, 47, 10) << 5 | fieldFromInstruction(Insn, 39, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

#else
static DecodeStatus decodeVenusNormalIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 17, 1) << 10 | fieldFromInstruction(Insn, 52, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 18, 7) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 17, 1) << 10 | fieldFromInstruction(Insn, 52, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 18, 7) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusShuffleIVV(MCInst &Inst, uint64_t Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 17, 1) << 10 | fieldFromInstruction(Insn, 52, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);
  Inst.addOperand(Op);
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 18, 7) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusShuffleIVX(MCInst &Inst, uint64_t Insn,
                                         uint64_t Address,
                                         const MCDisassembler *Decoder) {
  unsigned dest = fieldFromInstruction(Insn, 17, 1) << 10 | fieldFromInstruction(Insn, 52, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);
  Inst.addOperand(Op);
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 18, 7) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalMVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src2 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src1 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 17, 8) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalMVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src1 = fieldFromInstruction(Insn, 16, 1) << 10 | fieldFromInstruction(Insn, 42, 10) ;
  unsigned src2 = fieldFromInstruction(Insn, 15, 1) << 10 | fieldFromInstruction(Insn, 32, 10) ;
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src2, Address, Decoder);      
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 17, 8) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusCmxmul(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned vd2 = fieldFromInstruction(Insn, 56, 8) ;                                          
  unsigned vd1 = fieldFromInstruction(Insn, 48, 8) ;
  unsigned vs2 = fieldFromInstruction(Insn, 40, 8) ;
  unsigned vs1 = fieldFromInstruction(Insn, 32, 8) ;
  DecodeVenusRegisterClass(Inst, vd2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, vd1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, vd2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, vd1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, vs2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, vs1, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}


static DecodeStatus decodeVenusNormalSrc3OpsIVV(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src3 = fieldFromInstruction(Insn, 56, 8) ;                                          
  unsigned dest = fieldFromInstruction(Insn, 48, 8) ;
  unsigned src1 = fieldFromInstruction(Insn, 40, 8) ;
  unsigned src2 = fieldFromInstruction(Insn, 32, 8) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src3, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusNormalSrc3OpsIVX(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {
  unsigned src3 = fieldFromInstruction(Insn, 56, 8) ;                                          
  unsigned dest = fieldFromInstruction(Insn, 48, 8) ;
  unsigned src1 = fieldFromInstruction(Insn, 40, 8) ;
  unsigned src2 = fieldFromInstruction(Insn, 32, 8) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);
  DecodeGPRRegisterClass(Inst, src1, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src2, Address, Decoder);
  DecodeVenusRegisterClass(Inst, src3, Address, Decoder);        
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusShuffleClobberMove(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {                                         
  unsigned dest = fieldFromInstruction(Insn, 52, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);    
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  Inst.addOperand(Op);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}

static DecodeStatus decodeVenusRange(MCInst &Inst, uint64_t Insn,
                                           uint64_t Address,
                                           const MCDisassembler *Decoder) {                                         
  unsigned dest = fieldFromInstruction(Insn, 52, 10) ;
  DecodeVenusRegisterClass(Inst, dest, Address, Decoder);   
  unsigned uimm1_vew = fieldFromInstruction(Insn, 26, 1);
  decodeUImmOperand<1>(Inst, uimm1_vew, Address, Decoder);
  unsigned uimm1_vmask_read = fieldFromInstruction(Insn, 25, 1);
  DecodeStatus Result = decodeUImmOperand<1>(Inst, uimm1_vmask_read, Address, Decoder);
  (void)Result;
  uint32_t Immediate = 0;
  MCOperand Op = MCOperand::createImm(Immediate);
  Inst.addOperand(Op);
  unsigned gpr_avl = fieldFromInstruction(Insn, 15, 10) << 5 | fieldFromInstruction(Insn, 7, 5);
  DecodeGPRRegisterClass(Inst, gpr_avl, Address, Decoder);
  Inst.addOperand(Op);
  assert(Result == MCDisassembler::Success && "Invalid immediate");
  return MCDisassembler::Success;
}
#endif

DecodeStatus RISCVDisassembler::getInstruction(MCInst &MI, uint64_t &Size,
                                               ArrayRef<uint8_t> Bytes,
                                               uint64_t Address,
                                               raw_ostream &CS) const {
  // TODO: This will need modification when supporting instruction set
  // extensions with instructions > 32-bits (up to 176 bits wide).
  uint32_t Insn;
  uint64_t Insn64;
  DecodeStatus Result;
  bool Opc_Custom2 = 0;
  bool Opc_Custom1 = 0;
  bool Mnot = 0;
  bool barrier = 0;
  bool setcsr = 0;
  bool setcsrimm = 0;
  bool seti32ori64 = 0;


  if (Bytes.size() >= 4 ) {
    Opc_Custom2 = (Bytes[0] & 0x7F) == 0b1011011;
    Opc_Custom1 = (Bytes[0] & 0x7F) == 0b0101011;
    Mnot = (Bytes[3] & 0xF8) == 0 && (Bytes[1] & 0x70) == 0x70;
    // barrier_num = 0x0800705b;
    barrier = Bytes[3] == 0x08 && Bytes[2] == 0x00 && Bytes[1] == 0x70 && Bytes[0] == 0x5b;
    setcsr = (Bytes[3] & 0xF8) == 0x20 && (Bytes[1] & 0x70) == 0x70;
    setcsrimm = (Bytes[3] & 0xF8) == 0x28 && (Bytes[1] & 0x70) == 0x70;
    seti32ori64 = (Bytes[3] & 0xF8) == 0xE8 && (Bytes[1] & 0x70) == 0x70;
  }
  
  if ((Opc_Custom2 && !Mnot && !barrier && !setcsr && !setcsrimm && !seti32ori64)||Opc_Custom1) {
    if (STI.getFeatureBits()[RISCV::FeatureExtVenus] ){
      Insn64 = support::endian::read64le(Bytes.data());
      LLVM_DEBUG(dbgs() << "Trying Zvenus table \n");
      Result = decodeInstruction(DecoderTableZvenus_64, MI, Insn64, Address, this, STI);
      Size = 8;
    }
  }

  // It's a 32 bit instruction if bit 0 and 1 are 1.
  else if ((Bytes[0] & 0x3) == 0x3) {
    if (Bytes.size() < 4) {
      Size = 0;
      return MCDisassembler::Fail;
    }
    Insn = support::endian::read32le(Bytes.data());
    if (STI.getFeatureBits()[RISCV::FeatureStdExtZdinx] &&
        !STI.getFeatureBits()[RISCV::Feature64Bit]) {
      LLVM_DEBUG(dbgs() << "Trying RV32Zdinx table (Double in Integer and"
                           "rv32)\n");
      Result = decodeInstruction(DecoderTableRV32Zdinx32, MI, Insn, Address,
                                 this, STI);
      if (Result != MCDisassembler::Fail) {
        Size = 4;
        return Result;
      }
    }

    if (STI.getFeatureBits()[RISCV::FeatureStdExtZfinx]) {
      LLVM_DEBUG(dbgs() << "Trying RVZfinx table (Float in Integer):\n");
      Result = decodeInstruction(DecoderTableRVZfinx32, MI, Insn, Address, this,
                                 STI);
      if (Result != MCDisassembler::Fail) {
        Size = 4;
        return Result;
      }
    }
    
    LLVM_DEBUG(dbgs() << "Trying RISCV32 table :\n");
    Result = decodeInstruction(DecoderTable32, MI, Insn, Address, this, STI);
    Size = 4;
  } else {
    if (Bytes.size() < 2) {
      Size = 0;
      return MCDisassembler::Fail;
    }
    Insn = support::endian::read16le(Bytes.data());

    if (!STI.getFeatureBits()[RISCV::Feature64Bit]) {
      LLVM_DEBUG(
          dbgs() << "Trying RISCV32Only_16 table (16-bit Instruction):\n");
      // Calling the auto-generated decoder function.
      Result = decodeInstruction(DecoderTableRISCV32Only_16, MI, Insn, Address,
                                 this, STI);
      if (Result != MCDisassembler::Fail) {
        Size = 2;
        return Result;
      }
    }

    LLVM_DEBUG(dbgs() << "Trying RISCV_C table (16-bit Instruction):\n");
    // Calling the auto-generated decoder function.
    Result = decodeInstruction(DecoderTable16, MI, Insn, Address, this, STI);
    Size = 2;
  }

  return Result;
}
