//===-- RISCVInstrFormatsVenus.td - RISCV V Instruction Formats --*- tablegen -*-=//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
//  This file describes the RISC-V Venus extension instruction formats.
//
//===----------------------------------------------------------------------===//

class VenusFunct3<string name, bits<3> val> {
  string Name = name;
  bits<3> Value = val;
}
def Venus_OPIVV   : VenusFunct3<"OPIVV",  0b000>;
def Venus_OPIVX   : VenusFunct3<"OPIVX",  0b001>;
def Venus_OPIVI   : VenusFunct3<"OPIVI",  0b010>;
def Venus_OPMVV   : VenusFunct3<"OPMVV",  0b011>;
def Venus_OPMVX   : VenusFunct3<"OPMVX",  0b100>;
def Venus_OPMVI   : VenusFunct3<"OPMVI",  0b101>;
def Venus_OPMISC  : VenusFunct3<"OPMISC", 0b111>;

class VenusFunct5<string name, bits<5> val, RISCVOpcode opcode> {
  string Name = name;
  bits<5> Value = val;
  RISCVOpcode Opcode = opcode;
}
def Venus_vand    : VenusFunct5<"VAND",    0b00000, OPC_CUSTOM2>;
def Venus_vor     : VenusFunct5<"VOR",     0b00001, OPC_CUSTOM2>;
def Venus_vbrdcst : VenusFunct5<"VBRDCST", 0b00010, OPC_CUSTOM2>;
def Venus_vsll    : VenusFunct5<"VSLL",    0b00011, OPC_CUSTOM2>;
def Venus_vsrl    : VenusFunct5<"VSRL",    0b00100, OPC_CUSTOM2>;
def Venus_vsra    : VenusFunct5<"VSRA",    0b00101, OPC_CUSTOM2>;
def Venus_vxor    : VenusFunct5<"VXOR",    0b00110, OPC_CUSTOM2>;
def Venus_vseq    : VenusFunct5<"VSEQ",    0b10000, OPC_CUSTOM2>;
def Venus_vsne    : VenusFunct5<"VSNE",    0b10001, OPC_CUSTOM2>;
def Venus_vsltu   : VenusFunct5<"VSLTU",   0b10010, OPC_CUSTOM2>;
def Venus_vslt    : VenusFunct5<"VSLT",    0b10011, OPC_CUSTOM2>;
def Venus_vsleu   : VenusFunct5<"VSLEU",   0b10100, OPC_CUSTOM2>;
def Venus_vsle    : VenusFunct5<"VSLE",    0b10101, OPC_CUSTOM2>;
def Venus_vsgtu   : VenusFunct5<"VSGTU",   0b10110, OPC_CUSTOM2>;
def Venus_vsgt    : VenusFunct5<"VSGT",    0b10111, OPC_CUSTOM2>;
def Venus_vmnot   : VenusFunct5<"VMNOT",   0b00000, OPC_CUSTOM2>;
def Venus_range   : VenusFunct5<"RANGE",   0b00001, OPC_CUSTOM2>;
def Venus_shuffle : VenusFunct5<"SHUFFLE", 0b01000, OPC_CUSTOM2>;

def Venus_shuffle_clobber_move : VenusFunct5<"SHUFFLE_CLOBBER_MOVE", 0b00010, OPC_CUSTOM2>;

def Venus_vadd    : VenusFunct5<"VADD",    0b00000, OPC_CUSTOM1>;
def Venus_vrsub   : VenusFunct5<"VRSUB",   0b00001, OPC_CUSTOM1>;
def Venus_vsub    : VenusFunct5<"VSUB",    0b00010, OPC_CUSTOM1>;
def Venus_vmul    : VenusFunct5<"VMUL",    0b00011, OPC_CUSTOM1>;
def Venus_vmulh   : VenusFunct5<"VMULH",   0b00100, OPC_CUSTOM1>;
def Venus_vmulhu  : VenusFunct5<"VMULHU",  0b00101, OPC_CUSTOM1>;
def Venus_vmulhsu : VenusFunct5<"VMULHSU", 0b00110, OPC_CUSTOM1>;
def Venus_vmuladd : VenusFunct5<"VMULADD", 0b00111, OPC_CUSTOM1>;
def Venus_vmulsub : VenusFunct5<"VMULSUB", 0b01000, OPC_CUSTOM1>;
def Venus_vaddmul : VenusFunct5<"VADDMUL", 0b01001, OPC_CUSTOM1>;
def Venus_vsubmul : VenusFunct5<"VSUBMUL", 0b01010, OPC_CUSTOM1>;
def Venus_vcmxmul : VenusFunct5<"VCMXMUL", 0b01011, OPC_CUSTOM1>;
def Venus_vdiv    : VenusFunct5<"VDIV",    0b01100, OPC_CUSTOM1>;
def Venus_vrem    : VenusFunct5<"VREM",    0b01101, OPC_CUSTOM1>;
def Venus_vdivu   : VenusFunct5<"VDIVU",   0b01110, OPC_CUSTOM1>;
def Venus_vremu   : VenusFunct5<"VREMU",   0b01111, OPC_CUSTOM1>;
def Venus_vsadd   : VenusFunct5<"VSADD",   0b10000, OPC_CUSTOM1>;
def Venus_vsaddu  : VenusFunct5<"VSADDU",  0b10001, OPC_CUSTOM1>;
def Venus_vssub   : VenusFunct5<"VSSUB",   0b10010, OPC_CUSTOM1>;
def Venus_vssubu  : VenusFunct5<"VSSUBU",  0b10011, OPC_CUSTOM1>;
def Venus_vredand : VenusFunct5<"VREDAND", 0b11000, OPC_CUSTOM1>;
def Venus_vredor  : VenusFunct5<"VREDOR",  0b11001, OPC_CUSTOM1>;
def Venus_vredxor : VenusFunct5<"VREDXOR", 0b11010, OPC_CUSTOM1>;
def Venus_vredmax : VenusFunct5<"VREDMAX", 0b11011, OPC_CUSTOM1>;
def Venus_vredmaxu: VenusFunct5<"VREDMAXU",0b11100, OPC_CUSTOM1>;
def Venus_vredmin : VenusFunct5<"VREDMIN", 0b11101, OPC_CUSTOM1>;
def Venus_vredminu: VenusFunct5<"VREDMINU",0b11110, OPC_CUSTOM1>;
def Venus_vredsum : VenusFunct5<"VREDSUM", 0b11111, OPC_CUSTOM1>;

def Venus_setcsr  : VenusFunct5<"VSETCSR", 0b00100, OPC_CUSTOM2>;
def Venus_setcsrimm  : VenusFunct5<"VSETCSRIMM", 0b00101, OPC_CUSTOM2>;

def Venus_I32orI64  : VenusFunct5<"I32orI64", 0b11101, OPC_CUSTOM2>;
def Venus_shuffle_test  : VenusFunct5<"VSHUFFLE_TEST", 0b11111, OPC_CUSTOM2>;

// big-endian & small-endian switch macro
#define VENUS_LITTLE

class RVInstVenusNormal<dag outs, dag ins, string opcodestr, string argstr,
                        list<dag> pattern, InstFormat format> : Instruction {
  field bits<64> Inst;
  field bits<64> SoftFail = 0;  
  let Size = 8;

  bits<7> Opcode = 0;
#ifndef VENUS_LITTLE
  let Inst{38-32} = Opcode;
#else
  let Inst{6-0} = Opcode;
#endif

  let Namespace = "RISCV";

  let OutOperandList = outs;
  let InOperandList = ins;
  let AsmString = opcodestr # "\t" # argstr;
  let Pattern = pattern;

  // TODO: Target Specific Flags to be declared at below
  let TSFlags{4-0} = format.Value;
}

class RVInstVenusNormalI<VenusFunct5 funct5, VenusFunct3 funct3, 
                          string opcodestr, string argstr, dag ins, dag outs>
    : RVInstVenusNormal<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bit uimm1_vew;
    bit uimm1_vmask_read;
    bits<15> gpr_avl;
    bits<11> dest;
    bits<11> src1;
    bits<11> src2;
    
#ifndef VENUS_LITTLE
    let Inst{63-59} = funct5.Value;
    let Inst{58} = uimm1_vew;
    let Inst{57} = uimm1_vmask_read;
    let Inst{56-50} = gpr_avl{11-5};
    let Inst{49} = dest{10};
    let Inst{48} = src2{10};
    let Inst{47} = src1{10};
    let Inst{46-44} = funct3.Value;
    let Inst{43-39} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
    let Inst{31-30} = 0;
    let Inst{29-20} = dest{9-0};
    let Inst{19-10} = src2{9-0};
    let Inst{9-0} = src1{9-0};
#else
    let Inst{63-62} = 0;
    let Inst{61-52} = dest{9-0};
    let Inst{51-42} = src2{9-0};
    let Inst{41-32} = src1{9-0};
    let Inst{31-27} = funct5.Value;
    let Inst{26} = uimm1_vew;
    let Inst{25} = uimm1_vmask_read;
    let Inst{24-18} = gpr_avl{11-5};
    let Inst{17} = dest{10};
    let Inst{16} = src2{10};
    let Inst{15} = src1{10};
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
    
#endif
}

class RVInstVenusNormalM<VenusFunct5 funct5, VenusFunct3 funct3, 
                          string opcodestr, string argstr, dag ins, dag outs>
    : RVInstVenusNormal<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bit uimm1_vew;
    bit uimm1_vmask_read;
    bits<15> gpr_avl;
    bits<11> src1;
    bits<11> src2;

#ifndef VENUS_LITTLE    
    let Inst{63-59} = funct5.Value;
    let Inst{58} = uimm1_vew;
    let Inst{57} = uimm1_vmask_read;
    let Inst{56-49} = gpr_avl{12-5};
    let Inst{48} = src2{10};
    let Inst{47} = src1{10};
    let Inst{46-44} = funct3.Value;
    let Inst{43-39} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
    let Inst{31-30} = 0;
    let Inst{29-20} = 0;
    let Inst{19-10} = src2{9-0};
    let Inst{9-0} = src1{9-0};
#else
    let Inst{63-62} = 0;
    let Inst{61-52} = 0;
    let Inst{51-42} = src2{9-0};
    let Inst{41-32} = src1{9-0};
    let Inst{31-27} = funct5.Value;
    let Inst{26} = uimm1_vew;
    let Inst{25} = uimm1_vmask_read;
    let Inst{24-17} = gpr_avl{12-5};
    let Inst{16} = src2{10};
    let Inst{15} = src1{10};
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
#endif
}

class RVInstVenusCmxmul<VenusFunct5 funct5, VenusFunct3 funct3, 
                          string opcodestr, string argstr, dag ins, dag outs>
    : RVInstVenusNormal<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bit uimm1_vew;
    bit uimm1_vmask_read;
    bits<12> gpr_avl;
    bits<11> vd2;
    bits<11> vd1;
    bits<11> vs2;
    bits<11> vs1;

    let Inst{63-56} = vd2{7-0};
    let Inst{55-48} = vd1{7-0};
    let Inst{47-40} = vs2{7-0};
    let Inst{39-32} = vs1{7-0};
    let Inst{31-27} = funct5.Value;
    let Inst{26} = uimm1_vew;
    let Inst{25} = uimm1_vmask_read;
    let Inst{24-18} = gpr_avl{11-5};
    let Inst{17} = vd1{10};
    let Inst{16} = vs2{10};
    let Inst{15} = vs1{10};
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;

}

class RVInstVenusNormalSrc3Ops<VenusFunct5 funct5, VenusFunct3 funct3, 
                          string opcodestr, string argstr, dag ins, dag outs>
    : RVInstVenusNormal<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bit uimm1_vew;
    bit uimm1_vmask_read;
    bits<15> gpr_avl;
    bits<11> dest;
    bits<11> src1;
    bits<11> src2;
    bits<11> src3;

#ifndef VENUS_LITTLE
    let Inst{63-59} = funct5.Value;
    let Inst{58} = uimm1_vew;
    let Inst{57} = uimm1_vmask_read;
    let Inst{56-50} = gpr_avl{11-5};
    let Inst{49} = dest{10};
    let Inst{48} = src2{10};
    let Inst{47} = src1{10};
    let Inst{46-44} = funct3.Value;
    let Inst{43-39} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
    let Inst{31-24} = src3{7-0};
    let Inst{23-16} = dest{7-0};
    let Inst{15-8} = src2{7-0};
    let Inst{7-0} = src1{7-0};
#else
    let Inst{63-56} = src3{7-0};
    let Inst{55-48} = dest{7-0};
    let Inst{47-40} = src2{7-0};
    let Inst{39-32} = src1{7-0};
    let Inst{31-27} = funct5.Value;
    let Inst{26} = uimm1_vew;
    let Inst{25} = uimm1_vmask_read;
    let Inst{24-18} = gpr_avl{11-5};
    let Inst{17} = dest{10};
    let Inst{16} = src2{10};
    let Inst{15} = src1{10};
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
#endif    
}

class RVInstVenus32bit<VenusFunct5 funct5, VenusFunct3 funct3,
                       string opcodestr, string argstr, dag ins, dag outs>
    : RVInst<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bit uimm1_vew;
    bit uimm1_vmask_read;
    bits<15> gpr_avl;

    let Inst{31-27} = funct5.Value;
    let Inst{26} = uimm1_vew;
    let Inst{25} = uimm1_vmask_read;
    let Inst{24-15} = gpr_avl{14-5};
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = gpr_avl{4-0};
    let Opcode = funct5.Opcode.Value;
}

class RVInstVenusSetCSR<VenusFunct5 funct5, VenusFunct3 funct3,
                        string opcodestr, string argstr, dag ins, dag outs>
    : RVInst<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bits<12> csrreg;
    bits<5> scalarreg;

    let Inst{31-27} = funct5.Value;
    let Inst{26-15} = csrreg;
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = scalarreg;
    let Opcode = funct5.Opcode.Value;
    let DecoderMethod = "decodeVenusSetCSR";
}

class RVInstVenusSetCSRIMM<VenusFunct5 funct5, VenusFunct3 funct3,
                        string opcodestr, string argstr, dag ins, dag outs>
    : RVInst<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bits<12> csrimm;
    bits<5> csraddr;

    let Inst{31-27} = funct5.Value;
    let Inst{26-15} = csrimm;
    let Inst{14-12} = funct3.Value;
    let Inst{11-7} = csraddr;
    let Opcode = funct5.Opcode.Value;
    let DecoderMethod = "decodeVenusSetCSRIMM";
}

class RVInstVenusSetI32I64<VenusFunct5 funct5, VenusFunct3 funct3,
                        string opcodestr, string argstr, dag ins, dag outs>
    : RVInst<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    bits<1> I32orI64;
    let Inst{31-27} = funct5.Value;
    let Inst{26-15} = 0;
    let Inst{14-12} = funct3.Value;
    let Inst{11-8} = 0;
    let Inst{7} = I32orI64;
    let Opcode = funct5.Opcode.Value;
    let DecoderMethod = "decodeVenusi32ori64";
}

class RVInstVenusBarrier<string opcodestr, string argstr, dag ins, dag outs>
    : RVInst<outs, ins, opcodestr, argstr, [], InstFormatVenus>
{
    let Inst{31-0} = 0x0800705b;
    let DecoderMethod = "decodeVenusBarrier";
}


foreach INSTNAME = VenusParam.NormInstrSet in {

    def venus_ # INSTNAME # _ivv : 
                                RVInstVenusNormalI<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPIVV, 
                                    "vns_" # INSTNAME # ".ivv", 
                                    "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         !cast<VenusReg>("VNSR"):$src1,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows),
                                    (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusNormalIVV";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
    def venus_ # INSTNAME # _ivx : 
                                RVInstVenusNormalI<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPIVX, 
                                    "vns_" # INSTNAME # ".ivx", 
                                    "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         GPR:$src1,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows),
                                    (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusNormalIVX";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
} // VenusParam.NormInstrSet

foreach INSTNAME = VenusParam.MWInstrSet in {

    def venus_ # INSTNAME # _mvv :
                                RVInstVenusNormalM<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPMVV,
                                    "vns_" # INSTNAME # ".mvv", 
                                    "$src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         !cast<VenusReg>("VNSR"):$src1,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows), (outs)>{
        let DecoderMethod = "decodeVenusNormalMVV";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
    def venus_ # INSTNAME # _mvx :
                                RVInstVenusNormalM<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPMVX,
                                    "vns_" # INSTNAME # ".mvx", 
                                    "$src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         GPR:$src1,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows), (outs)>{
        let DecoderMethod = "decodeVenusNormalMVX";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
} // VenusParam.MWInstrSet

foreach INSTNAME = VenusParam.Src3OpsInstrSet in {

    def venus_ # INSTNAME # _ivv : 
                                RVInstVenusNormalSrc3Ops<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPIVV, 
                                    "vns_" # INSTNAME # ".ivv", 
                                    "$dest, $src1, $src2, $src3, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         !cast<VenusReg>("VNSR"):$src1,
                                         !cast<VenusReg>("VNSR"):$src3,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows),
                                    (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusNormalSrc3OpsIVV";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
    def venus_ # INSTNAME # _ivx : 
                                RVInstVenusNormalSrc3Ops<!cast<VenusFunct5>("Venus_v" # INSTNAME), Venus_OPIVX, 
                                    "vns_" # INSTNAME # ".ivx", 
                                    "$dest, $src1, $src2, $src3, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                         GPR:$src1,
                                         !cast<VenusReg>("VNSR"):$src3,
                                         uimm1:$uimm1_vew,
                                         uimm1:$uimm1_vmask_read,
                                         GPR:$gpr_avl,
                                         uimm16:$rows),
                                    (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusNormalSrc3OpsIVX";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }
} // VenusParam.Src3OpsInstrSet


def venus_cmxmul_ivv : RVInstVenusCmxmul<Venus_vcmxmul, Venus_OPIVV,
                                "vns_cmxmul.ivv", "$dest1, $dest2, $vd2, $vd1, $vs2, $vs1, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                (ins !cast<VenusReg>("VNSR"):$vd2,
                                     !cast<VenusReg>("VNSR"):$vd1,
                                     !cast<VenusReg>("VNSR"):$vs2,
                                     !cast<VenusReg>("VNSR"):$vs1,
                                     uimm1:$uimm1_vew,
                                     uimm1:$uimm1_vmask_read,
                                     GPR:$gpr_avl,
                                     uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest1,
                                      !cast<VenusReg>("VNSR"):$dest2)>{
        let Constraints = "$dest1 = $vd2,$dest2 = $vd1";
        let DecoderMethod = "decodeVenusCmxmul";
        let DecoderNamespace = !cast<string>("Zvenus_");
        let hasSideEffects = 1;
        let mayLoad = 0;
        let mayStore = 0;
    }


def venus_range_misc : RVInstVenusNormalI<Venus_range, Venus_OPMISC,
                                "vns_range.misc", "$dest, $uimm1_vew, $uimm1_vmask_read, $subid, $gpr_avl, $rows",
                                (ins uimm1:$uimm1_vew,
                                     uimm1:$uimm1_vmask_read,
                                     uimm16:$subid,
                                     GPR:$gpr_avl,
                                     uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusRange";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }

def venus_claim_misc : Pseudo<(outs !cast<VenusReg>("VNSR"):$dest),
                             (ins uimm1:$uimm1_vew,
                                  uimm1:$uimm1_vmask_read,
                                  uimm16:$subid,
                                  GPR:$gpr_avl,
                                  uimm16:$rows),
                             [], "vns_claim.misc", "$dest, $uimm1_vew, $uimm1_vmask_read, $subid, $gpr_avl, $rows">;

def venus_brdcst_ivv : RVInstVenusNormalI<Venus_vbrdcst, Venus_OPIVV, 
                                "vns_brdcst.ivv", 
                                "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                (ins !cast<VenusReg>("VNSR"):$src2,
                                        !cast<VenusReg>("VNSR"):$src1,
                                        uimm1:$uimm1_vew,
                                        uimm1:$uimm1_vmask_read,
                                        GPR:$gpr_avl,
                                        uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
    let Constraints = "$dest = $src2";
    let DecoderMethod = "decodeVenusNormalIVV";
    let DecoderNamespace = !cast<string>("Zvenus_");
}
def venus_brdcst_ivx : RVInstVenusNormalI<Venus_vbrdcst, Venus_OPIVX, 
                                "vns_brdcst.ivx", 
                                "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                (ins !cast<VenusReg>("VNSR"):$src2,
                                        GPR:$src1,
                                        uimm1:$uimm1_vew,
                                        uimm1:$uimm1_vmask_read,
                                        GPR:$gpr_avl,
                                        uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
    let Constraints = "$dest = $src2";
    let DecoderMethod = "decodeVenusNormalIVX";
    let DecoderNamespace = !cast<string>("Zvenus_");
}
let hasSideEffects = 1, mayLoad = 0, mayStore = 0 in {
    def venus_shuffle_test : RVInstVenusNormalI<Venus_shuffle_test, Venus_OPIVV,
                                    "vns_shuffle.test", "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                    (ins !cast<VenusReg>("VNSR"):$src2,
                                        !cast<VenusReg>("VNSR"):$src1,
                                        !cast<VenusReg>("VNSR"):$src3,
                                        uimm1:$uimm1_vew,
                                        uimm1:$uimm1_vmask_read,
                                        GPR:$gpr_avl,
                                        uimm16:$rows),
                                    (outs !cast<VenusReg>("VNSR"):$dest)>{
            // let Constraints = "@earlyclobber $dest, $dest = $src3";
            // let Constraints = "@earlyclobber $dest";
            let Constraints = "$dest = $src3";
            let DecoderMethod = "decodeVenusShuffleIVV";
            let DecoderNamespace = !cast<string>("Zvenus_");
        }
}

def venus_shuffle_ivv : RVInstVenusNormalI<Venus_shuffle, Venus_OPIVV,
                                "vns_shuffle.ivv", "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                (ins !cast<VenusReg>("VNSR"):$src2,
                                     !cast<VenusReg>("VNSR"):$src1,
                                     !cast<VenusReg>("VNSR"):$src3,
                                     uimm1:$uimm1_vew,
                                     uimm1:$uimm1_vmask_read,
                                     GPR:$gpr_avl,
                                     uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
        let Constraints = "$dest = $src3";
        let DecoderMethod = "decodeVenusShuffleIVV";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }

def venus_shuffle_ivx : RVInstVenusNormalI<Venus_shuffle, Venus_OPIVX,
                                "vns_shuffle.ivx", "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                                (ins !cast<VenusReg>("VNSR"):$src2,
                                     GPR:$src1,
                                     !cast<VenusReg>("VNSR"):$src3,
                                     uimm1:$uimm1_vew,
                                     uimm1:$uimm1_vmask_read,
                                     GPR:$gpr_avl,
                                     uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
        let Constraints = "$dest = $src3";
        let DecoderMethod = "decodeVenusShuffleIVX";
        let DecoderNamespace = !cast<string>("Zvenus_");
    }

def venus_shuffle_clobber_move : RVInstVenusNormalI<Venus_shuffle_clobber_move, Venus_OPMISC,
                                "vns_shuffle.clobbermove", "$dest, $uimm1_vew, $uimm1_vmask_read, $subid, $gpr_avl, $rows",
                                (ins uimm1:$uimm1_vew,
                                     uimm1:$uimm1_vmask_read,
                                     uimm16:$subid,
                                     GPR:$gpr_avl,
                                     uimm16:$rows),
                                (outs !cast<VenusReg>("VNSR"):$dest)>{
        let DecoderMethod = "decodeVenusShuffleClobberMove";
        let DecoderNamespace = !cast<string>("Zvenus_");
}

def venus_pseudo_ivv : Pseudo<(outs !cast<VenusReg>("VNSR"):$dest), 
                              (ins !cast<VenusReg>("VNSR"):$src2,
                                   !cast<VenusReg>("VNSR"):$src1,
                                   uimm1:$uimm1_vew,
                                   uimm1:$uimm1_vmask_read,
                                   GPR:$gpr_avl,
                                   uimm16:$rows), [], 
                                   "vns_pseudo.ivv", "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows">
                                   {
        let Constraints = "$dest = $src2";
    }

def venus_pseudo_ivx : Pseudo<(outs !cast<VenusReg>("VNSR"):$dest), 
                              (ins !cast<VenusReg>("VNSR"):$src2,
                                   GPR:$src1,
                                   uimm1:$uimm1_vew,
                                   uimm1:$uimm1_vmask_read,
                                   GPR:$gpr_avl,
                                   uimm16:$rows), [], 
                                   "vns_pseudo.ivx", "$dest, $src1, $src2, $uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows">
                                   {
        let Constraints = "$dest = $src2";
    }

def venus_scatter : Pseudo<(outs !cast<VenusReg>("VNSR"):$dest), 
                           (ins uimm1:$uimm1_vew,
                                uimm1:$uimm1_vmask_read,
                                uimm16:$subid,
                                GPR:$gpr_avl,
                                uimm16:$rows), [], "vns_scatter", "$dest, $uimm1_vew, $uimm1_vmask_read, $subid, $gpr_avl, $rows">;

def venus_vaddr_misc : Pseudo<(outs GPR:$vaddr), (ins !cast<VenusReg>("VNSR"):$src), [], 
                                     "vns_vaddr.misc", "$vaddr, $src">;


def venus_mnot_misc : RVInstVenus32bit<Venus_vmnot, Venus_OPMISC,
                        "vns_mnot.misc", "$uimm1_vew, $uimm1_vmask_read, $gpr_avl, $rows",
                        (ins uimm1:$uimm1_vew, uimm1:$uimm1_vmask_read, GPR:$gpr_avl, uimm16:$rows), 
                        (outs)>{
        let DecoderMethod = "decodeVenusMnot";
    }

def venus_delimit : Pseudo<(outs), (ins GPR:$id), [], "vns_delimit", "$id">;

let isCodeGenOnly = 0, isAsmParserOnly = 1, hasSideEffects = 1, mayLoad = 0,mayStore = 0 in {
    def venus_bind : Pseudo<(outs !cast<VenusReg>("VNSR"):$vd),
                                (ins GPR:$ptr, uimm16:$uimm16),
                                [], "vns_bind", "$vd, $ptr, $uimm16">;
}

let isCodeGenOnly = 0, isAsmParserOnly = 1, hasSideEffects = 1, mayLoad = 0, mayStore = 0 in {
    def venus_scalar_bind : Pseudo<(outs),
                                (ins GPR:$ptr),
                                [], "vns_scalar_bind", "$ptr">;
}

def venus_sink : Pseudo<(outs), (ins !cast<VenusReg>("VNSR"):$src), [],
                                     "vns_sink", "$src">;

def venus_barrier : RVInstVenusBarrier<"vns_barrier", " ",(ins), (outs)>;

def venus_setcsr : RVInstVenusSetCSR<Venus_setcsr, Venus_OPMISC,
                                     "vns_setcsr", "$scalarreg, $csrreg", (ins GPR:$csrreg, GPR:$scalarreg), (outs)>;

let hasSideEffects = 1, mayLoad = 0, mayStore = 0 in {
def venus_setcsrimm : RVInstVenusSetCSRIMM<Venus_setcsrimm, Venus_OPMISC,
                                     "vns_setcsrimm", "$csrimm, $csraddr", (ins uimm9:$csrimm,uimm5:$csraddr), (outs)>;
}

let hasSideEffects = 1 in {
def venus_i32ori64 : RVInstVenusSetI32I64<Venus_I32orI64, Venus_OPMISC,
                                     "vns_i32ori64", "$i32ori64", (ins uimm1:$i32ori64), (outs)>;
}