//===-- RISCVInstrInfoVenus.td - RISC-V 'V' instructions -------*- tablegen -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
///
/// This file describes the RISC-V instructions from the Venus Vector
/// extension
///
//===----------------------------------------------------------------------===//

include "RISCVInstrFormatsVenus.td"

//===----------------------------------------------------------------------===//
// Vector register and vector group type information.
//===----------------------------------------------------------------------===//
class VTypeInfoVenus<ValueType Vec,           VenusReg VectorReg, 
                     ValueType Scal = XLenVT, RegisterClass ScalarReg = GPR,
                     ValueType MW = isVoid>
{
    ValueType Vector = Vec;
    ValueType Scalar = Scal;
    ValueType MaskWrite = MW;
    VenusReg RegClass = VectorReg;
    RegisterClass ScalarRegClass = ScalarReg;
}

defset list<VTypeInfoVenus> AllVectorsVenus = {
    foreach VECTYPE = VenusParam.VenusValueTypeSet in {
        def VNSInfo_ # VECTYPE : VTypeInfoVenus<VECTYPE, !cast<VenusReg>("VNSR")>;
    }
} // AllVectors

//===----------------------------------------------------------------------===//
// Define classes & specify DAG instances in class Pat
//===----------------------------------------------------------------------===//

class VenusPatBinary<string intrinsic_name,
                        string instruction_name,
                        ValueType result_type,
                        ValueType op1_type,
                        ValueType op2_type,
                        VenusReg op_reg_class,
                        DAGOperand op2_kind> :
    Pat<(result_type (!cast<Intrinsic>(intrinsic_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))),
        (!cast<Instruction>(instruction_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinarySrc3Ops<string intrinsic_name,
                        string instruction_name,
                        ValueType result_type,
                        ValueType op1_type,
                        ValueType op2_type,
                        ValueType op3_type,
                        VenusReg op_reg_class,
                        DAGOperand op2_kind> :
    Pat<(result_type (!cast<Intrinsic>(intrinsic_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (op3_type op_reg_class:$vs3),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))),
        (!cast<Instruction>(instruction_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (op3_type op_reg_class:$vs3),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinarySrc3Ops_Shuffle<string intrinsic_name,
                        string instruction_name,
                        ValueType result_type,
                        ValueType op1_type,
                        ValueType op2_type,
                        ValueType op3_type,
                        VenusReg op_reg1_class,
                        DAGOperand op2_kind,
                        VenusReg op_reg3_class> :
    Pat<(result_type (!cast<Intrinsic>(intrinsic_name)
                     (op1_type op_reg1_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (op3_type op_reg3_class:$vs3),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))),
        (!cast<Instruction>(instruction_name)
                     (op1_type op_reg1_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (op3_type op_reg3_class:$vs3),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinaryVoidRet<string intrinsic_name,
                        string instruction_name,
                        ValueType op1_type,
                        ValueType op2_type,
                        VenusReg op_reg_class,
                        DAGOperand op2_kind> :
    Pat<(!cast<Intrinsic>(intrinsic_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows)),
        (!cast<Instruction>(instruction_name)
                     (op1_type op_reg_class:$vs1),
                     (op2_type op2_kind:$vs2),
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinary_Bind<string intrinsic_name,
                          string instruction_name,
                          ValueType result_type> :
    Pat<(result_type (!cast<VenusBindIntrinsic>(intrinsic_name)
                     (XLenVT GPR:$ptr),
                     (XLenVT timm:$rs1))),
        (!cast<Instruction>(instruction_name)
                     (XLenVT GPR:$ptr),
                     (XLenVT timm:$rs1))>;

class VenusPatBinary_ScalarBind<string intrinsic_name,
                                string instruction_name> :
    Pat<(!cast<VenusScalarBindIntrinsic>(intrinsic_name)
                     (XLenVT GPR:$ptr)),
        (!cast<Instruction>(instruction_name)
                     (XLenVT GPR:$ptr))>;

class VenusPatBinary_Range<string intrinsic_name,
                           string instruction_name,
                           ValueType result_type> :
    Pat<(result_type (!cast<VenusRangeIntrinsic>(intrinsic_name)
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT timm:$subid),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))),
        (!cast<Instruction>(instruction_name)
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT timm:$subid),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinary_Vaddr<string intrinsic_name,
                           string instruction_name,
                           ValueType result_type,
                           ValueType op1_type,
                           VenusReg op_reg_class> :
    Pat<(result_type (!cast<VenusVaddrIntrinsic>(intrinsic_name)
                     (op1_type op_reg_class:$vs1))),
        (!cast<Instruction>(instruction_name)
                     (op1_type op_reg_class:$vs1))>;

class VenusPatBinary_Mnot<string intrinsic_name,
                           string instruction_name> :
    Pat<(!cast<VenusMnotIntrinsic>(intrinsic_name)
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows)),
        (!cast<Instruction>(instruction_name)
                     (XLenVT timm:$vew),
                     (XLenVT timm:$vmask_read),
                     (XLenVT GPR:$avl),
                     (XLenVT timm:$rows))>;

class VenusPatBinary_Sink<string intrinsic_name,
                          string instruction_name,
                          ValueType op1_type,
                          VenusReg op_reg_class> :
    Pat<(!cast<VenusSinkIntrinsic>(intrinsic_name)
                    (op1_type op_reg_class:$src)),
        (!cast<Instruction>(instruction_name)
                    (op1_type op_reg_class:$src))>;

class VenusPatBinary_Delimit<string intrinsic_name,
                             string instruction_name> :
    Pat<(!cast<VenusDelimitIntrinsic>(intrinsic_name)
                    (XLenVT GPR:$id)),
        (!cast<Instruction>(instruction_name)
                    (XLenVT GPR:$id))>;

class VenusPatBinary_SetCSR<string intrinsic_name,
                            string instruction_name> :
    Pat<(!cast<VenusSetCSRIntrinsic>(intrinsic_name)
                    (XLenVT GPR:$addr),
                    (XLenVT GPR:$value)),
        (!cast<Instruction>(instruction_name)
                    (XLenVT GPR:$addr),
                    (XLenVT GPR:$value))>;

class VenusPatBinary_SetI32orI64<string intrinsic_name,
                            string instruction_name> :
    Pat<(!cast<VenusI32orI64Intrinsic>(intrinsic_name)
                    (XLenVT timm:$i32ori64)),
        (!cast<Instruction>(instruction_name)
                    (XLenVT timm:$i32ori64))>;

class VenusPatBinary_Barrier<string intrinsic_name,
                             string instruction_name> :
    Pat<(!cast<VenusBarrierIntrinsic>(intrinsic_name)),
        (!cast<Instruction>(instruction_name))>;

//===----------------------------------------------------------------------===//
// Define intrinsics
//===----------------------------------------------------------------------===//
foreach INFO = AllVectorsVenus in { 
foreach INSTNAME = VenusParam.NormInstrSet in {
    def : VenusPatBinary<"int_venus_" # INSTNAME # "_ivv", "venus_" # INSTNAME # "_ivv",
                            INFO.Vector, INFO.Vector, INFO.Vector, INFO.RegClass, INFO.RegClass>;
    def : VenusPatBinary<"int_venus_" # INSTNAME # "_ivx", "venus_"  #INSTNAME # "_ivx",
                            INFO.Vector, INFO.Vector, INFO.Scalar, INFO.RegClass, INFO.ScalarRegClass>;                      
} // foreach NormInstrSet

foreach INSTNAME = VenusParam.MWInstrSet in {
    def : VenusPatBinaryVoidRet<"int_venus_" # INSTNAME # "_mvv", "venus_" # INSTNAME # "_mvv",
                                INFO.Vector, INFO.Vector, INFO.RegClass, INFO.RegClass>;
    def : VenusPatBinaryVoidRet<"int_venus_" # INSTNAME # "_mvx", "venus_" # INSTNAME # "_mvx",
                                INFO.Vector, INFO.Scalar, INFO.RegClass, INFO.ScalarRegClass>;
} // foreach MWInstrSet

foreach INSTNAME = VenusParam.Src3OpsInstrSet in {
    def : VenusPatBinarySrc3Ops<"int_venus_" # INSTNAME # "_ivv", "venus_" # INSTNAME # "_ivv",
                            INFO.Vector, INFO.Vector, INFO.Vector, INFO.Vector, INFO.RegClass, INFO.RegClass>;
    def : VenusPatBinarySrc3Ops<"int_venus_" # INSTNAME # "_ivx", "venus_" # INSTNAME # "_ivx",
                            INFO.Vector, INFO.Vector, INFO.Scalar, INFO.Vector, INFO.RegClass, INFO.ScalarRegClass>;
} // foreach Src3OpsInstrSet

// def : VenusPatBinaryComplex<"int_venus_cmxmul_ivv", "venus_cmxmul_ivv", INFO.Vector, INFO.RegClass>;
def : VenusPatBinary<"int_venus_pseudo_ivv", "venus_pseudo_ivv", INFO.Vector, INFO.Vector, INFO.Vector, INFO.RegClass, INFO.RegClass>;
def : VenusPatBinary<"int_venus_pseudo_ivx", "venus_pseudo_ivx", INFO.Vector, INFO.Vector, INFO.Scalar, INFO.RegClass, INFO.ScalarRegClass>;
def : VenusPatBinary<"int_venus_brdcst_ivv", "venus_brdcst_ivv", INFO.Vector, INFO.Vector, INFO.Vector, INFO.RegClass, INFO.RegClass>;
def : VenusPatBinary<"int_venus_brdcst_ivx", "venus_brdcst_ivx", INFO.Vector, INFO.Vector, INFO.Scalar, INFO.RegClass, INFO.ScalarRegClass>;
def : VenusPatBinary_Range<"int_venus_range_misc", "venus_range_misc", INFO.Vector>;
def : VenusPatBinary_Range<"int_venus_claim_misc", "venus_claim_misc", INFO.Vector>;
def : VenusPatBinary_Range<"int_venus_scatter", "venus_scatter", INFO.Vector>;
def : VenusPatBinary_Range<"int_venus_shuffle_clobber_move", "venus_shuffle_clobber_move", INFO.Vector>;
def : VenusPatBinary_Vaddr<"int_venus_vaddr_misc", "venus_vaddr_misc", INFO.Scalar, INFO.Vector, INFO.RegClass>;

def : VenusPatBinary_Bind<"int_venus_bind", "venus_bind", INFO.Vector>;
def : VenusPatBinary_ScalarBind<"int_venus_scalar_bind", "venus_scalar_bind">;
def : VenusPatBinary_Sink<"int_venus_sink", "venus_sink", INFO.Vector, INFO.RegClass>;
def : VenusPatBinary_Mnot<"int_venus_mnot_misc", "venus_mnot_misc">;
def : VenusPatBinary_Delimit<"int_venus_delimit", "venus_delimit">;
def : VenusPatBinary_SetCSR<"int_venus_setcsr", "venus_setcsr">;
def : VenusPatBinary_Barrier<"int_venus_barrier", "venus_barrier">;

} // foreach INFO

def : VenusPatBinary_SetI32orI64<"int_venus_i32ori64", "venus_i32ori64">;

foreach RETVEC = AllVectorsVenus in {
    foreach TGTVEC = AllVectorsVenus in {
        foreach IDXVEC = AllVectorsVenus in {
            def : VenusPatBinarySrc3Ops<"int_venus_shuffle_ivv", "venus_shuffle_ivv", RETVEC.Vector, TGTVEC.Vector, IDXVEC.Vector, RETVEC.Vector, TGTVEC.RegClass, IDXVEC.RegClass>;
            def : VenusPatBinarySrc3Ops<"int_venus_shuffle_ivx", "venus_shuffle_ivx", RETVEC.Vector, TGTVEC.Vector, IDXVEC.Scalar, RETVEC.Vector, TGTVEC.RegClass, IDXVEC.ScalarRegClass>;
            def : VenusPatBinarySrc3Ops<"int_venus_shuffle_test", "venus_shuffle_test", RETVEC.Vector, TGTVEC.Vector, IDXVEC.Vector, RETVEC.Vector, TGTVEC.RegClass, IDXVEC.RegClass>;
        }
    }
}

//foreach RETVEC = AllVectorsVenus in {
//    foreach TGTVEC = AllVectorsVenus in {
//        foreach IDXVEC = AllVectorsVenus in {
//            def : VenusPatBinarySrc3Ops_Shuffle<"int_venus_shuffle_test", "venus_shuffle_test", RETVEC.Vector, RETVEC.Vector, TGTVEC.Vector, IDXVEC.Vector, RETVEC.RegClass, TGTVEC.RegClass, IDXVEC.RegClass>;
//        }
//    }
//}


