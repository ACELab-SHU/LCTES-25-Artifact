// Venus WiPU extension

include "llvm/IR/IntrinsicsVenusParam.td"

class VenusCommonIntrinsicsIVV: Intrinsic<[llvm_anyvector_ty],
                                       [LLVMMatchType<0>, LLVMMatchType<0>,
                                        llvm_anyint_ty, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<2>>, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<5>>]>;

class VenusCommonIntrinsicsIVVSrc3Ops: Intrinsic<[llvm_anyvector_ty],
                                       [LLVMMatchType<0>, LLVMMatchType<0>, LLVMMatchType<0>,
                                        llvm_anyint_ty, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<4>>, ImmArg<ArgIndex<6>>]>;

class VenusCommonIntrinsicsIVXSrc3Ops: Intrinsic<[llvm_anyvector_ty],
                                       [LLVMMatchType<0>, llvm_anyint_ty, LLVMMatchType<0>,
                                        LLVMMatchType<1>, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<4>>, ImmArg<ArgIndex<6>>]>;

class VenusComplexIntrinsic: Intrinsic< [llvm_anyvector_ty, LLVMMatchType<0>],
                                        [LLVMMatchType<0>, LLVMMatchType<0>, LLVMMatchType<0>, LLVMMatchType<0>,
                                         llvm_anyint_ty, LLVMMatchType<1>, LLVMMatchType<1>, LLVMMatchType<1>],
                                        [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<4>>, ImmArg<ArgIndex<5>>,ImmArg<ArgIndex<7>>]>;

class VenusCommonIntrinsicsIVX: Intrinsic<[llvm_anyvector_ty],
                                       [LLVMMatchType<0>, llvm_anyint_ty,
                                        LLVMMatchType<1>, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<2>>, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<5>>]>;

class VenusCommonIntrinsicsMVV: Intrinsic<[],
                                       [llvm_anyvector_ty, LLVMMatchType<0>,
                                        llvm_anyint_ty, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrWriteMem, IntrHasSideEffects, ImmArg<ArgIndex<2>>, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<5>>]>;

class VenusCommonIntrinsicsMVX: Intrinsic<[],
                                       [llvm_anyvector_ty, llvm_anyint_ty,
                                        LLVMMatchType<1>, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                       [IntrWriteMem, IntrHasSideEffects, ImmArg<ArgIndex<2>>, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<5>>]>;

// range
class VenusRangeIntrinsic: Intrinsic<[llvm_anyvector_ty], 
                                     [llvm_anyint_ty, LLVMMatchType<1>, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                     [IntrReadMem, IntrHasSideEffects, ImmArg<ArgIndex<0>>, ImmArg<ArgIndex<1>>, ImmArg<ArgIndex<2>>, ImmArg<ArgIndex<4>>]>;

// shuffle
class VenusShuffleIntrinsicIVV: Intrinsic<[llvm_anyvector_ty],
                                       [llvm_anyvector_ty, llvm_anyvector_ty, LLVMMatchType<0>, llvm_anyint_ty,
                                        LLVMMatchType<3>, llvm_anyint_ty, LLVMMatchType<3>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<4>>, ImmArg<ArgIndex<6>>]>;
class VenusShuffleIntrinsicIVX: Intrinsic<[llvm_anyvector_ty],
                                       [llvm_anyvector_ty, llvm_anyint_ty, LLVMMatchType<0>, llvm_anyint_ty,  
                                        LLVMMatchType<3>, llvm_anyint_ty, LLVMMatchType<3>],
                                       [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<4>>, ImmArg<ArgIndex<5>>]>;

class VenusShuffleTestIntrinsic: Intrinsic< [llvm_anyvector_ty],
                                        [llvm_anyvector_ty, llvm_anyvector_ty, LLVMMatchType<0>, 
                                         llvm_anyint_ty, LLVMMatchType<3>, llvm_anyint_ty, LLVMMatchType<3>],
                                        [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<3>>, ImmArg<ArgIndex<4>>,ImmArg<ArgIndex<6>>]>;

// venus vector addr, in essence we need an extra scalar register
// this returns the base address of venus VRF, as specified in `llc` command line
class VenusVaddrIntrinsic: Intrinsic<[llvm_i32_ty], [llvm_anyvector_ty], [IntrNoMem, IntrHasSideEffects]>;

// bind
class VenusBindIntrinsic: Intrinsic<[llvm_anyvector_ty], [llvm_ptr_ty, llvm_anyint_ty],
                                    [IntrNoMem, IntrHasSideEffects, ImmArg<ArgIndex<1>>]>;
class VenusScalarBindIntrinsic: Intrinsic<[], [llvm_ptr_ty],
                                    [IntrNoMem, IntrHasSideEffects]>;
// sink
class VenusSinkIntrinsic: Intrinsic<[], [llvm_anyvector_ty], [IntrWriteMem, IntrHasSideEffects]>;

// delimit intrinsic
// input: 0 - start of the intrinsic splitting
//        1 - end of the intrinsic splitting
class VenusDelimitIntrinsic: Intrinsic<[], [llvm_i32_ty], [IntrWriteMem, IntrHasSideEffects, ImmArg<ArgIndex<0>>]>;

// barrier
class VenusBarrierIntrinsic: Intrinsic<[], [], [IntrHasSideEffects]>;

// setcsr
class VenusSetCSRIntrinsic: Intrinsic<[], [llvm_i32_ty, llvm_i32_ty], [IntrHasSideEffects]>;

//seti32ori64
class VenusI32orI64Intrinsic: Intrinsic<[], [llvm_i32_ty], [IntrWriteMem, IntrHasSideEffects, ImmArg<ArgIndex<0>>]>;

foreach INSTNAME = VenusParam.NormInstrSet in {
    def int_venus_ # INSTNAME # _ivv  : VenusCommonIntrinsicsIVV;
    def int_venus_ # INSTNAME # _ivx  : VenusCommonIntrinsicsIVX;
}

foreach INSTNAME = VenusParam.MWInstrSet in {
    def int_venus_ # INSTNAME # _mvv  : VenusCommonIntrinsicsMVV;
    def int_venus_ # INSTNAME # _mvx  : VenusCommonIntrinsicsMVX;
}

foreach INSTNAME = VenusParam.Src3OpsInstrSet in {
    def int_venus_ # INSTNAME # _ivv  : VenusCommonIntrinsicsIVVSrc3Ops;
    def int_venus_ # INSTNAME # _ivx  : VenusCommonIntrinsicsIVXSrc3Ops;
}

def int_venus_cmxmul_ivv : VenusComplexIntrinsic;
def int_venus_range_misc : VenusRangeIntrinsic;
def int_venus_pseudo_ivv : VenusCommonIntrinsicsIVV;
def int_venus_pseudo_ivx : VenusCommonIntrinsicsIVX;
def int_venus_brdcst_ivv : VenusCommonIntrinsicsIVV;
def int_venus_brdcst_ivx : VenusCommonIntrinsicsIVX;
def int_venus_shuffle_ivv : VenusShuffleIntrinsicIVV;
def int_venus_shuffle_ivx : VenusShuffleIntrinsicIVX;
// Shuffle instruction can not be early-clobbered easily, so we allocate
// some venus registers as reserved. These swap spaces are for early-clobbered
// instructions. The venus_shuffle write result in the swap area if needed, and
// then use this intrinsic to move back to original register.
def int_venus_shuffle_clobber_move : VenusRangeIntrinsic;
// This intrinsic is called by middle-end only when there is a scatter shuffle
// i.e.:  `opd_element` > `op1_element`
def int_venus_scatter : VenusRangeIntrinsic;
def int_venus_claim_misc : VenusRangeIntrinsic;
def int_venus_vaddr_misc : VenusVaddrIntrinsic;

// mask reverse
class VenusMnotIntrinsic: Intrinsic<[], [llvm_anyint_ty, LLVMMatchType<1>, llvm_anyint_ty, LLVMMatchType<1>],
                                    [IntrWriteMem, IntrHasSideEffects, ImmArg<ArgIndex<0>>, ImmArg<ArgIndex<1>>, ImmArg<ArgIndex<3>>]>;

def int_venus_mnot_misc : VenusMnotIntrinsic;
def int_venus_bind : VenusBindIntrinsic;
def int_venus_scalar_bind : VenusScalarBindIntrinsic;
def int_venus_sink : VenusSinkIntrinsic;
def int_venus_delimit : VenusDelimitIntrinsic;
def int_venus_setcsr : VenusSetCSRIntrinsic;
def int_venus_barrier : VenusBarrierIntrinsic;
def int_venus_i32ori64 : VenusI32orI64Intrinsic;

def int_venus_shuffle_test : VenusShuffleTestIntrinsic;