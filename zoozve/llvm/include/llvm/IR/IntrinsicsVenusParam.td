def VenusParam {
    int NrRow = 256;
    int MaxReg = 2048;
    // Value Type Set per row
    list<ValueType> VenusValueTypeSet = [v4i64,  v8i64,  v16i64, v32i64,  v64i64,  v128i64,
                                         v8i32,  v16i32, v32i32, v64i32,  v128i32, v256i32,
                                         v16i16, v32i16, v64i16, v128i16, v256i16, v512i16,
                                         v32i8,  v64i8,  v128i8, v256i8,  v512i8,  v1024i8];
    // Scalar Value Type
    list<ValueType> VenusScalarValueTypeSet = [i8, i16, i32];
    // Mask Write Instruction Set
    list<string> MWInstrSet = ["seq", "sne", "sltu", "slt", "sleu", "sle", "sgtu", "sgt"];
    // Three Source Operands Instruction Set
    list<string> Src3OpsInstrSet = ["muladd", "mulsub", "addmul", "submul"];
    // Miscellanious Instruction Set
    list<string> MiscInstrSet = ["brdcst", "mnot", "range", "claim", "sink", "scalar_sink", "delimit", 
                                 "vaddr", "shuffle", "shuffle_clobber_move", "scatter", "barrier", "setcsr",
                                 "pseudo"];
    // Normal Instruction Set
    // broadcast intrinsic is IVX only, and checking is handled @ frontend
    list<string> NormInstrSet = !listconcat(["and", "or", "xor", "sll", "srl", "sra",
                                             "add", "sadd", "saddu", "rsub", "sub", "ssub", "ssubu",
                                             "mul", "mulh", "mulhu", "mulhsu",
                                             "div", "rem", "divu", "remu", 
                                             "redand", "redor", "redxor", "redmax", "redmaxu", "redmin", "redminu", "redsum"], MWInstrSet);
}

