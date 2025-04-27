#include "riscv_printf.h"
#include "venus.h"

/**
 *  OFDM Demodulation —— Version 1.5
 *  Data: 2025/01/16
 *  Author: shenyihao
 *  New Feature:
 *    1. Accelerated by using vcmxmul()
 *    2. Better support for 30kHz
 *    3. More accurate cos and sin values. Users are advised to change the cos
 * and sin value stored in x.bas
 *
 *  Copyright (c) 2025 by ACE_Lab, All Rights Reserved.
 */

typedef short __v2048i16 __attribute__((ext_vector_type(2048)));
typedef char  __v2500i8 __attribute__((ext_vector_type(2500)));

// ------------------------- Set fixed point -------------------------
short fft_fixed_vec[11] = {8, 7, 7, 7, 8, 7, 7, 8, 8, 8, 8};
short fft_shift_vec[11] = {1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1};

// short fft_fixed_vec[11] = {8, 7, 8, 7, 8, 7, 7, 8, 8, 8, 7};
// short fft_shift_vec[11] = {1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0};

// short fft1024_fixed_vec[10] = {8, 7, 7, 8, 7, 8, 7, 8, 8, 7};
// short fft1024_shift_vec[10] = {1, 0, 0, 1, 0, 1, 0, 1, 1, 0};

short fft1024_fixed_vec[10] = {7, 8, 7, 8, 7, 7, 8, 8, 8, 7};
short fft1024_shift_vec[10] = {0, 1, 0, 1, 0, 0, 1, 1, 1, 0};

// fixed point complex multiplication
VENUS_INLINE __v2500i8 MUL2048_8_FIXED(__v2500i8 a, __v2500i8 b, int fix_point, int length) {
  __v2500i8 high8;
  __v2500i8 low8;
  __v2500i8 result;
  short     high_shift = 8 - fix_point;

  low8   = vmul(a, b, MASKREAD_OFF, length);
  high8  = vmulh(a, b, MASKREAD_OFF, length);
  low8   = vsrl(low8, fix_point, MASKREAD_OFF, length);
  high8  = vsll(high8, high_shift, MASKREAD_OFF, length);
  result = vor(low8, high8, MASKREAD_OFF, length);
  return result;
};

typedef struct {
  short data;
} __attribute__((aligned(64))) short_struct;

int Task_nrOFDMDemodulation(__v2500i8 vin_real, __v2500i8 vin_imag, short_struct subCarrierSpace,
                            short_struct symbolNum, __v2500i8 cos_stage0, __v2500i8 cos_stage1, __v2500i8 cos_stage2,
                            __v2500i8 cos_stage3, __v2500i8 cos_stage4, __v2500i8 cos_stage5, __v2500i8 cos_stage6,
                            __v2500i8 cos_stage7, __v2500i8 cos_stage8, __v2500i8 cos_stage9, __v2500i8 cos_stage10,
                            __v2500i8 sin_stage0, __v2500i8 sin_stage1, __v2500i8 sin_stage2, __v2500i8 sin_stage3,
                            __v2500i8 sin_stage4, __v2500i8 sin_stage5, __v2500i8 sin_stage6, __v2500i8 sin_stage7,
                            __v2500i8 sin_stage8, __v2500i8 sin_stage9, __v2500i8 sin_stage10,
                            __v2048i16 shuffle_add_stage0, __v2048i16 shuffle_add_stage1, __v2048i16 shuffle_add_stage2,
                            __v2048i16 shuffle_add_stage3, __v2048i16 shuffle_add_stage4, __v2048i16 shuffle_add_stage5,
                            __v2048i16 shuffle_add_stage6, __v2048i16 shuffle_add_stage7, __v2048i16 shuffle_add_stage8,
                            __v2048i16 shuffle_add_stage9, __v2048i16 shuffle_add_stage10, __v2048i16 shuffle_wn_stage0,
                            __v2048i16 shuffle_wn_stage1, __v2048i16 shuffle_wn_stage2, __v2048i16 shuffle_wn_stage3,
                            __v2048i16 shuffle_wn_stage4, __v2048i16 shuffle_wn_stage5, __v2048i16 shuffle_wn_stage6,
                            __v2048i16 shuffle_wn_stage7, __v2048i16 shuffle_wn_stage8, __v2048i16 shuffle_wn_stage9,
                            __v2048i16 shuffle_wn_stage10, __v2048i16 targetIndices) {

  // input parameter
  int scs        = subCarrierSpace.data;
  int symbol_num = symbolNum.data;

  // -----------------------------Start OFDM Demod-------------------------------------

  int N_FFT         = 0;
  int cp_length     = 0;
  int symbol_offset = 0;

  if (scs == 15) {
    N_FFT = 2048;
    if (symbol_num == 0 || symbol_num == 7) {
      cp_length     = 160;
      symbol_offset = 80;
    } else {
      cp_length     = 144;
      symbol_offset = 72;
    }
  } else if (scs == 30) {
    N_FFT = 1024;
    // BUG
    // if (symbol_num == 0 || symbol_num == 7) {
    if (symbol_num == 0) {
      cp_length     = 88;
      symbol_offset = 44;
    } else {
      cp_length     = 72;
      symbol_offset = 36;
    }
  }

  // --------STEP 1 : Remove CP

  __v2048i16 Remove_CP_Index;
  vrange(Remove_CP_Index, N_FFT);
  Remove_CP_Index = vsadd(Remove_CP_Index, cp_length, MASKREAD_OFF, N_FFT);
  __v2048i16 shift_CP;
  vbrdcst(shift_CP, N_FFT, MASKREAD_OFF, N_FFT);
  vbrdcst(shift_CP, 0, MASKREAD_OFF, N_FFT - symbol_offset);
  Remove_CP_Index = vrsub(Remove_CP_Index, shift_CP, MASKREAD_OFF, N_FFT);

  __v2500i8 Data_without_CP_real;
  __v2500i8 Data_without_CP_imag;
  vclaim(Data_without_CP_real);
  vclaim(Data_without_CP_imag);
  vshuffle(Data_without_CP_real, Remove_CP_Index, vin_real, SHUFFLE_GATHER, N_FFT);
  vshuffle(Data_without_CP_imag, Remove_CP_Index, vin_imag, SHUFFLE_GATHER, N_FFT);

  //   Data_without_CP_real = vsra(Data_without_CP_real, 1, MASKREAD_OFF,
  //   N_FFT); Data_without_CP_imag = vsra(Data_without_CP_imag, 1,
  //   MASKREAD_OFF, N_FFT);

  //  STEP 1 'END'---------

  //  --------STEP 2 : FFT

  __v2500i8 OFDM_OutReal;
  __v2500i8 OFDM_OutImag;
  vclaim(OFDM_OutReal);
  vclaim(OFDM_OutImag);
  short stage_init       = 0;
  short calculate_length = 0;
  if (N_FFT == 2048) {
    stage_init       = 0;
    calculate_length = 1024;
  } else if (N_FFT == 1024) {
    stage_init       = 1;
    calculate_length = 512;
  } else {
    stage_init       = 0;
    calculate_length = 1024;
  }

  __v2048i16 shuffle_for_1024_Wn;
  vclaim(shuffle_for_1024_Wn);
  vrange(shuffle_for_1024_Wn, 2048);
  shuffle_for_1024_Wn = vsll(shuffle_for_1024_Wn, 1, MASKREAD_OFF, 2048);
  __v2500i8 Wn_cos;
  vclaim(Wn_cos);
  __v2500i8 Wn_sin;
  vclaim(Wn_sin);

  __v2048i16 copy_2048; // copy_2048 = [0 1 2 ... 2047]
  vclaim(copy_2048);
  vrange(copy_2048, 2048);
  __v2048i16 move_2048to1024; // move_2048to1024 = [1024 1025 1026 ... 3071]
  move_2048to1024 = vsadd(copy_2048, 1024, MASKREAD_OFF, 2048);
  __v2048i16 move_1024to512; // move_1024to512 = [512 513 514 ... 2559]
  move_1024to512 = vsadd(copy_2048, 512, MASKREAD_OFF, 2048);

  __v2500i8 data_real_down;
  __v2500i8 data_imag_down;
  vclaim(data_real_down);
  vclaim(data_imag_down);

  __v2500i8 tempAddResult_real;
  vclaim(tempAddResult_real);
  __v2500i8 tempAddResult_imag;
  vclaim(tempAddResult_imag);
  __v2500i8 tempWnResult_real;
  vclaim(tempWnResult_real);
  __v2500i8 tempWnResult_imag;
  vclaim(tempWnResult_imag);

  __v2500i8 *cmxreal_part = &tempWnResult_real;
  __v2500i8 *cmximag_part = &tempWnResult_imag;

  __v2500i8 a_sub_b_real;
  vclaim(a_sub_b_real);
  __v2500i8 a_sub_b_imag;
  vclaim(a_sub_b_imag);
  __v2500i8 cos_tempWnResult;
  __v2500i8 sin_tempWnResult;
  vclaim(cos_tempWnResult);
  vclaim(sin_tempWnResult);

  static int fraction = 0;

  if (N_FFT == 2048) {
    // butterfly
    // a ------- (a + b)
    //      |
    // b ------- (a - b)Wn

    //  Stage 0------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[0], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[0], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[0];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage0, cos_stage0, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    // reorder
    vshuffle(Data_without_CP_real, shuffle_add_stage0, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage0, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage0, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage0, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 1------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[1], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[1], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[1];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage1, cos_stage1, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage1, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage1, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage1, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage1, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 2------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[2], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[2], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[2];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage2, cos_stage2, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage2, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage2, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage2, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage2, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 3------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[3], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[3], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[3];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage3, cos_stage3, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage3, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage3, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage3, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage3, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 4------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[4], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[4], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[4];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage4, cos_stage4, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage4, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage4, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage4, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage4, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 5------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[5], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[5], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[5];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage5, cos_stage5, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage5, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage5, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage5, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage5, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 6------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[6], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[6], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[6];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage6, cos_stage6, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage6, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage6, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage6, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage6, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 7------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[7], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[7], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[7];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage7, cos_stage7, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage7, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage7, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage7, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage7, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 8------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[8], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[8], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[8];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage8, cos_stage8, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage8, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage8, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage8, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage8, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 9------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[9], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[9], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[9];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage9, cos_stage9, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage9, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_real, shuffle_wn_stage9, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_add_stage9, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage9, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //---------------------------------------------

    //  Stage 10------------------------------------
    vshuffle(data_real_down, move_2048to1024, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_2048to1024, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft_shift_vec[10], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft_shift_vec[10], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft_fixed_vec[10];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, sin_stage10, cos_stage10, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    //  STEP 3: fft shift
    vshuffle(OFDM_OutReal, shuffle_wn_stage10, tempAddResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(OFDM_OutReal, shuffle_add_stage10, tempWnResult_real, SHUFFLE_SCATTER, calculate_length);
    vshuffle(OFDM_OutImag, shuffle_wn_stage10, tempAddResult_imag, SHUFFLE_SCATTER, calculate_length);
    vshuffle(OFDM_OutImag, shuffle_add_stage10, tempWnResult_imag, SHUFFLE_SCATTER, calculate_length);
    //  STEP 3 'END'
    //---------------------------------------------
  } else if (N_FFT == 1024) {
    //  Stage 0------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage1, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage1, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[0], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[0], MASKREAD_OFF, calculate_length);

    // (a - b)Wn 
    fraction          = fft1024_fixed_vec[0];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage0, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage0, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage0, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage0, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 1------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage2, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage2, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[1], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[1], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[1];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage1, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage1, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage1, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage1, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 2------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage3, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage3, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[2], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[2], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[2];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage2, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage2, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage2, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage2, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 3------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage4, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage4, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[3], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[3], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[3];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage3, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage3, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage3, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage3, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 4------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage5, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage5, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[4], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[4], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[4];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage4, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage4, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage4, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage4, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 5------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage6, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage6, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[5], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[5], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[5];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage5, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage5, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage5, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage5, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 6------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage7, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage7, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[6], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[6], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[6];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage6, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage6, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage6, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage6, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 7------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage8, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage8, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[7], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[7], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[7];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage7, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage7, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage7, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage7, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 8------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage9, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage9, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[8], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[8], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[8];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    vshuffle(Data_without_CP_real, shuffle_add_stage8, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_real, shuffle_wn_stage8, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_add_stage8, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(Data_without_CP_imag, shuffle_wn_stage8, tempWnResult_imag, SHUFFLE_SCATTER, 1024);
    //---------------------------------------------

    //  Stage 9------------------------------------
    vshuffle(data_real_down, move_1024to512, Data_without_CP_real, SHUFFLE_GATHER, calculate_length);
    vshuffle(data_imag_down, move_1024to512, Data_without_CP_imag, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_cos, shuffle_for_1024_Wn, cos_stage10, SHUFFLE_GATHER, calculate_length);
    vshuffle(Wn_sin, shuffle_for_1024_Wn, sin_stage10, SHUFFLE_GATHER, calculate_length);
    //  a + b
    tempAddResult_real = vsadd(Data_without_CP_real, data_real_down, MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsadd(Data_without_CP_imag, data_imag_down, MASKREAD_OFF, calculate_length);
    tempAddResult_real = vsra(tempAddResult_real, fft1024_shift_vec[9], MASKREAD_OFF, calculate_length);
    tempAddResult_imag = vsra(tempAddResult_imag, fft1024_shift_vec[9], MASKREAD_OFF, calculate_length);

    // (a - b)Wn
    fraction          = fft1024_fixed_vec[9];
    tempWnResult_real = vssub(data_real_down, Data_without_CP_real, MASKREAD_OFF, calculate_length);
    tempWnResult_imag = vssub(data_imag_down, Data_without_CP_imag, MASKREAD_OFF, calculate_length);
    vsetshamt(fraction);
    vcmxmul(cmximag_part, cmxreal_part, tempWnResult_real, tempWnResult_imag, Wn_sin, Wn_cos, MASKREAD_OFF,
            calculate_length);
    vsetshamt(0);

    //  STEP 3 : fft shift
    vshuffle(OFDM_OutReal, shuffle_wn_stage9, tempAddResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(OFDM_OutReal, shuffle_add_stage9, tempWnResult_real, SHUFFLE_SCATTER, 1024);
    vshuffle(OFDM_OutImag, shuffle_wn_stage9, tempAddResult_imag, SHUFFLE_SCATTER, 1024);
    vshuffle(OFDM_OutImag, shuffle_add_stage9, tempWnResult_imag, SHUFFLE_SCATTER, 1024);

    //  STEP 3 'END'
    //---------------------------------------------
  }

  __v2048i16 targetIndices1024;
  vclaim(targetIndices1024);
  __v2500i8 OFDM_OutReal2;
  vclaim(OFDM_OutReal2);
  __v2500i8 OFDM_OutImag2;
  vclaim(OFDM_OutImag2);
  __v2048i16 temp512;
  vclaim(temp512);
  vbrdcst(temp512, 512, MASKREAD_OFF, 240);
  if (N_FFT == 2048) {
    vshuffle(OFDM_OutReal2, targetIndices, OFDM_OutReal, SHUFFLE_GATHER, 240);
    vshuffle(OFDM_OutImag2, targetIndices, OFDM_OutImag, SHUFFLE_GATHER, 240);
  } else if (N_FFT == 1024) {
    targetIndices1024 = vssub(temp512, targetIndices, MASKREAD_OFF, 240);
    vshuffle(OFDM_OutReal2, targetIndices1024, OFDM_OutReal, SHUFFLE_GATHER, 240);
    vshuffle(OFDM_OutImag2, targetIndices1024, OFDM_OutImag, SHUFFLE_GATHER, 240);
  }

  //   printf("\n\nFFT_Length:%hd\n\n", &N_FFT);

  vreturn(OFDM_OutReal2, 256, OFDM_OutImag2, 256);
}