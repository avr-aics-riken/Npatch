
/*
 * Npatch - Nagata Patch Library
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ パラメータ生成 関数 (Fortran)
///
////////////////////////////////////////////////////////////////////////////


#include "CalcGeo.h"
#include "FNpt.h"
#include <stdlib.h>


// 長田パッチパラメータ生成
void
fnpt_param_crt_ (
        NPT_REAL  p1[3],
        NPT_REAL  norm1[3],
        NPT_REAL  p2[3],
        NPT_REAL  norm2[3],
        NPT_REAL  p3[3],
        NPT_REAL  norm3[3],
        NPT_REAL  cp_side1_1[3],
        NPT_REAL  cp_side1_2[3],
        NPT_REAL  cp_side2_1[3],
        NPT_REAL  cp_side2_2[3],
        NPT_REAL  cp_side3_1[3],
        NPT_REAL  cp_side3_2[3],
        NPT_REAL  cp_center [3],
        int*      ret
   )
{
    *ret = npt_param_crt(
                 p1,norm1,p2,norm2,p3,norm3,
                 cp_side1_1,cp_side1_2,
                 cp_side2_1,cp_side2_2,
                 cp_side3_1,cp_side3_2,
                 cp_center
              );
}


// 長田パッチ（３次多項式） η、ξパラメータ取得
//    入力座標よりη、ξパラメータを求める
void
fnpt_cvt_pos_to_eta_xi_ (
        NPT_REAL  pos[3],
        NPT_REAL  p1[3],
        NPT_REAL  p2[3],
        NPT_REAL  p3[3],
        NPT_REAL* eta,
        NPT_REAL* xi
     )
{
    npt_cvt_pos_to_eta_xi(
                           pos,p1,p2,p3,
                           eta,xi
                         );
}


// 長田パッチ 近似曲面補正
//    入力：η、ξパラメータ
void
fnpt_correct_pnt_ (
        NPT_REAL*  eta,
        NPT_REAL*  xi,
        NPT_REAL   p1[3],
        NPT_REAL   p2[3],
        NPT_REAL   p3[3],
        NPT_REAL   cp_side1_1[3],
        NPT_REAL   cp_side1_2[3],
        NPT_REAL   cp_side2_1[3],
        NPT_REAL   cp_side2_2[3],
        NPT_REAL   cp_side3_1[3],
        NPT_REAL   cp_side3_2[3],
        NPT_REAL   cp_center [3],
        NPT_REAL   pos_o[3]
    )
{
    npt_correct_pnt (
                *eta, *xi,
                p1, p2, p3,
                cp_side1_1, cp_side1_2,
                cp_side2_1, cp_side2_2,
                cp_side3_1, cp_side3_2,
                cp_center,
                pos_o
            );
}


// 長田パッチ 近似曲面補正
void
fnpt_correct_pnt2_ (
        NPT_REAL  pos[3],
        NPT_REAL  p1[3],
        NPT_REAL  p2[3],
        NPT_REAL  p3[3],
        NPT_REAL  cp_side1_1[3],
        NPT_REAL  cp_side1_2[3],
        NPT_REAL  cp_side2_1[3],
        NPT_REAL  cp_side2_2[3],
        NPT_REAL  cp_side3_1[3],
        NPT_REAL  cp_side3_2[3],
        NPT_REAL  cp_center [3],
        NPT_REAL  pos_o[3]
    )
{
    npt_correct_pnt2 (
                pos,
                p1, p2, p3,
                cp_side1_1, cp_side1_2,
                cp_side2_1, cp_side2_2,
                cp_side3_1, cp_side3_2,
                cp_center,
                pos_o
            );
}

// 長田パッチ 頂点移動に伴う長田パッチパラメータ更新
void
fnpt_move_vertex_ (
        NPT_REAL  p1[3],
        NPT_REAL  p2[3],
        NPT_REAL  p3[3],
        NPT_REAL  cp_side1_1[3],
        NPT_REAL  cp_side1_2[3],
        NPT_REAL  cp_side2_1[3],
        NPT_REAL  cp_side2_2[3],
        NPT_REAL  cp_side3_1[3],
        NPT_REAL  cp_side3_2[3],
        NPT_REAL  cp_center [3],
        NPT_REAL  p1_n[3],
        NPT_REAL  p2_n[3],
        NPT_REAL  p3_n[3],
        NPT_REAL  cp_side1_1_n[3],
        NPT_REAL  cp_side1_2_n[3],
        NPT_REAL  cp_side2_1_n[3],
        NPT_REAL  cp_side2_2_n[3],
        NPT_REAL  cp_side3_1_n[3],
        NPT_REAL  cp_side3_2_n[3],
        NPT_REAL  cp_center_n [3]
    )
{
    npt_move_vertex (
            p1, p2, p3,
            cp_side1_1, cp_side1_2, cp_side2_1, cp_side2_2,
            cp_side3_1, cp_side3_2, cp_center,
            p1_n, p2_n, p3_n,
            cp_side1_1_n, cp_side1_2_n, cp_side2_1_n, cp_side2_2_n,
            cp_side3_1_n, cp_side3_2_n, cp_center_n
        );
}
