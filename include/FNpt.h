#ifndef _FNPT_H_
#define _FNPT_H_

/*
 * Nagata Patch Function
 *
 * Copyright (c) 2015-2016 RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ  関数 (Fortranインターフェース）
///
////////////////////////////////////////////////////////////////////////////

/** 実数型の指定
 * - デフォルトでは、NPT_REAL=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   NPT_REAL=doubleになる
 */
#ifdef _REAL_IS_DOUBLE_
#define NPT_REAL double
#else
#define NPT_REAL float
#endif


#ifndef INLINE
#define INLINE inline
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "Npt.h"

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif


////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ パラメータ生成 関数
///
////////////////////////////////////////////////////////////////////////////

/// 長田パッチパラメータ生成
///
/// @param [in]    p1           長田パッチ 頂点１座標
/// @param [in]    norm1        長田パッチ 頂点１法線ベクトル（単位ベクトル）
/// @param [in]    p2           長田パッチ 頂点２座標
/// @param [in]    norm2        長田パッチ 頂点２法線ベクトル（単位ベクトル）
/// @param [in]    p3           長田パッチ 頂点３座標
/// @param [in]    norm3        長田パッチ 頂点３法線ベクトル（単位ベクトル）
/// @param [out]   cp_side1_1   長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [out]   cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [out]   cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [out]   cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [out]   cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [out]   cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [out]   cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [out]   ret          リターンコード   =1 正常  =0 異常
/// @return 戻り値なし
/// @attention
///     エラー：許容誤差

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
   );


////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ 曲面補間 関数
///
////////////////////////////////////////////////////////////////////////////

///
/// 長田パッチ（３次多項式） η、ξパラメータ取得
///    入力座標よりη、ξパラメータを求める
///
/// @param [in]    pos    	入力点座標（p1p2p3三角形上の座標）
/// @param [in]    p1     	長田パッチ 頂点１座標
/// @param [in]    p2     	長田パッチ 頂点２座標
/// @param [in]    p3     	長田パッチ 頂点３座標
/// @param [out]   eta     	長田パッチ ηパラメータ
/// @param [out]   xi     	長田パッチ ξパラメータ
/// @return なし
/// @attention
///     入力点座標がp1p2p3三角形平面内でない場合、η、ξのパラメータは保証されない
///     (例)
///         頂点1         eta  = 0.0; xi = 0.0;
///         頂点2         eta  = 1.0; xi = 0.0;
///         頂点3         eta  = 1.0; xi = 1.0;
///         辺1の中点     eta  = 0.5; xi = 0.0;
///         辺2の中点     eta  = 1.0; xi = 0.5;
///         辺3の中点     eta  = 0.5; xi = 0.5;
///         3角形の重心   eta  = 2.0/3.0, xi = 0.5*2.0/3.0;
/// @attention
///     pos,p1,p2,p3 は同一座標系とする。
///
void
fnpt_cvt_pos_to_eta_xi_ (
        NPT_REAL  pos[3],
        NPT_REAL  p1[3],
        NPT_REAL  p2[3],
        NPT_REAL  p3[3],
        NPT_REAL* eta,
        NPT_REAL* xi
     );


///
/// 長田パッチ 近似曲面補正
///    入力：η、ξパラメータ
///
/// @param [in]    eta     		入力点座標 長田パッチ ηパラメータ
/// @param [in]    xi     		入力点座標 長田パッチ ξパラメータ
/// @param [in]    p1     		長田パッチ 頂点１座標
/// @param [in]    p2     		長田パッチ 頂点２座標
/// @param [in]    p3     		長田パッチ 頂点３座標
/// @param [in]    cp_side1_1	長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [out]   pos_o    	出力点座標（曲面補正後の点）
/// @return なし
/// @attention
///     ηとξのパラメータで３角形上の座標が決まる
///     (例)
///         頂点1         eta  = 0.0; xi = 0.0;
///         頂点2         eta  = 1.0; xi = 0.0;
///         頂点3         eta  = 1.0; xi = 1.0;
///         辺1の中点     eta  = 0.5; xi = 0.0;
///         辺2の中点     eta  = 1.0; xi = 0.5;
///         辺3の中点     eta  = 0.5; xi = 0.5;
///         3角形の重心   eta  = 2.0/3.0, xi = 0.5*2.0/3.0;
/// @attention
///     p1,p2,p3,cp_*,pos_o は同一座標系とする。(cp_*系は相対座標ではない）
///
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
    );


///
/// 長田パッチ 近似曲面補正
///
/// @param [in]    pos     		入力点座標（p1p2p3三角形上の座標）
/// @param [in]    p1     		長田パッチ 頂点１座標
/// @param [in]    p2     		長田パッチ 頂点２座標
/// @param [in]    p3     		長田パッチ 頂点３座標
/// @param [in]    cp_side1_1	長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [out]   pos_o    	出力点座標（曲面補正後の点）
/// @return なし
/// @attention
///     入力点座標がp1p2p3三角形平面内でない場合、出力座標は保証されない
/// @attention
///     pos,p1,p2,p3,cp_*,pos_o は同一座標系とする。(cp_*系は相対座標ではない）
///
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
    );


///
/// 長田パッチ 頂点移動に伴う長田パッチパラメータ更新
///     長田パッチパラメータの実体は制御点である。
///     頂点座標が移動するとそれに応じてパラメータを更新する必要がある
///
/// @param [in]    p1     	長田パッチ 頂点１座標
/// @param [in]    p2     	長田パッチ 頂点２座標
/// @param [in]    p3     	長田パッチ 頂点３座標
/// @param [in]    cp_side1_1	長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [in]    p1_n     	長田パッチ 移動後 頂点１座標
/// @param [in]    p2_n     	長田パッチ 移動後 頂点２座標
/// @param [in]    p3_n     	長田パッチ 移動後 頂点３座標
/// @param [out]   cp_side1_1_n	長田パッチ 移動後 p1p2辺の３次ベジェ制御点1
/// @param [out]   cp_side1_2_n 長田パッチ 移動後 p1p2辺の３次ベジェ制御点2
/// @param [out]   cp_side2_1_n 長田パッチ 移動後 p2p3辺の３次ベジェ制御点1
/// @param [out]   cp_side2_2_n 長田パッチ 移動後 p2p3辺の３次ベジェ制御点2
/// @param [out]   cp_side3_1_n 長田パッチ 移動後 p3p1辺の３次ベジェ制御点1
/// @param [out]   cp_side3_2_n 長田パッチ 移動後 p3p1辺の３次ベジェ制御点2
/// @param [out]   cp_center_n  長田パッチ 移動後 三角形中央の３次ベジェ制御点
/// @return なし
/// @attention 長田パッチパラメータ更新は速くないので、可能であれば
///      頂点移動時に同時に長田パッチパラメータ（制御点）を更新することを推奨する
///
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
    );

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif


#endif // _FNPT_H_
