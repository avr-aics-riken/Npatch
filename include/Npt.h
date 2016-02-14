#ifndef _NPT_H_
#define _NPT_H_

/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ  関数 (C++/C)
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

// 長田パッチ許容誤差
//    並行判定等に使用する
//    長田パッチの制御点を計算するかどうかに影響する
//    GEO_ALW_* の許容誤差は倍精度：e-5-e-8、単精度：e-3-e-5レベルのため
//    曲面補間するには厳しすぎる
//          0.01;    // 0.57 degree
//          0.005;   // 0.28 degree
//          0.003;   // 0.18degree
//          0.001;   // 0.057 degree
//          0.0001;  / 0.0057 degree
//    変更したい場合は、上位で先にNPT_ALW_Vを定義する
#ifndef NPT_ALW_V
#ifdef _REAL_IS_DOUBLE_
#define NPT_ALW_V  0.001
#else
#define NPT_ALW_V  0.005
#endif
#endif

#ifndef INLINE
#define INLINE inline
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 幾何演算系
#include "CalcGeo.h"

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif


////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ パラメータ生成 関数（インライン展開なし）
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
/// @return リターンコード   =1 正常  =0 異常
/// @attention
///     エラー：許容誤差

int
npt_param_crt(
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
        NPT_REAL  cp_center [3]
   );


////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ 曲面補間 関数（インライン展開あり）
///
////////////////////////////////////////////////////////////////////////////

///
/// 長田パッチ（３次多項式） η、ξパラメータ取得
///    入力座標よりη、ξパラメータを求める
///
/// @param [in]    pos      入力点座標（p1p2p3三角形上の座標）
/// @param [in]    p1       長田パッチ 頂点１座標
/// @param [in]    p2       長田パッチ 頂点２座標
/// @param [in]    p3       長田パッチ 頂点３座標
/// @param [out]   eta      長田パッチ ηパラメータ
/// @param [out]   xi       長田パッチ ξパラメータ
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
INLINE void
npt_cvt_pos_to_eta_xi(
        NPT_REAL  pos[3],
        NPT_REAL  p1[3],
        NPT_REAL  p2[3],
        NPT_REAL  p3[3],
        NPT_REAL* eta,
        NPT_REAL* xi
     )
{
    bool   bRet;
    NPT_REAL vec12[3], vec23[3];
    NPT_REAL pos12[3], pos23[3];
    NPT_REAL pos_x12[3], pos_x23[3], pos_x_dummy[3];

    vec12[0] = p2[0]-p1[0];
    vec12[1] = p2[1]-p1[1];
    vec12[2] = p2[2]-p1[2];
    vec23[0] = p3[0]-p2[0];
    vec23[1] = p3[1]-p2[1];
    vec23[2] = p3[2]-p2[2];

    // p1p2線分とpos点を通るp2p3に平行な線分との交点を求める
    //      pos点を通るp2p3に平行な線分の仮終点
    pos23[0] = pos[0] + vec23[0];
    pos23[1] = pos[1] + vec23[1];
    pos23[2] = pos[2] + vec23[2];
    //      線分の交点
    bRet = CalcCrossPointLine(
                p1, p2,         // 線分１
                pos, pos23,     // 線分２
                pos_x12,        // 線分１上の交点
                pos_x_dummy     // 線分２上の交点
            );
    if( !bRet ) {
        printf( "##### ERROR npt_cvt_pos_to_eta_xi\n" );
        printf( "     p1=%g %g %g\n",p1[0],p1[1],p1[2] );
        printf( "     p2=%g %g %g\n",p2[0],p2[1],p2[2] );
        printf( "     pos=%g %g %g\n",pos[0],pos[1],pos[2] );
        printf( "     pos23=%g %g %g\n",pos23[0],pos23[1],pos23[2] );
        exit;
    }

    // pos点を通るp1p2に平行な線分とp2p3線分との交点を求める
    //      pos点を通るp1p2に平行な線分の仮終点
    pos12[0] = pos[0] + vec12[0];
    pos12[1] = pos[1] + vec12[1];
    pos12[2] = pos[2] + vec12[2];
    //      線分の交点
    bRet = CalcCrossPointLine(
                pos, pos12,     // 線分１
                p2, p3,         // 線分２
                pos_x_dummy,    // 線分１上の交点
                pos_x23         // 線分２上の交点
            );
    if( !bRet ) {
        printf( "##### ERROR npt_cvt_pos_to_eta_xi\n" );
        printf( "     pos=%g %g %g\n",pos[0],pos[1],pos[2] );
        printf( "     pos12=%g %g %g\n",pos12[0],pos12[1],pos12[2] );
        printf( "     p2=%g %g %g\n",p2[0],p2[1],p2[2] );
        printf( "     p3=%g %g %g\n",p3[0],p3[1],p3[2] );
        exit;
    }

    // eta = p1とpos_x12の長さ / p1p2の長さ
    *eta = CalcLineSize(p1,pos_x12) / CalcVecSize( vec12 );

    // xi  = p2とpos_x23の長さ / p2p3の長さ
    *xi  = CalcLineSize(p2,pos_x23) / CalcVecSize( vec23 );

}


///
/// 長田パッチ 近似曲面補正
///    入力：η、ξパラメータ
///
/// @param [in]    eta          入力点座標 長田パッチ ηパラメータ
/// @param [in]    xi           入力点座標 長田パッチ ξパラメータ
/// @param [in]    p1           長田パッチ 頂点１座標
/// @param [in]    p2           長田パッチ 頂点２座標
/// @param [in]    p3           長田パッチ 頂点３座標
/// @param [in]    cp_side1_1   長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [out]   pos_o        出力点座標（曲面補正後の点）
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
INLINE void
npt_correct_pnt(
        NPT_REAL  eta,
        NPT_REAL  xi,
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
    NPT_REAL u,v,w;

    // x(u,v,w) =    p1*w*w*w + cp_side1_1*3*u*w*w + cp_side1_2*3*u*u*w
    //            +  p2*u*u*u + cp_side2_1*3*u*u*v + cp_side2_2*3*u*v*v
    //            +  p3*v*v*v + cp_side3_1*3*v*v*w + cp_dide3_2*3*v*w*w
    //            +  cp_center*6*u*v*w
    //
    //     u = eta - xi
    //     v = xi
    //     w = 1 - eta
    //     u + v + w = 1
    //

    u = eta - xi;
    v = xi;
    w = 1.0 - eta;

    pos_o[0] =   p1[0]*w*w*w + cp_side1_1[0]*3.0*u*w*w + cp_side1_2[0]*3.0*u*u*w
               + p2[0]*u*u*u + cp_side2_1[0]*3.0*u*u*v + cp_side2_2[0]*3.0*u*v*v
               + p3[0]*v*v*v + cp_side3_1[0]*3.0*v*v*w + cp_side3_2[0]*3.0*v*w*w
               + cp_center[0]*6.0*u*v*w;

    pos_o[1] =   p1[1]*w*w*w + cp_side1_1[1]*3.0*u*w*w + cp_side1_2[1]*3.0*u*u*w
               + p2[1]*u*u*u + cp_side2_1[1]*3.0*u*u*v + cp_side2_2[1]*3.0*u*v*v
               + p3[1]*v*v*v + cp_side3_1[1]*3.0*v*v*w + cp_side3_2[1]*3.0*v*w*w
               + cp_center[1]*6.0*u*v*w;

    pos_o[2] =   p1[2]*w*w*w + cp_side1_1[2]*3.0*u*w*w + cp_side1_2[2]*3.0*u*u*w
               + p2[2]*u*u*u + cp_side2_1[2]*3.0*u*u*v + cp_side2_2[2]*3.0*u*v*v
               + p3[2]*v*v*v + cp_side3_1[2]*3.0*v*v*w + cp_side3_2[2]*3.0*v*w*w
               + cp_center[2]*6.0*u*v*w;
}



///
/// 長田パッチ 近似曲面補正
///
/// @param [in]    pos          入力点座標（p1p2p3三角形上の座標）
/// @param [in]    p1           長田パッチ 頂点１座標
/// @param [in]    p2           長田パッチ 頂点２座標
/// @param [in]    p3           長田パッチ 頂点３座標
/// @param [in]    cp_side1_1   長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [out]   pos_o        出力点座標（曲面補正後の点）
/// @return なし
/// @attention
///     入力点座標がp1p2p3三角形平面内でない場合、出力座標は保証されない
/// @attention
///     pos,p1,p2,p3,cp_*,pos_o は同一座標系とする。(cp_*系は相対座標ではない）
///
INLINE void
npt_correct_pnt2(
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
    NPT_REAL eta, xi;

    // 入力点座標posのη、ξパラメータを求める
    npt_cvt_pos_to_eta_xi(
                pos,
                p1, p2, p3,
                &eta, &xi
            );

    // 曲面補正点を求める
    npt_correct_pnt(
            eta, xi,
            p1, p2, p3,
            cp_side1_1, cp_side1_2,
            cp_side2_1, cp_side2_2,
            cp_side3_1, cp_side3_2,
            cp_center,
            pos_o
    );

}



///
/// 長田パッチ 頂点移動に伴う長田パッチパラメータ更新
///     長田パッチパラメータの実体は制御点である。
///     頂点座標が移動するとそれに応じてパラメータを更新する必要がある
///
/// @param [in]    p1           長田パッチ 頂点１座標
/// @param [in]    p2           長田パッチ 頂点２座標
/// @param [in]    p3           長田パッチ 頂点３座標
/// @param [in]    cp_side1_1   長田パッチ p1p2辺の３次ベジェ制御点1
/// @param [in]    cp_side1_2   長田パッチ p1p2辺の３次ベジェ制御点2
/// @param [in]    cp_side2_1   長田パッチ p2p3辺の３次ベジェ制御点1
/// @param [in]    cp_side2_2   長田パッチ p2p3辺の３次ベジェ制御点2
/// @param [in]    cp_side3_1   長田パッチ p3p1辺の３次ベジェ制御点1
/// @param [in]    cp_side3_2   長田パッチ p3p1辺の３次ベジェ制御点2
/// @param [in]    cp_center    長田パッチ 三角形中央の３次ベジェ制御点
/// @param [in]    p1_n         長田パッチ 移動後 頂点１座標
/// @param [in]    p2_n         長田パッチ 移動後 頂点２座標
/// @param [in]    p3_n         長田パッチ 移動後 頂点３座標
/// @param [out]   cp_side1_1_n 長田パッチ 移動後 p1p2辺の３次ベジェ制御点1
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
INLINE void
npt_move_vertex(
        NPT_REAL    p1[3],
        NPT_REAL    p2[3],
        NPT_REAL    p3[3],
        NPT_REAL    cp_side1_1[3],
        NPT_REAL    cp_side1_2[3],
        NPT_REAL    cp_side2_1[3],
        NPT_REAL    cp_side2_2[3],
        NPT_REAL    cp_side3_1[3],
        NPT_REAL    cp_side3_2[3],
        NPT_REAL    cp_center [3],
        NPT_REAL    p1_n[3],
        NPT_REAL    p2_n[3],
        NPT_REAL    p3_n[3],
        NPT_REAL    cp_side1_1_n[3],
        NPT_REAL    cp_side1_2_n[3],
        NPT_REAL    cp_side2_1_n[3],
        NPT_REAL    cp_side2_2_n[3],
        NPT_REAL    cp_side3_1_n[3],
        NPT_REAL    cp_side3_2_n[3],
        NPT_REAL    cp_center_n [3]
    )
{
    NPT_REAL len12, vec13[3];
    NPT_REAL x_axis[3], y_axis[3], z_axis[3];
    NPT_REAL x_axis_n[3], y_axis_n[3], z_axis_n[3];
    int      i;

    // ローカル座標軸方向の決定
    //  x_axis
    CalcLineVec( p1, p2, x_axis, &len12 );
    //  z_axis
    CalcVec( p1, p3, vec13 );
    CalcOutProduct( x_axis, vec13, z_axis );
    CalcNormalize( z_axis );
    //  y_axis
    CalcOutProduct( z_axis, x_axis, y_axis );

    // 移動後のワールド座標軸方向の決定
    //  x_axis
    CalcLineVec( p1_n, p2_n, x_axis_n, &len12 );
    //  z_axis
    CalcVec( p1_n, p3_n, vec13 );
    CalcOutProduct( x_axis_n, vec13, z_axis_n );
    CalcNormalize( z_axis_n );
    //  y_axis
    CalcOutProduct( z_axis_n, x_axis_n, y_axis_n );

    NPT_REAL  change_mat[4][4];       // 変換マトリクス
    NPT_REAL  change_mat_n[4][4];     // 変換マトリクス
    NPT_REAL  change_mat_all[4][4];   // 最終変換マトリクス

    // ワールド --> ローカル変換マトリクス
    Calc_3dMat4TranAxis( p1, x_axis, y_axis, z_axis, change_mat );

    // ローカル --> ワールド変換マトリクス
    //Calc_3dMat4TranAxisInv( p1, x_axis_n, y_axis_n, z_axis_n, change_mat_n );
    Calc_3dMat4TranAxisInv( p1_n, x_axis_n, y_axis_n, z_axis_n, change_mat_n );

    // 最終変換マトリクス
    Calc_3dMat4Multi44( change_mat, change_mat_n, change_mat_all );
    //Calc_3dMat4Multi44( change_mat_n, change_mat, change_mat_all );

    NPT_REAL pnts_in[7][4], pnts_out[7][4];

    // 配列に詰める
    pnts_in[0][0]=cp_side1_1[0];  pnts_in[0][1]=cp_side1_1[1];  pnts_in[0][2]=cp_side1_1[2];  pnts_in[0][3]=1.0;
    pnts_in[1][0]=cp_side1_2[0];  pnts_in[1][1]=cp_side1_2[1];  pnts_in[1][2]=cp_side1_2[2];  pnts_in[1][3]=1.0;
    pnts_in[2][0]=cp_side2_1[0];  pnts_in[2][1]=cp_side2_1[1];  pnts_in[2][2]=cp_side2_1[2];  pnts_in[2][3]=1.0;
    pnts_in[3][0]=cp_side2_2[0];  pnts_in[3][1]=cp_side2_2[1];  pnts_in[3][2]=cp_side2_2[2];  pnts_in[3][3]=1.0;
    pnts_in[4][0]=cp_side3_1[0];  pnts_in[4][1]=cp_side3_1[1];  pnts_in[4][2]=cp_side3_1[2];  pnts_in[4][3]=1.0;
    pnts_in[5][0]=cp_side3_2[0];  pnts_in[5][1]=cp_side3_2[1];  pnts_in[5][2]=cp_side3_2[2];  pnts_in[5][3]=1.0;
    pnts_in[6][0]=cp_center[0];   pnts_in[6][1]=cp_center[1];   pnts_in[6][2]=cp_center[2];   pnts_in[6][3]=1.0;

    // 座標変換
    for ( i=0; i<7; i++ ) {
        Calc_3dMat4Multi14( pnts_in[i], change_mat_all, pnts_out[i] );
    }

    // 配列から戻す
    cp_side1_1_n[0]=pnts_out[0][0];  cp_side1_1_n[1]=pnts_out[0][1];  cp_side1_1_n[2]=pnts_out[0][2];
    cp_side1_2_n[0]=pnts_out[1][0];  cp_side1_2_n[1]=pnts_out[1][1];  cp_side1_2_n[2]=pnts_out[1][2];
    cp_side2_1_n[0]=pnts_out[2][0];  cp_side2_1_n[1]=pnts_out[2][1];  cp_side2_1_n[2]=pnts_out[2][2];
    cp_side2_2_n[0]=pnts_out[3][0];  cp_side2_2_n[1]=pnts_out[3][1];  cp_side2_2_n[2]=pnts_out[3][2];
    cp_side3_1_n[0]=pnts_out[4][0];  cp_side3_1_n[1]=pnts_out[4][1];  cp_side3_1_n[2]=pnts_out[4][2];
    cp_side3_2_n[0]=pnts_out[5][0];  cp_side3_2_n[1]=pnts_out[5][1];  cp_side3_2_n[2]=pnts_out[5][2];
    cp_center_n [0]=pnts_out[6][0];  cp_center_n [1]=pnts_out[6][1];  cp_center_n [2]=pnts_out[6][2];
}

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif


#endif // _NPT_H_
