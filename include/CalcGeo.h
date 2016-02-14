#ifndef _CALC_GEO_H_
#define _CALC_GEO_H_

/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 幾何演算系
///
////////////////////////////////////////////////////////////////////////////

/** 実数型の指定
 * - デフォルトでは、GEO_REAL=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   GEO_REAL=doubleになる
 */
#ifdef _REAL_IS_DOUBLE_
#define GEO_REAL double
#else
#define GEO_REAL float
#endif

#define INLINE inline

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// マトリック演算系
#include "CalcGeo_Matrix.h"

// 許容誤差
#ifndef GEO_ALW_L

#ifdef _REAL_IS_DOUBLE_
#define GEO_ALW_L  (1.0e-5)
#define GEO_ALW_V  (1.0e-8)
#else
#define GEO_ALW_L  (1.0e-3)
#define GEO_ALW_V  (1.0e-5)
#endif

#endif

#ifndef PAI
#define PAI 3.14159265358979323846       /* πの値 */
#endif

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif

#include <stdbool.h>

// acos(x),asin(x) で x が
//      -1.0 ～ 1.0 の間にないと誤動作するので調整する
#define CAL_STDVAL(x)               \
    {                               \
        if ( (x) > 1.0 ) {          \
            (x) = 1.0;              \
        } else if ( (x) < -1.0 ) {  \
            (x) = -1.0;             \
        }                           \
    }


//====================================================================

///
/// ベクトルのサイズ（長さ）を求める
/// @param [in]  vec     方向ベクトル
/// @return ベクトルのサイズ（長さ）
///
INLINE GEO_REAL
CalcVecSize(
       GEO_REAL    vec[3]  // [in]  ベクトル
     )
{
    return sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}

///
/// 線分の長さ取得
/// @param [in]  pp     始点
/// @param [in]  lp     終点
/// @return ベクトルのサイズ（長さ）
///
INLINE GEO_REAL
CalcLineSize(
       GEO_REAL    pp[3],      // [in]  始点
       GEO_REAL    lp[3]       // [in]  終点
     )
{
    GEO_REAL    vec[3];
    vec[0] = lp[0] - pp[0];
    vec[1] = lp[1] - pp[1];
    vec[2] = lp[2] - pp[2];

    return sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}

///
/// ベクトルを正規化する（長さ１のベクトルとする）
/// @param [inout]  vec     方向ベクトル   out:長さ１に変更
/// @return なし
///
INLINE void
CalcNormalize(
       GEO_REAL    vec[3]  // [inout]  in :ベクトル  (長さ任意）
                         //          out:ベクトル  (長さ1）
     )
{
    GEO_REAL len = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
    vec[0] /= len;
    vec[1] /= len;
    vec[2] /= len;
}


///
/// ベクトルを正規化する（長さ１のベクトルとする）
/// @param [inout]  vec     方向ベクトル   out:長さ１に変更
/// @return true 正常 false 同一点
///
INLINE bool
CalcNormalize2(
       GEO_REAL    vec[3],    // [in]  ベクトル  (長さ任意）
       GEO_REAL    vec_o[3]   // [out] ベクトル  (長さ1）
     )
{
    GEO_REAL len = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
    if( len > GEO_ALW_L ) {
        vec_o[0] = vec[0] / len;
        vec_o[1] = vec[1] / len;
        vec_o[2] = vec[2] / len;
        return true;
    } else {
        vec_o[0] = 0.0;
        vec_o[1] = 0.0;
        vec_o[2] = 0.0;
        return false;
   }
}


///
/// 線分のベクトル作成（単位ベクトルではない）
/// @param [inout]  vec     ベクトル  （単位ベクトルではない）
/// @return なし
///
INLINE void
CalcVec(
       GEO_REAL    pp[3],      // [in]  始点
       GEO_REAL    lp[3],      // [in]  終点
       GEO_REAL    vec[3]      // [out] ベクトル
     )
{
    vec[0] = lp[0] - pp[0];
    vec[1] = lp[1] - pp[1];
    vec[2] = lp[2] - pp[2];
}


///
/// 線分のベクトルと長さ取得
/// @param [inout]  vec     方向ベクトル   out:長さ１に変更
/// @return true 正常 false 同一点
///
INLINE bool
CalcLineVec(
       GEO_REAL    pp[3],      // [in]  始点
       GEO_REAL    lp[3],      // [in]  終点
       GEO_REAL    vec[3],     // [out] ベクトル
       //GEO_REAL&   length      // [out] 長さ
       GEO_REAL*   length      // [out] 長さ
     )
{
    vec[0] = lp[0] - pp[0];
    vec[1] = lp[1] - pp[1];
    vec[2] = lp[2] - pp[2];

    *length = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );

    if( *length > GEO_ALW_L ) {
        vec[0] /= *length;
        vec[1] /= *length;
        vec[2] /= *length;
        return true;
    } else {
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
        *length = 0.0;
        return false;
    }
}


///
/// ベクトルの内積
/// @param [in]   vec1     ベクトル１
/// @param [in]   vec1     ベクトル２
/// @return 内積値
///
INLINE GEO_REAL
CalcInProduct(
       GEO_REAL    vec1[3], // [in]  ベクトル１
       GEO_REAL    vec2[3]  // [in]  ベクトル２
   )
{
    return ( vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2] );
}


///
/// ベクトルの外積
/// @param [in]   vec1     ベクトル１
/// @param [in]   vec1     ベクトル２
/// @param [out]  vec_o    外積ベクトル
/// @return なし
///
INLINE void
CalcOutProduct(
       GEO_REAL    vec1 [3],  // [in]  ベクトル１
       GEO_REAL    vec2 [3],  // [in]  ベクトル２
       GEO_REAL    vec_o[3]   // [out] 外積ベクトル
   )
{
    vec_o[ 0 ] = vec1[ 1 ] * vec2[ 2 ] - vec1[ 2 ] * vec2[ 1 ];
    vec_o[ 1 ] = vec1[ 2 ] * vec2[ 0 ] - vec1[ 0 ] * vec2[ 2 ];
    vec_o[ 2 ] = vec1[ 0 ] * vec2[ 1 ] - vec1[ 1 ] * vec2[ 0 ];
}


///
/// ベクトル間の角度（ラジアン）を求める
/// @param [in]  vec1     ベクトル１（単位ベクトル）
/// @param [in]  vec2     ベクトル２（単位ベクトル）
/// @return ベクトル間の角度（０～πラジアン）
///                       3次元のため、裏側からも見れるため角度は（０～π）で返す
///                       >πの角度を求めるためには見る方向も指定する必要あり
///
INLINE GEO_REAL
CalcVecAngle(
       GEO_REAL    vec1[3],  // [in]  :ベクトル１（単位ベクトル）
       GEO_REAL    vec2[3]   // [in]  :ベクトル２（単位ベクトル）
     )
{
    GEO_REAL asw;
    double asw_double;

    asw = CalcInProduct( vec1, vec2 );
    asw_double = asw;

    // asw が -1.0 ～ 1.0 の間にないとacosが誤動作するので調整する
    CAL_STDVAL(asw_double);

    double rad_double = acos( asw_double );  // ベクトル間の角度（０～π）

    return (GEO_REAL)(rad_double);
}


///
/// ベクトル間の角度（度数）を求める
/// @param [in]  vec1     ベクトル１（単位ベクトル）
/// @param [in]  vec2     ベクトル２（単位ベクトル）
/// @return ベクトル間の角度（０～１８０度）
///                       3次元のため、裏側からも見れるため角度は（０～１８０）で返す
///                       >１８０度の角度を求めるためには見る方向も指定する必要あり
///
INLINE GEO_REAL
CalcVecAngleDegree(
       GEO_REAL    vec1[3],  // [in]  :ベクトル１（単位ベクトル）
       GEO_REAL    vec2[3]   // [in]  :ベクトル２（単位ベクトル）
     )
{
    GEO_REAL rad = CalcVecAngle( vec1, vec2 );

    return (180.0*rad / PAI);
}


///
/// 点と方向ベクトルより平面の方程式のDを求める
///            平面の方程式  Ax + By + Cz = D
/// @param [in]  pos        点座標
/// @param [in]  vec        方向ベクトル（正規化済）
/// @return 平面の方程式のD   (原点からの距離+-あり）
///
INLINE GEO_REAL
CalcPlaneD(
       GEO_REAL    pos[3], // [in]  平面上の点座標
       GEO_REAL    vec[3]  // [in]  面の法線ベクトル（正規化済）
                         //           A:vec[0], B:vec[1], C:vec[2] にあたる
   )
{
    return CalcInProduct( pos, vec );
}


///
/// ２平面の交線（無限線分）取得
/// @param [in]  vec1        平面１の法線ベクトル(正規化済）
/// @param [in]  d1          平面１の原点からの距離
/// @param [in]  vec2        平面２の法線ベクトル(正規化済）
/// @param [in]  d2          平面２の原点からの距離
/// @param [out] pos         面の交線の通過点（原点からの最短距離）
/// @param [out] vec         面の交線のベクトル
/// @return 終了コード   true:成功  false:失敗 (２平面が並行など）
///
INLINE bool
CalcIntersectionLine(
       GEO_REAL    vec1[3],   // [in]  平面1  法線ベクトル(単位ベクトル）
                            //                平面の方程式  Ax + By + Cz = D
                            //                A:vec[0], B:vec[1], C:vec[2]
       GEO_REAL    d1,        // [in]
       GEO_REAL    vec2[3],   // [in]  平面2   法線ベクトル(単位ベクトル）
                            //                平面の方程式  Ax + By + Cz = D
                            //                A:vec[0], B:vec[1], C:vec[2]
       GEO_REAL    d2,        // [in]
       GEO_REAL    pos[3],    // [out] 面の交線の通過点（原点からの最短距離）
       GEO_REAL    vec[3]     // [out] 面の交線のベクトル
   )
{
    GEO_REAL d3;

    // 交線のベクトルを求める
    CalcOutProduct( vec1, vec2, vec );

    GEO_REAL len = CalcVecSize( vec );
#if 1
    if( len < GEO_ALW_V )  {
#else
    if( len < GEO_ALW_L )  {
#endif
       return false;   // ２面が並行
    }

    vec[0] /= len;
    vec[1] /= len;
    vec[2] /= len;

    d3 = 0.0;  // 原点を通る面を想定する

    // 3平面の交点を求める（以下の式を解けば良い）
    //     A1x + B1y + C1z = D1     vec1,d1
    //     A2x + B2y + C2z = D2     vec2,d2
    //     A3x + B3y + C3z = D3     vec,d3=0.0
    {
        GEO_REAL  b_c;
        GEO_REAL  a_c;
        GEO_REAL  a_b;
        GEO_REAL  size_det;
        GEO_REAL  d_c;
        GEO_REAL  d_b;
        GEO_REAL  a_d;
        GEO_REAL  detinv;

        /*  */
        b_c = vec2[1]*vec[2] - vec[1]*vec2[2];
        a_c = vec2[0]*vec[2] - vec[0]*vec2[2];
        a_b = vec2[0]*vec[1] - vec[0]*vec2[1];

        size_det = vec1[0]*b_c - vec1[1]*a_c + vec1[2]*a_b;
        if ( fabs( size_det ) < GEO_ALW_V ) {      /* 許容誤差以下 */
            return false;                 /* 少なくとも２平面が平行 */
        }

        d_c = d3*vec2[2] - d2*vec[2];
        d_b = d3*vec2[1] - d2*vec[1];
        a_d = vec[0]*d2  - vec2[0]*d3;
        detinv = (GEO_REAL)1.0/size_det;

        pos[0] = (vec1[1]*d_c + d1*b_c      - vec1[2]*d_b )*detinv;
        pos[1] = (-d1*a_c     - vec1[0]*d_c - vec1[2]*a_d )*detinv;
        pos[2] = (vec1[1]*a_d + vec1[0]*d_b + d1*a_b      )*detinv;
    }

    return true;
}


///
/// ２平面の交線（無限線分）取得
///           CalcIntersectionLineと求め方と通過点の位置が違うのみ
/// @param [in]  vec1        平面１の法線ベクトル(正規化済）
/// @param [in]  d1          平面１の原点からの距離
/// @param [in]  vec2        平面２の法線ベクトル(正規化済）
/// @param [in]  d2          平面２の原点からの距離
/// @param [out] pos         面の交線の通過点（x-y平面上、x-z平面上、y-z平面上のいずれか）
/// @param [out] vec         面の交線のベクトル
/// @return 終了コード   true:成功  false:失敗 (２平面が並行など）
///

INLINE bool
CalcIntersectionLine2(
       GEO_REAL    vec1[3],   // [in]  平面1  法線ベクトル(単位ベクトル）
                            //                平面の方程式  Ax + By + Cz = D
                            //                A:vec[0], B:vec[1], C:vec[2]
       GEO_REAL    d1,        // [in]         原点からの距離(+-あり）
       GEO_REAL    vec2[3],   // [in]  平面2   法線ベクトル(単位ベクトル）
                            //                平面の方程式  Ax + By + Cz = D
                            //                A:vec[0], B:vec[1], C:vec[2]
       GEO_REAL    d2,        // [in]         原点からの距離(+-あり）
       GEO_REAL    pos[3],    // [out] 面の交線の通過点
                            //          内部処理上は、X-Y平面、X-Z平面、Y-Z平面のいずれか
       GEO_REAL    vec[3]     // [out] 面の交線のベクトル(単位ベクトル）
   )
{
    GEO_REAL    vec_wk[3];

    // 交線のベクトルを求める
    CalcOutProduct( vec1, vec2, vec_wk );

    GEO_REAL len = CalcVecSize( vec_wk );
    if( len < GEO_ALW_V )  {
       return false;   // ２面が並行
    }

    // ベクトルの正規化
    vec[0] = vec_wk[0]/len;
    vec[1] = vec_wk[1]/len;
    vec[2] = vec_wk[2]/len;

    //
    // 交線を通過する点を求める
    //

    // 交線ベクトルの成分のうち、z成分が大きいならば
    //     0.5**2 + 0.5**2 + 0.5**2 = 0.75
    //     0.6**2 + 0.6**2 + 0.6**2 = 1.08  0.6を閾値にすると成り立たない場合がある
    //
    if ( fabs(vec[2]) > 0.5 )
    {
        pos[0] = ( d1 * vec2[1] - d2 * vec1[1] ) / vec_wk[2];
        pos[1] = ( d1 * vec2[0] - d2 * vec1[0] ) / ( - vec_wk[2] );
        pos[2] = 0.0;
    }
    // 交線ベクトルの成分のうち、y成分が大きいならば
    else if ( fabs(vec[1]) > 0.5 )
    {
        pos[0] = ( d1 * vec2[2] - d2 * vec1[2] ) / ( - vec_wk[1] );
        pos[1] = 0.0;
        pos[2] = ( d1 * vec2[0] - d2 * vec1[0] ) / vec_wk[1];
    }
    // 交線ベクトルの成分のうち、x成分が大きいならば
    else if ( fabs(vec[0]) > 0.5 )
    {
        pos[0] = 0.0;
        pos[1] = ( d1 * vec2[2] - d2 * vec1[2] ) / vec_wk[0];
        pos[2] = ( d1 * vec2[1] - d2 * vec1[1] ) / ( - vec_wk[0] );
    }
    else
    {
        return false;  // エラー
    }

    return true;
}


///
/// 点から線分上に垂線を下した点を求める
/// @param [in]   pnt       点座標
/// @param [in]   pos       線分の通過点
/// @param [in]   vec       線分の方向ベクトル（単位ベクトル）
/// @param [out]  pos_x     点から線分に下した垂線との交点
/// @return 終了コード   なし
///

INLINE void
CalcNearPosOnLine(
       GEO_REAL    pnt[3],   // [in]  点座標
       GEO_REAL    pos[3],   // [in]  線分の通過点
       GEO_REAL    vec[3],   // [in]  線分の方向ベクトル（単位ベクトル）
       GEO_REAL    pos_x[3]  // [out] 点から線分に下した垂線との交点
   )
{
    GEO_REAL vec_wk[3];

    vec_wk[0] = pnt[0] - pos[0];
    vec_wk[1] = pnt[1] - pos[1];
    vec_wk[2] = pnt[2] - pos[2];

    GEO_REAL dist_wk = CalcInProduct( vec, vec_wk );

    // 点から線分に下した垂線との交点
    pos_x[0] = pos[0] + vec[0]*dist_wk;
    pos_x[1] = pos[1] + vec[1]*dist_wk;
    pos_x[2] = pos[2] + vec[2]*dist_wk;

}


///
/// ベクトルの回転
/// @param [in]   rot_vec   回転軸のベクトル
/// @param [in]   rad       回転角 (rad)
/// @param [in]   veci      入力ベクトル
/// @param [out]  veco      回転後のベクトル
/// @return 終了コード   なし
///
INLINE void
CalcVecRotate2(
         GEO_REAL rot_vec[3], // [in]  回転軸のベクトル
         GEO_REAL rad,        // [in]  回転角 (rad)
         GEO_REAL veci[3],    // [in]  入力ベクトル
         GEO_REAL veco[3]     // [out] 回転後のベクトル
      )
{
    GEO_REAL org[ 3 ];
    GEO_REAL mat4[ 4 ][ 4 ];
    GEO_REAL matA[ 4 ];
    GEO_REAL matB[ 4 ];


    org[ 0 ] = 0.0;
    org[ 1 ] = 0.0;
    org[ 2 ] = 0.0;

    Calc_3dMat4Rot2( org, rot_vec, rad, mat4 );

    matA[ 0 ] = veci[ 0 ];
    matA[ 1 ] = veci[ 1 ];
    matA[ 2 ] = veci[ 2 ];
    matA[ 3 ] = 0;

    Calc_3dMat4Multi14( matA, mat4, matB );

    veco[ 0 ] = matB[ 0 ];
    veco[ 1 ] = matB[ 1 ];
    veco[ 2 ] = matB[ 2 ];
}


///
/// ベクトルのミラー
///           線対称のベクトル取得
///           ベクトルの回転の180度固定版
///           単位ベクトルである必要はありません
/// @param [in]   rot_vec   回転軸のベクトル
/// @param [in]   veci      入力ベクトル
/// @param [out]  veco      線対称にミラーしたベクトル
/// @return 終了コード   なし
///

INLINE void
CalcVecMirror(
         GEO_REAL rot_vec[3], // [in]  回転軸のベクトル
         GEO_REAL veci[3],    // [in]  入力ベクトル
         GEO_REAL veco[3]     // [out] 線対称にミラーしたベクトル
      )
{
    GEO_REAL pos_x[ 3 ];

    // ベクトルを線分とみなして、回転軸ベクトル上の
    // 交点を求める
    GEO_REAL dist_wk = CalcInProduct( rot_vec, veci );

    // 点から線分に下した垂線との交点
    pos_x[0] = rot_vec[0]*dist_wk;
    pos_x[1] = rot_vec[1]*dist_wk;
    pos_x[2] = rot_vec[2]*dist_wk;

    // ミラーベクトル設定
    veco[0] = pos_x[0] + ( pos_x[0] - veci[0] );
    veco[1] = pos_x[1] + ( pos_x[1] - veci[1] );
    veco[2] = pos_x[2] + ( pos_x[2] - veci[2] );
}


///
/// ２線分の交点（最近点）を求める
/// @param [in]   pp1       線分１始点
/// @param [in]   lp1       線分１終点
/// @param [in]   pp2       線分２始点
/// @param [in]   lp2       線分２終点
/// @param [out]  pos_x1    線分１上の最近点
/// @param [out]  pos_x2    線分２上の最近点
/// @return 終了コード   true:成功  false:失敗 (２線分が並行など）
///
INLINE bool
CalcCrossPointLine(
         GEO_REAL pp1[3],     // [in]  線分１始点
         GEO_REAL lp1[3],     // [in]  線分１終点
         GEO_REAL pp2[3],     // [in]  線分２始点
         GEO_REAL lp2[3],     // [in]  線分２終点
         GEO_REAL pos_x1[3],  // [out] 線分１上の最近点
         GEO_REAL pos_x2[3]   // [out] 線分２上の最近点
      )
{
    bool   bRet;
    GEO_REAL vec1[3], uni_vec1[3];
    GEO_REAL vec2[3], uni_vec2[3];
    GEO_REAL wk1, wk2;
    GEO_REAL vec_pp1_pp2[3];

    CalcVec( pp1, lp1, vec1 );
    bRet = CalcNormalize2( vec1, uni_vec1 );
    if( !bRet ) {
        return false;
    }

    CalcVec( pp2, lp2, vec2 );
    bRet = CalcNormalize2( vec2, uni_vec2 );
    if( !bRet ) {
        return false;
    }

    wk1 = CalcInProduct( uni_vec1, uni_vec2 );
    wk2 = 1.0 - wk1*wk1;
    if( wk2 < GEO_ALW_V ) {
        return false;  // 線分が並行
    }

    CalcVec( pp1, pp2, vec_pp1_pp2 );

    GEO_REAL d1 = (   CalcInProduct(vec_pp1_pp2,uni_vec1)
                  - wk1*CalcInProduct(vec_pp1_pp2,uni_vec2) ) / wk2;

    GEO_REAL d2 = (   wk1*CalcInProduct(vec_pp1_pp2,uni_vec1)
                  - CalcInProduct(vec_pp1_pp2,uni_vec2)     ) / wk2;

    // 線分１上の最近点
    pos_x1[0] = pp1[0] + d1*uni_vec1[0];
    pos_x1[1] = pp1[1] + d1*uni_vec1[1];
    pos_x1[2] = pp1[2] + d1*uni_vec1[2];
    // 線分２上の最近点
    pos_x2[0] = pp2[0] + d2*uni_vec2[0];
    pos_x2[1] = pp2[1] + d2*uni_vec2[1];
    pos_x2[2] = pp2[2] + d2*uni_vec2[2];

    return true;
}


///
/// 点列座標を指定された座標系に変換する（ワールド-->ローカル）
///
/// @param [in]   orig      指定座標系の原点となる点
/// @param [in]   x_axis    指定座標系のX軸となるベクトル
/// @param [in]   y_axis    指定座標系のY軸となるベクトル
/// @param [in]   z_axis    指定座標系のZ軸となるベクトル
/// @param [in]   n		    点列座標数
/// @param [in]   pnti      点列座標
/// @param [in]   pnto      指定座標系での点列座標
/// @return なし
///
INLINE void
Calc_3dTransAxisPntn(
		GEO_REAL orig[3],
		GEO_REAL x_axis[3],
		GEO_REAL y_axis[3],
		GEO_REAL z_axis[3],
        int n,
        GEO_REAL pntni[][3],
        GEO_REAL pntno[][3]
    )
{
    GEO_REAL  change_mat[4][4];       /* 変換マトリクス   */
    int     loop_cnt;               /* 座標変換カウンタ */
    GEO_REAL  src_point[4];           /* 変換前座標       */
    GEO_REAL  dsp_point[4];           /* 変換後座標       */

    /* ワールド --> ローカル変換マトリクス */
    Calc_3dMat4TranAxis( orig, x_axis, y_axis, z_axis, change_mat );

    /* 点列座標の変換 */
    for ( loop_cnt = 0; loop_cnt < n; loop_cnt++ ) {
        src_point[0] = pntni[loop_cnt][0];
        src_point[1] = pntni[loop_cnt][1];
        src_point[2] = pntni[loop_cnt][2];
        src_point[3] = 1.0;

        Calc_3dMat4Multi14( src_point, change_mat, dsp_point );

        /* 変換後の座標セーブ */
        pntno[loop_cnt][0] = dsp_point[0];
        pntno[loop_cnt][1] = dsp_point[1];
        pntno[loop_cnt][2] = dsp_point[2];
    }
}


///
/// 点列座標を指定された座標系に逆変換する（ローカル-->ワールド）
///
/// @param [in]   orig      指定座標系の原点となる点（ワールド座標系から見た座標）
/// @param [in]   x_axis    指定座標系のX軸となるベクトル（ワールド座標系から見たベクトル）
/// @param [in]   y_axis    指定座標系のY軸となるベクトル（ワールド座標系から見たベクトル）
/// @param [in]   z_axis    指定座標系のZ軸となるベクトル（ワールド座標系から見たベクトル）
/// @param [in]   n		    点列座標数
/// @param [in]   pnti      点列座標（入力）
/// @param [in]   pnto      点列座標（逆変換後）
/// @return なし
///
INLINE void
Calc_3dTransAxisPntnInv(
		GEO_REAL orig[3],
		GEO_REAL x_axis[3],
		GEO_REAL y_axis[3],
		GEO_REAL z_axis[3],
        int n,
        GEO_REAL pntni[][3],
        GEO_REAL pntno[][3]
    )
{
    GEO_REAL  change_mat[4][4];       /* 変換マトリクス   */
    int       loop_cnt;               /* 座標変換カウンタ */
    GEO_REAL  src_point[4];           /* 変換前座標       */
    GEO_REAL  dsp_point[4];           /* 変換後座標       */

    /* ローカル --> ワールド変換マトリクス */
    Calc_3dMat4TranAxisInv( orig, x_axis, y_axis, z_axis, change_mat );

    /* 点列座標の変換 */
    for ( loop_cnt = 0; loop_cnt < n; loop_cnt++ ) {
        src_point[0] = pntni[loop_cnt][0];
        src_point[1] = pntni[loop_cnt][1];
        src_point[2] = pntni[loop_cnt][2];
        src_point[3] = 1.0;

        Calc_3dMat4Multi14( src_point, change_mat, dsp_point );

        /* 変換後の座標セーブ */
        pntno[loop_cnt][0] = dsp_point[0];
        pntno[loop_cnt][1] = dsp_point[1];
        pntno[loop_cnt][2] = dsp_point[2];
    }
}

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif


#endif // _CALC_GEO_H_
