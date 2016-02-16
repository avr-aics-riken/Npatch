
/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチ パラメータ生成 関数
///
////////////////////////////////////////////////////////////////////////////


#include "CalcGeo.h"
#include "Npt.h"
#include <stdlib.h>

//------------------------------------------------------------------
//  プロトタイプ宣言： プライベート関数
//------------------------------------------------------------------
void npt_param_calcControlPointEdge( NPT_REAL p1[3], NPT_REAL norm1[3], NPT_REAL d1, NPT_REAL p2[3], NPT_REAL norm2[3], NPT_REAL d2,
           NPT_REAL cp1_e[3], NPT_REAL cp2_e[3] );
void npt_param_calcControlPointCenter( NPT_REAL p1[3], NPT_REAL p2[3], NPT_REAL p3[3],
           NPT_REAL cp1_p1p2[3], NPT_REAL cp2_p1p2[3], NPT_REAL cp1_p2p3[3], NPT_REAL cp2_p2p3[3], NPT_REAL cp1_p3p1[3], NPT_REAL cp2_p3p1[3],
           NPT_REAL cp_center[3] );
void npt_param_calcP11( NPT_REAL p1[3], NPT_REAL norm1[3], NPT_REAL d1, NPT_REAL p2[3], NPT_REAL norm2[3], NPT_REAL d2, NPT_REAL norm_base[3],
           NPT_REAL p11[3] );
void npt_param_correctP11( NPT_REAL p11[3], NPT_REAL p1[3], NPT_REAL p2[3], NPT_REAL norm_base[3],
           NPT_REAL p11_0[3], NPT_REAL p11_1[3] );

// #################################################################
//    公開関数
// #################################################################

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
/// @return リターンコード   =0 正常  !=0 異常
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
   )
{
    NPT_REAL d1, d2, d3;

    // 接平面：平面の方程式 D 計算
    //    平面の方程式  Ax + By + Cz = D
    d1 = CalcPlaneD( p1, norm1 );
    d2 = CalcPlaneD( p2, norm2 );
    d3 = CalcPlaneD( p3, norm3 );


    //-------------------
    //  p1->p2辺 制御点
    //-------------------
    npt_param_calcControlPointEdge(
           p1,          // [in]  頂点１座標
           norm1,       // [in]  頂点１ベクトル
           d1,          // [in]  頂点１ 原点からの距離(+-)
           p2,          // [in]  頂点２座標
           norm2,       // [in]  頂点２ベクトル
           d2,          // [in]  頂点２原点からの距離(+-)
           cp_side1_1,  // [out] p1p2辺の制御点1
           cp_side1_2   // [out] p1p2辺の制御点2
       );

    //-------------------
    //  p2->p3辺 制御点
    //-------------------
    npt_param_calcControlPointEdge(
           p2,          // [in]  頂点２座標
           norm2,       // [in]  頂点２ベクトル
           d2,          // [in]  頂点２原点からの距離(+-)
           p3,          // [in]  頂点３座標
           norm3,       // [in]  頂点３ベクトル
           d3,          // [in]  頂点３原点からの距離(+-)
           cp_side2_1,  // [out] p2p3辺の制御点1
           cp_side2_2   // [out] p2p3辺の制御点2
       );

    //-------------------
    //  p3->p1辺 制御点
    //-------------------
    npt_param_calcControlPointEdge(
           p3,          // [in]  頂点３座標
           norm3,       // [in]  頂点３ベクトル
           d3,          // [in]  頂点３原点からの距離(+-)
           p1,          // [in]  頂点１座標
           norm1,       // [in]  頂点１ベクトル
           d1,          // [in]  頂点１ 原点からの距離(+-)
           cp_side3_1,  // [out] p3p1辺の制御点1
           cp_side3_2   // [out] p3p1辺の制御点2
       );

    //-------------------
    //  中央制御点
    //-------------------
    npt_param_calcControlPointCenter(
           p1,          // [in]  頂点１座標
           p2,          // [in]  頂点２座標
           p3,          // [in]  頂点２座標
           cp_side1_1,  // [in]  p1p2辺の制御点1
           cp_side1_2,  // [in]  p1p2辺の制御点2
           cp_side2_1,  // [in]  p2p3辺の制御点1
           cp_side2_2,  // [in]  p2p3辺の制御点2
           cp_side3_1,  // [in]  p3p1辺の制御点1
           cp_side3_2,  // [in]  p3p1辺の制御点2
           cp_center    // [out] 中央制御点
        );

    return 0;
}


// #################################################################
//    非公開（プライベート）関数
// #################################################################

/// 長田パッチ各辺の制御点取得
void
npt_param_calcControlPointEdge(
           NPT_REAL        p1[3],        // [in]  頂点１座標
           NPT_REAL        norm1[3],     // [in]  頂点１ベクトル
           NPT_REAL        d1,           // [in]  頂点１ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
           NPT_REAL        p2[3],        // [in]  頂点２座標
           NPT_REAL        norm2[3],     // [in]  頂点２ベクトル
           NPT_REAL        d2,           // [in]  頂点２ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
           NPT_REAL        cp1_e[3],     // [out] 辺から求める制御点１
           NPT_REAL        cp2_e[3]      // [out] 辺から求める制御点２
       )
{
    NPT_REAL vec_wk[3];
    NPT_REAL norm_wk1[3], norm_wk2[3];
    NPT_REAL norm_base[3]; // 曲線平面II（基準面）の法線ベクトル
    NPT_REAL p11[3];       // 制御点p11
    NPT_REAL p11_0[3];     // 制御点座標（補正後 1番目の制御点用）
    NPT_REAL p11_1[3];     // 制御点座標（補正後 2番目の制御点用）

    //-------------------------------------------
    // 制御点を置く、基準面を求める
    //     2等分した平面 曲線平面II（基準面）とする
    //     基準面： p1, norm_base で定義される
    //-------------------------------------------

    //  p1->p2辺
    vec_wk[0] = p2[0] - p1[0];
    vec_wk[1] = p2[1] - p1[1];
    vec_wk[2] = p2[2] - p1[2];

    //     法線ベクトル
    CalcOutProduct( vec_wk, norm1, norm_wk1 ); // at p1
    CalcOutProduct( vec_wk, norm2, norm_wk2 ); // at p2
    //     2等分して基準となる法線ベクトルとする
    norm_base[0] = norm_wk1[0] + norm_wk2[0];
    norm_base[1] = norm_wk1[1] + norm_wk2[1];
    norm_base[2] = norm_wk1[2] + norm_wk2[2];
    CalcNormalize( norm_base );  // 単位ベクトル化

    //-------------------------------------------
    // 制御点決定
    //-------------------------------------------

    // 制御点p11取得
    //     2次多項式と同様
    npt_param_calcP11(
              p1,        // [in]  頂点１座標
              norm1,     // [in]  頂点１ベクトル
              d1,        // [in]  頂点１ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
              p2,        // [in]  頂点２座標
              norm2,     // [in]  頂点２ベクトル
              d2,        // [in]  頂点２ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
              norm_base, // [in]  曲線平面II（基準面）の法線ベクトル
              p11        // [out] 制御点座標
         );

    // 制御点p11逆行補正対応
    npt_param_correctP11(
              p11,       // [in]  制御点座標
              p1,        // [in]  頂点１座標
              p2,        // [in]  頂点２座標
              norm_base, // [in]  曲線平面II（基準面）の法線ベクトル
              p11_0,     // [out] 制御点座標（補正後 1番目の制御点用）
              p11_1      // [out] 制御点座標（補正後 2番目の制御点用）
         );

    // 制御点(3次多項式用) 設定
    cp1_e[0] = ( p1[0] + 2.0*p11_0[0] ) / 3.0;
    cp1_e[1] = ( p1[1] + 2.0*p11_0[1] ) / 3.0;
    cp1_e[2] = ( p1[2] + 2.0*p11_0[2] ) / 3.0;

    cp2_e[0] = ( p2[0] + 2.0*p11_1[0] ) / 3.0;
    cp2_e[1] = ( p2[1] + 2.0*p11_1[1] ) / 3.0;
    cp2_e[2] = ( p2[2] + 2.0*p11_1[2] ) / 3.0;

#if 0
    // デバッグ用 cp1_eとcp2_eの順番が正当か確認する
    NPT_REAL vec_wk1[3], vec_wk2[3];
    CalcVec( p1,    p2, vec_wk1 );
    CalcVec( cp1_e, cp2_e, vec_wk2 );
    NPT_REAL asw = CalcInProduct( vec_wk1, vec_wk2 );
    if( asw < 0.0 ) {
        printf("#### ERROR calcControlPointEdge cp1_e,cp2_e order\n");
        printf("  p1 = %lg %lg %lg\n",p1[0],p1[1],p1[2]);
        printf("  p2 = %lg %lg %lg\n",p2[0],p2[1],p2[2]);
        printf("  cp1_e = %lg %lg %lg\n",cp1_e[0],cp1_e[1],cp1_e[2]);
        printf("  cp2_e = %lg %lg %lg\n",cp2_e[0],cp2_e[1],cp2_e[2]);
    }
#endif
}

// 中央の制御点取得（３次多項式用）
void
npt_param_calcControlPointCenter(
           NPT_REAL        p1[3],        // [in]  頂点１座標
           NPT_REAL        p2[3],        // [in]  頂点２座標
           NPT_REAL        p3[3],        // [in]  頂点３座標
           NPT_REAL        cp1_p1p2[3],  // [in]  p1p2辺の制御点1
           NPT_REAL        cp2_p1p2[3],  // [in]  p1p2辺の制御点2
           NPT_REAL        cp1_p2p3[3],  // [in]  p2p3辺の制御点1
           NPT_REAL        cp2_p2p3[3],  // [in]  p2p3辺の制御点2
           NPT_REAL        cp1_p3p1[3],  // [in]  p3p1辺の制御点1
           NPT_REAL        cp2_p3p1[3],  // [in]  p3p1辺の制御点2
           NPT_REAL        cp_center[3]  // [out] 中央制御点
      )
{
    cp_center[0] =   (   cp1_p1p2[0] + cp2_p1p2[0]
                       + cp1_p2p3[0] + cp2_p2p3[0]
                       + cp1_p3p1[0] + cp2_p3p1[0] ) / 4.0
                   - (   p1[0] + p2[0] + p3[0]     ) / 6.0;

    cp_center[1] =   (   cp1_p1p2[1] + cp2_p1p2[1]
                       + cp1_p2p3[1] + cp2_p2p3[1]
                       + cp1_p3p1[1] + cp2_p3p1[1] ) / 4.0
                   - (   p1[1] + p2[1] + p3[1]     ) / 6.0;

    cp_center[2] =   (   cp1_p1p2[2] + cp2_p1p2[2]
                       + cp1_p2p3[2] + cp2_p2p3[2]
                       + cp1_p3p1[2] + cp2_p3p1[2] ) / 4.0
                   - (   p1[2] + p2[2] + p3[2]     ) / 6.0;

}


// 制御点p11取得
void
npt_param_calcP11(
           NPT_REAL        p1[3],        // [in]  頂点１座標
           NPT_REAL        norm1[3],     // [in]  頂点１ベクトル（単位ベクトル）
           NPT_REAL        d1,           // [in]  頂点１ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
           NPT_REAL        p2[3],        // [in]  頂点２座標
           NPT_REAL        norm2[3],     // [in]  頂点２ベクトル（単位ベクトル）
           NPT_REAL        d2,           // [in]  頂点２ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
           NPT_REAL        norm_base[3], // [in]  曲線平面II（基準面）の法線ベクトル
           NPT_REAL        p11[3]        // [out] 制御点座標
   )
{
    bool bRet;
    NPT_REAL pos_line[3], vec_line[3];  // 面の交線の通過点, 面の交線のベクトル
    NPT_REAL pos_m[3];  // p1-p2,p2-p3,p3-p1辺の中点
    NPT_REAL pos_x[3];  // 辺の中点から垂線に下した交点
    NPT_REAL d1_base;   // 基準面 p1での原点からの距離(+-)
    NPT_REAL d_pos_x_base;   // 基準面 pos_xでの原点からの距離(+-)


    //------------------------------------------------
    // 平面の並行判定
    //------------------------------------------------
    NPT_REAL asw = CalcInProduct( norm1, norm2 );
    if( (1.0-fabs(asw)) < NPT_ALW_V ) { // 平面が同方向
        // 制御点座標をp1,p2の中点とする
        p11[0] = ( p1[0] + p2[0] ) / 2.0;
        p11[1] = ( p1[1] + p2[1] ) / 2.0;
        p11[2] = ( p1[2] + p2[2] ) / 2.0;
        return;
    }

    //------------------------------------------------
    // p1,p2の接平面とp1->p2辺が重なると補間出来ないための処理
    //------------------------------------------------
    NPT_REAL vec_p1p2[3]; //  p1->p2辺のベクトル
    CalcVec( p1, p2, vec_p1p2 );
    CalcNormalize( vec_p1p2 );
    asw = CalcInProduct( vec_p1p2, norm1 );
    if( fabs(asw) < NPT_ALW_V ) {
        // 制御点座標をp1,p2の中点とする
        p11[0] = ( p1[0] + p2[0] ) / 2.0;
        p11[1] = ( p1[1] + p2[1] ) / 2.0;
        p11[2] = ( p1[2] + p2[2] ) / 2.0;
        return;
    }
    asw = CalcInProduct( vec_p1p2, norm2 );
    if( fabs(asw) < NPT_ALW_V ) {
        // 制御点座標をp1,p2の中点とする
        p11[0] = ( p1[0] + p2[0] ) / 2.0;
        p11[1] = ( p1[1] + p2[1] ) / 2.0;
        p11[2] = ( p1[2] + p2[2] ) / 2.0;
        return;
    }

    //------------------------------------------------
    // p1-p2 辺の処理
    //    p1とp2情報より制御点座標を求める
    //------------------------------------------------

    // p1とp2の接平面より交線を求める
    bRet = CalcIntersectionLine(
                  norm1, d1,    // [in]  p1 接平面情報
                  norm2, d2,    // [in]  p2 接平面情報
                  pos_line,     // [out] 面の交線の通過点
                  vec_line      // [out] 面の交線のベクトル
              );

    if( bRet )   {
    } else {
        // 交線なし
        // 制御点座標をp1,p2の中点とする
        p11[0] = ( p1[0] + p2[0] ) / 2.0;
        p11[1] = ( p1[1] + p2[1] ) / 2.0;
        p11[2] = ( p1[2] + p2[2] ) / 2.0;
        return;
    }

    // 辺の中点
    pos_m[0] = ( p1[0] + p2[0] ) / 2.0;
    pos_m[1] = ( p1[1] + p2[1] ) / 2.0;
    pos_m[2] = ( p1[2] + p2[2] ) / 2.0;

    // 辺の中点から垂線に下した交点
    CalcNearPosOnLine(
                 pos_m,      // [in]  点座標
                 pos_line,   // [in]  線分の通過点
                 vec_line,   // [in]  線分の方向ベクトル（単位ベクトル）
                 pos_x       // [out] 点から線分に下した垂線との交点
           );

    //------------------------------------------------
    // 曲線平面II（基準面）に投影する
    //------------------------------------------------

    //      曲線平面II（基準面） p1の原点からの距離を求める
    d1_base = CalcPlaneD( p1, norm_base );

    //      曲線平面II（基準面） pos_xの原点からの距離を求める
    d_pos_x_base = CalcPlaneD( pos_x, norm_base );

    //      基準面に投影する
    p11[0] = pos_x[0] + norm_base[0]*( d1_base - d_pos_x_base );
    p11[1] = pos_x[1] + norm_base[1]*( d1_base - d_pos_x_base );
    p11[2] = pos_x[2] + norm_base[2]*( d1_base - d_pos_x_base );
}


// 制御点p11逆行対応補正（３次多項式用）
void
npt_param_correctP11(
           NPT_REAL        p11[3],       // [in]  制御点座標
           NPT_REAL        p1[3],        // [in]  頂点１座標
           NPT_REAL        p2[3],        // [in]  頂点２座標
           NPT_REAL        norm_base[3], // [in]  曲線平面II（基準面）の法線ベクトル
           NPT_REAL        p11_0[3],     // [out] 制御点座標（補正後 1番目の制御点用）
           NPT_REAL        p11_1[3]      // [out] 制御点座標（補正後 2番目の制御点用）
   )
{
    int  mode;   // -1: 始点側逆行  0: p1-p2内補正不要  1:終点側逆行
    NPT_REAL vec_p1_p2[3];    // p1->p2ベクトル
    NPT_REAL vec_p1_p11[3];   // p1->p11ベクトル
    NPT_REAL vec_p2_p11[3];   // p2->p11ベクトル
    NPT_REAL vec_p1_p11_mirror[3];   // p1->p11ベクトルのミラー
    NPT_REAL vec_p2_p11_mirror[3];   // p2->p11ベクトルのミラー
    NPT_REAL pos_x2[3];
    NPT_REAL pos_wk[3];
    NPT_REAL len_wk;
    NPT_REAL asw;
    bool   bRet;

    // デフォルト
    p11_0[0] = p11[0]; p11_0[1] = p11[1]; p11_0[2] = p11[2];
    p11_1[0] = p11[0]; p11_1[1] = p11[1]; p11_1[2] = p11[2];

    // p1->p2ベクトル
    bRet = CalcLineVec( p1, p2, vec_p1_p2, &len_wk );
    if( !bRet ) return;   // 同一点

    // p1->p11ベクトル
    bRet = CalcLineVec( p1, p11, vec_p1_p11, &len_wk );
    if( !bRet ) return;   // 同一点

    // p2->p11ベクトル
    bRet = CalcLineVec( p2, p11, vec_p2_p11, &len_wk );
    if( !bRet ) return;   // 同一点

    // modeの決定
    asw = CalcInProduct( vec_p1_p2, vec_p1_p11 );
    if( asw < 0.0 ) {
        mode = -1;   // -1: 始点側逆行
    } else {
        asw = CalcInProduct( vec_p1_p2, vec_p2_p11 );
        if( asw > 0.0 ) {
            mode = 1;   // 1:終点側逆行
        } else {
            mode = 0;   // 補正なし
            p11_0[0] = p11[0]; p11_0[1] = p11[1]; p11_0[2] = p11[2];
            p11_1[0] = p11[0]; p11_1[1] = p11[1]; p11_1[2] = p11[2];
            return;  // 逆行していないため終了
        }
    }

    // p1->p11ベクトル のミラー
    CalcVecMirror(
                  vec_p1_p2,    // 回転軸のベクトル
                  vec_p1_p11,   //
                  vec_p1_p11_mirror // ミラーしたベクトル
               );

    // p2->p11ベクトル のミラー
    CalcVecMirror(
                  vec_p1_p2,    // 回転軸のベクトル
                  vec_p2_p11,   //
                  vec_p2_p11_mirror // ミラーしたベクトル
               );

    // 始点側逆行補正
    if( mode == -1 ) {
        // ---- p11_0 制御点座標（補正後 1番目の制御点用）設定 ------

        // サブルーチンの引数に合わせるため適当に終点位置を決める
        pos_wk[0] = p2[0] + vec_p2_p11_mirror[0]*1.0;
        pos_wk[1] = p2[1] + vec_p2_p11_mirror[1]*1.0;
        pos_wk[2] = p2[2] + vec_p2_p11_mirror[2]*1.0;

        // p11_0 制御点座標（補正後 1番目の制御点用）を求める
        //    線分と線分の交点
        bRet = CalcCrossPointLine(
                     p1, p11,       // 線分１
                     p2, pos_wk,    // 線分２
                     p11_0,         // 線分１上の交点
                     pos_x2         // 線分２上の交点
                 );
        if( !bRet ) {
            printf("#### WARNING correctP11:CalcCrossPointLine() 1-1\n");
            printf("  p1 = %lg %lg %lg\n",p1[0],p1[1],p1[2]);
            printf("  p11= %lg %lg %lg\n",p11[0],p11[1],p11[2]);
            printf("  p2 = %lg %lg %lg\n",p2[0],p2[1],p2[2]);
            printf("  pos_wk = %lg %lg %lg\n",pos_wk[0],pos_wk[1],pos_wk[2]);

            // 安全のため制御点を中点にする
            p11_0[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_0[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_0[2] = ( p1[2] + p2[2] ) / 2.0;
        }

        // ---- p11_1 制御点座標（補正後 2番目の制御点用）設定 ------

        // サブルーチンの引数に合わせるため適当に終点位置を決める
        pos_wk[0] = p1[0] + vec_p1_p11_mirror[0]*1.0;
        pos_wk[1] = p1[1] + vec_p1_p11_mirror[1]*1.0;
        pos_wk[2] = p1[2] + vec_p1_p11_mirror[2]*1.0;

        // p11_1 制御点座標（補正後 1番目の制御点用）を求める
        //    線分と線分の交点
        bRet = CalcCrossPointLine(
                     p2, p11,       // 線分１
                     p1, pos_wk,    // 線分２
                     p11_1,         // 線分１上の交点
                     pos_x2         // 線分２上の交点
                 );
        if( !bRet ) {
            printf("#### WARNING correctP11:CalcCrossPointLine() 1-2\n");
            printf("  p2 = %lg %lg %lg\n",p2[0],p2[1],p2[2]);
            printf("  p11= %lg %lg %lg\n",p11[0],p11[1],p11[2]);
            printf("  p1 = %lg %lg %lg\n",p1[0],p1[1],p1[2]);
            printf("  pos_wk = %lg %lg %lg\n",pos_wk[0],pos_wk[1],pos_wk[2]);

            // 安全のため制御点を中点にする
            p11_1[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_1[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_1[2] = ( p1[2] + p2[2] ) / 2.0;
        }
    }
    else  // mode==1
    // 終点側逆行補正
    //     始点側逆行補正と同じにしてある。問題あるなら変更する
    if ( mode == 1 ) {
        // ---- p11_0 制御点座標（補正後 1番目の制御点用）設定 ------

        // サブルーチンの引数に合わせるため適当に終点位置を決める
        pos_wk[0] = p2[0] + vec_p2_p11_mirror[0]*1.0;
        pos_wk[1] = p2[1] + vec_p2_p11_mirror[1]*1.0;
        pos_wk[2] = p2[2] + vec_p2_p11_mirror[2]*1.0;

        // p11_0 制御点座標（補正後 1番目の制御点用）を求める
        //    線分と線分の交点
        bRet = CalcCrossPointLine(
                     p1, p11,       // 線分１
                     p2, pos_wk,    // 線分２
                     p11_0,         // 線分１上の交点
                     pos_x2         // 線分２上の交点
                 );
        if( !bRet ) {
            printf("#### WARNING correctP11:CalcCrossPointLine() 2-1\n");
            printf("  p1 = %lg %lg %lg\n",p1[0],p1[1],p1[2]);
            printf("  p11= %lg %lg %lg\n",p11[0],p11[1],p11[2]);
            printf("  p2 = %lg %lg %lg\n",p2[0],p2[1],p2[2]);
            printf("  pos_wk = %lg %lg %lg\n",pos_wk[0],pos_wk[1],pos_wk[2]);

            // 安全のため制御点を中点にする
            p11_0[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_0[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_0[2] = ( p1[2] + p2[2] ) / 2.0;
        }

        // ---- p11_1 制御点座標（補正後 2番目の制御点用）設定 ------

        // サブルーチンの引数に合わせるため適当に終点位置を決める
        pos_wk[0] = p1[0] + vec_p1_p11_mirror[0]*1.0;
        pos_wk[1] = p1[1] + vec_p1_p11_mirror[1]*1.0;
        pos_wk[2] = p1[2] + vec_p1_p11_mirror[2]*1.0;

        // p11_1 制御点座標（補正後 1番目の制御点用）を求める
        //    線分と線分の交点
        bRet = CalcCrossPointLine(
                     p2, p11,       // 線分１
                     p1, pos_wk,    // 線分２
                     p11_1,         // 線分１上の交点
                     pos_x2         // 線分２上の交点
                 );
        if( !bRet ) {
            printf("#### WARNING correctP11:CalcCrossPointLine() 2-2\n");
            printf("  p2 = %lg %lg %lg\n",p2[0],p2[1],p2[2]);
            printf("  p11= %lg %lg %lg\n",p11[0],p11[1],p11[2]);
            printf("  p1 = %lg %lg %lg\n",p1[0],p1[1],p1[2]);
            printf("  pos_wk = %lg %lg %lg\n",pos_wk[0],pos_wk[1],pos_wk[2]);

            // 安全のため制御点を中点にする
            p11_1[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_1[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_1[2] = ( p1[2] + p2[2] ) / 2.0;
        }
    }

}
