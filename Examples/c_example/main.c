
/*
 * Nagata Patch Library Test
 *
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// Cインターフェース サンプル
///   三角形の頂点列データを入力として、以下の処理を行う
///     ・頂点の法線ベクトルを求める
///     ・長田パッチの生成
///     ・三角形の各辺の中点の長田パッチ上の曲面補間点を求める
///     ・三角形の頂点と上記の曲面補間点より、三角形数を４倍としたデータを生成する
///   - 当サンプルは頂点の法線ベクトルを求めるのに単純にループさせているが
///     三角形数が多い場合には処理時間がかかるようになるので注意
///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Npt.h"

#define NMAX 20

// 面の法線ベクトル
void get_tri_normal( NPT_REAL p1[3], NPT_REAL p2[3], NPT_REAL p3[3], NPT_REAL norm[3] ) 
{
    GEO_REAL vec12[3], vec13[3];
    vec12[0]=p2[0]-p1[0];  vec12[1]=p2[1]-p1[1]; vec12[2]=p2[2]-p1[2];
    vec13[0]=p3[0]-p1[0];  vec13[1]=p3[1]-p1[1]; vec13[2]=p3[2]-p1[2];
    CalcOutProduct( vec12, vec13, norm );
    CalcNormalize( norm );
}

// STLファイル出力
void output_stl_file( char* file_name, int num, 
                      NPT_REAL tri[][3][3], NPT_REAL norm[][3] ) 
{
    FILE* fp;
    int i,j;

    fp = fopen( file_name, "w");
    if( fp == NULL )  {
        printf("#### Error: stl file open error file_name=%s\n",file_name);
        exit(1);
    }
    fprintf(fp,"solid model1\n");
    for(i=0; i<num; i++ ) {
        fprintf(fp,"    facet normal %f %f %f\n",norm[i][0],norm[i][1],norm[i][2]);
        fprintf(fp,"    outer loop\n");
        for(j=0; j<3; j++ ) {
            fprintf(fp,"        vertex %f %f %f\n",tri[i][j][0],tri[i][j][1],tri[i][j][2]);
        }
        fprintf(fp,"    endloop\n");
        fprintf(fp,"    endfacet\n");
    }
    fprintf(fp,"endsolid model1\n");

    fclose(fp);
}

// NPTファイル出力
void output_npt_file( char* file_name, int num, 
                      NPT_REAL tri[][3][3], NPT_REAL npatch[][7][3] ) 
{
    FILE* fp;
    int i,j;

    fp = fopen( file_name, "w");
    if( fp == NULL )  {
        printf("#### Error: npt file open error file_name=%s\n",file_name);
        exit(1);
    }
    fprintf(fp,"%d\n",num);
    for(i=0; i<num; i++ ) {
        fprintf(fp,"facet\n");
        for(j=0; j<3; j++ ) {
            fprintf(fp,"    vertex %f %f %f\n",tri[i][j][0],tri[i][j][1],tri[i][j][2]);
        }
        for(j=0; j<7; j++ ) {
            fprintf(fp,"    coef%d %f %f %f\n",(j+1),npatch[i][j][0],npatch[i][j][1],npatch[i][j][2]);
        }
    }

    fclose(fp);
}


//----------------------------------------------------
//  メインルーチン
//----------------------------------------------------

int main(int argc, char** argv )
{
    NPT_REAL  p0[3];
    NPT_REAL  p1[3];
    NPT_REAL  p2[3];
    NPT_REAL  p3[3];
    NPT_REAL  p4[3];
    NPT_REAL  p5[3];
    NPT_REAL  p6[3];
    NPT_REAL  p7[3];
#if 0
    NPT_REAL  midpnt[3]; // 中間点
#endif

    int       num_tri;
    NPT_REAL  tri       [NMAX][3][3];  // 三角形の座標
    NPT_REAL  plane_norm[NMAX][3];     // 三角形の法線ベクトル
    NPT_REAL  vtx_norm  [NMAX][3][3];  // 三角形の頂点の法線ベクトル
    NPT_REAL  npatch    [NMAX][7][3];  // 長田パッチパラメータ

    int i,j, iret;
    //int ip,jp,iflg;
    int ip,jp;
    NPT_REAL  norm_tmp[3];
    NPT_REAL  eps = 0.01;
    int num_4;
    NPT_REAL  eta, xi;
    NPT_REAL  p12[3],p23[3],p31[3];
    NPT_REAL  tri_4      [4*NMAX][3][3];  // 長田パッチ補間により４倍したポリゴン
    NPT_REAL  plane_norm4[4*NMAX][3];     // ４倍したポリゴンの法線ベクトル

    char* file_name_stl_in   ="stl_in.stl";
    char* file_name_stl4_out ="stl4_out.stl";
    char* file_name_npt_out  ="npt_out.npt";

    printf( "#### Program Start\n");

    // 頂点の設定
#if 0
    // 三角形３個（開放端あり）
    num_tri = 3;
    p0[0]=  0.0;  p0[1]=  0.0;  p0[2]=  0.0;
    p1[0]=100.0;  p1[1]=  0.0;  p1[2]=  0.0;
    p2[0]=  0.0;  p2[1]=100.0;  p2[2]=  0.0;
    p3[0]=  0.0;  p3[1]=  0.0;  p3[2]=100.0;

    tri[0][0][0]=p0[0];  tri[0][0][1]=p0[1];  tri[0][0][2]=p0[2];
    tri[0][1][0]=p2[0];  tri[0][1][1]=p2[1];  tri[0][1][2]=p2[2];
    tri[0][2][0]=p1[0];  tri[0][2][1]=p1[1];  tri[0][2][2]=p1[2];

    tri[1][0][0]=p0[0];  tri[1][0][1]=p0[1];  tri[1][0][2]=p0[2];
    tri[1][1][0]=p1[0];  tri[1][1][1]=p1[1];  tri[1][1][2]=p1[2];
    tri[1][2][0]=p3[0];  tri[1][2][1]=p3[1];  tri[1][2][2]=p3[2];

    tri[2][0][0]=p0[0];  tri[2][0][1]=p0[1];  tri[2][0][2]=p0[2];
    tri[2][1][0]=p3[0];  tri[2][1][1]=p3[1];  tri[2][1][2]=p3[2];
    tri[2][2][0]=p2[0];  tri[2][2][1]=p2[1];  tri[2][2][2]=p2[2];
#endif

#if 1
    // 直方体
    num_tri = 12;
    p0[0]=  0.0;  p0[1]=  0.0;  p0[2]=  0.0;
    p1[0]=100.0;  p1[1]=  0.0;  p1[2]=  0.0;
    p2[0]=100.0;  p2[1]=100.0;  p2[2]=  0.0;
    p3[0]=  0.0;  p3[1]=100.0;  p3[2]=  0.0;
    p4[0]=  0.0;  p4[1]=  0.0;  p4[2]=100.0;
    p5[0]=100.0;  p5[1]=  0.0;  p5[2]=100.0;
    p6[0]=100.0;  p6[1]=100.0;  p6[2]=100.0;
    p7[0]=  0.0;  p7[1]=100.0;  p7[2]=100.0;
    // 三角形の設定
    //     底面 p0->p3->p2->p1
    tri[0][0][0]=p0[0];  tri[0][0][1]=p0[1];  tri[0][0][2]=p0[2];
    tri[0][1][0]=p3[0];  tri[0][1][1]=p3[1];  tri[0][1][2]=p3[2];
    tri[0][2][0]=p1[0];  tri[0][2][1]=p1[1];  tri[0][2][2]=p1[2];

    tri[1][0][0]=p2[0];  tri[1][0][1]=p2[1];  tri[1][0][2]=p2[2];
    tri[1][1][0]=p1[0];  tri[1][1][1]=p1[1];  tri[1][1][2]=p1[2];
    tri[1][2][0]=p3[0];  tri[1][2][1]=p3[1];  tri[1][2][2]=p3[2];
    //     上面 p4->p5->p6->p7
    tri[2][0][0]=p4[0];  tri[2][0][1]=p4[1];  tri[2][0][2]=p4[2];
    tri[2][1][0]=p5[0];  tri[2][1][1]=p5[1];  tri[2][1][2]=p5[2];
    tri[2][2][0]=p7[0];  tri[2][2][1]=p7[1];  tri[2][2][2]=p7[2];

    tri[3][0][0]=p6[0];  tri[3][0][1]=p6[1];  tri[3][0][2]=p6[2];
    tri[3][1][0]=p7[0];  tri[3][1][1]=p7[1];  tri[3][1][2]=p7[2];
    tri[3][2][0]=p5[0];  tri[3][2][1]=p5[1];  tri[3][2][2]=p5[2];
    //     Front  p0->p1->p5->p4
    tri[4][0][0]=p0[0];  tri[4][0][1]=p0[1];  tri[4][0][2]=p0[2];
    tri[4][1][0]=p1[0];  tri[4][1][1]=p1[1];  tri[4][1][2]=p1[2];
    tri[4][2][0]=p4[0];  tri[4][2][1]=p4[1];  tri[4][2][2]=p4[2];

    tri[5][0][0]=p5[0];  tri[5][0][1]=p5[1];  tri[5][0][2]=p5[2];
    tri[5][1][0]=p4[0];  tri[5][1][1]=p4[1];  tri[5][1][2]=p4[2];
    tri[5][2][0]=p1[0];  tri[5][2][1]=p1[1];  tri[5][2][2]=p1[2];
    //     Right  p1->p2->p6->p5
    tri[6][0][0]=p1[0];  tri[6][0][1]=p1[1];  tri[6][0][2]=p1[2];
    tri[6][1][0]=p2[0];  tri[6][1][1]=p2[1];  tri[6][1][2]=p2[2];
    tri[6][2][0]=p5[0];  tri[6][2][1]=p5[1];  tri[6][2][2]=p5[2];

    tri[7][0][0]=p6[0];  tri[7][0][1]=p6[1];  tri[7][0][2]=p6[2];
    tri[7][1][0]=p5[0];  tri[7][1][1]=p5[1];  tri[7][1][2]=p5[2];
    tri[7][2][0]=p2[0];  tri[7][2][1]=p2[1];  tri[7][2][2]=p2[2];
    //     Back  p2->p3->p7->p6
    tri[8][0][0]=p2[0];  tri[8][0][1]=p2[1];  tri[8][0][2]=p2[2];
    tri[8][1][0]=p3[0];  tri[8][1][1]=p3[1];  tri[8][1][2]=p3[2];
    tri[8][2][0]=p6[0];  tri[8][2][1]=p6[1];  tri[8][2][2]=p6[2];

    tri[9][0][0]=p7[0];  tri[9][0][1]=p7[1];  tri[9][0][2]=p7[2];
    tri[9][1][0]=p6[0];  tri[9][1][1]=p6[1];  tri[9][1][2]=p6[2];
    tri[9][2][0]=p3[0];  tri[9][2][1]=p3[1];  tri[9][2][2]=p3[2];
    //     Left  p3->p0->p4->p7
    tri[10][0][0]=p3[0];  tri[10][0][1]=p3[1];  tri[10][0][2]=p3[2];
    tri[10][1][0]=p0[0];  tri[10][1][1]=p0[1];  tri[10][1][2]=p0[2];
    tri[10][2][0]=p7[0];  tri[10][2][1]=p7[1];  tri[10][2][2]=p7[2];

    tri[11][0][0]=p4[0];  tri[11][0][1]=p4[1];  tri[11][0][2]=p4[2];
    tri[11][1][0]=p7[0];  tri[11][1][1]=p7[1];  tri[11][1][2]=p7[2];
    tri[11][2][0]=p0[0];  tri[11][2][1]=p0[1];  tri[11][2][2]=p0[2];
#endif


    // 面の法線ベクトル設定
    for( i=0; i<num_tri; i++ ) {
        get_tri_normal( tri[i][0], tri[i][1], tri[i][2], plane_norm[i] ); 
    }
 
    // 確認のためSTLファイルに出力
    output_stl_file( file_name_stl_in, num_tri, tri, plane_norm );


    // 頂点ベクトルの設定
    //    点数が少ないので単純ループ
    for(i=0; i<num_tri; i++ ) {  // 3角形のloop
        for(ip=0; ip<3; ip++ ) { // 頂点のloop
           //printf("---- vtx_norm set i=%d ip=%d\n",i,ip);
           norm_tmp[0] = plane_norm[i][0]; // 自身の法線ベクトル設定
           norm_tmp[1] = plane_norm[i][1];
           norm_tmp[2] = plane_norm[i][2];

           for(j=0; j<num_tri; j++ ) { // 相手の3角形のloop
               if( i == j ) continue;
               //iflg = 0; 
               for(jp=0; jp<3; jp++ ) { // 相手の頂点のloop

                   if( fabs(tri[i][ip][0]-tri[j][jp][0]) < eps  && 
                       fabs(tri[i][ip][1]-tri[j][jp][1]) < eps  && 
                       fabs(tri[i][ip][2]-tri[j][jp][2]) < eps       ) {

                       //printf("   match j=%d jp=%d\n",j,jp);

                       norm_tmp[0] += plane_norm[j][0]; // 接している面の法線ベクトル
                       norm_tmp[1] += plane_norm[j][1];
                       norm_tmp[2] += plane_norm[j][2];
                       break;
                   }
               }
           }
           CalcNormalize2( norm_tmp, vtx_norm[i][ip] );
           // Debug : 頂点の法線ベクトル出力
           //printf("---- vtx_norm[%d][%d]= %f %f %f\n",i,ip,
           //              vtx_norm[i][ip][0],vtx_norm[i][ip][1],vtx_norm[i][ip][2] );
        }
    } 

    // 長田パッチ変換
    for(i=0; i<num_tri; i++ ) {  // 3角形のloop
       //printf("***** npt_param_crt() i=%d\n",i);
       iret = npt_param_crt (
                  tri[i][0], vtx_norm[i][0], 
                  tri[i][1], vtx_norm[i][1], 
                  tri[i][2], vtx_norm[i][2], 
                  npatch[i][0], npatch[i][1], npatch[i][2], npatch[i][3], 
                  npatch[i][4], npatch[i][5], npatch[i][6]
                );
       if( iret != 1 ) {
           printf("#### Error npt_param_crt() ret=%d i=%d\n",iret,i);
           exit(1);
       }
    }

    // 確認のため長田パッチファイルに出力
    output_npt_file( file_name_npt_out, num_tri, tri, npatch );

#if 0
    // eta,xi取得テスト
    //    三角形１の辺１の中点
    midpnt[0]=( tri[0][0][0]+tri[0][1][0])/2.0;
    midpnt[1]=( tri[0][0][1]+tri[0][1][1])/2.0;
    midpnt[2]=( tri[0][0][2]+tri[0][1][2])/2.0;
    npt_cvt_pos_to_eta_xi( midpnt, tri[0][0], tri[0][1], tri[0][2], &eta, &xi );
    printf("---- edge1 mid point eta=%f xi=%f\n",eta,xi);
    //    三角形１の辺２の中点
    midpnt[0]=( tri[0][1][0]+tri[0][2][0])/2.0;
    midpnt[1]=( tri[0][1][1]+tri[0][2][1])/2.0;
    midpnt[2]=( tri[0][1][2]+tri[0][2][2])/2.0;
    npt_cvt_pos_to_eta_xi( midpnt, tri[0][0], tri[0][1], tri[0][2], &eta, &xi );
    printf("---- edge2 mid point eta=%f xi=%f\n",eta,xi);
    //    三角形１の辺３の中点
    midpnt[0]=( tri[0][2][0]+tri[0][0][0])/2.0;
    midpnt[1]=( tri[0][2][1]+tri[0][0][1])/2.0;
    midpnt[2]=( tri[0][2][2]+tri[0][0][2])/2.0;
    npt_cvt_pos_to_eta_xi( midpnt, tri[0][0], tri[0][1], tri[0][2], &eta, &xi );
    printf("---- edge3 mid point eta=%f xi=%f\n",eta,xi);
#endif

    // 各辺の中点を補間して点を接点として追加
    //     ポリゴン数が４倍になる
    num_4 = 0;
    for(i=0; i<num_tri; i++ ) {  // 3角形のloop
        // 辺１の中点の曲面補間点
        eta=0.5; xi=0.0;
        npt_correct_pnt(
                eta, xi,
                tri[i][0], tri[i][1], tri[i][2],
                npatch[i][0], npatch[i][1], npatch[i][2], npatch[i][3],
                npatch[i][4], npatch[i][5], npatch[i][6],
                p12
             );
        //printf("---- correct i=%d p12= %f %f %f\n",i, p12[0],p12[1],p12[2]);

        // 辺２の中点の曲面補間点
        eta=1.0; xi=0.5;
        npt_correct_pnt(
                eta, xi,
                tri[i][0], tri[i][1], tri[i][2],
                npatch[i][0], npatch[i][1], npatch[i][2], npatch[i][3],
                npatch[i][4], npatch[i][5], npatch[i][6],
                p23
             );
        //printf("---- correct i=%d p23= %f %f %f\n",i, p23[0],p23[1],p23[2]);

        // 辺３の中点の曲面補間点
        eta=0.5; xi=0.5;
        npt_correct_pnt(
                eta, xi,
                tri[i][0], tri[i][1], tri[i][2],
                npatch[i][0], npatch[i][1], npatch[i][2], npatch[i][3],
                npatch[i][4], npatch[i][5], npatch[i][6],
                p31
             );
        //printf("---- correct i=%d p31= %f %f %f\n",i, p31[0],p31[1],p31[2]);


        // ポリゴン１個目
        tri_4[num_4][0][0]=tri[i][0][0];
        tri_4[num_4][0][1]=tri[i][0][1];
        tri_4[num_4][0][2]=tri[i][0][2]; 
        tri_4[num_4][1][0]=p12[0];
        tri_4[num_4][1][1]=p12[1];
        tri_4[num_4][1][2]=p12[2];
        tri_4[num_4][2][0]=p31[0];
        tri_4[num_4][2][1]=p31[1];
        tri_4[num_4][2][2]=p31[2];
        num_4 += 1;

        //ポリゴン２個目
        tri_4[num_4][0][0]=tri[i][1][0];
        tri_4[num_4][0][1]=tri[i][1][1];
        tri_4[num_4][0][2]=tri[i][1][2]; 
        tri_4[num_4][1][0]=p23[0];
        tri_4[num_4][1][1]=p23[1];
        tri_4[num_4][1][2]=p23[2];
        tri_4[num_4][2][0]=p12[0];
        tri_4[num_4][2][1]=p12[1];
        tri_4[num_4][2][2]=p12[2];
        num_4 += 1;

        //ポリゴン３個目
        tri_4[num_4][0][0]=tri[i][2][0];
        tri_4[num_4][0][1]=tri[i][2][1];
        tri_4[num_4][0][2]=tri[i][2][2]; 
        tri_4[num_4][1][0]=p31[0];
        tri_4[num_4][1][1]=p31[1];
        tri_4[num_4][1][2]=p31[2];
        tri_4[num_4][2][0]=p23[0];
        tri_4[num_4][2][1]=p23[1];
        tri_4[num_4][2][2]=p23[2];
        num_4 += 1;

        //ポリゴン４個目
        tri_4[num_4][0][0]=p12[0];
        tri_4[num_4][0][1]=p12[1];
        tri_4[num_4][0][2]=p12[2]; 
        tri_4[num_4][1][0]=p23[0];
        tri_4[num_4][1][1]=p23[1];
        tri_4[num_4][1][2]=p23[2];
        tri_4[num_4][2][0]=p31[0];
        tri_4[num_4][2][1]=p31[1];
        tri_4[num_4][2][2]=p31[2];
        num_4 += 1;
    }

    // 面の法線ベクトル設定
    for( i=0; i<num_4; i++ ) {
        get_tri_normal( tri_4[i][0], tri_4[i][1], tri_4[i][2], plane_norm4[i] ); 
    }

    // 曲面補間した４倍のポリゴンSTLファイルに出力
    output_stl_file( file_name_stl4_out, num_4, tri_4, plane_norm4 );



    printf( "#### Program End\n");

    return 0;
}
