/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_MeshPatch.h
 * @brief  Npatch MeshPatch Class Header
 * @author aics
 */

//********************************************
// class MeshPatch
//********************************************
#ifndef _NPT_MESH_PATCH_
#define _NPT_MESH_PATCH_

#include <math.h>


#define CALC_EPSILON_L  1.0e-5
#define CALC_EPSILON_V  1.0e-8


class MeshPatch
{

public:

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	MeshPatch(void);
	~MeshPatch(void);


	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	double m_alwAngle;

private:

	double	m_p1[3];        // 頂点１座標
	double	m_p2[3];        // 頂点２座標
	double	m_p3[3];        // 頂点３座標
	double	m_norm1[3];     // 頂点１法線ベクトル（単位ベクトル）
	double	m_norm2[3];     // 頂点２法線ベクトル（単位ベクトル）
	double	m_norm3[3];     // 頂点３法線ベクトル（単位ベクトル）
	double	m_d1;           // 頂点１接平面の原点からの距離
	double	m_d2;           // 頂点２接平面の原点からの距離
	double	m_d3;           // 頂点３接平面の原点からの距離

	// ++++++++++++ 3次多項式用
	double	m_cp1_p1p2[3];  // p1p2辺の３次ベジェ制御点1
	double	m_cp2_p1p2[3];  // p1p2辺の３次ベジェ制御点2
	double	m_cp1_p2p3[3];  // p2p3辺の３次ベジェ制御点1
	double	m_cp2_p2p3[3];  // p2p3辺の３次ベジェ制御点2
	double	m_cp1_p3p1[3];  // p3p1辺の３次ベジェ制御点1
	double	m_cp2_p3p1[3];  // p3p1辺の３次ベジェ制御点2
	double	m_cp_center[3]; // 中央の３次ベジェ制御点

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------

public:

	// ++++++++長田パッチ係数を求めるメイン関数　++++++++
	int CalcGetCoef2Deg(  double norm[3], double p1[3], double norm1[3], double p2[3], double norm2[3], 
						  double p3[3], double norm3[3], double C1[3], double C2[3], double C3[3] );

	// +++++++ ζ、ηを入力し、評価点を求める ++++++++++++
	int Evaluate2Deg( double eta, double zeta, double p1[3], double p2[3], double p3[3], 
					  double c1[3], double c2[3], double c3[3], double pos[3] );

private:

	// ++++++++++++++ 係数を求める ++++++++++++++++++
	void CalcPatchCoef(	double norm[3], double p1[3], double norm1[3], double d1, 
										double p2[3], double norm2[3], double d2, double cvec[3] );

	// ++++++++++++++ 交線を求める ++++++++++++++++++
	bool CalcIntersectionLine(	double vec1[3], double d1,  
								double vec2[3], double d2, 
								double pos[3], double vec[3] );

	// ++++++++++++++ Xmを求める ++++++++++++++++++
	void CalcNearPosOnLine( double pnt[3], double pos[3], double vec[3], double pos_x[3] );


	// ++++++++++++++++++++++++++++++++++++++++++++++
	//    　　　　三次ベジェ曲面
	// ++++++++++++++++++++++++++++++++++++++++++++++
	
	// ++++++++++++++ Xmを求める ++++++++++++++++++
public:

	int	CalcGetCoef3Deg( double norm[3],  double p1[3], double norm1[3], 
										  double p2[3], double norm2[3], 
										  double p3[3], double norm3[3], 
										  double C1[3], double C2[3], double C3[3], double C4[3],
										  double C5[3], double C6[3], double C7[3] );

	int Evaluate3Deg( double eta, double zeta, double p1[3], double p2[3], double p3[3], 
					  double c1[3], double c2[3], double c3[3], double c4[3], double c5[3], 
					  double c6[3], double c7[3], double pos_o[3] );

private:

	int		calcControlPoint_3deg( void );

	void	calcControlPointEdge( double p1[3], double norm1[3], double d1,
								  double p2[3], double norm2[3], double d2,        
								  double cp1_e[3], double cp2_e[3]      );


	void	calcControlPointCenter( double p1[3], double p2[3], double p3[3], 
									double cp1_p1p2[3], double cp2_p1p2[3], double cp1_p2p3[3],
						            double cp2_p2p3[3], double cp1_p3p1[3], double cp2_p3p1[3],
									double cp_center[3] );

	void	calcP11( double p1[3], double norm1[3], double d1,          
					 double p2[3], double norm2[3], double d2,        
					 double norm_base[3], double p11[3]      );


	void	correctP11( double p11[3], double p1[3], double p2[3],       
						double norm_base[3], double p11_0[3], double p11_1[3] );



	// ----------------------------------------------
	//         inline 関数
	// ----------------------------------------------
public:
	double	CalcVecSize( double vec[3] );
	double	CalcInProduct( double vec1[3], double vec2[3] );

private:
	void	CalcOutProduct( double vec1[3], double vec2[3], double vec_0[3] );
	double	CalcPlaneD( double pos[3], double vec[3] );

	void	CalcNormalize( double vec[3] );
	bool	CalcNormalize( double vec[3], double vec_o[3] );
	bool	CalcLineVec( double pp[3], double lp[3], double vec[3], double& length );
	void	CalcVecMirror( double rot_vec[3], double veci[3], double veco[3] );
	void	CalcVec( double pp[3], double lp[3], double vec[3] );
	bool	CalcCrossPointLine( double pp1[3], double lp1[3], double pp2[3], double lp2[3], 
								double pos_x1[3], double pos_x2[3] );


};

// ++++++++++++++ ベクトルのサイズ（長さ）を求める ++++++++++++++++++
inline double MeshPatch::CalcVecSize( double vec[3] )
{
	return sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}

// ++++++++ 点と方向ベクトルより、平面方程式Dを求める +++++++++++++++
inline double MeshPatch::CalcPlaneD( double pos[3], double vec[3] )
{
	return CalcInProduct( pos, vec );
}

// +++++++++++++++++++++++ ベクトルの内積 ++++++++++++++++++++++++++++
inline double MeshPatch::CalcInProduct( double vec1[3], double vec2[3] )
{
	return ( vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2] );
}

// +++++++++++++++++++++++ ベクトルの外積 ++++++++++++++++++++++++++++
inline void MeshPatch::CalcOutProduct( double vec1[3], double vec2[3], double vec_0[3] )
{
	vec_0[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
	vec_0[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
	vec_0[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];	
}

// +++++++++++++++++++++++ CalcNormalize ++++++++++++++++++++++++++++
inline void MeshPatch::CalcNormalize( double vec[3] )
{
    double len = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
    vec[0] /= len;
    vec[1] /= len;
    vec[2] /= len;
}

// +++++++++++++++++++++++ CalcNormalize ++++++++++++++++++++++++++++
inline bool MeshPatch::CalcNormalize( double vec[3], double vec_o[3] )
{
    double len = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
    if( len > CALC_EPSILON_L ) {
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

// +++++++++++++++++++++++ CalcLineVec ++++++++++++++++++++++++++++
inline bool MeshPatch::CalcLineVec( double pp[3], double lp[3], double vec[3], double& length )
{
    vec[0] = lp[0] - pp[0];
    vec[1] = lp[1] - pp[1];
    vec[2] = lp[2] - pp[2];

    length = sqrt ( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );

    if( length > CALC_EPSILON_L ) {
        vec[0] /= length;
        vec[1] /= length;
        vec[2] /= length;
        return true;
    } else {
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
        length = 0.0;
        return false;
    }
}

// +++++++++++++++++++++++ CalcVecMirror ++++++++++++++++++++++++++++
inline void MeshPatch::CalcVecMirror( double rot_vec[3], double veci[3], double veco[3] )
{
    double pos_x[ 3 ];

    // ベクトルを線分とみなして、回転軸ベクトル上の交点を求める
    double dist_wk = CalcInProduct( rot_vec, veci );
    
    // 点から線分に下した垂線との交点
    pos_x[0] = rot_vec[0]*dist_wk;
    pos_x[1] = rot_vec[1]*dist_wk;
    pos_x[2] = rot_vec[2]*dist_wk;

    // ミラーベクトル設定
    veco[0] = pos_x[0] + ( pos_x[0] - veci[0] );
    veco[1] = pos_x[1] + ( pos_x[1] - veci[1] );
    veco[2] = pos_x[2] + ( pos_x[2] - veci[2] );
}

// +++++++++++++++++++++++ CalcVec ++++++++++++++++++++++++++++
inline void MeshPatch::CalcVec( double pp[3], double lp[3], double vec[3] )
{
    vec[0] = lp[0] - pp[0];
    vec[1] = lp[1] - pp[1];
    vec[2] = lp[2] - pp[2];
}


// +++++++++++++++++++++++ CalcCrossPointLine ++++++++++++++++++++++++++++
inline bool MeshPatch::CalcCrossPointLine( double pp1[3],	 double lp1[3], 
										   double pp2[3],	 double lp2[3],   
										   double pos_x1[3], double pos_x2[3] )
{
    bool   bRet;
    double vec1[3], uni_vec1[3];
    double vec2[3], uni_vec2[3];
    double wk1, wk2;
    double vec_pp1_pp2[3];
    
    CalcVec( pp1, lp1, vec1 );
    bRet = CalcNormalize( vec1, uni_vec1 );
    if( !bRet ) {
        return false;
    }

    CalcVec( pp2, lp2, vec2 );
    bRet = CalcNormalize( vec2, uni_vec2 );
    if( !bRet ) {
        return false;
    }

    wk1 = CalcInProduct( uni_vec1, uni_vec2 );
    wk2 = 1.0 - wk1*wk1;
    if( wk2 < CALC_EPSILON_V ) {
        return false;  // 線分が並行
    }
    
    CalcVec( pp1, pp2, vec_pp1_pp2 );
    
    double d1 = (   CalcInProduct(vec_pp1_pp2,uni_vec1) 
                  - wk1*CalcInProduct(vec_pp1_pp2,uni_vec2) ) / wk2;
    
    double d2 = (   wk1*CalcInProduct(vec_pp1_pp2,uni_vec1) 
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

#endif // _NPT_MESH_PATCH_