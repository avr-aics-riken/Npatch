/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_MeshPatch.cxx
 * @brief  Npatch MeshPatch Class
 * @author aics
 */


#include "npt_MeshPatch.h"


/**
 * @brief           コンストラクタ
 */
MeshPatch::MeshPatch(void)
{

	m_alwAngle = 0.0031;
}

/**
 * @brief           デストラクタ
 */
MeshPatch::~MeshPatch(void)
{
}


/**
 * @brief           長田パッチの係数を求めるメイン関数
 * @param[in]       norm    面の法線ベクトル（単位ベクトル）
 * @param[in]       p1      頂点１の座標
 * @param[in]       norm1   頂点１の法線ベクトル（単位ベクトル）
 * @param[in]       p2      頂点２の座標
 * @param[in]       norm2   頂点２の法線ベクトル（単位ベクトル）
 * @param[in]       p3      頂点３の座標
 * @param[in]       norm3   頂点３の法線ベクトル（単位ベクトル）
 * @param[in,out]   C1      長田パッチの係数C１
 * @param[in,out]   C2      長田パッチの係数C２
 * @param[in,out]   C3      長田パッチの係数C３
 * @return          1
 */
int MeshPatch::CalcGetCoef2Deg(	double norm[3], 
								double p1[3], double norm1[3], 
								double p2[3], double norm2[3], 
								double p3[3], double norm3[3],
								double C1[3], double C2[3],	double C3[3] )
{

	double d1, d2, d3;

	d1 = CalcPlaneD( p1, norm1 );
	d2 = CalcPlaneD( p2, norm2 );
	d3 = CalcPlaneD( p3, norm3 );

	// ++++++++++ 2次曲面　+++++++++++++++++
	CalcPatchCoef( norm, p1, norm1, d1, p2, norm2, d2, C1 );
	CalcPatchCoef( norm, p2, norm2, d2, p3, norm3, d3, C2 );
	CalcPatchCoef( norm, p3, norm3, d3, p1, norm1, d1, C3 );

	return 1;
}

/**
 * @brief           長田パッチの係数を求めるメイン関数
 * @param[in]       norm    法線
 * @param[in]       p1      頂点1
 * @param[in]       norm1   頂点1の法線ベクトル
 * @param[in]       p2      頂点2
 * @param[in]       norm2   頂点2の法線ベクトル
 * @param[in]       p3      頂点3
 * @param[in]       norm3   頂点3の法線ベクトル
 * @param[out]      C1      3次曲面頂点1
 * @param[out]      C2      3次曲面頂点2
 * @param[out]      C3      3次曲面頂点3
 * @param[out]      C4      3次曲面頂点4
 * @param[out]      C5      3次曲面頂点5
 * @param[out]      C6      3次曲面頂点6
 * @param[out]      C7      3次曲面頂点7
 * @return          1
 */
int MeshPatch::CalcGetCoef3Deg(	double norm[3], 
							double p1[3], double norm1[3], 
							double p2[3], double norm2[3], 
							double p3[3], double norm3[3],
							double C1[3], double C2[3], double C3[3], double C4[3],
							double C5[3], double C6[3], double C7[3] )
{

   // メンバー変数に転送
    m_p1[0] = p1[0];   m_p1[1] = p1[1];  m_p1[2] = p1[2];
    m_p2[0] = p2[0];   m_p2[1] = p2[1];  m_p2[2] = p2[2];
    m_p3[0] = p3[0];   m_p3[1] = p3[1];  m_p3[2] = p3[2];

    m_norm1[0] = norm1[0]; m_norm1[1] = norm1[1]; m_norm1[2] = norm1[2];
    m_norm2[0] = norm2[0]; m_norm2[1] = norm2[1]; m_norm2[2] = norm2[2];
    m_norm3[0] = norm3[0]; m_norm3[1] = norm3[1]; m_norm3[2] = norm3[2];
    
    // 接平面：平面の方程式 D 計算
    //    平面の方程式  Ax + By + Cz = D
    m_d1 = CalcPlaneD( p1, norm1 );
    m_d2 = CalcPlaneD( p2, norm2 );
    m_d3 = CalcPlaneD( p3, norm3 );

	// ++++++++++ 3次曲面　+++++++++++++++++
	calcControlPoint_3deg();

	C1[0] = m_cp1_p1p2[0];
	C1[1] = m_cp1_p1p2[1];
	C1[2] = m_cp1_p1p2[2];

	C2[0] = m_cp2_p1p2[0];
	C2[1] = m_cp2_p1p2[1];
	C2[2] = m_cp2_p1p2[2];

	C3[0] = m_cp1_p2p3[0];
	C3[1] = m_cp1_p2p3[1];
	C3[2] = m_cp1_p2p3[2];

	C4[0] = m_cp2_p2p3[0];
	C4[1] = m_cp2_p2p3[1];
	C4[2] = m_cp2_p2p3[2];

	C5[0] = m_cp1_p3p1[0];
	C5[1] = m_cp1_p3p1[1];
	C5[2] = m_cp1_p3p1[2];

	C6[0] = m_cp2_p3p1[0];
	C6[1] = m_cp2_p3p1[1];
	C6[2] = m_cp2_p3p1[2];

	C7[0] = m_cp_center[0];
	C7[1] = m_cp_center[1];
	C7[2] = m_cp_center[2];

	return 1;

}


/**
 * @brief           長田パッチの係数を求める
 * @param[in]       norm    面の法線ベクトル（単位ベクトル）
 * @param[in]       p1      頂点１の座標
 * @param[in]       norm1   頂点１の法線ベクトル（単位ベクトル）
 * @param[in]       d1      頂点１原点からの距離
 * @param[in]       p2      頂点２の座標
 * @param[in]       norm2   頂点２の法線ベクトル（単位ベクトル）
 * @param[in]       d2      頂点２原点からの距離
 * @param[out]      cvec    長田パッチの係数
 */
void MeshPatch::CalcPatchCoef( double norm[3], double p1[3], double norm1[3], double d1, 
											   double p2[3], double norm2[3], double d2, double cvec[3] )
{

	bool	bRet;
	double	pos_line[3], vec_line[3];
	double	pos_m[3];
	double	pos_x[3];
	
	double alw1_norm = m_alwAngle;

	double alw2_norm = 0.0001;	//頂点ベクトル並行チェック用　　0.0057° 
	double alw_plane = 0.003;	//面と頂点ベクトル並行チェック用　　0.18° 

	double asw_norm = 1.0 - fabs( CalcInProduct( norm1, norm2 ) );

	if( asw_norm < alw1_norm ){						// 頂点ベクトルが同方向（穏やかなチェック）

		double asw_plane = 1.0 - fabs( CalcInProduct( norm, norm1 ) );
		
		if( asw_plane < alw_plane ){				// 面と頂点ベクトルが同方向
			if( asw_norm < alw2_norm ){				// 頂点ベクトルが同方向			
				cvec[0] = 0.0;
				cvec[1] = 0.0;
				cvec[2] = 0.0;
				return;
			}
		}
		else{
			cvec[0] = 0.0;
			cvec[1] = 0.0;
			cvec[2] = 0.0;
			return;
		}
	}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++
	//   P1-P2辺の処理、P1とP2の情報より、係数ベクトルを求める
	//+++++++++++++++++++++++++++++++++++++++++++++++++++

	// +++++++++ P1とP2の切平面より交線を求める
	bRet = CalcIntersectionLine( norm1, d1, norm2, d2, pos_line, vec_line );
	
	if( bRet ){									//交線あり	
		pos_m[0] = ( p1[0] + p2[0] ) / 2.0;
		pos_m[1] = ( p1[1] + p2[1] ) / 2.0;
		pos_m[2] = ( p1[2] + p2[2] ) / 2.0;
		
		// +++++++++++ 辺の中点から垂線に下した交点
		CalcNearPosOnLine( pos_m, pos_line, vec_line, pos_x );

		// ++++++++++++ 係数ベクトルを求める
		cvec[0] = ( pos_m[0] - pos_x[0] ) * 2.0;
		cvec[1] = ( pos_m[1] - pos_x[1] ) * 2.0;
		cvec[2] = ( pos_m[2] - pos_x[2] ) * 2.0;
	}
	else{										//交線なし
		cvec[0] = 0.0;
		cvec[1] = 0.0;
		cvec[2] = 0.0;
	}

	return;

}

/**
 * @brief           長田パッチ係数より、曲面上の点を求める
 * @param[in]       eta     η
 * @param[in]       zeta    ζ
 * @param[in]       p1      頂点１の座標
 * @param[in]       p2      頂点２の座標
 * @param[in]       p3      頂点３の座標
 * @param[in]       c1      長田パッチ係数C１
 * @param[in]       c2      長田パッチ係数C２
 * @param[in]       c3      長田パッチ係数C３
 * @param[out]      pos     求められた点の座標
 * @return          0
 */
int MeshPatch::Evaluate2Deg( double eta, double zeta, double p1[3], double p2[3], double p3[3], 
							 double c1[3], double c2[3], double c3[3], double pos[3] )
{

	double A1 = ( 1.0 - eta  );
	double A2 = ( eta - zeta );

	pos[0] = p1[0]*A1 + p2[0]*A2 + p3[0]*zeta - c1[0]*A1*A2 - c2[0]*A2*zeta - c3[0]*A1*zeta;
	pos[1] = p1[1]*A1 + p2[1]*A2 + p3[1]*zeta - c1[1]*A1*A2 - c2[1]*A2*zeta - c3[1]*A1*zeta;
	pos[2] = p1[2]*A1 + p2[2]*A2 + p3[2]*zeta - c1[2]*A1*A2 - c2[2]*A2*zeta - c3[2]*A1*zeta;

	return 0;
}

/**
 * @brief           ２平面の交線（無現線分）
 * @param[in]       vec1    平面１の法線ベクトル（単位ベクトル）
 * @param[in]       d1      頂点１原点からの距離
 * @param[in]       vec2    平面２の法線ベクトル（単位ベクトル）
 * @param[in]       d2      頂点２原点からの距離
 * @param[out]      pos     面の交線の通過点（原点からの最短距離）
 * @param[out]      vec     面の交線ベクトル
 * @return          true
 */
bool MeshPatch::CalcIntersectionLine( double vec1[3], double d1, 
											 double vec2[3], double d2, 
											 double pos[3],  double vec[3] )
{

	double d3;

	//++++++++++++　交線のベクトルを求める　
	CalcOutProduct( vec1, vec2, vec );
	double len = CalcVecSize( vec );
	
	// +++++++ 許容誤差が認めない場合、両面が並行
	//if( len < CALC_EPSILON_L ) return false;
	if( len < CALC_EPSILON_V ) return false;

	vec[0] /= len;
	vec[1] /= len;
	vec[2] /= len;		
	d3 = 0.0;			// 原点を通過する面を想定する

	// +++++++ ３平面の交点を求める（以下の方程式を解く）
	//
	//    A1x + B1y + C1z = D1   vec1, d1
	//    A2x + B2y + C2z = D2   vec2, d2
	//    A3x + B3y + C3z = D3   vec,  d3 = 0.0
	//
	// +++++++++++++++++++++++++++++++++++++++++++++++
	{
		double b_c, a_c, a_b, size_det;
		double d_c, d_b, a_d, detinv;
	
		b_c = vec2[1]*vec[2] - vec[1]*vec2[2];
		a_c = vec2[0]*vec[2] - vec[0]*vec2[2];	
		a_b = vec2[0]*vec[1] - vec[0]*vec2[1];

		size_det = vec1[0]*b_c - vec1[1]*a_c + vec1[2]*a_b;

		if( fabs(size_det) < CALC_EPSILON_V ){
			return false;			// 少なくとも両平面並行
		}

		d_c = d3*vec2[2] - d2*vec[2];
		d_b = d3*vec2[1] - d2*vec[1];
		a_d = d2*vec[0]  - d3*vec2[0];
		detinv = (double)1.0/size_det;

		pos[0] = ( vec1[1]*d_c +      d1*b_c - vec1[2]*d_b )*detinv;
		pos[1] = (   (-d1)*a_c - vec1[0]*d_c - vec1[2]*a_d )*detinv;
		pos[2] = ( vec1[1]*a_d + vec1[0]*d_b +      d1*a_b )*detinv;
	}

	return true;
}

/**
 * @brief           点から線分上に垂線を下した点を求める
 * @param[in]       pnt     点座標
 * @param[in]       pos     線分の通過点
 * @param[in]       vec     線分の方向ベクトル（単位ベクトル）
 * @param[out]      pos_x   点から線分上に垂線を下した点
 */
void MeshPatch::CalcNearPosOnLine( double pnt[3], double pos[3], 
								   double vec[3], double pos_x[3] )
{
	double vec_wk[3];
	
	vec_wk[0] = pnt[0] - pos[0];
	vec_wk[1] = pnt[1] - pos[1];
	vec_wk[2] = pnt[2] - pos[2];
	
	double dis_wk = CalcInProduct( vec, vec_wk );

	// +++++++ 点から線分上に垂線を下した点
	pos_x[0] = pos[0] + vec[0]*dis_wk;
	pos_x[1] = pos[1] + vec[1]*dis_wk;
	pos_x[2] = pos[2] + vec[2]*dis_wk;		

}

//++++++++++++++++++++++++++++++++++++++++++++++++++
//   　　　　　　　　三次ベジェ曲面
//++++++++++++++++++++++++++++++++++++++++++++++++++

/**
 * @brief           長田パッチ係数より、曲面上の点を求める
 * @param[in]       eta     η
 * @param[in]       zeta    ζ
 * @param[in]       p1      頂点1の座標
 * @param[in]       p2      頂点2の座標
 * @param[in]       p3      頂点3の座標
 * @param[in]       c1      長田パッチ生成頂点座標
 * @param[in]       c2      長田パッチ生成頂点座標
 * @param[in]       c3      長田パッチ生成頂点座標
 * @param[in]       c4      長田パッチ生成頂点座標
 * @param[in]       c5      長田パッチ生成頂点座標
 * @param[in]       c6      長田パッチ生成頂点座標
 * @param[in]       c7      長田パッチ生成頂点座標
 * @param[out]      pos_o   求められた点の座標
 * @return          1
 */
int MeshPatch::Evaluate3Deg( double eta, double zeta, double p1[3], double p2[3], double p3[3], 
							 double c1[3], double c2[3], double c3[3], double c4[3], double c5[3], 
							 double c6[3], double c7[3], double pos_o[3] )
{
    double u,v,w;
    
    // x(u,v,w) =    p1*w*w*w + cp1_p1p2*3*u*w*w + cp2_p1p2*3*u*u*w 
    //            +  p2*u*u*u + cp1_p2p3*3*u*u*v + cp2_p2p3*3*u*v*v
    //            +  p3*v*v*v + cp1_p3p1*3*v*v*w + cp2_p3p1*3*v*w*w
    //            +  cp_center*6*u*v*w
    //
    //     u = eta - zeta
    //     v = zeta
    //     w = 1 - eta
    //     u + v + w = 1
    //
    
    u = eta - zeta;
    v = zeta;
    w = 1.0 - eta;

    pos_o[0] =   p1[0]*w*w*w + c1[0]*3.0*u*w*w + c2[0]*3.0*u*u*w
               + p2[0]*u*u*u + c3[0]*3.0*u*u*v + c4[0]*3.0*u*v*v
               + p3[0]*v*v*v + c5[0]*3.0*v*v*w + c6[0]*3.0*v*w*w
			   + c7[0]*6.0*u*v*w;

    pos_o[1] =   p1[1]*w*w*w + c1[1]*3.0*u*w*w + c2[1]*3.0*u*u*w
               + p2[1]*u*u*u + c3[1]*3.0*u*u*v + c4[1]*3.0*u*v*v
               + p3[1]*v*v*v + c5[1]*3.0*v*v*w + c6[1]*3.0*v*w*w
			   + c7[1]*6.0*u*v*w;

    pos_o[2] =   p1[2]*w*w*w + c1[2]*3.0*u*w*w + c2[2]*3.0*u*u*w
               + p2[2]*u*u*u + c3[2]*3.0*u*u*v + c4[2]*3.0*u*v*v
               + p3[2]*v*v*v + c5[2]*3.0*v*v*w + c6[2]*3.0*v*w*w
			   + c7[2]*6.0*u*v*w;

    return 1;
}


/**
 * @brief           メッシュのベジェ制御点取得（３次多項式用）
 */
int MeshPatch::calcControlPoint_3deg()
{
    // ++++++++ p1->p2辺 制御点
    calcControlPointEdge( m_p1, m_norm1, m_d1, m_p2, m_norm2, m_d2, m_cp1_p1p2, m_cp2_p1p2 );

    //  p2->p3辺 制御点
    calcControlPointEdge( m_p2, m_norm2, m_d2, m_p3, m_norm3, m_d3, m_cp1_p2p3, m_cp2_p2p3 );
 
    //  p3->p1辺 制御点
    calcControlPointEdge( m_p3, m_norm3, m_d3, m_p1, m_norm1, m_d1, m_cp1_p3p1, m_cp2_p3p1 );

    //  中央制御点
    calcControlPointCenter( m_p1, m_p2, m_p3, m_cp1_p1p2, m_cp2_p1p2, m_cp1_p2p3, 
								  m_cp2_p2p3, m_cp1_p3p1, m_cp2_p3p1, m_cp_center );

    return 1;
}


/**
 * @brief           各辺の制御点取得（３次多項式用）
 * @param[in]       p1      頂点１座標
 * @param[in]       norm1   頂点１ベクトル
 * @param[in]       d1      頂点１ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
 * @param[in]       p2      頂点２座標
 * @param[in]       norm2   頂点２ベクトル
 * @param[in]       d2      頂点２原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
 * @param[out]      cp1_e   辺から求める制御点１
 * @param[out]      cp2_e   辺から求める制御点２
 */
void MeshPatch::calcControlPointEdge( double p1[3], double norm1[3], double d1, 
									  double p2[3], double norm2[3], double d2,  
									  double cp1_e[3], double cp2_e[3] )
{
    double vec_wk[3];
    double norm_wk1[3], norm_wk2[3];
    double norm_base[3]; // 曲線平面II（基準面）の法線ベクトル
    double p11[3];       // 制御点p11
    double p11_0[3];     // 制御点座標（補正後 1番目の制御点用）
    double p11_1[3];     // 制御点座標（補正後 2番目の制御点用）

    //-------------------------------------------
    // 制御点を置く、基準面を求める
    //     2等分した平面 曲線平面II（基準面）とする
    //     基準面： p1, norm_base で定義される
    //-------------------------------------------
        
    //  p1->p2辺
    vec_wk[0] = p2[0] - p1[0];
    vec_wk[1] = p2[1] - p1[1];
    vec_wk[2] = p2[2] - p1[2];
    
    // ++++++++ 法線ベクトル
    CalcOutProduct( vec_wk, norm1, norm_wk1 ); // at p1
    CalcOutProduct( vec_wk, norm2, norm_wk2 ); // at p2
    
	// +++++++  2等分して基準となる法線ベクトルとする
    norm_base[0] = norm_wk1[0] + norm_wk2[0];
    norm_base[1] = norm_wk1[1] + norm_wk2[1];
    norm_base[2] = norm_wk1[2] + norm_wk2[2];
    CalcNormalize( norm_base );  // 単位ベクトル化
    
    //-------------------------------------------
    // 制御点決定
    //-------------------------------------------

    // 制御点p11取得 ( 2次多項式と同様 )
    calcP11( p1, norm1, d1, p2, norm2, d2, norm_base, p11 );

    // 制御点p11逆行補正対応
    correctP11( p11, p1, p2, norm_base, p11_0, p11_1 );
    
    // 制御点(3次多項式用) 設定(式5.14,5.15)
    cp1_e[0] = ( p1[0] + 2.0*p11_0[0] ) / 3.0;
    cp1_e[1] = ( p1[1] + 2.0*p11_0[1] ) / 3.0;
    cp1_e[2] = ( p1[2] + 2.0*p11_0[2] ) / 3.0;
    
    cp2_e[0] = ( p2[0] + 2.0*p11_1[0] ) / 3.0;
    cp2_e[1] = ( p2[1] + 2.0*p11_1[1] ) / 3.0;
    cp2_e[2] = ( p2[2] + 2.0*p11_1[2] ) / 3.0;

}

/**
 * @brief           中央の制御点取得（３次多項式用）
 * @param[in]       p1          頂点１座標
 * @param[in]       p2          頂点２座標
 * @param[in]       p3          頂点３座標
 * @param[in]       cp1_p1p2    p1p2辺の制御点1
 * @param[in]       cp2_p1p2    p1p2辺の制御点2
 * @param[in]       cp1_p2p3    p2p3辺の制御点1
 * @param[in]       cp2_p2p3    p2p3辺の制御点2
 * @param[in]       cp1_p3p1    p3p1辺の制御点1
 * @param[in]       cp2_p3p1    p3p1辺の制御点2
 * @param[out]      cp_center   中央制御点
 */
void MeshPatch::calcControlPointCenter( double p1[3], double p2[3], double p3[3],      
										double cp1_p1p2[3], double cp2_p1p2[3],
										double cp1_p2p3[3], double cp2_p2p3[3],
										double cp1_p3p1[3], double cp2_p3p1[3],
										double cp_center[3]
      )
{
    cp_center[0] =   (   cp1_p1p2[0] + cp2_p1p2[0]  + cp1_p2p3[0] + cp2_p2p3[0] + cp1_p3p1[0] + cp2_p3p1[0] ) / 4.0
                   - (   p1[0] + p2[0] + p3[0]     ) / 6.0;

    cp_center[1] =   (   cp1_p1p2[1] + cp2_p1p2[1]  + cp1_p2p3[1] + cp2_p2p3[1] + cp1_p3p1[1] + cp2_p3p1[1] ) / 4.0
                   - (   p1[1] + p2[1] + p3[1]     ) / 6.0;

    cp_center[2] =   (   cp1_p1p2[2] + cp2_p1p2[2]  + cp1_p2p3[2] + cp2_p2p3[2] + cp1_p3p1[2] + cp2_p3p1[2] ) / 4.0
                   - (   p1[2] + p2[2] + p3[2]     ) / 6.0;

}



/**
 * @brief           制御点p11取得
 * @param[in]       p1          頂点１座標
 * @param[in]       norm1       頂点１ベクトル
 * @param[in]       d1          頂点１ 原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
 * @param[in]       p2          頂点２座標
 * @param[in]       norm2       頂点２ベクトル
 * @param[in]       d2          頂点２原点からの距離(+-)  平面の方程式  Ax + By + Cz = D
 * @param[in]       norm_base   曲線平面II（基準面）の法線ベクトル
 * @param[out]      p11         制御点座標
 */
void MeshPatch::calcP11( double p1[3], double norm1[3], double d1,  
						 double p2[3], double norm2[3], double d2, 
						 double norm_base[3], double p11[3] )
{
    bool bRet;
    double pos_line[3], vec_line[3];	// 面の交線の通過点, 面の交線のベクトル
    double pos_m[3];					// p1-p2,p2-p3,p3-p1辺の中点
    double pos_x[3];					// 辺の中点から垂線に下した交点
    double d1_base;						// 基準面 p1での原点からの距離(+-) 
    double d_pos_x_base;				// 基準面 pos_xでの原点からの距離(+-) 
    
    
    //------------------------------------------------
    // 平面の並行判定
    //    CalcIntersectionLine内の許容誤差はe-5-e-8レベルのため
    //    曲面補間するには厳しすぎる
    //------------------------------------------------
    double alw = 0.001;   // 0.057 degree

    double asw = CalcInProduct( norm1, norm2 );
    if( (1.0-fabs(asw)) < alw ) { // 平面が同方向
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
    bRet = CalcIntersectionLine( norm1, d1, norm2, d2, pos_line, vec_line );
    
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
    CalcNearPosOnLine( pos_m, pos_line, vec_line, pos_x );
    
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

/**
 * @brief           制御点p11逆行対応補正（３次多項式用）
 * @param[in]       p11         制御点座標
 * @param[in]       p1          頂点１座標
 * @param[in]       p2          頂点２座標
 * @param[in]       p3          頂点３座標
 * @param[in]       norm_base   曲線平面II（基準面）の法線ベクトル1
 * @param[out]      p11_0       制御点座標（補正後 1番目の制御点用）
 * @param[out]      p11_1       制御点座標（補正後 2番目の制御点用）
 */
void MeshPatch::correctP11( double p11[3], double p1[3], double p2[3],  
						    double norm_base[3], double p11_0[3], double p11_1[3] )
{
    int  mode;   // -1: 始点側逆行  0: p1-p2内補正不要  1:終点側逆行
    double vec_p1_p2[3];			// p1->p2ベクトル
    double vec_p1_p11[3];			// p1->p11ベクトル
    double vec_p2_p11[3];			// p2->p11ベクトル
    double vec_p1_p11_mirror[3];	// p1->p11ベクトルのミラー
    double vec_p2_p11_mirror[3];	// p2->p11ベクトルのミラー
    double pos_x2[3];
    double pos_wk[3];
    double len_wk;
    double asw;
    bool   bRet;

    // デフォルト
    p11_0[0] = p11[0]; p11_0[1] = p11[1]; p11_0[2] = p11[2];
    p11_1[0] = p11[0]; p11_1[1] = p11[1]; p11_1[2] = p11[2];
    
    // p1->p2ベクトル
    bRet = CalcLineVec( p1, p2, vec_p1_p2, len_wk );
    if( !bRet ) return;

    // p1->p11ベクトル
    bRet = CalcLineVec( p1, p11, vec_p1_p11, len_wk );
    if( !bRet ) return;

    // p2->p11ベクトル
    bRet = CalcLineVec( p2, p11, vec_p2_p11, len_wk );
    if( !bRet ) return;

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
            return;
        }
    }
    
    // p1->p11ベクトル のミラー
    CalcVecMirror( vec_p1_p2, vec_p1_p11, vec_p1_p11_mirror );

    // p2->p11ベクトル のミラー
    CalcVecMirror( vec_p1_p2, vec_p2_p11, vec_p2_p11_mirror );

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

            // 安全のため制御点を中点にする
            p11_1[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_1[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_1[2] = ( p1[2] + p2[2] ) / 2.0;
        }
    }
    else 
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

            // 安全のため制御点を中点にする
            p11_1[0] = ( p1[0] + p2[0] ) / 2.0;
            p11_1[1] = ( p1[1] + p2[1] ) / 2.0;
            p11_1[2] = ( p1[2] + p2[2] ) / 2.0;
        }
    }
}
