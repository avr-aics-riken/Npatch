/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Face.cxx
 * @brief  Npatch CFace Class
 * @author aics
 */

#include <stdio.h>
#include "npt_Face.h"


//********************************************
// class CFace
//********************************************


// ----------------------------------------------
//                  コンストラクタ
// ----------------------------------------------


/**
 * @brief           コンストラクタ　（初期化） 
 */
CFace::CFace()
{
	for(int i=0; i<3; i++) m_pFace[i]   = NULL;
	for(int i=0; i<3; i++) m_pVertex[i] = NULL;
}

/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       pVertex1    第一頂点の座標
 * @param[in]       pVertex2    第二頂点の座標
 * @param[in]       pVertex3    第三頂点の座標
 */
CFace::CFace( CVertex *pVertex1, 
				  CVertex *pVertex2, 
				  CVertex *pVertex3  )
{
	// +++++++++ Face ++++++++++++
	for(int i=0; i<3; i++ )  m_pFace[i] = NULL;
	
	// +++++++++ Vertices ++++++++
	Set( pVertex1, pVertex2, pVertex3);

	// +++++++++ Normal +++++++++
	m_Normal.Set(0.0f,0.0f,0.0f);
}

// ---------------------------------------------
//             クラス操作                     //
// ---------------------------------------------

/**
 * @brief           ARRAY情報のクリア
 */
inline void CFace::Clear()
{
	int i;
	for( i=0; i<3; i++ ) m_pFace[i]		= NULL;
	for( i=0; i<3; i++ ) m_pVertex[i]	= NULL;
}

/**
 * @brief           頂点の代入
 * @param[in]       pVertex1    第一頂点の座標
 * @param[in]       pVertex2    第二頂点の座標
 * @param[in]       pVertex3    第三頂点の座標
 */
inline void CFace::Set( CVertex *pVertex1, 
						  CVertex *pVertex2, 
						  CVertex *pVertex3  )
{
	m_pVertex[0] = pVertex1;
	m_pVertex[1] = pVertex2;
	m_pVertex[2] = pVertex3;
}

// ---------------------------------------------
//             法線の計算                     //
// ---------------------------------------------

/**
 * @brief           
 * @param[in]       面の法線計算
 */
void CFace::CalculateNormal()
{
	CVector u( m_pVertex[0], m_pVertex[1] );
	CVector v( m_pVertex[0], m_pVertex[2] );

	u.Cross(v);
	m_Normal.Set(u);
	m_Normal.NormalizeL2();
}

/**
 * @brief           面上の節点の法線計算
 * @param[in]       vector  法線ベクトル
 * @param[in]       index   面上の頂点番号（0、または1、または2）
 * @param[in]       Num     出来た法線の総数
 */
void CFace::SetNormal_New( CVector &vector, 
							 int index, int Num )
{
	m_NewNormal[index].x = vector.x()/Num;
	m_NewNormal[index].y = vector.y()/Num;
	m_NewNormal[index].z = vector.z()/Num;
}

/**
 * @brief           長田係数を取得
 * @param[in]       c_aws   長田パッチ係数
 * @param[in]       index   係数の番号（C1、C2、C3）
 */
void CFace::SetCoef( double c_aws[3], int index )
{
	m_Coef[index].x = c_aws[0];
	m_Coef[index].y = c_aws[1];
	m_Coef[index].z = c_aws[2];
}



