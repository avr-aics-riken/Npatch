/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

//********************************************
// class CFace
//********************************************
#include <stdio.h>
#include "npt_Face.h"


// ----------------------------------------------
//                  コンストラクタ
// ----------------------------------------------


// **************************************************************
//  関数名： CFace
//  機能  ： コンストラクタ　（初期化） 
//   in   ： なし
//  out   ： なし
// **************************************************************
CFace::CFace()
{
	for(int i=0; i<3; i++) m_pFace[i]   = NULL;
	for(int i=0; i<3; i++) m_pVertex[i] = NULL;
}

// **************************************************************
//  関数名： CFace
//  機能  ： コンストラクタ　（初期化） 
//   in   ： なし
//  out   ： なし
// **************************************************************
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

// **************************************************************
//  関数名： Clear
//  機能  ： ARRAY情報のクリア
//   in   ： なし
//  out   ： なし
// **************************************************************
inline void CFace::Clear()
{
	int i;
	for( i=0; i<3; i++ ) m_pFace[i]		= NULL;
	for( i=0; i<3; i++ ) m_pVertex[i]	= NULL;
}

// **************************************************************
//  関数名： Set
//  機能  ： 頂点の代入
//   in   ： pVertex1  第一頂点の座標
//   in   ： pVertex2  第二頂点の座標
//   in   ： pVertex3  第三頂点の座標
//  out   ： なし
// **************************************************************
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

// **************************************************************
//  関数名： CalculateNormal
//  機能  ： 面の法線計算
//   in   ： なし
//  out   ： なし
// **************************************************************
void CFace::CalculateNormal()
{
	CVector u( m_pVertex[0], m_pVertex[1] );
	CVector v( m_pVertex[0], m_pVertex[2] );

	u.Cross(v);
	m_Normal.Set(u);
	m_Normal.NormalizeL2();
}

// **************************************************************
//  関数名： SetNormal_New
//  機能  ： 面上の節点の法線計算
//   in   ： &vector  		法線ベクトル
//   in   ： index  		面上の頂点番号（0、または1、または2）
//   in   ： Num  			出来た法線の総数
//  out   ： なし
// **************************************************************
void CFace::SetNormal_New( CVector &vector, 
							 int index, int Num )
{
	m_NewNormal[index].x = vector.x()/Num;
	m_NewNormal[index].y = vector.y()/Num;
	m_NewNormal[index].z = vector.z()/Num;
}

// **************************************************************
//  関数名： SetCoef
//  機能  ： 長田係数を取得
//   in   ： c_aws[3]		長田パッチ係数
//   in   ： index  		int index		係数の番号（C1、C2、C3）
//  out   ： なし
// **************************************************************
void CFace::SetCoef( double c_aws[3], int index )
{
	m_Coef[index].x = c_aws[0];
	m_Coef[index].y = c_aws[1];
	m_Coef[index].z = c_aws[2];
}



