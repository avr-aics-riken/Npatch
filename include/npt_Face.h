/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Face.h
 * @brief  Npatch CFace Class Header
 * @author aics
 */

//********************************************
// class CFace
//********************************************
#ifndef _NPT_FACE_
#define _NPT_FACE_

#include "npt_Vector.h"


typedef struct{
	float x;
	float y;
	float z;	
} st_ReNormal;

typedef struct{
	double x;
	double y;
	double z;	
} nPatch;

class CFace
{

public :

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	CFace();
	CFace( CVertex *pVertex1, CVertex *pVertex2, CVertex *pVertex3 );
	virtual ~CFace() {}

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------

	// ++++++++++++++++ 面の操作 ++++++++++++++++++++
	void Clear();
	void Set( CVertex *pVertex1, CVertex *pVertex2, CVertex *pVertex3 );

	CVertex *v(int index) { return m_pVertex[index%6]; }

	// ++++++++++++++++ 法線計算 +++++++++++++++++++
	CVector * GetNormal() { return &m_Normal; }
	void		SetNormal( CVector &vector ) { m_Normal.Set(vector); }
	void		CalculateNormal();
	void		SetNormal_New( CVector &vector, int index, int Num );


	void		SetCoef( double c_aws[3], int index );


	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------
	st_ReNormal m_NewNormal[3];	// 頂点の法線
	nPatch		m_Coef[7];      // 長田パッチの係数

private :

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------
	CVertex	*m_pVertex[3];      // 頂点
	CFace		*m_pFace[3];        // 隣接面
	CVector	m_Normal;           // 法線



};

#endif // _NPT_FACE_

