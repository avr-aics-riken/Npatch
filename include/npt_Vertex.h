/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

//********************************************
// class CVertex
//********************************************
#ifndef _NPT_VERTEX_
#define _NPT_VERTEX_

#include "npt_Vector.h"
#include "npt_Array.h"

class CFace;


class CVertex
{

public :

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	CVertex( const float x, const float y, const float z); 
	virtual ~CVertex() {}


	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------

	// ++++++++++++++++ 設定 ++++++++++++++++++++
	void Set( CVertex *pVertex );

	float x(void) { return m_Coord[0]; }
	float y(void) { return m_Coord[1]; }
	float z(void) { return m_Coord[2]; }


	// ++++++++++++++++ 法線関連 ++++++++++++++++++++
	void		SetNormal(CVector &vector) { m_Normal.Set(vector); }
	CVector * GetNormal(void)				 { return &m_Normal;     }

	// ++++++++++++++++ 隣接面関連 ++++++++++++++++++
	int		GetNumOfFaceNeighbor(void);
	void	AddNeighbor( CFace *pFace );
	void	RemoveAllFaceNeighbor(void);
	CFace *GetFaceNeighbor( const unsigned int index) 
								{ return m_ArrayFaceNeighbor[index]; }


private :

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------

	float				m_Coord[3];				// 座標
	CVector			m_Normal;				// 法線
	CArray<CFace>   m_ArrayFaceNeighbor;	// 隣接面

};



inline CVertex::CVertex( const float x, const float y, const float z) 
{
	m_Coord[0] = x; 
	m_Coord[1] = y; 
	m_Coord[2] = z; 
}

inline void CVertex::Set( CVertex *pVertex )
{
	m_Coord[0] = pVertex->x();
	m_Coord[1] = pVertex->y();
	m_Coord[2] = pVertex->z();
}

#endif // _NPT_VERTEX_

