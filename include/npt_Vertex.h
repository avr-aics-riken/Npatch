/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   npt_Vertex.h
 * @brief  Npatch CVertex Class Header
 * @author aics    
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
/**
 * @brief           デストラクタ　（バッファを削除） 
 */
	virtual ~CVertex() {}


	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------

	// ++++++++++++++++ 設定 ++++++++++++++++++++
	void Set( CVertex *pVertex );

/**
 * @brief           X座標の取得
 * @return          X座標
 */
	float x(void) { return m_Coord[0]; }
/**
 * @brief           Y座標の取得
 * @return          Y座標
 */
	float y(void) { return m_Coord[1]; }
/**
 * @brief           Z座標の取得
 * @return          Z座標
 */
	float z(void) { return m_Coord[2]; }


	// ++++++++++++++++ 法線関連 ++++++++++++++++++++
/**
 * @brief           法線ベクトルの登録
 * @param[in]       vector	法線ベクトル
 */
	void		SetNormal(CVector &vector) { m_Normal.Set(vector); }
/**
 * @brief           法線ベクトルの取得
 * @return          法線ベクトル
 */
	CVector * GetNormal(void)				 { return &m_Normal;     }

	// ++++++++++++++++ 隣接面関連 ++++++++++++++++++
	int		GetNumOfFaceNeighbor(void);
	void	AddNeighbor( CFace *pFace );
	void	RemoveAllFaceNeighbor(void);
/**
 * @brief           隣接面の取得
 * @param[in]       index   インデックス
 * @return          隣接面
 */
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



/**
 * @brief           頂点座標の登録
 * @param[in]       x	x座標
 * @param[in]       y	y座標
 * @param[in]       z	z座標
 */
inline CVertex::CVertex( const float x, const float y, const float z) 
{
	m_Coord[0] = x; 
	m_Coord[1] = y; 
	m_Coord[2] = z; 
}

/**
 * @brief           頂点座標の登録
 * @param[in]       pVertex	頂点座標
 */
inline void CVertex::Set( CVertex *pVertex )
{
	m_Coord[0] = pVertex->x();
	m_Coord[1] = pVertex->y();
	m_Coord[2] = pVertex->z();
}

#endif // _NPT_VERTEX_

