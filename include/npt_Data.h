/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Data.h
 * @brief  Npatch CData Class Header
 * @author aics
 */


//********************************************
// class CData
//********************************************
#ifndef _NPT_DATA_
#define _NPT_DATA_

#include "npt_Vertex.h"
#include "npt_Vector.h"
#include "npt_Face.h"
#include "npt_IOforStl.h"

class CData  
{

public:

// *******************************************
//   データ構造
// *******************************************
typedef struct
{
	double	x;
	double	y;
	double	z;
} VertexD;

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	CData();
	virtual ~CData();


	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	void	MainProcess(  char* infile, int degeree);
	void	SaveNPT( int OutType );
	
	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------

	// ***************************************
	//     　　　入力用
	// ***************************************
	char	m_InFileName[256];
	char	m_OutFileName[256];
	int		m_FileType;
	int		m_ExType;

	// ***************************************
	//     　　　アイテム総数
	// ***************************************
	int		m_FaceNum;
	int		m_VertexNum;

	// ***************************************
	//     　　　ポインタ
	// ***************************************
	IOForStl	readStl;	//STL読込クラスの継承

	CArray< CVertex >	m_ArrayVertex;
	CArray< CFace   >   m_ArrayFace;
	VertexD					*bAsw[3];	//長田バッチの計算点

protected:

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	int		StoreStl();
	int		Store_AsciiStlAsMesh( CArray<CFace>*   pArrayFace, CArray<CVertex>* pArrayVertex );
	int		Store_BinaryStlAsMesh( CArray<CFace>* pArrayFace, CArray<CVertex>* pArrayVertex );

	int		Store_AsciiNptAsMesh( CArray<CFace>*   pArrayFace, CArray<CVertex>* pArrayVertex );
	int		Store_BinaryNptAsMesh( CArray<CFace>* pArrayFace, CArray<CVertex>* pArrayVertex );

	// ++++++++++++++ 共有頂点の平均化　+++++++++++++++++++++++

public :

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	void	Vertex_FreeArray();  //描画用、将来protectedに回帰


protected:

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	int		Vertex_BuildAdjacency();
	int		Vertex_CalcNormalPerVertex(void);
	int		Vertex_CalcNormalPerFace(void);

//public :
	// ++++++++++++++ 長田ﾊﾞｯﾁ　+++++++++++++++++++++++
	int		MeshPatch_GetCoef2Deg();
	int		MeshPatch_GetValue2Deg();		// 一辺中の一点（η、ζ）の評価値を求める
	int		MeshPatch_GetValue3Point2Deg();// 三辺中の中点（η、ζ）の評価値を求める
	
	int		MeshPatch_GetCoef3Deg();
	int		MeshPatch_GetValue3Deg();		// 一辺中の一点（η、ζ）の評価値を求める
	int		MeshPatch_GetValue3Point3Deg();	// 三辺中の中点（η、ζ）の評価値を求める

	void	MeshPatch_freeCoef();
	
	double	m_CheckAngle;
	int		m_degree;		// 多項式の次数
	double	m_zeta, m_eta;

	// ++++++++++ inline ++++++++++++++++++++++++++++++
	void	SetValue( double p[3], double x, double y, double z );
	void	SetValueSub( double vec[3], double x1, double y1, double z1, 
									    double x2, double y2, double z2 );

};

inline void CData::SetValue( double p[3], double x, double y, double z )
{
	p[0] = x; 
	p[1] = y;
	p[2] = z;
}

inline void CData::SetValueSub( double vec[3], double x1, double y1, double z1, 
											   double x2, double y2, double z2 )
{
		vec[0] = x2 - x1;
		vec[1] = y2 - y1;
		vec[2] = z2 - z1;
}

#endif // _NPT_DATA_
