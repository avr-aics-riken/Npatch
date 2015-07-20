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
#include "npt_Vertex.h"

//----------------------------------------------------------
//               隣接面の処理                             //
//----------------------------------------------------------

// *********************************************************
//  関数名： AddNeighbor()
//  機能  ： 隣の面を追加
//   in   ： なし   
//  out   ： なし
// *********************************************************
void CVertex::AddNeighbor(CFace *pFace)
{
	m_ArrayFaceNeighbor.Add( pFace );
}

// *********************************************************
//  関数名： GetNumOfFaceNeighbor()
//  機能  ： 隣の面の総数を取得
//   in   ： なし    
//  out   ： 面の数       
// *********************************************************
int CVertex::GetNumOfFaceNeighbor(void)
{
	return m_ArrayFaceNeighbor.GetSize();
}

// *********************************************************
//  関数名： RemoveAllFaceNeighbor()
//  機能  ： 隣の面を削除
//   in   ： なし
//  out   ： なし       
// *********************************************************
void CVertex::RemoveAllFaceNeighbor(void)
{
	m_ArrayFaceNeighbor.RemoveAll();
}
