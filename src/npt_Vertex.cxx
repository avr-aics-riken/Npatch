/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Vertex.cxx
 * @brief  Npatch CVertex Class
 * @author aics
 */

#include "npt_Vertex.h"


//********************************************
// class CVertex
//********************************************

//----------------------------------------------------------
//               隣接面の処理                             //
//----------------------------------------------------------

// *********************************************************
//  関数名： AddNeighbor()
//  機能  ： 隣の面を追加
//   in   ： なし   
//  out   ： なし
// *********************************************************
/**
 * @brief           隣の面を追加
 * @param[in]       pFace   面
 */
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
/**
 * @brief           隣の面の総数を取得
 * @return          面の数
 */
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
/**
 * @brief           隣の面を削除
 */
void CVertex::RemoveAllFaceNeighbor(void)
{
	m_ArrayFaceNeighbor.RemoveAll();
}
