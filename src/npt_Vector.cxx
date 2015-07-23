/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   npt_Vector.cxx
 * @brief  Npatch CVector Class
 * @author aics    
 */

#include <math.h>

#include "npt_Vertex.h"
#include "npt_Vector.h"


//********************************************
// class CVector
//********************************************


// ----------------------------------------------
//                  コンストラクタ
// ----------------------------------------------


// ********************************************
//  関数名： CFace()
//  機能  ： コンストラクタ
// ********************************************
/**
 * @brief           デストラクタ （バッファを削除） 
 */
CVector::~CVector()
{
}

// **************************************************************
//  関数名： CVector
//  機能  ： コンストラクタ　（初期化） 
//   in   ： vector 　ベクトル　例；B[3]
//  out   ： なし
// **************************************************************
/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       vector 　ベクトル　例；B[3]
 */
CVector::CVector( CVector &vector )
{
	Set( &vector );
}

// **************************************************************
//  関数名： CVector
//  機能  ： コンストラクタ　（初期化） 
//   in   ： *pVector 　ベクトルのポインタ　
//  out   ： なし
// **************************************************************
/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       pVector ベクトル
 */
CVector::CVector( CVector *pVector )
{
	Set( pVector );
}

// **************************************************************
//  関数名： CVector
//  機能  ： コンストラクタ　（初期化） 
//   in   ： x　　X座標
//   in   ： ｙ　 X座標　
//   in   ： ｚ   Z座標　　
//  out   ： なし
// **************************************************************
/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       x   X座標
 * @param[in]       y   Y座標
 * @param[in]       z   Z座標
 */
CVector::CVector( const float x, const float y, const float z)
{
	Set( x, y, z);
}

// **************************************************************
//  関数名： CVector
//  機能  ： コンストラクタ　（初期化） 
//   in   ： x　　X座標
//   in   ： ｙ　 Y座標　
//   in   ： ｚ   Z座標　　
//  out   ： なし
// **************************************************************
/**
 * @brief           コンストラクタ　（初期化）
 * @param[in]       pVertex1    始点
 * @param[in]       pVertex2    終点
 */
CVector::CVector( CVertex *pVertex1, CVertex *pVertex2 )
{
	m_x = pVertex2->x() - pVertex1->x();
	m_y = pVertex2->y() - pVertex1->y();
	m_z = pVertex2->z() - pVertex1->z();
}
	
// ********************************************
//  関数名： CVector()
//  機能  ： コンストラクタ（ベクトルの長さ）
//   in   ： ベクトル
// ********************************************
/**
 * @brief           コンストラクタ（ベクトルの長さ）
 * @param[in]       vertex1     始点
 * @param[in]       vertex2     終点
 */
CVector::CVector( CVertex &vertex1, CVertex &vertex2 )
{
	m_x = vertex2.x() - vertex1.x();
	m_y = vertex2.y() - vertex1.y();
	m_z = vertex2.z() - vertex1.z();
}


// ----------------------------------------------
//                  その他
// ----------------------------------------------

// ********************************************
//  関数名： Clear()
//  機能  ： ベクトルの初期化
// ********************************************
/**
 * @brief           ベクトルの初期化
 */
void CVector::Clear(void)
{
	Set( 0.0f, 0.0f, 0.0f);
}

// ********************************************
//  関数名： Set()
//  機能  ： ベクトル格納
//   in   ： 座標など
// ********************************************
void CVector::Set( const float x, const float y, const float z )
/**
 * @brief           ベクトル格納
 * @param[in]       x   X座標
 * @param[in]       y   Y座標
 * @param[in]       z   Z座標
 */
{
	m_x = x;
	m_y = y;
	m_z = z;
}

// ********************************************
//  関数名： Set()
//  機能  ： ベクトル格納
//   in   ： ベクトル
// ********************************************
/**
 * @brief           ベクトル格納
 * @param[in]       pVector ベクトル
 */
void CVector::Set( CVector *pVector )
{
	Set( pVector->x(), pVector->y(), pVector->z());
}


// ********************************************
//  関数名： Set()
//  機能  ： ベクトル格納
//   in   ： ベクトル
// ********************************************
/**
 * @brief           ベクトル格納
 * @param[in]       pVector ベクトル
 */
void CVector::Set( CVector *pVector, int index )
{
	Set( pVector->x(), pVector->y(), pVector->z());
}

// ********************************************
//  関数名： Set()
//  機能  ： ベクトル格納
//   in   ： ベクトル
// ********************************************
/**
 * @brief           ベクトル格納
 * @param[in]       vector  ベクトル
 */
void CVector::Set(CVector &vector)
{
	Set( vector.x(), vector.y(), vector.z());
}


// ********************************************
//  関数名： Set()
//  機能  ： ベクトル格納
//   in   ： ベクトル
// ********************************************
/**
 * @brief           ベクトル格納
 * @param[in]       pVertex1    始点
 * @param[in]       pVertex2    終点
 */
void CVector::Set( CVertex *pVertex1, CVertex *pVertex2 )
{
	m_x = pVertex2->x() - pVertex1->x();
	m_y = pVertex2->y() - pVertex1->y();
	m_z = pVertex2->z() - pVertex1->z();
}

//********************************************
// Operator + 
//********************************************
/**
 * @brief           演算子 +
 * @param[in]       pVector ベクトル
 */
void CVector::operator+=( CVector* pVector )
{
	m_x += pVector->x();
	m_y += pVector->y();
	m_z += pVector->z();
}

// ********************************************
//  関数名： Cross()
//  機能  ： 外積（本ベクトルと他のベクトル）
//   in   ： ベクトル
// ********************************************
/**
 * @brief           外積（本ベクトルと他のベクトル）
 * @param[in]       vector ベクトル
 */
void CVector::Cross( CVector& vector )
{
	double x = (double)m_y*(double)vector.z() - (double)m_z*(double)vector.y();
	double y = (double)m_z*(double)vector.x() - (double)m_x*(double)vector.z();
	double z = (double)m_x*(double)vector.y() - (double)m_y*(double)vector.x();
	Set( (float)x, (float)y, (float)z );
}

// ********************************************
//  関数名： NormalizeL2()
//  機能  ： ベクトルの正規化
// ********************************************
/**
 * @brief           ベクトルの正規化
 */
void CVector::NormalizeL2(void)
{
	double norm = GetNormL2();
	if( norm != 0.0f )
	{
		m_x = (float)(m_x / norm);
		m_y = (float)(m_y / norm);
		m_z = (float)(m_z / norm);
	}
}
	
// ********************************************
//  関数名： GetNormL2()
//  機能  ： ベクトルの長さ
// ********************************************
/**
 * @brief           ベクトルの長さ
 * @return          ベクトルの長さ
 */
double CVector::GetNormL2(void)
{
	return sqrt((double)m_x*(double)m_x + 
	            (double)m_y*(double)m_y + 
				(double)m_z*(double)m_z  );
}



// ** EOF **



