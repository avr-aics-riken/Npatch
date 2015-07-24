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


/**
 * @brief           デストラクタ （バッファを削除） 
 */
CVector::~CVector()
{
}

/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       vector 　ベクトル　例；B[3]
 */
CVector::CVector( CVector &vector )
{
	Set( &vector );
}

/**
 * @brief           コンストラクタ　（初期化） 
 * @param[in]       pVector ベクトル
 */
CVector::CVector( CVector *pVector )
{
	Set( pVector );
}

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

/**
 * @brief           ベクトルの初期化
 */
void CVector::Clear(void)
{
	Set( 0.0f, 0.0f, 0.0f);
}

/**
 * @brief           ベクトル格納
 * @param[in]       x   X座標
 * @param[in]       y   Y座標
 * @param[in]       z   Z座標
 */
void CVector::Set( const float x, const float y, const float z )
{
	m_x = x;
	m_y = y;
	m_z = z;
}

/**
 * @brief           ベクトル格納
 * @param[in]       pVector ベクトル
 */
void CVector::Set( CVector *pVector )
{
	Set( pVector->x(), pVector->y(), pVector->z());
}


/**
 * @brief           ベクトル格納
 * @param[in]       pVector ベクトル
 * @param[in]       index   インデックス
 */
void CVector::Set( CVector *pVector, int index )
{
	Set( pVector->x(), pVector->y(), pVector->z());
}

/**
 * @brief           ベクトル格納
 * @param[in]       vector  ベクトル
 */
void CVector::Set(CVector &vector)
{
	Set( vector.x(), vector.y(), vector.z());
}


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



