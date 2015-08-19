/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Array.h
 * @brief  Npatch CData Class Header
 * @author aics
 */

//********************************************
// class CArray
//********************************************
#ifndef _NPT_ARRAY_
#define _NPT_ARRAY_

//typedef unsigned char BYTE;

#include <string.h>
#include <stdio.h>

#include <iostream>
#include <cassert>

typedef unsigned char BYTE;

template<class T>
class CArray
{

// Implementation
private:
	T** m_pData;     // the actual array of data
	int m_nSize;     // # of elements (upperBound - 1)
	int m_nMaxSize;  // max allocated
	int m_nGrowBy;   // grow amount

public:

	// Construction
/**
 * @brief           コンストラクタ
 */
CArray(){
	m_pData = 0;
	m_nSize = 0;
	m_nMaxSize = 0;
	m_nGrowBy = 0;
}

// Destruction
/**
 * @brief           デストラクタ
 */
~CArray()	{
	delete [] (BYTE*)m_pData;
	}

float m_scale;

// Attributes
/**
 * @brief           配列サイズの取得
 * @return          配列サイズ
 */
int GetSize()		{ return m_nSize; }
/**
 * @brief           配列上限取得
 * @return          配列上限
 */
int GetUpperBound() { return m_nSize-1; }

// Operations

// Really clean up
/**
 * @brief           領域の解放
 */
void Free(){
	for(int i=0;i<m_nSize;i++){
		if(m_pData[i] != NULL){
			delete m_pData[i];
		  m_pData[i] = NULL;
		}
	}
	RemoveAll();
}

// ++++++++ Clean up pointers array 
/**
 * @brief           全削除
 */
void RemoveAll() { SetSize(0); }

// ++++++++ Accessing elements
/**
 * @brief           値の取得
 * @param[in]       nIndex	取得するインデックス
 * @return          値
 */
T* GetAt(int nIndex){ 
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex]; 
}

// +++++++++ SetAt
/**
 * @brief           値の登録
 * @param[in]       nIndex		登録するインデックス
 * @param[in]       newElement	登録する値
 */
void SetAt(int nIndex, T* newElement){ 
	assert(nIndex >= 0 && nIndex < m_nSize);
	m_pData[nIndex] = newElement; 
}

// ++++++++++ Add
/**
 * @brief           値を末尾に追加
 * @param[in]       newElement	追加する値
 * @return          末尾のインデックス
 */
int Add(T* newElement){ 
	int nIndex = m_nSize;
	SetAtGrow(nIndex, newElement);
	return nIndex; 
}

// ++++++++++ overloaded operator helpers
/**
 * @brief           配列アクセス
 * @param[in]       nIndex	取得するインデックス
 * @return          値
 */
T* operator[](int nIndex) { return GetAt(nIndex); }

// ++++++++++ Potentially growing the array
/**
 * @brief           値の登録（配列上限を指定されたインデックス以上に保つ）
 * @param[in]       nIndex		登録するインデックス
 * @param[in]       newElement	登録する値
 */
void SetAtGrow(int nIndex, T* newElement){

	assert(nIndex >= 0);
	if (nIndex >= m_nSize)	SetSize(nIndex+1);
	m_pData[nIndex] = newElement;
}

// +++++++++ SetSize
/**
 * @brief           配列サイズの設定
 * @param[in]       nNewSize	配列サイズ
 * @param[in]       nGrowBy		
 */
void SetSize(int nNewSize, int nGrowBy = -1){

	assert(nNewSize >= 0);

	if(nGrowBy != -1) m_nGrowBy = nGrowBy;
		
	if(nNewSize == 0){
		delete[] (BYTE*)m_pData;
		m_pData = NULL;
		m_nSize = 0;
		m_nMaxSize = 0;
	}
	else 
		if(m_pData == NULL){
			m_pData = (T**) new BYTE[nNewSize * sizeof(T*)];
			memset(m_pData, 0, nNewSize * sizeof(T*));  // zero fill
			m_nSize = nNewSize;
			m_nMaxSize = nNewSize;
		}
		else 
			if(nNewSize <= m_nMaxSize){
				if (nNewSize > m_nSize)	memset(&m_pData[m_nSize], 0, (nNewSize-m_nSize) * sizeof(T*));
				m_nSize = nNewSize;
			}
			else{
				int nGrowBy = m_nGrowBy;
				if (nGrowBy == 0) nGrowBy = std::min(1024, std::max(4, m_nSize / 8));
				//if (nGrowBy == 0) nGrowBy = Amin(1024, Amax(4, m_nSize / 8));

				int nNewMax;
				if (nNewSize < m_nMaxSize + nGrowBy) nNewMax = m_nMaxSize + nGrowBy;
				else								 nNewMax = nNewSize;
					
				assert(nNewMax >= m_nMaxSize);  // no wrap around

				T** pNewData = (T**) new BYTE[nNewMax * sizeof(T*)];
					
				// copy new data from old
				memcpy(pNewData, m_pData, m_nSize * sizeof(T*));
				
				// construct remaining elements
				assert(nNewSize > m_nSize);
				memset(&pNewData[m_nSize], 0, (nNewSize-m_nSize) * sizeof(T*));
					
				// get rid of old stuff (note: no destructors called)
				delete[] (BYTE*)m_pData;
				m_pData = pNewData;
				m_nSize = nNewSize;
				m_nMaxSize = nNewMax;
		}
}


};


/**
 * @brief           最大値計算
 * @param[in]       a	比較値1
 * @param[in]       b	比較値2
 * @return          最大値
 */
inline const int Amax( const int a, const int b )
{
	return a < b ? b : a;
}

/**
 * @brief           最小値計算
 * @param[in]       a	比較値1
 * @param[in]       b	比較値2
 * @return          最小値
 */
inline const int Amin( const int a, const int b )
{
	return b < a ? b : a;
}



#endif // _NPT_ARRAY_
