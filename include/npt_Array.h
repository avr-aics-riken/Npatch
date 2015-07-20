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
CArray(){
	m_pData = 0;
	m_nSize = 0;
	m_nMaxSize = 0;
	m_nGrowBy = 0;
}

// Destruction
~CArray()	{
	delete [] (BYTE*)m_pData;
	}

float m_scale;

// Attributes
int GetSize()		{ return m_nSize; }
int GetUpperBound() { return m_nSize-1; }

// Operations

// Really clean up
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
void RemoveAll() { SetSize(0); }

// ++++++++ Accessing elements
T* GetAt(int nIndex){ 
	assert(nIndex >= 0 && nIndex < m_nSize);
	return m_pData[nIndex]; 
}

// +++++++++ SetAt
void SetAt(int nIndex, T* newElement){ 
	assert(nIndex >= 0 && nIndex < m_nSize);
	m_pData[nIndex] = newElement; 
}

// ++++++++++ Add
int Add(T* newElement){ 
	int nIndex = m_nSize;
	SetAtGrow(nIndex, newElement);
	return nIndex; 
}

// ++++++++++ overloaded operator helpers
T* operator[](int nIndex) { return GetAt(nIndex); }

// ++++++++++ Potentially growing the array
void SetAtGrow(int nIndex, T* newElement){

	assert(nIndex >= 0);
	if (nIndex >= m_nSize)	SetSize(nIndex+1);
	m_pData[nIndex] = newElement;
}

// +++++++++ SetSize
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


inline const int Amax( const int a, const int b )
{
	return a < b ? b : a;
}

inline const int Amin( const int a, const int b )
{
	return b < a ? b : a;
}



#endif // _NPT_ARRAY_
