/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

//********************************************
// class CVector
//********************************************
#ifndef _NPT_VECTOR_
#define _NPT_VECTOR_

#define DBL_EPS 1e-6

class CVertex;

class CVector
{



public :

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	CVector() { m_x = m_y = m_z = 0.0f; } 
	virtual ~CVector();
	
	CVector( const float x, const float y, const float z );
	CVector( CVector &vector  );
	CVector( CVector *pVector );
	CVector( CVertex *pVertex1, CVertex *pVertex2 );
	CVector( CVertex &vertex1,  CVertex &vertex2  );


	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------

	// ++++++++++++++++ 設定 ++++++++++++++++++++
	void Clear(void);
	void Set( CVector *pVector );
	void Set( CVector *pVector, int index );

	void Set( CVector &vector );
	void Set( CVertex *pVertex1,CVertex *pVertex2 );
	void Set( const float x, const float y, const float z );


	void x(float x) {	m_x = x; }
	void y(float y) {	m_y = y; }
	void z(float z) {	m_z = z; } 

	float x(void) { return m_x; }
	float y(void) { return m_y; }
	float z(void) { return m_z; }

	// ++++++++++++++++ その他 ++++++++++++++++++++
	void operator+=( CVector* pVector );

	void	Cross( CVector& vector );
	void	NormalizeL2(void);
	double	GetNormL2(void);

private :

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------
	float m_x;
	float m_y;
	float m_z;

};


#endif // _NPT_VECTOR_
