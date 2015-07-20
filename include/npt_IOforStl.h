/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_IOforStl.h
 * @brief  Npatch IOForStl Class Header
 * @author aics
 */

//********************************************
// class IOForStl
//********************************************
#ifndef _NPT_IO_FOR_STL_
#define _NPT_IO_FOR_STL_

#include "npt_Face.h"
#include "npt_Vertex.h"

// ****************************************
//     　　　定数定義
// ****************************************
#define MAX_LINE_STL	100000
#define MAX_WORD_STL	1000
#define MAX_HEAD_SIZE	80


class IOForStl
{
public:

	// ----------------------------------------------
	//                  コンストラクタ
	// ----------------------------------------------
	IOForStl(void);
	~IOForStl(void);

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	void	LoadStl( char *InFile,int ExType, int deg );

	int		ReadWord();
	int		OffsetToString( char *string );

	int		SaveNPT_Ascii(char *name, CArray<CFace>& m_Face, int FaceNum, int degree );
	int		SaveNPT_Binary(char *name, CArray<CFace>& m_Face, int FaceNum, int degree );

	int		LoadNPT_Binary( char *filename,  int mod  );
	void	Delet_Buff();

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------

	// ***************************************
	//     　　　入力用
	// ***************************************
	int		m_FileType;

	char*	m_pBuffer;
	int		m_IndexBuffer;
	char    m_pBufferWord[MAX_WORD_STL];

	// ***************************************
	//     　　　アイテム総数
	// ***************************************
	int		FaceNum;
	int		VertexNum;
	
	int		Stl_Type;      // ( 0=ASCII, 1=binary )
	int		Npt_Type;      // ( 0=ASCII, 1=binary )
	int		Degree;
	
protected:

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	int		LoadStlFile( char *InFile );
	int		GetStlNum_Ascii();
	int		GetStlNum_Binary (char* InFile );
	int		GetNptNum_Ascii( char* InFile );
	int		GetNptNum_Binary( char* InFile );
	
	int		ReadLine();

	// ----------------------------------------------
	//                  メンバ変数 
	// ----------------------------------------------
	int		m_SizeFile;
	char    m_pBufferLine[MAX_LINE_STL];

public:

	// ----------------------------------------------
	//                  メンバ関数 
	// ----------------------------------------------
	int		CheckStlType( char* InFile );

	float	ReadFloat( FILE *fp );
	int		ReadInt( FILE* fp, int* n, int swap );

};

#endif // _NPT_IO_FOR_STL_