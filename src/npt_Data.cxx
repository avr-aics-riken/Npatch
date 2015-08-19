/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_Data.cxx
 * @brief  Npatch CData Class
 * @author aics
 */

#include <string.h>
#include <stdio.h>

#include "npt_Data.h"
#include "npt_MeshPatch.h"

/**
 * @brief コンストラクタ （初期化）
 */
CData::CData()
{

	m_FaceNum	= 0;
	m_VertexNum	= 0;

	m_CheckAngle = 0.18;
}

/**
 * @brief デストラクタ （バッファを削除）
 */
CData::~CData()
{
	//	+++++++++ メモリ開放 ++++++++++++++++++++
	Vertex_FreeArray();

}

/**
 * @brief           メイン関数
 * @param[in]       infile    入力ファイル名
 * @param[in]       degeree   曲面の次数
 */
void CData::MainProcess( char* infile, int degeree )
{

	strcpy( m_InFileName, infile );
	m_degree = degeree;

	char Ex[256]; 
	int expen = strlen(m_InFileName)-3;
	strncpy( Ex, &m_InFileName[expen], 3 );


	m_ExType = 0;
	if( strcmp( Ex, "npt" ) == 0 )  m_ExType = 1;

	strncpy( m_OutFileName, m_InFileName, strlen( m_InFileName )-4  );
	strcat( m_OutFileName, "_lib.npt" );

	std::cout<< "InFileName = " << m_InFileName << std::endl;
	std::cout<< "OutFileName = " << m_OutFileName << std::endl;
	
	std::cout<< "Ex = " << Ex << std::endl;	
	std::cout<< "ExType = " << m_ExType << std::endl;


	// +++++++++++ LoadStl ++++++++++++++++++
	readStl.LoadStl( m_InFileName, m_ExType, m_degree );

	m_FaceNum	= readStl.FaceNum;
	m_VertexNum = readStl.VertexNum;	

	// ++++++ データ格納 +++++++++++++++++
	StoreStl();

	// ++++++++++++++++++ 長田ﾊﾞｯﾁの係数を求める ++++++++++++++++++++++++++++
	if( m_degree == 2 )			MeshPatch_GetCoef2Deg();
	else if( m_degree == 3 )	MeshPatch_GetCoef3Deg();

	// ++++++ 長田係数より、曲面上の点を求める、将来は不要となる
	//SaveNPT();

}

/**
 * @brief NPTファイルの出力
 * @param[in] OutType 出力形式　（0: Ascii   1: Binary）
 */
void CData::SaveNPT( int OutType )
{
	if( OutType == 0 )		readStl.SaveNPT_Ascii( m_OutFileName, m_ArrayFace, m_FaceNum, m_degree );
	else if( OutType == 1 )	readStl.SaveNPT_Binary( m_OutFileName, m_ArrayFace, m_FaceNum, m_degree );
}

/**
 * @brief           ASCII形式のデータを読み込むためのプロセス
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::StoreStl()
{

	m_ArrayVertex.SetSize(	m_VertexNum );
	m_ArrayFace.SetSize(	m_FaceNum	);

	if( readStl.Stl_Type == 0 ){
		if( m_ExType == 0  )	 Store_AsciiStlAsMesh( &m_ArrayFace, &m_ArrayVertex );
		else if( m_ExType == 1 ) Store_AsciiNptAsMesh( &m_ArrayFace, &m_ArrayVertex );
	}
	else if( readStl.Stl_Type == 1 ){
		if( m_ExType == 0  )	 Store_BinaryStlAsMesh( &m_ArrayFace, &m_ArrayVertex );
		else if( m_ExType == 1 ) Store_BinaryNptAsMesh( &m_ArrayFace, &m_ArrayVertex );
	}
	
	// +++++++++++++ 頂点法線ベクトルの平均化　+++++++++++++
	Vertex_BuildAdjacency();
	Vertex_CalcNormalPerFace();
	Vertex_CalcNormalPerVertex();

	return 1;
}

/**
 * @brief           ASCII形式のデータを読み込みと格納　　（STL）
 * @param[in]       pArrayFace      面のバファ
 * @param[in]       pArrayVertex    頂点のバファ
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Store_AsciiStlAsMesh( CArray<CFace> *pArrayFace, CArray<CVertex> *pArrayVertex )
{
	int   tmpNum;
	float *tmp_x, *tmp_y, *tmp_z; 

	// +++++++++ TMP－ARRAYの領域確保 +++++++++++++++
	tmp_x = new float[m_VertexNum];
	tmp_y = new float[m_VertexNum];
	tmp_z = new float[m_VertexNum];

	int tmp = readStl.m_IndexBuffer;

	int		i,  k, kk;
	int		seNum1, seNum2, seNum3;
	bool	bFlag1, bFlag2, bFlag3;
	float	x1,y1,z1,x2,y2,z2,x3,y3,z3;

	int		v1, v2, v3;

	int		iVertex = 0;
	int		iFace	= 0;

	kk = 0;
	
	// ***********************************************
	//        面毎にデータ処理                      //
	// ***********************************************
	for( i=0; i<m_FaceNum; i++){

		// +++++ outer loop文字列を探し、その直後に処理 ++++
		if( readStl.OffsetToString("outer loop") ){
			readStl.m_IndexBuffer += 11;
			
			readStl.ReadWord();
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &x1 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &y1 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &z1 );

			readStl.ReadWord();
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &x2 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &y2 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &z2 );

			readStl.ReadWord();
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &x3 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &y3 );
			readStl.ReadWord();	sscanf( readStl.m_pBufferWord, "%f", &z3 );

			seNum1 = seNum2 = seNum3 = iVertex;
			bFlag1 = bFlag2 = bFlag3 = false;
			v1 = v2 = v3 = 0;

			if( i > 0 ){

				for( k=0; k<tmpNum; k++ ){

					if( bFlag1 && bFlag2 && bFlag3 ) break; 

					if( !bFlag1 ){
						if( tmp_x[k] == x1 && tmp_y[k] == y1 && tmp_z[k] == z1 ){
							seNum1 = k;
							bFlag1 = true;
						}
					}

					if( !bFlag2 ){
						if( tmp_x[k] == x2 && tmp_y[k] == y2 && tmp_z[k] == z2 ){
							seNum2 = k;
							bFlag2 = true;
						}
					}

					if( !bFlag3 ){
						if( tmp_x[k] == x3 && tmp_y[k] == y3 && tmp_z[k] == z3 ){
							seNum3 = k;
							bFlag3 = true;
						}
					}
				}
			}

			// +++++ 上記ではない場合、節点をARRAYに格納、次のTEMP座標に代入
			if( !bFlag1 ){
				pArrayVertex->SetAt( iVertex++, new CVertex( x1, y1, z1 ) );
				tmp_x[kk] = x1;	tmp_y[kk] = y1;	tmp_z[kk] = z1;
				seNum1 = kk;
				kk++;				
			}
			
			if( !bFlag2 ){
				pArrayVertex->SetAt( iVertex++, new CVertex( x2, y2, z2 ) );
				tmp_x[kk] = x2; tmp_y[kk] = y2;	tmp_z[kk] = z2;
				seNum2 = kk;
				kk++;				
			}
			if( !bFlag3 ){
				pArrayVertex->SetAt( iVertex++, new CVertex( x3, y3, z3 ) );
				tmp_x[kk] = x3;	tmp_y[kk] = y3; tmp_z[kk] = z3;
				seNum3 = kk;
				kk++;				
			}
			
			v1 = seNum1;
			v2 = seNum2;
			v3 = seNum3;

			tmpNum = kk;
			
			// ++++++++  要素と節点の関係を格納 ++++++++++++++++++++
			CFace *pFace = new CFace( pArrayVertex->GetAt(v1), pArrayVertex->GetAt(v2), pArrayVertex->GetAt(v3) );
			
			// ++++++++++++++++++  要素を格納 ++++++++++++++++++++++
			pArrayFace->SetAt( iFace++, pFace );
		}
	}

	readStl.m_IndexBuffer = tmp;
	
	// +++++++++++ 領域開放 ++++++++++++++++

	delete [] tmp_x;
	delete [] tmp_y;
	delete [] tmp_z;

	// +++++++ 新しい節点総数を格納 ++++++++
	pArrayVertex->SetSize( iVertex );

	return 1;
}

/**
 * @brief           BINARY形式のデータを読み込みと格納　（STL）
 * @param[in]       pArrayFace      面のバファ
 * @param[in]       pArrayVertex    頂点のバファ
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Store_BinaryStlAsMesh( CArray<CFace>* pArrayFace, 
							      CArray<CVertex>* pArrayVertex )
{
	int i, j, k;
	int nf = 0, nb = 0;
	char s[256];
	float x, y,z, dumy;
	FILE *fp;

	// +++++++++++++ バイナリファイルのＯＰＥＮ ++++++++ 
	if ((fp = fopen( m_InFileName, "rb" )) == NULL)				 return -1;

	// +++++++++++++ モデル名（ヘット）を読み込み ++++++
	if (fread(s,sizeof(char),MAX_HEAD_SIZE,fp) != MAX_HEAD_SIZE) return -1;

	// ++++++++++ 面の総数を読み込む( nf：整数 ）+++++++
	s[MAX_HEAD_SIZE] = '\0';
	readStl.ReadInt(fp, &nf, 0 );

	// +++++ 節点座標を格納するための領域を確保 ++++++
	float *bin_x = NULL;
	float *bin_y = NULL;
	float *bin_z = NULL;

	nb = nf*3;
	bin_x = new float[nb];
	bin_y = new float[nb];
	bin_z = new float[nb];

    // +++ 節点情報などを読み、一時Arrayに代入 +++
	k = 0;
	for (i=0;i<nf;i++){
		
	   // +++++++++++  法線を読み込む ++++++++++
		dumy = readStl.ReadFloat( fp );
		dumy = readStl.ReadFloat( fp );
		dumy = readStl.ReadFloat( fp );

		// ++++++++ 座標を読み込む（代入） ++++++
		for (j=0; j<3; j++){
			x = readStl.ReadFloat( fp );
			y = readStl.ReadFloat( fp );
			z = readStl.ReadFloat( fp );
			
			bin_x[k] = x;
			bin_y[k] = y;
			bin_z[k] = z;			
			k++;
		}

	   // ++++++++++++ ギャップを読み込む ++++++
		fgetc(fp);
		fgetc(fp);
	}

	fclose(fp);

	// ++++++++++++++++++++++++++++ やっと、上部おわり　++++++++++++++++++++++++++++++++++

	// ***********************************************************
	//   頂点構成を計算し、節点、コンビティブティをArrayに代入  //
	// ***********************************************************	
	int   tmpNum;
	float *tmp_x, *tmp_y, *tmp_z; 

	// +++++++++ TMP－ARRAYの領域確保 +++++++++++++++
	tmp_x = new float[m_VertexNum];
	tmp_y = new float[m_VertexNum];
	tmp_z = new float[m_VertexNum];

	int tmp = readStl.m_IndexBuffer;

	int		seNum1, seNum2, seNum3;
	bool	bFlag1, bFlag2, bFlag3;
	float	x1, y1, z1, x2, y2, z2, x3, y3, z3;
	int		v1, v2, v3;

	int		iVertex = 0;
	int		iFace	= 0;
	int		vCount	= 0;
		
	// ***********************************************
	//        面毎にデータ処理                      //
	// ***********************************************
	int kk = 0;
	for( i=0; i<m_FaceNum; i++){

		x1 = bin_x[vCount  ];  y1 = bin_y[vCount  ];  z1 = bin_z[vCount  ];
		x2 = bin_x[vCount+1];  y2 = bin_y[vCount+1];  z2 = bin_z[vCount+1];
		x3 = bin_x[vCount+2];  y3 = bin_y[vCount+2];  z3 = bin_z[vCount+2];

		seNum1 = seNum2 = seNum3 = iVertex;
		bFlag1 = bFlag2 = bFlag3 = false;
		v1 = v2 = v3 = 0;

		if( i > 0 ){

			for( k=0; k<tmpNum; k++ ){

				if( bFlag1 && bFlag2 && bFlag3 ) break; 

				if( !bFlag1 ){
					if( tmp_x[k] == x1 && tmp_y[k] == y1 && tmp_z[k] == z1 ){
						seNum1 = k;
						bFlag1 = true;
					}
				}

				if( !bFlag2 ){
					if( tmp_x[k] == x2 && tmp_y[k] == y2 && tmp_z[k] == z2 ){
						seNum2 = k;
						bFlag2 = true;
					}
				}

				if( !bFlag1 ){
					if( tmp_x[k] == x3 && tmp_y[k] == y3 && tmp_z[k] == z3 ){
						seNum3 = k;
						bFlag3 = true;
					}
				}
			}
		}

		// +++++ 上記ではない場合、節点をARRAYに格納、次のTEMP座標に代入
		if( !bFlag1 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x1, y1, z1 ) );
			tmp_x[kk]	= x1;	tmp_y[kk]	= y1;	tmp_z[kk]	= z1;
			seNum1 = kk;
			kk++;				
		}
			
		if( !bFlag2 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x2, y2, z2 ) );
			tmp_x[kk] = x2;	 tmp_y[kk] = y2;	tmp_z[kk] = z2;
			seNum2 = kk;
			kk++;				
		}
		if( !bFlag3 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x3, y3, z3 ) );
			tmp_x[kk] = x3;	tmp_y[kk] = y3;  tmp_z[kk] = z3;
			seNum3 = kk;
			kk++;				
		}

		v1 = seNum1;
		v2 = seNum2;
		v3 = seNum3;

		tmpNum	= kk;
		vCount	= vCount + 3;

		// ++++++++  要素と節点の関係を格納 ++++++++++++++++++++
		CFace *pFace = new CFace( pArrayVertex->GetAt(v1), pArrayVertex->GetAt(v2), pArrayVertex->GetAt(v3) );
			
		// ++++++++++++++++++  要素を格納 ++++++++++++++++++++++
		pArrayFace->SetAt( iFace++, pFace );
	}

	// +++++++++++ 領域開放 ++++++++++++++++
	delete [] bin_x;
	delete [] bin_y; 
	delete [] bin_z; 

	delete [] tmp_x;
	delete [] tmp_y;
	delete [] tmp_z;

	// +++++++ 新しい節点総数を格納 ++++++++
	pArrayVertex->SetSize( iVertex );

	return 1;
}

/**
 * @brief           ASCII形式のデータを読み込みと格納　　（NPT）
 * @param[in]       pArrayFace      面のバファ
 * @param[in]       pArrayVertex    頂点のバファ
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Store_AsciiNptAsMesh( CArray<CFace>* pArrayFace, 
							      CArray<CVertex>* pArrayVertex )
{


	FILE *fp;
	char	dumy1[128], dumy2[128], dumy3[128], dumy4[128];

	if( ( fp = fopen( m_InFileName, "r" ) ) == NULL) return 0;

	int deg;
	if( m_degree == 2 ) deg = 3;
	else if( m_degree == 3 ) deg = 7;


	int   tmpNum;
	float *tmp_x, *tmp_y, *tmp_z; 

	// +++++++++ TMP－ARRAYの領域確保 +++++++++++++++
	tmp_x = new float[m_VertexNum];
	tmp_y = new float[m_VertexNum];
	tmp_z = new float[m_VertexNum];

	int		i,  k, kk;
	int		seNum1, seNum2, seNum3;
	bool	bFlag1, bFlag2, bFlag3;
	float	x1,y1,z1,x2,y2,z2,x3,y3,z3;

	int		v1, v2, v3;
	int		iVertex = 0;
	int		iFace	= 0;

	fgets( dumy1, sizeof(dumy1), fp);
	fgets( dumy1, sizeof(dumy1), fp);
	fgets( dumy1, sizeof(dumy1), fp);

	kk = 0;
	// ***********************************************
	//        面毎にデータ処理                      //
	// ***********************************************
	for( i=0; i<m_FaceNum; i++){

		//fgets( dumy, sizeof(dumy), fp);

		fscanf( fp, "%s%s%s", dumy1,dumy2,dumy3 );

		fscanf( fp, "%s %f %f %f", dumy1, &x1, &y1, &z1 );
		fscanf( fp, "%s %f %f %f", dumy1, &x2, &y2, &z2 );
		fscanf( fp, "%s %f %f %f", dumy1, &x3, &y3, &z3 );

		for ( int j=0; j<deg; j++ ){
			fscanf( fp, "%s%s%s%s", dumy1,dumy2,dumy3, dumy4 );
		}


		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// ２番目要素以後の場合、座標をTMP座標に代入し新の座標を比較して、
		// 同じならば、節点の番号を取得するのみ、実座標を省略
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		seNum1 = seNum2 = seNum3 = iVertex;
		bFlag1 = bFlag2 = bFlag3 = false;
		v1 = v2 = v3 = 0;

		if( i > 0 ){

			for( k=0; k<tmpNum; k++ ){
				if( bFlag1 && bFlag2 && bFlag3 ) break; 

				if( !bFlag1 ){
					if( tmp_x[k] == x1 && tmp_y[k] == y1 && tmp_z[k] == z1 ){
						seNum1 = k;
						bFlag1 = true;
					}
				}

				if( !bFlag2 ){
					if( tmp_x[k] == x2 && tmp_y[k] == y2 && tmp_z[k] == z2 ){
						seNum2 = k;
						bFlag2 = true;
					}
				}

				if( !bFlag3 ){
					if( tmp_x[k] == x3 && tmp_y[k] == y3 && tmp_z[k] == z3 ){
						seNum3 = k;
						bFlag3 = true;
					}
				}
			}
		}


		// +++++ 上記ではない場合、節点をARRAYに格納、次のTEMP座標に代入
		if( !bFlag1 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x1, y1, z1 ) );
			tmp_x[kk] = x1;	tmp_y[kk] = y1;	tmp_z[kk] = z1;
			seNum1 = kk;
			kk++;				
		}
			
		if( !bFlag2 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x2, y2, z2 ) );
			tmp_x[kk] = x2; tmp_y[kk] = y2;	tmp_z[kk] = z2;
			seNum2 = kk;
			kk++;				
		}
		if( !bFlag3 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x3, y3, z3 ) );
			tmp_x[kk] = x3;	tmp_y[kk] = y3; tmp_z[kk] = z3;
			seNum3 = kk;
			kk++;				
		}

		v1 = seNum1;
		v2 = seNum2;
		v3 = seNum3;

		tmpNum = kk;

		// ++++++++  要素と節点の関係を格納 ++++++++++++++++++++
		CFace *pFace = new CFace( pArrayVertex->GetAt(v1), pArrayVertex->GetAt(v2), pArrayVertex->GetAt(v3) );
			
		// ++++++++++++++++++  要素を格納 ++++++++++++++++++++++
		pArrayFace->SetAt( iFace++, pFace );

	}
	
	fclose(fp);

	// +++++++++++ 領域開放 ++++++++++++++++

	delete [] tmp_x;
	delete [] tmp_y;
	delete [] tmp_z;

	// +++++++ 新しい節点総数を格納 ++++++++
	pArrayVertex->SetSize( iVertex );

	return 1;
}

/**
 * @brief           BINARY形式のデータを読み込みと格納　　（NPT）
 * @param[in]       pArrayFace      面のバファ
 * @param[in]       pArrayVertex    頂点のバファ
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Store_BinaryNptAsMesh( CArray<CFace>* pArrayFace, 
							      CArray<CVertex>* pArrayVertex )
{
	int i, j, k;
	int nf = 0, nb = 0;
	float x, y,z, dumy;
	FILE *fp;

	// +++++++++++++ バイナリファイルのＯＰＥＮ ++++++++ 
	if ((fp = fopen( m_InFileName, "rb" )) == NULL)				 return -1;

	int deg;
	if( m_degree == 2 ) deg = 3;
	else if( m_degree == 3 ) deg = 7;

	// ++++++++++ 面の総数を読み込む( nf：整数 ）+++++++
	readStl.ReadInt(fp, &nf, 0 );

	// +++++ 節点座標を格納するための領域を確保 ++++++
	float *bin_x = NULL;
	float *bin_y = NULL;
	float *bin_z = NULL;

	nb = nf*3;
	bin_x = new float[nb];
	bin_y = new float[nb];
	bin_z = new float[nb];

    // +++ 節点情報などを読み、一時Arrayに代入 +++
	k = 0;
	for (i=0;i<nf;i++){

		// ++++++++ 座標を読み込む（代入） ++++++
		for (j=0; j<3; j++){
			x = readStl.ReadFloat( fp );
			y = readStl.ReadFloat( fp );
			z = readStl.ReadFloat( fp );
			
			bin_x[k] = x;
			bin_y[k] = y;
			bin_z[k] = z;			
			k++;
		}

		for (j=0; j<deg; j++){
			dumy = readStl.ReadFloat( fp );
			dumy = readStl.ReadFloat( fp );
			dumy = readStl.ReadFloat( fp );
		}
	}
	fclose(fp);

	// ++++++++++++++++++++++++++++ やっと、上部おわり　++++++++++++++++++++++++++++++++++

	// ***********************************************************
	//   頂点構成を計算し、節点、コンビティブティをArrayに代入  //
	// ***********************************************************	
	int   tmpNum;
	float *tmp_x, *tmp_y, *tmp_z; 

	// +++++++++ TMP－ARRAYの領域確保 +++++++++++++++
	tmp_x = new float[m_VertexNum];
	tmp_y = new float[m_VertexNum];
	tmp_z = new float[m_VertexNum];

	int tmp = readStl.m_IndexBuffer;

	int		seNum1, seNum2, seNum3;
	bool	bFlag1, bFlag2, bFlag3;
	float	x1, y1, z1, x2, y2, z2, x3, y3, z3;
	int		v1, v2, v3;

	int		iVertex = 0;
	int		iFace	= 0;
	int		vCount	= 0;
		
	// ***********************************************
	//        面毎にデータ処理                      //
	// ***********************************************
	int kk = 0;
	for( i=0; i<nf; i++){

		x1 = bin_x[vCount  ];  y1 = bin_y[vCount  ];  z1 = bin_z[vCount  ];
		x2 = bin_x[vCount+1];  y2 = bin_y[vCount+1];  z2 = bin_z[vCount+1];
		x3 = bin_x[vCount+2];  y3 = bin_y[vCount+2];  z3 = bin_z[vCount+2];

		seNum1 = seNum2 = seNum3 = iVertex;
		bFlag1 = bFlag2 = bFlag3 = false;
		v1 = v2 = v3 = 0;

		if( i > 0 ){

			for( k=0; k<tmpNum; k++ ){

				if( bFlag1 && bFlag2 && bFlag3 ) break; 

				if( !bFlag1 ){
					if( tmp_x[k] == x1 && tmp_y[k] == y1 && tmp_z[k] == z1 ){
						seNum1 = k;
						bFlag1 = true;
					}
				}

				if( !bFlag2 ){
					if( tmp_x[k] == x2 && tmp_y[k] == y2 && tmp_z[k] == z2 ){
						seNum2 = k;
						bFlag2 = true;
					}
				}

				if( !bFlag1 ){
					if( tmp_x[k] == x3 && tmp_y[k] == y3 && tmp_z[k] == z3 ){
						seNum3 = k;
						bFlag3 = true;
					}
				}
			}
		}

		// +++++ 上記ではない場合、節点をARRAYに格納、次のTEMP座標に代入
		if( !bFlag1 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x1, y1, z1 ) );
			tmp_x[kk]	= x1;	tmp_y[kk]	= y1;	tmp_z[kk]	= z1;
			seNum1 = kk;
			kk++;				
		}
			
		if( !bFlag2 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x2, y2, z2 ) );
			tmp_x[kk] = x2;	 tmp_y[kk] = y2;	tmp_z[kk] = z2;
			seNum2 = kk;
			kk++;				
		}
		if( !bFlag3 ){
			pArrayVertex->SetAt( iVertex++, new CVertex( x3, y3, z3 ) );
			tmp_x[kk] = x3;	tmp_y[kk] = y3;  tmp_z[kk] = z3;
			seNum3 = kk;
			kk++;				
		}

		v1 = seNum1;
		v2 = seNum2;
		v3 = seNum3;

		tmpNum	= kk;
		vCount	= vCount + 3;

		// ++++++++  要素と節点の関係を格納 ++++++++++++++++++++
		CFace *pFace = new CFace( pArrayVertex->GetAt(v1), pArrayVertex->GetAt(v2), pArrayVertex->GetAt(v3) );
			
		// ++++++++++++++++++  要素を格納 ++++++++++++++++++++++
		pArrayFace->SetAt( iFace++, pFace );
	}

	// +++++++++++ 領域開放 ++++++++++++++++
	delete [] bin_x;
	delete [] bin_y; 
	delete [] bin_z; 

	delete [] tmp_x;
	delete [] tmp_y;
	delete [] tmp_z;

	// +++++++ 新しい節点総数を格納 ++++++++
	pArrayVertex->SetSize( iVertex );

	return 1;
}

// ----------------------------------------------------------
//                      共有頂点の平均化 
// ----------------------------------------------------------

/**
 * @brief           共通エージ及び頂点によって、隣接面を作成
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Vertex_BuildAdjacency()
{
	int i;

	// ++++++++++ 要素と節点の存在をチェックする ++++++++
	int nVertex = m_ArrayVertex.GetSize();
	int nFace	= m_ArrayFace.GetSize();

	if( nVertex < 3 || nFace <= 0 )	return 0;

	// ++++ まず、隣接面を一回クリアする +++++++++++++++
	for( i=0; i<nVertex; i++){
		m_ArrayVertex[i]->RemoveAllFaceNeighbor();
	}

	//*********************************************
	//  面毎に、隣接頂点を設定する 
	//*********************************************
	for( i=0; i<nFace; i++ ){
		CFace *pFaceCurrent = m_ArrayFace[i];
		for( int j=0; j<3; j++){
			pFaceCurrent->v(j)->AddNeighbor( pFaceCurrent );
		}
	}
	return true;
}

//********************************************
//       頂点の法線計算
//********************************************
/**
 * @brief           頂点の法線計算
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Vertex_CalcNormalPerVertex(void)
{

	int nFace = m_ArrayFace.GetSize();
	int i, j, k;

	CFace   *mFace;
	CVector vector;

	// **************************************************
	//       面法線を頂点法線にする（初期化）
	// **************************************************
	for( i=0; i<nFace; i++ ){
		mFace = m_ArrayFace[i];
		vector = mFace->GetNormal();
		for( j=0; j<3; j++){
			mFace->SetNormal_New(vector, j, 1);
		}
	}

	// **************************************************
	//       頂点法線を計算する。
	// **************************************************
	for( i=0; i<nFace; i++ ){

		// ++++++++++ 当面を取得 +++++++++++++++++
		mFace = m_ArrayFace[i];
		
		for( j=0; j<3; j++){
			// +++++++++ 当面頂点を取得 +++++++
			CVertex *pVertex = mFace->v(j);

			int NbNormal = 0;
			float tt;

			CVector hvector;
			CVector fvector;
			CVector gvector;
			
			// ++++++++++ この頂点の隣接面の総数を取得 ++++++
			int NumOfFaceNeighbor = pVertex->GetNumOfFaceNeighbor();
			
			// ++++++++++ この面の法線を取得 ++++++++++++++++			
			fvector = mFace->GetNormal();
			
			for( k=0; k<NumOfFaceNeighbor; k++){
				//  ++++++++++++++++ 隣接面を取得 +++++++++++++
				CFace *pFace = pVertex->GetFaceNeighbor(k);
			
				//  ++++++++++++ 隣接面の法線を取得 +++++++++++
				hvector = pFace->GetNormal();

				// ++++++++++ 両法線の内積を求める +++++++++++++++
				tt = hvector.x()*fvector.x() + hvector.y()*fvector.y() 
					                         + hvector.z()*fvector.z(); 

				// +++ 面間の角度が10°以上の場合、法線を加算 ++++
				if ( tt >0.85 ){
					NbNormal++;
					gvector += pFace->GetNormal();
				}
			}

			// ++++  隣接面が10°以上ならば、新しい法線にする ++++
			if(NbNormal>=1)
			{
				// +++++++++ ベクトルの正規化 +++++++++
				vector.NormalizeL2();
			
				// ++++ 頂点に新しい法線を格納 ++++++++ 
				mFace->SetNormal_New(gvector, j, NbNormal);
			}
		}
	}

	return 1;
}

// *********************************************************
//  関数名： Vertex_CalcNormalPerFace
//  機能  ： 面の法線計算
//   in   ： なし
//  out   ： エラーコード　（0：エラー　1：正常）
// *********************************************************
/**
 * @brief           面の法線計算
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::Vertex_CalcNormalPerFace(void)
{

	for( int i=0; i<m_ArrayFace.GetSize(); i++){
		CFace *pFace = m_ArrayFace[i];
		pFace->CalculateNormal();
	}
	return 1;
}

// *********************************************************
//  関数名： Vertex_FreeArray
//  機能  ： 面と頂点のバファを削除
//   in   ： なし
//  out   ： なし
// *********************************************************

//********************************************
// Free
//********************************************
/**
 * @brief           面と頂点のバファを削除
 */
void CData::Vertex_FreeArray()
{
	m_ArrayVertex.Free();
	m_ArrayFace.Free();
}


// ----------------------------------------------------------
//                      長田ﾊﾞｯﾁ
// ----------------------------------------------------------

// *********************************************************
//  関数名： MeshPatch_GetCoef2Deg
//  機能  ： 長田パッチの係数を求める（二次曲面））
//   in   ： なし
//  out   ： エラーコード　（0：エラー　1：正常）
// *********************************************************
/**
 * @brief           長田パッチの係数を求める（二次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetCoef2Deg()
{


	double norm[3]; 
	double p1[3], norm1[3];
	double p2[3], norm2[3];
	double p3[3], norm3[3];
	double c1_aws[3], c2_aws[3], c3_aws[3];


	MeshPatch mesh_path;

	mesh_path.m_alwAngle = (3.1415926)*m_CheckAngle/180;

	CVector *pVector;

	for( int i=0; i<m_ArrayFace.GetSize(); i++){

		CFace *pFace = m_ArrayFace[i];
		pVector = pFace->GetNormal();

		// +++++++++++ Normal on Face  ++++++++++++++++
		SetValue( norm, pVector->x(), pVector->y(), pVector->z() );

		// +++++++++++ Vertex coordinate and normal vector on First Point ++++++++++++++++
		SetValue( p1,    pFace->v(0)->x(),		  pFace->v(0)->y(),		   pFace->v(0)->z() );
		SetValue( norm1, pFace->m_NewNormal[0].x, pFace->m_NewNormal[0].y, pFace->m_NewNormal[0].z );

		// +++++++++++ Vertex coordinate and normal vector on Second Point ++++++++++++++++
		SetValue( p2,    pFace->v(1)->x(),		  pFace->v(1)->y(),		   pFace->v(1)->z() );
		SetValue( norm2, pFace->m_NewNormal[1].x, pFace->m_NewNormal[1].y, pFace->m_NewNormal[1].z );

		// +++++++++++ Vertex coordinate and normal vector on Third Point ++++++++++++++++
		SetValue( p3,    pFace->v(2)->x(),		  pFace->v(2)->y(),		   pFace->v(2)->z() );
		SetValue( norm3, pFace->m_NewNormal[2].x, pFace->m_NewNormal[2].y, pFace->m_NewNormal[2].z );

		mesh_path.CalcGetCoef2Deg( norm, p1, norm1, p2, norm2, p3, norm3, c1_aws, c2_aws, c3_aws );

		pFace->SetCoef(c1_aws, 0);
		pFace->SetCoef(c2_aws, 1);
		pFace->SetCoef(c3_aws, 2);
	}

	return 1;

}

/**
 * @brief           長田パッチの係数を求める（三次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetCoef3Deg()
{


	double norm[3]; 
	double p1[3], norm1[3];
	double p2[3], norm2[3];
	double p3[3], norm3[3];
	double aswC1[3], aswC2[3], aswC3[3], aswC4[3], aswC5[3], aswC6[3], aswC7[3];


	MeshPatch mesh_path;

	mesh_path.m_alwAngle = (3.1415926)*m_CheckAngle/180;

	CVector *pVector;

	for( int i=0; i<m_ArrayFace.GetSize(); i++){

		CFace *pFace = m_ArrayFace[i];
		pVector = pFace->GetNormal();


		// +++++++++++ Normal on Face  ++++++++++++++++
		SetValue( norm, pVector->x(), pVector->y(), pVector->z() );

		// +++++++++++ Vertex coordinate and normal vector on First Point ++++++++++++++++
		SetValue( p1,    pFace->v(0)->x(),		  pFace->v(0)->y(),		   pFace->v(0)->z() );
		SetValue( norm1, pFace->m_NewNormal[0].x, pFace->m_NewNormal[0].y, pFace->m_NewNormal[0].z );

		// +++++++++++ Vertex coordinate and normal vector on Second Point ++++++++++++++++
		SetValue( p2,    pFace->v(1)->x(),		  pFace->v(1)->y(),		   pFace->v(1)->z() );
		SetValue( norm2, pFace->m_NewNormal[1].x, pFace->m_NewNormal[1].y, pFace->m_NewNormal[1].z );

		// +++++++++++ Vertex coordinate and normal vector on Third Point ++++++++++++++++
		SetValue( p3,    pFace->v(2)->x(),		  pFace->v(2)->y(),		   pFace->v(2)->z() );
		SetValue( norm3, pFace->m_NewNormal[2].x, pFace->m_NewNormal[2].y, pFace->m_NewNormal[2].z );

		mesh_path.CalcGetCoef3Deg( norm, p1, norm1, p2, norm2, p3, norm3, 
									aswC1, aswC2, aswC3, aswC4, aswC5, aswC6, aswC7 );

		pFace->SetCoef(aswC1, 0);
		pFace->SetCoef(aswC2, 1);
		pFace->SetCoef(aswC3, 2);
		pFace->SetCoef(aswC4, 3);
		pFace->SetCoef(aswC5, 4);
		pFace->SetCoef(aswC6, 5);
		pFace->SetCoef(aswC7, 6);
	}

	return 1;

}

/**
 * @brief           長田係数及びη、ζより、曲面上の点を求める（二次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetValue2Deg()
{

	double p1[3], c1[3];
	double p2[3], c2[3];
	double p3[3], c3[3];

	double		pos_aws[3];
	MeshPatch	mesh_path;

	bAsw[0]	= new VertexD[m_FaceNum];

	for( int i=0; i<m_FaceNum; i++){

		CFace *pFace = m_ArrayFace[i];
		
		SetValue( p1, pFace->v(0)->x(), pFace->v(0)->y(), pFace->v(0)->z() );
		SetValue( p2, pFace->v(1)->x(),	pFace->v(1)->y(), pFace->v(1)->z() );
		SetValue( p3, pFace->v(2)->x(),	pFace->v(2)->y(), pFace->v(2)->z() );

		SetValue( c1, pFace->m_Coef[0].x, pFace->m_Coef[0].y,  pFace->m_Coef[0].z );
		SetValue( c2, pFace->m_Coef[1].x, pFace->m_Coef[1].y,  pFace->m_Coef[1].z );
		SetValue( c3, pFace->m_Coef[2].x, pFace->m_Coef[2].y,  pFace->m_Coef[2].z );

		// +++++ 近似曲面補正ありの追加した頂点の座標
		mesh_path.Evaluate2Deg( m_eta, m_zeta, p1, p2, p3, c1, c2, c3, pos_aws ); 

		bAsw[0][i].x = (float) pos_aws[0];
		bAsw[0][i].y = (float) pos_aws[1];
		bAsw[0][i].z = (float) pos_aws[2];
	}

	return 1;
}

/**
 * @brief           長田係数及びη、ζより、曲面上の点を求める（三次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetValue3Deg()
{

	double p1[3], p2[3], p3[3];
	double c1[3], c2[3], c3[3], c4[3], c5[3], c6[3], c7[3];

	double		pos_aws[3];
	MeshPatch	mesh_path;

	bAsw[0]	= new VertexD[m_FaceNum];

	for( int i=0; i<m_FaceNum; i++){

		CFace *pFace = m_ArrayFace[i];

		SetValue( p1, pFace->v(0)->x(), pFace->v(0)->y(), pFace->v(0)->z() );
		SetValue( p2, pFace->v(1)->x(),	pFace->v(1)->y(), pFace->v(1)->z() );
		SetValue( p3, pFace->v(2)->x(),	pFace->v(2)->y(), pFace->v(2)->z() );

		SetValue( c1, pFace->m_Coef[0].x, pFace->m_Coef[0].y,  pFace->m_Coef[0].z );
		SetValue( c2, pFace->m_Coef[1].x, pFace->m_Coef[1].y,  pFace->m_Coef[1].z );
		SetValue( c3, pFace->m_Coef[2].x, pFace->m_Coef[2].y,  pFace->m_Coef[2].z );
		SetValue( c4, pFace->m_Coef[3].x, pFace->m_Coef[3].y,  pFace->m_Coef[3].z );
		SetValue( c5, pFace->m_Coef[4].x, pFace->m_Coef[4].y,  pFace->m_Coef[4].z );
		SetValue( c6, pFace->m_Coef[5].x, pFace->m_Coef[5].y,  pFace->m_Coef[5].z );
		SetValue( c7, pFace->m_Coef[6].x, pFace->m_Coef[6].y,  pFace->m_Coef[6].z );

		// +++++ 近似曲面補正ありの追加した頂点の座標
		mesh_path.Evaluate3Deg( m_eta, m_zeta, p1, p2, p3, 
						c1, c2, c3, c4, c5, c6, c7, pos_aws );

		bAsw[0][i].x = (float) pos_aws[0];
		bAsw[0][i].y = (float) pos_aws[1];
		bAsw[0][i].z = (float) pos_aws[2];
	}

	return 1;
}

/**
 * @brief           長田係数及びη、ζより、曲面上の各辺の中点を求める（二次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetValue3Point2Deg()
{

	double p1[3], c1[3];
	double p2[3], c2[3];
	double p3[3], c3[3];

	double		pos_aws[3];
	MeshPatch	mesh_path;

	bAsw[0]	= new VertexD[m_FaceNum];
	bAsw[1]	= new VertexD[m_FaceNum];
	bAsw[2]	= new VertexD[m_FaceNum];

	for( int i=0; i<m_FaceNum; i++){

		CFace *pFace = m_ArrayFace[i];


		SetValue( p1, pFace->v(0)->x(), pFace->v(0)->y(), pFace->v(0)->z() );
		SetValue( p2, pFace->v(1)->x(),	pFace->v(1)->y(), pFace->v(1)->z() );
		SetValue( p3, pFace->v(2)->x(),	pFace->v(2)->y(), pFace->v(2)->z() );

		SetValue( c1, pFace->m_Coef[0].x, pFace->m_Coef[0].y,  pFace->m_Coef[0].z );
		SetValue( c2, pFace->m_Coef[1].x, pFace->m_Coef[1].y,  pFace->m_Coef[1].z );
		SetValue( c3, pFace->m_Coef[2].x, pFace->m_Coef[2].y,  pFace->m_Coef[2].z );

		// +++++ 近似曲面補正ありの追加した頂点の座標
		mesh_path.Evaluate2Deg( 0.5, 0.0, p1, p2, p3, c1, c2, c3, pos_aws ); 
		bAsw[0][i].x = (float) pos_aws[0];
		bAsw[0][i].y = (float) pos_aws[1];
		bAsw[0][i].z = (float) pos_aws[2];

		mesh_path.Evaluate2Deg( 1.0, 0.5, p1, p2, p3, c1, c2, c3, pos_aws ); 
		bAsw[1][i].x = (float) pos_aws[0];
		bAsw[1][i].y = (float) pos_aws[1];
		bAsw[1][i].z = (float) pos_aws[2];

		mesh_path.Evaluate2Deg( 0.5, 0.5, p1, p2, p3, c1, c2, c3, pos_aws ); 
		bAsw[2][i].x = (float) pos_aws[0];
		bAsw[2][i].y = (float) pos_aws[1];
		bAsw[2][i].z = (float) pos_aws[2];

	}

	m_eta	= 0.5;
	m_zeta	= 0.0;

	return 1;

}

/**
 * @brief           長田係数及びη、ζより、曲面上の各辺の中点を求める（三次曲面）
 * @return          エラーコード　（0：エラー　1：正常）
 */
int CData::MeshPatch_GetValue3Point3Deg()
{

	double p1[3], p2[3], p3[3];
	double c1[3], c2[3], c3[3], c4[3], c5[3], c6[3], c7[3];

	double		pos_aws[3];
	MeshPatch	mesh_path;

	bAsw[0]	= new VertexD[m_FaceNum];
	bAsw[1]	= new VertexD[m_FaceNum];
	bAsw[2]	= new VertexD[m_FaceNum];

	for( int i=0; i<m_FaceNum; i++){

		CFace *pFace = m_ArrayFace[i];


		SetValue( p1, pFace->v(0)->x(), pFace->v(0)->y(), pFace->v(0)->z() );
		SetValue( p2, pFace->v(1)->x(),	pFace->v(1)->y(), pFace->v(1)->z() );
		SetValue( p3, pFace->v(2)->x(),	pFace->v(2)->y(), pFace->v(2)->z() );

		SetValue( c1, pFace->m_Coef[0].x, pFace->m_Coef[0].y,  pFace->m_Coef[0].z );
		SetValue( c2, pFace->m_Coef[1].x, pFace->m_Coef[1].y,  pFace->m_Coef[1].z );
		SetValue( c3, pFace->m_Coef[2].x, pFace->m_Coef[2].y,  pFace->m_Coef[2].z );
		SetValue( c4, pFace->m_Coef[3].x, pFace->m_Coef[3].y,  pFace->m_Coef[3].z );
		SetValue( c5, pFace->m_Coef[4].x, pFace->m_Coef[4].y,  pFace->m_Coef[4].z );
		SetValue( c6, pFace->m_Coef[5].x, pFace->m_Coef[5].y,  pFace->m_Coef[5].z );
		SetValue( c7, pFace->m_Coef[6].x, pFace->m_Coef[6].y,  pFace->m_Coef[6].z );

		// +++++ 近似曲面補正ありの追加した頂点の座標
		mesh_path.Evaluate3Deg( 0.5, 0.0, p1, p2, p3, c1, c2, c3, c4, c5, c6, c7, pos_aws ); 
		bAsw[0][i].x = (float) pos_aws[0];
		bAsw[0][i].y = (float) pos_aws[1];
		bAsw[0][i].z = (float) pos_aws[2];

		mesh_path.Evaluate3Deg( 1.0, 0.5, p1, p2, p3, c1, c2, c3, c4, c5, c6, c7, pos_aws ); 
		bAsw[1][i].x = (float) pos_aws[0];
		bAsw[1][i].y = (float) pos_aws[1];
		bAsw[1][i].z = (float) pos_aws[2];

		mesh_path.Evaluate3Deg( 0.5, 0.5, p1, p2, p3, c1, c2, c3, c4, c5, c6, c7, pos_aws ); 
		bAsw[2][i].x = (float) pos_aws[0];
		bAsw[2][i].y = (float) pos_aws[1];
		bAsw[2][i].z = (float) pos_aws[2];

	}

	m_eta	= 0.5;
	m_zeta	= 0.0;
	
	return 0;
}

/**
 * @brief           長田係数に使う計算結果を削除
 */
void CData::MeshPatch_freeCoef()
{
	
	if( bAsw[0] != NULL ){
		delete [] bAsw[0];
		bAsw[0] = NULL;
	}

	if( bAsw[1] != NULL ){
		delete [] bAsw[1];
		bAsw[1] = NULL;
	}

	if( bAsw[2] != NULL ){
		delete [] bAsw[2];
		bAsw[2] = NULL;
	}
	
}
