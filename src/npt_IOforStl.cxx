/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   npt_IOforStl.cxx
 * @brief  Npatch IOForStl Class
 * @author aics
 */

#include <string.h>
#include <stdio.h>

#include "npt_IOforStl.h"

/**
 * @brief           コンストラクタ　（初期化）
 */
IOForStl::IOForStl(void)
{
	m_pBuffer		 = NULL;
	m_pBufferLine[0] = '\0';
	m_IndexBuffer	 = 0;

	FaceNum		= 0;
	VertexNum	= 0;
}

/**
 * @brief           デストラクタ　（バッファを削除）
 */
IOForStl::~IOForStl(void)
{
	Delet_Buff();
}

/**
 * @brief           バッファ削除
 */
void IOForStl::Delet_Buff()
{
	if( m_pBuffer != NULL ) delete [] m_pBuffer;
	m_IndexBuffer	= 0;
	m_pBuffer		= NULL;
}

/**
 * @brief           IOクラスのメイン関数
 * @param[in]       InFile  入力ファイル名
 * @param[in]       ExType  拡張子　0：Stl　1：Npt
 * @param[in]       deg     曲面の次数
 */
void IOForStl::LoadStl( char *InFile, int ExType, int deg )
{
	
	LoadStlFile( InFile );
	Stl_Type = CheckStlType( InFile );

	if( ExType == 0 ){
		if( Stl_Type == 0 ){
			FaceNum		= GetStlNum_Ascii();
			VertexNum	= FaceNum*3;
		}
		else if( Stl_Type == 1 ){
			FaceNum		= GetStlNum_Binary( InFile );
			VertexNum	= FaceNum*3;
		}
	} else if( ExType == 1 ){
		if( Stl_Type == 0 ){
			FaceNum		= GetNptNum_Ascii( InFile );
			VertexNum	= FaceNum*3;
		}
		else if( Stl_Type == 1 ){
			FaceNum		= GetNptNum_Binary( InFile );
			VertexNum	= FaceNum*3;
		}		
	}
	
	std::cout<< "InFile = " << InFile << std::endl;
	std::cout<< "Stl_Type = " << Stl_Type << std::endl;
	std::cout<< "ExType = " << ExType << std::endl;
	std::cout<< "FaceNum = " << FaceNum << std::endl;
}

/**
 * @brief           STLファイルの読み込み
 * @param[in]       InFile  STLファイル名
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::LoadStlFile( char *InFile )
{

	FILE* fp;

	if( ( fp = fopen( InFile, "rb" ) ) == NULL) return 0;

	//ファイルの大きさセット
	fseek(fp, 0, SEEK_END);
	m_SizeFile = ftell(fp);

	// ++++++++++++++ データ用領域を確保 ++++++++++++++++++++++++++++++
	m_pBuffer = new char[m_SizeFile];
	if(m_pBuffer == NULL) return 0;

	//	読み込みます
	fseek(fp, 0, SEEK_SET);
	fread(m_pBuffer, 1, m_SizeFile, fp);
	fclose(fp);

	return 1;
}

/**
 * @brief           バッファから、行を抽出
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::ReadLine()
{
	m_pBufferLine[0] = '\0';
	int i=0;
	do
		m_pBufferLine[i++] = m_pBuffer[m_IndexBuffer++];
	while(m_pBuffer[m_IndexBuffer-1] != '\n' &&
		    //i < 120 &&
				m_IndexBuffer < m_SizeFile);

	m_pBufferLine[i-1] = '\0';

	return 1;
}


/**
 * @brief           バッファ行から、必要の文字を探し出す
 * @param[in]       string  文字列
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::OffsetToString(char *string)
{

	while(m_IndexBuffer < m_SizeFile)
	{
		ReadLine();
		char *adr = strstr(m_pBufferLine,string);
		if(strstr(m_pBufferLine,string) != NULL)
		{
			m_IndexBuffer = m_IndexBuffer - strlen(m_pBufferLine) - 1 + (adr-m_pBufferLine);
			assert(m_IndexBuffer >= 0);
			return 1;
		}
	}

	return 0;
}

/**
 * @brief           バッファ上に一WORDを読む
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::ReadWord()
{

	m_pBufferWord[0] = '\0';
	int i=0;

	// Jump to next valid character
	while((m_pBuffer[m_IndexBuffer] == '\n' || 
		     m_pBuffer[m_IndexBuffer] == '\t' || 
		     m_pBuffer[m_IndexBuffer] == '\r' || 
		     m_pBuffer[m_IndexBuffer] == ' ') &&
				 m_IndexBuffer < m_SizeFile)
		m_IndexBuffer++;

	// Check eof
	if(m_IndexBuffer >= m_SizeFile)
		return 0;

	do
		m_pBufferWord[i++] = m_pBuffer[m_IndexBuffer++];
	while(m_pBuffer[m_IndexBuffer-1] != '\n' && 
		    m_pBuffer[m_IndexBuffer-1] != '\t' && 
		    m_pBuffer[m_IndexBuffer-1] != '\r' && 
		    m_pBuffer[m_IndexBuffer-1] != ' ' && 
		    i < MAX_WORD_STL &&
				m_IndexBuffer < m_SizeFile);

	m_pBufferWord[i-1] = '\0';

	return 1;
}

/**
 * @brief           STLの三角形の総数を取得（STLのAsciiコード）
 * @return          要素の総数
 */
int IOForStl::GetStlNum_Ascii(void)
{

	int tmp = m_IndexBuffer;
	int nf = 0;
	
	// ++++++ データ行を読込、FACET NORMALの場合、総数を加算 +++++++ 
	while(m_IndexBuffer < m_SizeFile){
		ReadLine();
		if(strstr(m_pBufferLine,"facet normal") != NULL || 
		   strstr(m_pBufferLine,"FACET NORMAL") != NULL) nf++;
	}

	m_IndexBuffer = tmp;

	return nf;
}

/**
 * @brief           STLの三角形の総数を取得（NptのAsciiコード）
 * @param[in]       InFile  入力ファイル名、XXX.npt
 * @return          要素の総数
 */
int IOForStl::GetNptNum_Ascii( char* InFile ){

	int nf = 0;
	FILE *fp;

	if( ( fp = fopen( InFile, "r" ) ) == NULL) return 0;

	fscanf(fp,"%d",&nf );
	fclose(fp);

	return nf;
}

/**
 * @brief           STLの三角形の総数を取得（Nptの	Binaryコード）
 * @param[in]       InFile  入力ファイル名、XXX.npt
 * @return          要素の総数
 */
int IOForStl::GetNptNum_Binary( char* InFile ){

	int nf = 0, nb = 0;
	FILE *fp;

	// ***********************************************
	//     バイナリデータの読込                     //
	// ***********************************************
	if ((fp = fopen(InFile, "rb" )) == NULL) {
		return -1;
	}

	// ++++++++++++++ 要素総数の読み込み ++++++++++++++++++++++++++++
	ReadInt(fp, &nf, 0 );

	fclose(fp);

	return nf;
}

/**
 * @brief           長田係数及び頂点座標を出力(Asicc)
 * @param[in]       name    出力ファイル名
 * @param[in]       wFace   三角形情報
 * @param[in]       FaceNum 三角形の総数
 * @param[in]       degree  曲面の次数
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::SaveNPT_Ascii( char *name, CArray<CFace>& wFace, int FaceNum, int degree )
{


	char OutName[256];
	strcpy( OutName, name );
	
//	strcat( OutName, "_New.NPT" );

//	if( degree == 2 )	strcat( OutName, "_Ascii_2.NPT" );
//	else if( degree == 3 )	strcat( OutName, "_Ascii_3.NPT" );

	std::cout<< "InFile = " << OutName << std::endl;


	FILE *fp;
	if( ( fp = fopen( OutName, "w" ) ) == NULL) return 0;
	fprintf( fp, " %d \n", FaceNum );
	fprintf( fp, "  file name : %s \n\n", OutName );


	CFace   *pFace;
	for( int i=0; i<FaceNum; i++){
		
		pFace = wFace[i];
		assert(pFace != NULL);

		// ++++++++++++++++++++++++  面のIDの書出し ++++++++++++++++++++++++
		fprintf( fp, "  facet ID  %10d \n", i );

		// ++++++++++++++++++++++++  頂点の座標書出し ++++++++++++++++++++++
		for( int j=0; j<3; j++){
			fprintf( fp, "      vertex  % e  % e  % e\n", pFace->v(j)->x(), pFace->v(j)->y(), pFace->v(j)->z());
		}

		// ++++++++++++++++++++++++  長田パッチの係数書出し ++++++++++++++++++++++
		fprintf( fp, "      coef1   % e  % e  % e\n", pFace->m_Coef[0].x, pFace->m_Coef[0].y, pFace->m_Coef[0].z );
		fprintf( fp, "      coef2   % e  % e  % e\n", pFace->m_Coef[1].x, pFace->m_Coef[1].y, pFace->m_Coef[1].z );
		fprintf( fp, "      coef3   % e  % e  % e\n", pFace->m_Coef[2].x, pFace->m_Coef[2].y, pFace->m_Coef[2].z );
		if( degree == 3 ){
			fprintf( fp, "      coef4   % e  % e  % e\n", pFace->m_Coef[3].x, pFace->m_Coef[3].y, pFace->m_Coef[3].z );
			fprintf( fp, "      coef5   % e  % e  % e\n", pFace->m_Coef[4].x, pFace->m_Coef[4].y, pFace->m_Coef[4].z );
			fprintf( fp, "      coef6   % e  % e  % e\n", pFace->m_Coef[5].x, pFace->m_Coef[5].y, pFace->m_Coef[5].z );
			fprintf( fp, "      coef7   % e  % e  % e\n", pFace->m_Coef[6].x, pFace->m_Coef[6].y, pFace->m_Coef[6].z );
		
		}
	}

	// ********** 終了キーワットの書出し +++++++++++++++++++
	fprintf( fp, "end Face Informaion\n");

	fclose(fp);

	return 1;
}

/**
 * @brief           長田係数及び頂点座標を出力(Binary)
 * @param[in]       name    出力ファイル名
 * @param[in]       wFace   三角形情報
 * @param[in]       FaceNum 三角形の総数
 * @param[in]       degree  曲面の次数
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::SaveNPT_Binary( char *name, CArray<CFace>& wFace, int FaceNum, int degree )
{

	char OutName[256];
	strcpy( OutName, name );
	
	FILE *fpw;
	fpw = fopen( OutName, "wb" );
	if( fpw == NULL ){
		printf( "書込用 %sファイルが開けません\n", OutName );
		return -1;
	}

	int NumOfLoop;
	float   x,y,z;
	int     i, j;

	if( degree == 2 ) NumOfLoop = 3;
	else if( degree == 3 ) NumOfLoop = 7;

	CFace   *pFace;
		
	// ++++++++++++++ 要素総数の書出し +++++++++++++++++++++++++++++		

	int NbFace = FaceNum;
	fwrite( &NbFace, sizeof(int), 1, fpw );


	for( i=0; i<NbFace; i++){
		pFace = wFace[i];
		//ASSERT(pFace != NULL);

		// ++++++++++++++++++++++++  頂点ループの書出し ++++++++++++++++++++
		for( j=0; j<3; j++){
			x = pFace->v(j)->x();
			y = pFace->v(j)->y();
			z = pFace->v(j)->z();

			fwrite( &x, sizeof(float), 1, fpw );
			fwrite( &y, sizeof(float), 1, fpw );
			fwrite( &z, sizeof(float), 1, fpw );
		}

		for( j=0; j<NumOfLoop; j++){
			x = pFace->m_Coef[j].x;
			y = pFace->m_Coef[j].y;
			z = pFace->m_Coef[j].z;

			fwrite( &x, sizeof(float), 1, fpw );
			fwrite( &y, sizeof(float), 1, fpw );
			fwrite( &z, sizeof(float), 1, fpw );
		}
	}


	fclose(fpw);

	return 1;
}

/**
 * @brief           長田係数及び頂点座標を読込(Binary)
 * @param[in]       filename    出力ファイル名
 * @param[in]       mod         曲面の次数
 * @return          エラーコード　（0：エラー　1：正常）
 */
int	IOForStl::LoadNPT_Binary( char *filename, int mod ){

	int i, j, k;
	int nf = 0, nb = 0, nc = 0;
	float x, y, z;
	FILE *fp;
	int deg;

	// ***********************************************
	//     バイナリデータの読込                     //
	// ***********************************************
	if ((fp = fopen(filename, "rb" )) == NULL) {
		return -1;
	}

	// ++++++++++++++ 要素総数の読み込み ++++++++++++++++++++++++++++
	ReadInt(fp, &nf, 0 );

	// ++++++++++++++ 座標格納用TEMPARRAYの領域を確保 +++++++++++++++	
	float *bin_x, *bin_y, *bin_z;
	float *bin_cx, *bin_cy, *bin_cz;

	nb = nf*3;
	bin_x = new float[nb];
	bin_y = new float[nb];
	bin_z = new float[nb];

	if( mod == 2 ) deg = 3;
	if( mod == 3 ) deg = 7;

	nc = nf*deg;
	bin_cx = new float[nc];
	bin_cy = new float[nc];
	bin_cz = new float[nc];

	k = 0;
	// ***********************************************
	//              座標を読込                      //
	// ***********************************************
	for (i=0;i<nf;i++){

	   // +++++++座標を読み込み、ARRAYに格納 +++++++++
		for (j=0; j<3; j++){
			x = ReadFloat( fp );
			y = ReadFloat( fp );
			z = ReadFloat( fp );
			
			bin_x[k] = x;
			bin_y[k] = y;
			bin_z[k] = z;			
			k++;
		}
	}

	// ***********************************************
	//              係数を読込                      //
	// ***********************************************
	k = 0;
	for (i=0;i<nf;i++){

	   // +++++++座標を読み込み、ARRAYに格納 +++++++++
		for ( j=0; j<deg; j++ ){
			x = ReadFloat( fp );
			y = ReadFloat( fp );
			z = ReadFloat( fp );
			
			bin_cx[k] = x;
			bin_cy[k] = y;
			bin_cz[k] = z;			
			k++;
		}
	}
	fclose(fp);

	// +++++++++ TEM ARRAY 領域を開放 +++++++++++
	delete [] bin_x; 
	delete [] bin_y; 
	delete [] bin_z;

	delete [] bin_cx; 
	delete [] bin_cy; 
	delete [] bin_cz;

	FaceNum		= nf;
	VertexNum	= nb;

	return 1;

}

/**
 * @brief           Asciiか、Binaryかをチェックする
 * @param[in]       InFile  入力ファイル名
 * @return          判断結果　（0：Ascii 1：Binary ）
 */
int IOForStl::CheckStlType( char* InFile )
{

	// ++++++++++++ ﾃﾞｨﾌｫﾙﾄをASCIIにする +++++++++
	Stl_Type = 0;

	FILE*	pFILE;
	int 	c;

	if( NULL == InFile )
	{
		return -1;
	}

	pFILE = fopen( InFile, "rb" );
	if( NULL == pFILE )
	{
		return -1;
	}
	while( 0 == feof( pFILE ) )
	{
		c = fgetc( pFILE );
		if( EOF == c )
		{
			break;
		}
		else if( c <= 8 )
		{
			// アスキーコードの08H以下が検出されたので、おそらくバイナリファイルでしょう。
			fclose( pFILE );
			Stl_Type = 1;
			return Stl_Type;
		}
	}
	fclose( pFILE );


	return Stl_Type ;

}

/**
 * @brief           三角形パッチの総数を読込
 * @param[in]       InFile  入力ファイル名
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::GetStlNum_Binary (char* InFile )
{

	int nf = 0;
	char s[256];
	FILE *fp;

	// ***********************************************
	//     バイナリデータの読込                     //
	// ***********************************************
	if ((fp = fopen(InFile, "rb" )) == NULL) return -1;

	// ++++++++++++++ ヘッダーを読み込み ++++++++++++++++++++++++++++
	if (fread(s,sizeof(char),MAX_HEAD_SIZE,fp) != MAX_HEAD_SIZE ) return -1;

	// ++++++++++++++ 要素総数の読み込み ++++++++++++++++++++++++++++
	s[MAX_HEAD_SIZE] = '\0';
	ReadInt(fp, &nf, 0 );

	return nf;

}


/**
 * @brief           BINARY形式のINT型の読込
 * @param[in]       fp      入力ファイル名
 * @param[in]       n       データ
 * @param[in]       swap    SWAP
 * @return          エラーコード　（0：エラー　1：正常）
 */
int IOForStl::ReadInt(FILE *fp, int *n, int swap)
{
  unsigned char *cptr,tmp;

   if (fread(n, 4, 1, fp) != 1)
      return -1;
   if (swap) {
      cptr = (unsigned char *)n;
      tmp = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
   }
   return -1;
}

/**
 * @brief           BINARY形式のFloatT型の読込
 * @param[in]       fp  入力ファイル名
 * @return          Float型の数字結果
 */
float IOForStl::ReadFloat( FILE *fp )
{
  union 
    {
      int   int_value;
      float float_value;
    } value;
	value.int_value  =  fgetc(fp) & 0xFF;
	value.int_value |= (fgetc(fp) & 0xFF) << 0x08;
	value.int_value |= (fgetc(fp) & 0xFF) << 0x10;
	value.int_value |= (fgetc(fp) & 0xFF) << 0x18;
    return(value.float_value);

}
