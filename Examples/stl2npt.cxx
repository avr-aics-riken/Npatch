/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   main.cxx
 * @brief  Test program
 * @author aics
 */

#include <stdio.h>
#include <stdlib.h>

#include "npt_Data.h"

int main(int argc, char *argv[]) {

	int degree;
	char InFileName[256];
	int ExType;

	strcpy( InFileName,  argv[1] );

	degree = atoi( argv[2] );
	ExType = atoi( argv[3] );

	std::cout<<" " << std::endl;
	std::cout<<" " << std::endl;
	std::cout<< "Stl File = " << InFileName << std::endl;
	std::cout<<"degree = " << degree << std::endl;
	std::cout<<"Extype = " << ExType << std::endl;

   	CData *pData= NULL;
   	pData = new CData;

   	pData->MainProcess( InFileName, degree );
  
	std::cout<<" " << std::endl;
	std::cout<<"Face Size = " << pData->m_ArrayFace.GetSize() << std::endl;
	std::cout<<" " << std::endl;

	std::cout<< "++++++++ Show you Nagata Patch coef in flowing +++++" << std::endl;

	float x, y, z;
	for( int i=0; i<5; i++){
		CFace *pFace = pData->m_ArrayFace[i];

		std::cout<<" " << std::endl;
		x = pFace->m_Coef[0].x;
		y = pFace->m_Coef[0].y;
		z = pFace->m_Coef[0].z;
		std::cout<<"Face["<< i << "][Coef 1 ]" << x << y << z << std::endl;

		x = pFace->m_Coef[1].x;
		y = pFace->m_Coef[1].y;
		z = pFace->m_Coef[1].z;
		std::cout<<"Face["<< i << "][Coef 2 ]" << x << y << z << std::endl;

		x = pFace->m_Coef[2].x;
		y = pFace->m_Coef[2].y;
		z = pFace->m_Coef[2].z;
		std::cout<<"Face["<< i << "][Coef 3 ]" << x << y << z << std::endl;

		if( degree == 3 ){	

		x = pFace->m_Coef[3].x;
		y = pFace->m_Coef[3].y;
		z = pFace->m_Coef[3].z;
		std::cout<<"Face["<< i << "][Coef 4 ]" << x << y << z << std::endl;

		x = pFace->m_Coef[4].x;
		y = pFace->m_Coef[4].y;
		z = pFace->m_Coef[4].z;
		std::cout<<"Face["<< i << "][Coef 5 ]" << x << y << z << std::endl;

		x = pFace->m_Coef[5].x;
		y = pFace->m_Coef[5].y;
		z = pFace->m_Coef[5].z;
		std::cout<<"Face["<< i << "][Coef 6 ]" << x << y << z << std::endl;

		x = pFace->m_Coef[6].x;
		y = pFace->m_Coef[6].y;
		z = pFace->m_Coef[6].z;
		std::cout<<"Face["<< i << "][Coef 7 ]" << x << y << z << std::endl;

		}

	}

	pData->SaveNPT( ExType );

   	delete pData;

  	return 1;
}

