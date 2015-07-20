#ifndef _NPT_UTIL_H_
#define _NPT_UTIL_H_

/*
 * NpatchLib - Nagata Patch Library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   npt_Utility.h
 * @brief  Utility Class Header
 * @author aics    
 */

#include <string>
#include "npt_Version.h"


class npt_Utility {

public:


  /** コンストラクタ */
  npt_Utility();

  /** デストラクタ */
  ~npt_Utility();

  
  /** バージョンを出力する */
  static std::string getVersionInfo()
  {
    std::string str(NPT_VERSION_NO);
    return str;
  }


};


#endif // _NPT_UTIL_H_
