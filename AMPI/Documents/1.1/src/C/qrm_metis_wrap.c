/* ##############################################################################################
**
** Copyright 2012 CNRS, INPT
**  
** This file is part of qr_mumps.
**  
** qr_mumps is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as 
** published by the Free Software Foundation, either version 3 of 
** the License, or (at your option) any later version.
**  
** qr_mumps is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**  
** You can find a copy of the GNU Lesser General Public License
** in the qr_mumps/doc directory.
**
** ##############################################################################################*/


/*##############################################################################################*/
/** @file qrm_metis_wrap.c
 * FIXME: add comments
 *
 * $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
 * $Author: abuttari $
 * $Version: 1.1$
 * $Revision: 1980 $
 *
 **/
/*##############################################################################################*/


#if defined(have_metis)
#include "metis.h"
#endif
#include <stdio.h>

void qrm_metis(int *n, int *iptr, int *jcn, int *cperm, int *iperm){
#if defined(have_metis)


#if defined(METIS_VER_MAJOR)
  int options [METIS_NOPTIONS] ;	

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING]=1;

  METIS_NodeND(n, iptr, jcn, NULL, options, cperm, iperm);

#else
  int options [8] ;	
  int numflag;

  options[0]=0;

  numflag=1;
  METIS_NodeND(n, iptr, jcn, &numflag, &options[0], cperm, iperm);
#endif

  return;
#endif
}







