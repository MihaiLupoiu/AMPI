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
/** @file qrm_colamd_wrap.c
 * This filetains a wrapper for colamd routines
 *
 * $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
 * $Author: abuttari $
 * $Version: 1.1$
 * $Revision: 1980 $
 *
 **/
/*##############################################################################################*/


#if defined(have_colamd)
#include "colamd.h"
#endif
#include <stdio.h>

void qrm_colamd(int n_row, int n_col, int Alen, int *A, int *p, int *err){
#if defined(have_colamd)
  int stats [COLAMD_STATS] ;	
  double knobs [COLAMD_KNOBS] ;
  

  colamd_set_defaults (knobs) ;

  knobs[0]=-1;
  knobs[1]=10;
  knobs[2]=1;

  *err = colamd (n_row, n_col, Alen, A, p,
		 knobs, stats) ;
  if(!*err)
    printf("Error in COLAMD! %d\n",stats[3]);

  return;
#endif
}


void qrm_colamd_recommended(int *alen, int nnz, int n_row, int n_col){

#if defined(have_colamd)
  *alen = colamd_recommended(nnz, n_row, n_col);
  if(!alen)
    printf("Error in COLAMD_RECOMMENDED!\n");
#endif
}





