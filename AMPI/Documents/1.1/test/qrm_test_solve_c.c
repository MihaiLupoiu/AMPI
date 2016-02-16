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


/**##############################################################################################
* @file qrm_coverage.F90
* This file contains coverage tests for the C interface
*
* $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
* $Author: abuttari $
* $Version: 1.1$
* $Revision: 1980 $
*
##############################################################################################*/

#include "_qrm_mumps.h"
#include <string.h>
#include <stdio.h>

_qrm_real_c _qrm_test_solve_c(struct _qrm_spmat_type_c *qrm_spmat, _qrm_data_c *b, 
                  _qrm_data_c *x, _qrm_data_c *r){

  _qrm_real_c err1, err2;
  int i;

  for(i=0; i<qrm_spmat->m; i++)
    r[i] = b[i];

  if(qrm_spmat->m >= qrm_spmat->n) {
    _qrm_least_squares_c(qrm_spmat, b, x, 1);
  } else {
    _qrm_min_norm_c(qrm_spmat, b, x, 1);
  }

  _qrm_residual_norm_c(qrm_spmat, r, x, 1, &err1);
  _qrm_residual_orth_c(qrm_spmat, r, 1, &err2);

  if(err1 < err2)
    err2 = err1;

  return err2;

}
