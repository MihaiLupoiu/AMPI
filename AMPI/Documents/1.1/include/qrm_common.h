! /* ##############################################################################################
! **
! ** Copyright 2012 CNRS, INPT
! **  
! ** This file is part of qr_mumps.
! **  
! ** qr_mumps is free software: you can redistribute it and/or modify
! ** it under the terms of the GNU Lesser General Public License as 
! ** published by the Free Software Foundation, either version 3 of 
! ** the License, or (at your option) any later version.
! **  
! ** qr_mumps is distributed in the hope that it will be useful,
! ** but WITHOUT ANY WARRANTY; without even the implied warranty of
! ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! ** GNU Lesser General Public License for more details.
! **  
! ** You can find a copy of the GNU Lesser General Public License
! ** in the qr_mumps/doc directory.
! **
! ** ##############################################################################################*/
!  
!  
! /*##############################################################################################*/
! /** @file qrm_common.h
!  * Common header file
!  *
!  * $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!  * $Author: abuttari $
!  * $Version: 1.1$
!  * $Revision: 1980 $
!  *
!  */
! /* ############################################################################################## */

#define __QRM_PRNT_ERR(X)  if(qrm_eunit.gt.0) write(qrm_eunit,X)
#define __QRM_PRNT_MSG(X)  if(qrm_ounit.gt.0) write(qrm_ounit,X)
#define __QRM_PRNT_DBG(X)  if(qrm_dunit.gt.0) write(qrm_dunit,X)

#define __QRM_CHECK_RET(NAME,STR,GOTO) if(qrm_err_stack%nelem .gt. 0 ) then;\
call qrm_err_push(17,NAME,aed=STR);\
goto GOTO;\
endif

#define trace 1
