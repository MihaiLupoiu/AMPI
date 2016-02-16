!! ##############################################################################################
!!
!! Copyright 2012 CNRS, INPT
!!  
!! This file is part of qr_mumps.
!!  
!! qr_mumps is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as 
!! published by the Free Software Foundation, either version 3 of 
!! the License, or (at your option) any later version.
!!  
!! qr_mumps is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!  
!! You can find a copy of the GNU Lesser General Public License
!! in the qr_mumps/doc directory.
!!
!! ##############################################################################################


!! ##############################################################################################
!> @file qrm_matnrm.F90
!! this file contains a routine that computes the norm of a sparse matrix
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the matrix norm. The return value is
!! a real scalar
!!
!! @param[in] qrm_mat  the inpur A matrix
!! @param[in]   ntype  the norm type. It can be one of these
!!                     1  : 1-norm
!!                     i  : inf-norm
!!                     f  : Frobenius-norm
!! @param[out]    nrm  the output norm
!!
subroutine _qrm_matnrm(qrm_mat, ntype, nrm)

  use _qrm_spmat_mod
  use qrm_string_mod
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), intent(in) :: qrm_mat
  _qrm_real                         :: nrm
  character                         :: ntype

  _qrm_real, allocatable  :: tmp(:)
  integer                 :: r, c, i
  _qrm_real :: _rxnrm2

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_matnrm'

  call qrm_err_act_save(err_act)

  if(qrm_str_tolower(ntype) .eq. 'i') then
     call qrm_aalloc(tmp, qrm_mat%m)
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)
     tmp = _qrm_zero
     do i=1, qrm_mat%nz
        r = qrm_mat%irn(i)
        tmp(r) = tmp(r)+abs(qrm_mat%val(i))
     end do
     nrm = maxval(tmp)
  else if(qrm_str_tolower(ntype) .eq. '1') then
     call qrm_aalloc(tmp, qrm_mat%n)
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)
     tmp = _qrm_zero
     do i=1, qrm_mat%nz
        c = qrm_mat%jcn(i)
        tmp(c) = tmp(c)+abs(qrm_mat%val(i))
     end do
     nrm = maxval(tmp)
  else if(qrm_str_tolower(ntype) .eq. 'f') then
     nrm = _rxnrm2(qrm_mat%nz, qrm_mat%val, 1)
  else
     call qrm_err_push(15, 'qrm_matnrm')
     goto 9999
  end if

  call qrm_adealloc(tmp)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return

end subroutine _qrm_matnrm
