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
!> @file qrm_apply2d.F90
!! This file contains a routine that applies Q or Q' to a set of vectors
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function applies Q or Q^T to a set of vectors
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether Q or Q^T will be applied. Only the
!!                      first character is important.
!! 
!! @param[in,out] b     a 2d array containing the vectors to which Q will 
!!                      be applied. 
!!
subroutine _qrm_apply2d(qrm_mat, transp, b)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_solve_mod, protect => _qrm_apply2d
#if defined(_OPENMP)
  use omp_lib
#endif
  implicit none

  type(_qrm_spmat_type), intent(in) :: qrm_mat
  _qrm_data, intent(inout)          :: b(:,:)
  character(len=*), intent(in)       :: transp

  integer :: nb, rhs_nthreads, nrhs, i, keeph
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_apply'
  
  call qrm_err_act_save(err_act)
  
  ! immediately check if the facto was done. Otherwise push an error and return
  if(.not. qrm_mat%fdata%ok) then
     call qrm_err_push(14, 'qrm_apply')
     goto 9999
  end if

  call qrm_get(qrm_mat, 'qrm_keeph', keeph)
  if(keeph .eq. qrm_no_) then
     call qrm_err_push(30, 'qrm_apply')
     goto 9999
  end if
  
  ! blocking to deal with multiple rhs
  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  call qrm_get(qrm_mat, 'qrm_rhsnthreads', rhs_nthreads)
  nrhs = size(b,2)
  if(nb.le.0) nb = nrhs

#if defined(_OPENMP)
  call omp_set_nested(.true.)
#endif
  if(qrm_str_tolower(transp(1:1)) .eq. 't') then

     !$omp parallel do num_threads(rhs_nthreads) private(i)
     do i=1, nrhs, nb
        call qrm_apply_qt(qrm_mat, b(:,i:min(nrhs,i+nb-1)))
     end do
     !$omp end parallel do
     __QRM_CHECK_RET(name,'qrm_apply_qt',9999)
     
  else 
     !$omp parallel do num_threads(rhs_nthreads) private(i)
     do i=1, nrhs, nb
        call qrm_apply_q(qrm_mat, b(:,i:min(nrhs,i+nb-1)))
     end do
     !$omp end parallel do
     __QRM_CHECK_RET(name,'qrm_apply_q',9999)
  end if

#if defined(_OPENMP)
  call omp_set_nested(.false.)
#endif

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_apply2d


!> @brief This function applies Q or Q^T to a single vector
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether Q or Q^T will be applied. Only the
!!                      first character is important.
!! 
!! @param[in,out] b     a 1d array containing the vector to which Q will 
!!                      be applied. 
!!
subroutine _qrm_apply1d(qrm_mat, transp, b)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_string_mod
  use _qrm_solve_mod, protect => _qrm_apply1d
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type), intent(in) :: qrm_mat
  _qrm_data, intent(inout)          :: b(:)
  character(len=*), intent(in)       :: transp
  
  _qrm_data, pointer :: pnt(:,:)
  integer :: n, keeph
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_apply1d'
  
  call qrm_err_act_save(err_act)

  ! immediately check if the facto was done. Otherwise push an error and return
  if( qrm_mat%fdata%done .eq. 0) then
     call qrm_err_push(14, 'qrm_apply1d')
     goto 9999
  end if

  call qrm_get(qrm_mat, 'qrm_keeph', keeph)
  if(keeph .eq. qrm_no_) then
     call qrm_err_push(30, 'qrm_apply')
     goto 9999
  end if
  
  
  n = size(b,1)

  call _qrm_remap_pnt(b, pnt, n)

  if(qrm_str_tolower(transp(1:1)) .eq. 't') then
     call qrm_apply_qt(qrm_mat, pnt)
     __QRM_CHECK_RET(name,'qrm_apply_qt',9999)
  else 
     call qrm_apply_q(qrm_mat, pnt)
     __QRM_CHECK_RET(name,'qrm_apply_q',9999)
  end if

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
   
end subroutine _qrm_apply1d
