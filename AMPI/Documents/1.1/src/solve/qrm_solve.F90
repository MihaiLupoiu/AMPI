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
!> @file qrm_solve2d.F90
!! This file contains a routine that solves for R or R' against multiple vectors
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


! #include "prec.h"
#include "qrm_common.h"


!> @brief This function solves for R or R' against multiple vectors
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether R or R^T will be solved
!!                      for. Only the first character is important.
!! 
!! @param[in]     b     a 2d array containing the RHS vectors
!!
!! @param[out]    x     a 2d array containing the solution vectors
!!
subroutine _qrm_solve2d(qrm_mat, transp, b, x)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use qrm_string_mod
  use _qrm_utils_mod
  use _qrm_solve_mod, protect => _qrm_solve2d
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(inout)  :: b(:,:)
  _qrm_data, intent(out)    :: x(:,:)
  character(len=*) :: transp

  integer :: i, nb, nrhs, rhs_nthreads
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_solve'
  
  call qrm_err_act_save(err_act)

  if(.not. qrm_mat%fdata%ok) then
     call qrm_err_push(14, 'qrm_solve')
     return
  end if

  ! FIXME x should not be initialized to 0
  x = _qrm_zero
  
  ! blocking to deal with multiple rhs
  call qrm_get(qrm_mat, 'qrm_rhsnb', nb)
  call qrm_get(qrm_mat, 'qrm_rhsnthreads', rhs_nthreads)
  nrhs = size(b,2)
  if(nb.le.0) nb = nrhs

  !$ call omp_set_nested(.true.)
  if(qrm_str_tolower(transp(1:1)) .eq. 't') then
     !$omp parallel do num_threads(rhs_nthreads) private(i)
     do i=1, nrhs, nb
        call qrm_solve_rt(qrm_mat, b(:,i:min(nrhs,i+nb-1)), x(:,i:min(nrhs,i+nb-1)))
     end do
     !$omp end parallel do
     __QRM_CHECK_RET(name,'qrm_solve_rt',9999)
  else 
     !$omp parallel do num_threads(rhs_nthreads) private(i)
     do i=1, nrhs, nb
        call qrm_solve_r(qrm_mat, b(:,i:min(nrhs,i+nb-1)), x(:,i:min(nrhs,i+nb-1)))
     end do
     !$omp end parallel do
     __QRM_CHECK_RET(name,'qrm_solve_r',9999)
  end if

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return

end subroutine _qrm_solve2d



!> @brief This function solves for R or R' against a single vector
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in] transp    a string saying whether R or R^T will be solved
!!                      for. Only the first character is important.
!! 
!! @param[in]     b     a 1d array containing the RHS vector 
!!
!! @param[out]    x     a 1d array containing the solution vector
!!
subroutine _qrm_solve1d(qrm_mat, transp, b, x)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_fdata_mod
  use qrm_string_mod
  use _qrm_utils_mod
  use _qrm_solve_mod, protect => _qrm_solve1d
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(in)  :: b(:)
  _qrm_data, intent(out) :: x(:)
  character(len=*) :: transp

  _qrm_data, pointer :: pntb(:,:), pntx(:,:)
  integer :: n
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_solve1d'
  
  call qrm_err_act_save(err_act)

  if( qrm_mat%fdata%done .eq. 0) then
     call qrm_err_push(14, 'qrm_solve')
     goto 9999
  end if

  ! FIXME x should not be initialized to 0
  x = _qrm_zero
  
  n = size(b,1)
  call _qrm_remap_pnt(b, pntb, n)
  n = size(x,1)
  call _qrm_remap_pnt(x, pntx, n)

  if(qrm_str_tolower(transp(1:1)) .eq. 't') then
     call qrm_solve_rt(qrm_mat, pntb, pntx)
     __QRM_CHECK_RET(name,'qrm_solve_rt',9999)
  else 
     call qrm_solve_r(qrm_mat, pntb, pntx)
     __QRM_CHECK_RET(name,'qrm_solve_r',9999)
  end if

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return

end subroutine _qrm_solve1d

