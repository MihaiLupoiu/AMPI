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
!> @file qrm_residual_orth.F90
!! This file contains a subroutine for computing the scaled norm of the product A'*r
!! (orthogonality or residual to the image of A)
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This routine computes the scaled norm of the product A'*r for multiple residuals

!> This routine computes the norm of the scaled product A'*r
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!
!! @param[in]           r the residuals. A 2D array of leading dimension qrm_mat%m. 
!! 
!! @param[out]        nrm the output norms norm
!! 
subroutine _qrm_residual_orth2d(qrm_mat, r, nrm)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data             :: r(:,:)
  _qrm_real             :: nrm(:)

  _qrm_data, allocatable :: tmp(:,:)
  _qrm_real, allocatable :: rnrm(:)
  _qrm_real              :: anrm

  ! error management
  integer                      :: err_act
  character(len=*), parameter  :: name='qrm_residual_orth'

  call qrm_err_act_save(err_act)

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  call qrm_aalloc(tmp, qrm_mat%n, size(r,2))
  call qrm_aalloc(rnrm, size(r,2))

  ! compute A'*r
  call qrm_matmul(qrm_mat, 't', _qrm_one, r, _qrm_zero, tmp)
   
  call qrm_vecnrm(r  , qrm_mat%m, '2', rnrm)
  call qrm_vecnrm(tmp, qrm_mat%n, '2',  nrm)
  call qrm_matnrm(qrm_mat, 'f', anrm)
  nrm = nrm/(rnrm*anrm)

  call qrm_adealloc(tmp)
  call qrm_adealloc(rnrm)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_residual_orth2d

!> @brief This routine computes the scaled norm of the product A'*r

!> This routine computes the norm of the scaled product A'*r
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!
!! @param[in]           r the residual. A 1D array of leading dimension qrm_mat%m. 
!! 
!! @param[out]        nrm norm
!! 
subroutine _qrm_residual_orth1d(qrm_mat, r, nrm)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data             :: r(:)
  _qrm_real             :: nrm

  _qrm_data, allocatable :: tmp(:)
  _qrm_real              :: rnrm
  _qrm_real              :: anrm

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_residual_orth'

  call qrm_err_act_save(err_act)

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  call qrm_aalloc(tmp, qrm_mat%n)

  ! compute A'*r
  call qrm_matmul(qrm_mat, 't', _qrm_one, r, _qrm_zero, tmp)
   
  call qrm_vecnrm(r,   qrm_mat%m, '2', rnrm)
  call qrm_vecnrm(tmp, qrm_mat%n, '2',  nrm)
  call qrm_matnrm(qrm_mat, 'f', anrm)

  nrm = nrm/(rnrm*anrm)
  call qrm_adealloc(tmp)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_residual_orth1d





!> @brief This routine computes the scaled norm of the product A'*r

!> This routine computes the norm of the scaled product A'*r
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain the residuals
!! @param[out]          x the solutions
!! 
!! @param[out]        nrm norms
!! 
subroutine _qrm_residual_and_orth2d(qrm_mat, b, x, nrm)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data             :: b(:,:), x(:,:)
  _qrm_real             :: nrm(:)

  _qrm_data, allocatable :: tmp(:,:)
  _qrm_real, allocatable :: rnrm(:)

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_residual_orth'

  call qrm_err_act_save(err_act)

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  ! compute the residual
  call qrm_matmul(qrm_mat, 'n', -_qrm_one, x, _qrm_one, b)

  call qrm_aalloc(tmp , qrm_mat%n, size(x,2))
  call qrm_aalloc(rnrm, size(x,2))

  ! compute A'*r
  call qrm_matmul(qrm_mat, 't', _qrm_one, b, _qrm_zero, tmp)
  
  call qrm_vecnrm(b,   qrm_mat%m, '2', rnrm)
  call qrm_vecnrm(tmp, qrm_mat%n, '2',  nrm)
  nrm = nrm/rnrm

  call qrm_adealloc(tmp)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_residual_and_orth2d


!> @brief This routine computes the scaled norm of the product A'*r

!> This routine computes the norm of the scaled product A'*r
!! for a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain the residuals
!! @param[out]          x the solution
!! 
!! @param[out]        nrm norm
!! 
subroutine _qrm_residual_and_orth1d(qrm_mat, b, x, nrm)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_utils_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data             :: b(:), x(:)
  _qrm_real             :: nrm

  _qrm_data, allocatable :: tmp(:)
  _qrm_real              :: rnrm

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_residual_orth'

  call qrm_err_act_save(err_act)

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  ! compute the residual
  call qrm_matmul(qrm_mat, 'n', -_qrm_one, x, _qrm_one, b)

  call qrm_aalloc(tmp, qrm_mat%n)

  ! compute A'*r
  call qrm_matmul(qrm_mat, 't', _qrm_one, b, _qrm_zero, tmp)
  
  call qrm_vecnrm(b,   qrm_mat%m, '2', rnrm)
  call qrm_vecnrm(tmp, qrm_mat%n, '2',  nrm)
  nrm = nrm/rnrm

  call qrm_adealloc(tmp)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_residual_and_orth1d
