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
!> @file qrm_analysis_mod.F90
!! This file contains a subroutine for computing the least-squares solution of a problem
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This routine computes the least-squares solution of a problem

!> This routine computed the least-squares solution of an overdetermined system
!! Ax=b with multiple RHSs.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output the original data will be unchanged and
!!                        the result of the analysis and factorization phases will
!!                        be stored in the adata and fdata fields, respectively.
!! @param[in,out]       b the RHSs. A 2D array of leading dimension qrm_mat%m. On output
!!                        it will contain Q'*b
!! @param[out]          x the solution, i.e., R\Q'*b
!! 
subroutine _qrm_least_squares2d(qrm_mat, b, x)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use _qrm_analysis_mod
  use _qrm_factorization_mod
  use _qrm_solve_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data :: b(:,:), x(:,:)

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_least_squares'

  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Entering the least-squares driver")')

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  if(qrm_mat%m .lt. qrm_mat%n) then
     call qrm_err_push(30, name,ied=(/qrm_mat%m,qrm_mat%n,0,0,0/))
     goto 9999
  end if

  ! analysis
  call _qrm_analyse(qrm_mat, 'n')
  __QRM_CHECK_RET(name,'qrm_analyse',9999)
  ! factorization
  call _qrm_factorize(qrm_mat, 'n')
  __QRM_CHECK_RET(name,'qrm_factorize',9999)

  call _qrm_apply2d(qrm_mat, 't', b)
  __QRM_CHECK_RET(name,'qrm_apply',9999)
  call _qrm_solve2d(qrm_mat, 'n', b, x)
  __QRM_CHECK_RET(name,'qrm_solve',9999)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_least_squares2d




!> @brief This routine computes the least-squares solution of a problem

!> This routine computed the least-squares solution of an overdetermined system
!! Ax=b with a single RHS.
!!
!! @param[in]     qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output the original data will be unchanged and
!!                        the result of the analysis and factorization phases will
!!                        be stored in the adata and fdata fields, respectively.
!! @param[in,out]       b the RHSs. A 1D array of leading dimension qrm_mat%m. On output
!!                        it will contain Q'*b
!! @param[out]          x the solution, i.e., R\Q'*b
!! 
subroutine _qrm_least_squares1d(qrm_mat, b, x)
  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_error_mod
  use _qrm_analysis_mod
  use _qrm_factorization_mod
  use _qrm_solve_mod
  implicit none

  type(_qrm_spmat_type) :: qrm_mat
  _qrm_data :: b(:), x(:)

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_least_squares'

  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Entering the least-squares driver")')

  call _qrm_check_spmat(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  if(qrm_mat%m .lt. qrm_mat%n) then
     call qrm_err_push(30, name,ied=(/qrm_mat%m,qrm_mat%n,0,0,0/))
     goto 9999
  end if

  ! analysis
  call _qrm_analyse(qrm_mat, 'n')
  __QRM_CHECK_RET(name,'qrm_analyse',9999)

  ! factorization
  call _qrm_factorize(qrm_mat, 'n')
  __QRM_CHECK_RET(name,'qrm_factorize',9999)

  call _qrm_apply1d(qrm_mat, 't', b)
  __QRM_CHECK_RET(name,'qrm_apply',9999)
  call _qrm_solve1d(qrm_mat, 'n', b, x)
  __QRM_CHECK_RET(name,'qrm_solve',9999)

  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return
  
end subroutine _qrm_least_squares1d
