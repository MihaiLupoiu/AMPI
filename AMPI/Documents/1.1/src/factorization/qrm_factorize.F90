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
!> @file qrm_factorize.F90
!! This file contains the main factorization driver
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine is the main factorization driver
!! 
!! @param[in,out] qrm_mat the problem containing the matrix to be factorized.
!!
!! @param[in] transp whether to factorize the input matrix or its
!!            transpose. Accepted values are 't' or 'n'
!!
subroutine _qrm_factorize(qrm_mat, transp)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_factorization_mod, protect => _qrm_factorize
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_string_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  character, optional, intent(in) :: transp

  integer                         :: i, totnnz, info, h
  real(kind(1.d0))                :: t1, t2
  type(_qrm_front_type), pointer  :: front
  character                       :: itransp
  integer, pointer                :: tmp(:)
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_factorize'

  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Entering the factorization driver")')

  ! immediately check if the analysis was done. Otherwise push an error and return
  if(.not. qrm_mat%adata%ok) then
     call qrm_err_push(13, 'qrm_factorize')
     goto 9999
  end if

  call _qrm_check_spmat(qrm_mat, qrm_factorize_)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  

  if(present(transp)) then
     itransp = transp
  else
     itransp = 'n'
  end if

  ! in case transp=='t' switch temporarily the row and column indices as well as m and n
  if(qrm_str_tolower(itransp) .eq. 't') then
     tmp => qrm_mat%irn
     qrm_mat%irn => qrm_mat%jcn
     qrm_mat%jcn => tmp
     i = qrm_mat%m
     qrm_mat%m = qrm_mat%n
     qrm_mat%n = i
#if defined(zprec) || defined(cprec)
     qrm_mat%val = conjg(qrm_mat%val)
#endif
  end if

  ! initialize the data for the facto
  call _qrm_factorization_init(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_factorization_init',9998)

  !$ call omp_set_num_threads(1)

  call _qrm_factorization_core(qrm_mat)
  __QRM_CHECK_RET(name,'qrm_factorization_core',9998)

9998 continue
  
  qrm_mat%gstats(qrm_nnz_r_) = 0
  qrm_mat%gstats(qrm_nnz_h_) = 0
  do i=1, qrm_mat%adata%nnodes
     qrm_mat%gstats(qrm_nnz_r_) =   qrm_mat%gstats(qrm_nnz_r_) + &
          & qrm_mat%fdata%front_list(i)%rsize
  end do

  if(qrm_mat%icntl(qrm_keeph_) .eq. qrm_yes_) then
     do i=1, qrm_mat%adata%nnodes
     qrm_mat%gstats(qrm_nnz_h_) =   qrm_mat%gstats(qrm_nnz_h_) + &
          & qrm_mat%fdata%front_list(i)%hsize
     end do
  end if
     
     

  ! in case transp=='t' switch temporarily the row and column indices as well as m and n
  if(qrm_str_tolower(itransp) .eq. 't') then
     tmp => qrm_mat%irn
     qrm_mat%irn => qrm_mat%jcn
     qrm_mat%jcn => tmp
     i = qrm_mat%m
     qrm_mat%m = qrm_mat%n
     qrm_mat%n = i
#if defined(zprec) || defined(cprec)
     qrm_mat%val = conjg(qrm_mat%val)
#endif
  end if

  call qrm_err_get(info)
  if(info .gt. 0) goto 9999

  ! the factorization was succesfully performed
  qrm_mat%fdata%ok = .true.

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_factorize

