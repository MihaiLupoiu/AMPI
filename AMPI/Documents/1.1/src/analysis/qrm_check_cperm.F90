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
!> @file qrm_check_cperm.F90
!! This file contains the routine that check whether a provided permutation is good
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine simply checks whether a column permutation provided
!! by the user makes sens.
!! 
!! @param[in] cperm  the permutation to be checked
!! @param[in] n      the size of the permutation verctor, i.e. the number of
!!                   columns in the matrix
!!                   
subroutine qrm_check_cperm(cperm, n)
  use qrm_error_mod
  use qrm_mem_mod
  implicit none

  integer :: cperm(:)
  integer :: n

  integer, allocatable :: tmp(:)
  integer :: i
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_check_cperm'

  call qrm_err_act_save(err_act)

  call qrm_aalloc(tmp, n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  tmp = 0

  do i=1, n
     if (cperm(i) .gt. n .or. cperm(i) .lt. 1) then
        call qrm_err_push(8, name)
        call qrm_adealloc(tmp)
        goto 9999
     end if
     if(tmp(cperm(i)) .gt. 0) then
        call qrm_err_push(8, name)
        call qrm_adealloc(tmp)
        goto 9999
     else
        tmp(cperm(i)) = 1
     end if
  end do

  call qrm_adealloc(tmp)
     
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
  
end subroutine qrm_check_cperm
