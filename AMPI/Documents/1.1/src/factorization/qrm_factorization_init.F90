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
!> @file qrm_factorization_init.F90
!! This file contains a subroutine that initializes the factorization
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h" 
!> @brief This subroutine initializes the data structures needed for
!! the actual factorization.

!> The main task achieved by this routine is the creation of what we
!! call (in mumps terminology) the arrowheads. Basically it builds a
!! list of @link _qrm_fdata_mod::_qrm_front_type @endlink elements
!! (each one corresponding to one front) and associates to each of
!! them the related coefficients of the original matrix in CSR
!! format. This coefficients will be assembled into the front matrix
!! at the moment of its activation (this is done by the @link
!! _qrm_init_front @endlink routine).
!!
!! param[in] qrm_mat the usual blob associated to the problem
!!
subroutine _qrm_factorization_init(qrm_mat)
 
  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_error_mod
  implicit none
  
  type(_qrm_spmat_type), target :: qrm_mat
  
  integer :: f, nrowsf, nvalsf, i, j, roff, r, c, lrow, itmp, k, m, n
  integer, allocatable :: rmap(:), rcnt(:), vcnt(:), row_to_frow(:)
  type(_qrm_front_type), pointer :: front
  _qrm_data :: vtmp
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_analyse'

  call qrm_err_act_save(err_act)
  
  call _qrm_fdata_destroy(qrm_mat%fdata)
  __QRM_CHECK_RET(name,'qrm_fdata_destroy',9999)
  allocate(qrm_mat%fdata%front_list(qrm_mat%adata%nnodes))

  ! rmap contains a mapping of the rows onto fronts
  call qrm_aalloc(rmap, qrm_mat%m)
  call qrm_aalloc(row_to_frow, qrm_mat%m)
  call qrm_aalloc(rcnt, qrm_mat%m)
  call qrm_aalloc(vcnt, qrm_mat%adata%nnodes)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  rcnt = 0
  vcnt = 0
  m = qrm_mat%m
  n = qrm_mat%n

  qrm_mat%fdata%nfronts = qrm_mat%adata%nnodes
#if defined(_OPENMP)
  call omp_init_lock(qrm_mat%fdata%lock)
#endif

  ! build the row mapping
  ! TODO: stair(i) may be set to be the first row of front i instead
  ! of the last
  roff = 1
  do f = 1, qrm_mat%adata%nnodes
     ! for all the rows belonging to this front
     do i = roff, qrm_mat%adata%stair(f)
        rmap(qrm_mat%adata%rperm(i)) = f
        row_to_frow(qrm_mat%adata%rperm(i)) = i-roff+1
     end do
     roff = qrm_mat%adata%stair(f)+1
  end do

  ! first pass to count
  do k=1, qrm_mat%nz
     j = qrm_mat%jcn(k)
     i = qrm_mat%irn(k)
     if((j.le.0) .or. (j.gt.n) .or. (i.le.0) .or. (i.gt.m) ) cycle
     f = rmap(i)
     ! count the number of values per row
     rcnt(i) = rcnt(i)+1
     ! count the number of nonzeros from the original matrix in front f
     vcnt(f) = vcnt(f)+1
  end do

  ! pass to allocate
  roff = 0
  qrm_mat%fdata%front_list(:)%status=0
  do f = 1, qrm_mat%adata%nnodes
     front      => qrm_mat%fdata%front_list(f)
     nrowsf = qrm_mat%adata%stair(f) - roff
     front%anrows = nrowsf
     front%num = f
     do i=qrm_mat%adata%childptr(f), qrm_mat%adata%childptr(f+1)-1
        c = qrm_mat%adata%child(i)
        if(qrm_mat%adata%small(c) .eq. 0) front%status = front%status-1
     end do
#if defined(_OPENMP)
     call omp_init_lock(front%lock)
#endif
     call qrm_aalloc(front%aiptr, max(nrowsf+1,2))
     call qrm_aalloc(front%ajcn,  vcnt(f))
     call qrm_aalloc(front%aval,  vcnt(f))
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)

     front%aiptr(1:2) = 1
     do i = 1, nrowsf-1
        front%aiptr(i+2) = front%aiptr(i+1)+rcnt(qrm_mat%adata%rperm(roff+i))
     end do

     roff = qrm_mat%adata%stair(f)
  end do



  ! pass to fill
  do k=1, qrm_mat%nz
     j = qrm_mat%jcn(k)
     i = qrm_mat%irn(k)
     if((j.le.0) .or. (j.gt.n) .or. (i.le.0) .or. (i.gt.m) ) cycle
     r     = i
     c     = j
     f     = rmap(r)
     lrow  = row_to_frow(r)
     front => qrm_mat%fdata%front_list(f)
     front%ajcn(front%aiptr(lrow+1)) = c
     front%aval(front%aiptr(lrow+1)) = qrm_mat%val(k)
     front%aiptr(lrow+1) = front%aiptr(lrow+1)+1
  end do

  qrm_mat%fdata%done = 0

  if(qrm_mat%adata%ncsing .gt. 0) then
     qrm_mat%fdata%done = 1
     ! the first front may contain all the singletons. For this
     ! reason, each diagonal element must be the first of its
     ! row. This is done to make the solve easier and speed it -up.
     ! we just sweep over all the elements of the front and swap the
     ! diagonal element with the first along its row. It may be
     ! optimized TODO
     front => qrm_mat%fdata%front_list(1)
     do i=1, front%anrows
        do j=front%aiptr(i), front%aiptr(i+1)-1
           c = front%ajcn(j)
           if(qrm_mat%adata%cperm(i) .eq. c) then
              ! this is a diagonal element; swap
              vtmp = front%aval(front%aiptr(i))
              front%aval(front%aiptr(i)) = front%aval(j)
              front%aval(j) = vtmp
              itmp = front%ajcn(front%aiptr(i))
              front%ajcn(front%aiptr(i)) = c
              front%ajcn(j) = itmp
           end if
        end do
     end do
  end if
  


  call qrm_adealloc(rmap)
  call qrm_adealloc(row_to_frow)
  call qrm_adealloc(rcnt)
  call qrm_adealloc(vcnt)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
 
  return
 
end subroutine _qrm_factorization_init
