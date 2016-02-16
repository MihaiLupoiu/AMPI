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
!> @file qrm_activate_front.F90
!! This file contains the routines for activation and cleaning of a front
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine activates a front

!> The activation of a front implies the treatment of all the small
!! subtrees rooted at the front's children. Small trees are
!! recognizable since their root is flagged with a 1 in the
!! adata%small array
!!
!! @param[in,out] qrm_mat the whole problem. this obviously contains
!!                        the fornt to be activate
!!
!! @param[in] fnum the id of the front to be activated
!!
!! @param[in,out] flops a counter for flops. this is updated wrt to
!!                      the flops performed so far
!!
subroutine _qrm_activate_front(qrm_mat, fnum, flops)

  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_sort_mod
  use _qrm_factorization_mod, protect => _qrm_activate_front
  implicit none
  
  type(_qrm_spmat_type), target :: qrm_mat
  integer                        :: fnum
  real(kind(1.d0))    :: flops
    
  integer :: n_rows_orig, i, j, row, col, c, child, roff, crows, f
  integer :: m, n, np, rsize, hsize, mr, mc, mn, k, p, ne
  integer :: frow, fcol, nb, b, father, thn
  type(_qrm_front_type), pointer :: front
  type(qrm_adata_type), pointer :: adata

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_activate_front'

  call qrm_err_act_save(err_act)


  thn = 0
  !$ thn=omp_get_thread_num()

  front  => null()

  ! to make things easier
  front      => qrm_mat%fdata%front_list(fnum)
  father     =  qrm_mat%adata%parent(fnum)
  adata      => qrm_mat%adata

  ! sweep over the children to check whether there are small subtrees
  ! that have to be treated
  do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
     c = adata%child(p)
     if(adata%small(c) .eq. 1) call _qrm_do_subtree(qrm_mat, c, flops)
     __QRM_CHECK_RET(name,'qrm_do_subtree',9999)
  end do
  
  call _qrm_init_front(qrm_mat, fnum, .true.)
  __QRM_CHECK_RET(name,'qrm_init_front',9999)
  
  front%owner = thn

  ! clean
  do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
     c = adata%child(p)
     if(adata%small(c) .eq. 1) then
        call _qrm_clean_front(qrm_mat, c)
        __QRM_CHECK_RET(name,'qrm_clean_front',9999)
        
        !$ call omp_set_lock(qrm_mat%fdata%lock)
        qrm_mat%fdata%done = qrm_mat%fdata%done+1
        !$ call omp_unset_lock(qrm_mat%fdata%lock)
        
     end if
     
  end do
  
  call qrm_err_act_restore(err_act)
  return
  
9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  
  return
  
end subroutine _qrm_activate_front






!> @brief This routine performs the cleaning of a front.

!> Cleaning a front means saving the parts corresponding tot he R and
!! Q factors, and then freeing all the memory that is not needed
!! anymore
!!
!! @param[in,out] qrm_mat the whole problem. this obviously contains
!!                        the fornt to be cleaned
!!
!! @param[in] fnum the id of the front to be activated
!!
subroutine _qrm_clean_front(qrm_mat, fnum)

  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_factorization_mod, protect => _qrm_clean_front
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  integer                        :: fnum
  

  type(_qrm_front_type), pointer :: front
  integer :: i, j

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_clean_front'

  call qrm_err_act_save(err_act)

  front  => qrm_mat%fdata%front_list(fnum)
  ! if(front%num.eq.3) then
     ! open(4,file="ffr.m",action="write")
     ! ! open(5,file="ffi.m",action="write")
     ! write(4,'("Ffr=[")',advance='no') 
     ! ! write(5,'("Ffi=[")',advance='no')
     ! do i=1, front%m
        ! do j=1, front%n
           ! write(4,'(f20.15,x)',advance='no')real(front%front(i,j))
           ! ! write(5,'(f20.15,x)',advance='no')aimag(front%front(i,j))
        ! end do
        ! write(4,'(" ")')
        ! ! write(5,'(" ")')
     ! end do
     ! write(4,'("];")')
     ! ! write(5,'("];")')
     ! close(4)
     ! ! close(5)
  ! end if
  
  call qrm_adealloc(front%rowmap)
  call qrm_adealloc(front%ptable)
  call qrm_adealloc(front%colmap)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  if(qrm_mat%icntl(5) .gt. 0) then
     ! h has to be stored
     call _qrm_store_h(front)
     __QRM_CHECK_RET(name,'qrm_store_h',9999)
  end if
  call qrm_adealloc(front%t)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call _qrm_store_r(front)
  __QRM_CHECK_RET(name,'qrm_store_r',9999)

  call qrm_adealloc(front%front)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)
  
  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_clean_front


#if defined (RFPF)

! !==========================================================================================
! !==========================================================================================
! subroutine _qrm_store_h(front)
!   ! this subroutine is meant to store the Householder vectors
!   ! computed during the facto fron front%front to front%h
!   use qrm_mem_mod
!   use _qrm_fdata_mod
!   use _qrm_rfpf_mod
!   implicit none
  
!   type(_qrm_front_type) :: front
  
!   integer :: i, j, cnt, p, c, tsize, cs, hsize

!   ! error management
!   integer                         :: err_act
!   character(len=*), parameter     :: name='qrm_store_h'

!   call qrm_err_act_save(err_act)
  
!   hsize=0
!   j = 0
!   do
!      cs = min(front%nb, front%ne-j) 
!      if(cs .le. 0) exit
!      j  = j+cs
!      hsize = hsize+ cs*(cs+1)/2 + (max(front%stair(j),j) -j)*cs
!   end do
  
!   call qrm_aalloc(front%h, hsize)
!   __QRM_CHECK_RET(name,'qrm_aalloc',9999)
!   cnt=1
!   do c=1, front%ne, front%nb
!      tsize = min(front%nb, front%ne-c+1)
!      ! store V1 (the triangular part)
!      call _qrm_to_rfpf('l', 'u', tsize, front%front(c,c), front%m, front%h(cnt))
!      ! store V2 (the rectangular part)
!      cnt = cnt+(tsize*(tsize+1))/2
!      do j = c, c+tsize-1
!         do i = c+tsize, front%stair(c+tsize-1)
!            front%h(cnt) = front%front(i,j)
!            cnt = cnt+1
!         end do
!      end do
!   end do
  
!   call qrm_err_act_restore(err_act)
!   return

! 9999 continue ! error management
!   call qrm_err_act_restore(err_act)
!   if(err_act .eq. qrm_abort_) then
!      call qrm_err_check()
!   end if

!   return
! end subroutine _qrm_store_h



! !==========================================================================================
! !==========================================================================================
! subroutine _qrm_store_r(front)
!   ! this subroutine is meant to store the R factor
!   ! computed during the facto fron front%front to front%r
!   use qrm_mem_mod
!   use _qrm_fdata_mod
!   use _qrm_rfpf_mod
!   implicit none
  
!   type(_qrm_front_type) :: front
  
!   integer :: i, j, cnt, p, tsize, r, rsize

!   ! error management
!   integer                         :: err_act
!   character(len=*), parameter     :: name='qrm_store_r'

!   call qrm_err_act_save(err_act)
  
!   rsize = (front%npiv*(front%npiv+1)/2 + front%npiv*(front%n-front%npiv))
!   call qrm_aalloc(front%r, rsize)
!   __QRM_CHECK_RET(name,'qrm_aalloc',9999)

!   cnt=1
  
!   do r=1, front%npiv, front%nb
!      tsize = min(front%nb, front%npiv-r+1)
!      call _qrm_to_rfpf('u', 'n', tsize, front%front(r,r), front%m, front%r(cnt))
!      cnt = cnt+(tsize*(tsize+1))/2
!      do j=r+tsize, front%n
!         do i= r, r+tsize - 1 
!            front%r(cnt) = front%front(i,j)
!            cnt = cnt+1
!         end do
!      end do
!   end do
  
!   call qrm_err_act_restore(err_act)
!   return

! 9999 continue ! error management
!   call qrm_err_act_restore(err_act)
!   if(err_act .eq. qrm_abort_) then
!      call qrm_err_check()
!   end if

!   return
! end subroutine _qrm_store_r


#else
!==========================================================================================
!==========================================================================================
subroutine _qrm_store_h(front)
  ! this subroutine is meant to store the Householder vectors
  ! computed during the facto fron front%front to front%h
  use qrm_mem_mod
  use _qrm_fdata_mod
  use _qrm_rfpf_mod
  implicit none
  
  type(_qrm_front_type) :: front
  
  integer :: i, j, cnt, p, c, cs, hsize, m, pk, jp, k

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_store_h'

  call qrm_err_act_save(err_act)
  
  hsize=0
  front%hsize=0

  cnt=1
  do jp = 1, front%ne, front%nb
     pk = min(front%nb, front%ne-jp+1) 
     if(pk .le. 0) exit
     
     do j = jp, jp+pk-1, front%ib
        k = min(front%ib, jp+pk - j)
        if(k .le. 0) exit 
        m = max(front%stair(j+k-1),j+k-1) - j+1
        hsize = hsize+k*m
     end do
  end do

  call qrm_aalloc(front%h, hsize)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  cnt=1
  j = 0
  p = 1

  outer: do jp = 1, front%ne, front%nb
     pk = min(front%nb, front%ne-jp+1) 
     if(pk .le. 0) exit outer
     
     inner: do j = jp, jp+pk-1, front%ib
        k = min(front%ib, jp+pk - j)
        if(k .le. 0) exit inner
        m = max(front%stair(j+k-1),j+k-1) - j+1
        front%hsize = front%hsize+k*(k+1)/2 + k*(m-k)
        do c=1, k
           front%h(cnt:cnt+c-1) = front%t(1:c, j+c-1)
           front%h(cnt+c:cnt+m-1) = front%front(j+c:j+m-1,j+c-1)
           cnt = cnt+m
        end do
     end do inner
     p = p+1
  end do outer

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
end subroutine _qrm_store_h



!==========================================================================================
!==========================================================================================
subroutine _qrm_store_r(front)
  ! this subroutine is meant to store the R factor
  ! computed during the facto fron front%front to front%r
  use qrm_mem_mod
  use _qrm_fdata_mod
  use _qrm_rfpf_mod
  implicit none
  
  type(_qrm_front_type) :: front
  
  integer :: i, j, cnt, rsize, n, c, jp, k, pk

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_store_r'

  call qrm_err_act_save(err_act)

  rsize=0
  front%rsize=0

  cnt=1
  do jp = 1, front%npiv, front%nb
     pk = min(front%nb, front%npiv-jp+1) 
     if(pk .le. 0) exit
     
     do j = jp, jp+pk-1, front%ib
        k = min(front%ib, jp+pk - j)
        if(k .le. 0) exit 
        n = front%n-j+1
        rsize = rsize+k*n
     end do
  end do
  
  call qrm_aalloc(front%r, rsize)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  front%rsize = front%npiv*(front%npiv+1)/2 + &
       & front%npiv*(front%n-front%npiv)

  cnt=1
  outer: do jp = 1, front%npiv, front%nb
     pk = min(front%nb, front%npiv-jp+1) 
     if(pk .le. 0) exit outer
     
     inner: do j = jp, jp+pk-1, front%ib
        k = min(front%ib, jp+pk - j)
        if(k .le. 0) exit inner
        n = front%n-j+1

        do c=1, n
           front%r(cnt:cnt+k-1) = front%front(j:j+k-1,j+c-1)
           cnt = cnt+k
        end do
     end do inner
  end do outer
  
  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
end subroutine _qrm_store_r
#endif


