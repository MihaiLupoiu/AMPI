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
!> @file qrm_init_front.F90
!! This file contains a routine which performs the initialization of a front
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This routine initializes a front.

!> Front initialization consists in
!! creating the structure of the front and all the data structures
!! that allow the assembly of the front itself.
!!
!! Assume the column permutation is:
!! @verbatim
!! (/7, 3, 8, 1, 9, 5, 6, 2/)
!! @endverbatim
!!
!! Values from the original matrix A:
!! @verbatim
!!     7 3 1 5 6
!!   +----------
!!  3| a   a a a
!!  6| a a     a
!!  2|       a a
!! @endverbatim
!!
!!
!! Values from child X:
!! @verbatim
!!     7 3 8 9 2
!!   +----------
!!  5| x x x x x
!!  1|   x x x x
!!  7|     x x x
!! @endverbatim
!!
!!
!! Values from child Y:
!! @verbatim
!!     3 1 5
!!   +------ 
!!  4| y y y
!!  8|   y y
!!  9|     y
!! @endverbatim
!!
!!
!! Final result:
!! @verbatim
!!     7 3 8 1 9 5 6 2
!!   +---------------- 
!!  3| a     a   a a        
!!  6| a a         a
!!  5| x x x   x     x
!!  1|   x x   x     x
!!  4|   y   y   y
!!  7|     x   x     x
!!  8|       y   y
!!  2|           a a
!!  9|           y
!! @endverbatim
!!
!! For a description of the following fields see @link _qrm_fdata_mod::_qrm_front_type @endlink
!! @verbatim
!! front%stair=(/3, 5, 6, 7, 7, 9/)
!! front%cols=(/7, 3, 8, 1, 9, 5, 6, 2/)
!! front%rows=(/3, 6, 5, 1, 4, 7, 8, 2, 9/)
!! front%colmap(/4, 8, 2, 0, 6, 7, 1, 3, 5/)
!! X%rowmap=(/3, 4, 6/)
!! Y%rowmap=(/5, 7, 9/)
!! @endverbatim
!!
!! @param[in,out] qrm_mat the big problem data structure
!!
!! @param[in] fnum the id of the front to be activated
!! 
!! @param[in] par whether this front will be treated by multiple
!!                threads or not. this is needed because in the case
!!                where the front belongs to a sequential subtree, it
!!                can be assembled immediatebly; otherwise, references
!!                will be computed for later assembly
!!
!! @param[in] work optional workspace
!!

subroutine _qrm_init_front(qrm_mat, fnum, par, work)
  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_sort_mod
  implicit none
  
  type(_qrm_spmat_type), target :: qrm_mat
  integer                        :: fnum
  logical                        :: par
  integer, optional, target      :: work(:)

  integer :: n_rows_orig, i, j, row, col, c, child, roff, f
  integer :: m, n, npiv, k, p, ne, cbr, cbc
  integer :: nb, b, father

  integer, allocatable :: first(:)
  integer, pointer     :: gcolmap(:)

  type(_qrm_front_type), pointer :: front, cfront
  type(qrm_adata_type), pointer  :: adata
  type(_qrm_fdata_type), pointer :: fdata
  logical                        :: map, sfront

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_init_front'
  
  call qrm_err_act_save(err_act)

  front  => null()
  cfront => null()
  nullify(gcolmap)

  ! to make things easier
  adata      => qrm_mat%adata
  fdata      => qrm_mat%fdata
  front      => fdata%front_list(fnum)
  front%m    =  adata%nfrows(fnum)
  front%n    =  adata%rc(fnum)
  if(par) then
     front%nb   =  qrm_mat%icntl(qrm_nb_)
  else
     front%nb   =  qrm_mat%icntl(qrm_ib_)
  end if
  front%ib   =  qrm_mat%icntl(qrm_ib_)
  father     =  adata%parent(fnum)

  if( (front%n .le. 0) .or. (front%m .le. 0)) then
     ! nothing to do here. Mark the front as done
     front%status = qrm_done_
     goto 10
  end if

  if (fnum .eq. 1) then
     roff = 0
  else
     roff = adata%stair(fnum-1)
  end if

  ! set up a few variables that come handy in the activation and
  ! factorization of the front
  m          = front%m
  n          = front%n
  ! The number of eliminations to be performed on the front f
  front%ne   = min(m,n) 
  ! The number of panels in f
  front%np   = max((front%ne-1)/front%nb + 1,0)
  ! The number of block-columns in f
  front%nc   = (front%n-1)/front%nb + 1
  ! The number of pivots in f
  front%npiv = min(adata%cp_ptr(fnum+1)-adata%cp_ptr(fnum),front%ne)
  cbr        = front%m-front%npiv    ! The number of rows in R
  cbc        = front%n-front%npiv    ! The number of rows in the cb

  ! if the size of the front is smaller than ib, then there's no need
  ! to compute front%stair
  sfront = front%ne .lt. front%ib
  ! if(sfront) write(*,*)'====> sfront!!!'
  ! allocate front%cols and copy the list f column indices in it
  call qrm_aalloc(front%cols, front%n) 
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  front%cols(1:front%n) = adata%fcol(adata%fcol_ptr(fnum): &
       & adata%fcol_ptr(fnum+1)-1)

  if(par) then
     ! front%ptable is the progress table for keeping track of the work
     ! done on this front. ptable(i) is initialized to the number of
     ! contributions that have to be assembled into block-column i
     call qrm_aalloc(front%ptable, front%nc )
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)
     front%ptable = 0
  end if

  ! build the gcolmap for front f. gcolmap(k)=j means that global
  ! column k is column j inside front f
  if(present(work)) then
      gcolmap => work
   else
      call qrm_palloc(gcolmap, qrm_mat%n)
      __QRM_CHECK_RET(name,'qrm_palloc',9999)
   end if

  do j=1, front%n
     k = front%cols(j)
     gcolmap(k)=j
  end do

  call qrm_aalloc(front%stair, front%n+1)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  if(.not. sfront) then
     front%stair = 0
     call qrm_aalloc(first, front%anrows)
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)

     first = front%n+1
     ! count in the rows from the original matrix
     do i=1, front%anrows
        ! sweep this row and determine its first coefficient
        do p=front%aiptr(i), front%aiptr(i+1)-1
           j = gcolmap(front%ajcn(p))
           ! TODO: can be optimized id ajcn is sorted
           if (j .lt. first(i)) first(i)=j
        end do
        front%stair(first(i)+1) = front%stair(first(i)+1)+1
     end do

     ! count in rows coming from CBs
     do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
        c = adata%child(p)
        cfront => fdata%front_list(c)
        
        ! ne is the number of Householder vectors computed on the
        ! child c. npiv is the number of fully assembled pivots in c
        ne = cfront%ne
        npiv = cfront%npiv
        if(npiv.eq.ne) cycle
        
        ! count in all the rows on the CB of c
        do i=1, ne-npiv
           f = gcolmap(cfront%cols(npiv+i)) !cfront%colmap(i)
           front%stair(f+1) = front%stair(f+1)+1
        end do
     end do
     
     ! finalize stair
     do i=2, front%n+1
        front%stair(i) = front%stair(i)+front%stair(i-1)
     end do
  else
     front%stair = front%m
  end if

  call qrm_aalloc(front%rows, front%m)
  call qrm_aalloc(front%front, front%m, front%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  front%front = _qrm_zero
  row = 0
  ! at this point we're ready to assemble the rows from the original matrix
  ! FIXME: deal with the presence of duplicates
  do i=1, front%anrows
     if(sfront) then
        row = i
     else
        f = first(i) ! f is the front-local column index of the first
        ! coefficient in this row 
        front%stair(f) = front%stair(f)+1
        row = front%stair(f)
     end if

     front%rows(row) = adata%rperm(roff+i)
     !sweep this row and assemble it
     do p=front%aiptr(i), front%aiptr(i+1)-1
        col = gcolmap(front%ajcn(p))
        front%front(row,col) = front%front(row,col)+front%aval(p)
     end do
  end do


  ! now we can build the row-mappings for assemblying the children
  ! later or, if required, assemble the children directly
  do p = adata%childptr(fnum), adata%childptr(fnum+1)-1
     c = adata%child(p)
     cfront => fdata%front_list(c)

     map = par .and. (adata%small(c) .ne. 1)

     ! ne is the number of Householder vectors computed on the
     ! child c. npiv is the number of fully assembled pivots in c
     ne   = cfront%ne
     npiv = cfront%npiv
     if(npiv.eq.ne) cycle

     ! this is the row mapping on the child. cfront%rowmap(k)=i
     ! means that the k-th row of cfront will be assembled into the
     ! i-th row of front
     do i=npiv+1, ne
        if(sfront) then
           row = row+1
        else
           f = gcolmap(cfront%cols(i))
           front%stair(f) = front%stair(f)+1
           row = front%stair(f)
        end if
        front%rows(row) = cfront%rows(i)
        cfront%rowmap(i-npiv) = row
     end do

     if(.not. map) then
        do j=npiv+1, cfront%n
           ! this is the column mapping on the child. cfront%colmap(k)=j
           ! means that the k-th column of cfront will be assembled into the
           ! j-th column of front
           cfront%colmap(j-npiv) = gcolmap(cfront%cols(j))

           ! fill in the front with coefficients front chil front cfront
           do i=npiv+1, min(j,ne)
              row = cfront%rowmap(i-npiv)
              col = cfront%colmap(j-npiv)
              front%front(row, col) = cfront%front(i,j)
           end do
        end do
     else
        ! initialize ptable
        do j=1, cfront%n - npiv
           ! this is the column mapping on the child. cfront%colmap(k)=j
           ! means that the k-th column of cfront will be assembled into the
           ! j-th column of front
           cfront%colmap(j) = gcolmap(cfront%cols(npiv+j))

           ! count the number of contributions that have to be
           ! assembled inside the corresponding block-column
           b = (cfront%colmap(j)-1)/front%nb+1  ! j in cfront goes into b-column b of front
           front%ptable(b) = front%ptable(b)-1
        end do
     end if
  end do

  call qrm_adealloc(first)
  if(present(work)) then
     nullify(gcolmap)
  else
     call qrm_pdealloc(gcolmap)
     nullify(gcolmap)
  end if
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  ! allocate tau and t
  call qrm_aalloc(front%tau, front%ne)
  call qrm_aalloc(front%t, front%ib, front%ne)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)


  ! FIXME: check if this is really necessary
  front%tau = _qrm_zero
  front%t   = _qrm_zero

  if (father .ne. 0) then
     call qrm_aalloc(front%rowmap, cbr)
     call qrm_aalloc(front%colmap, cbc)
     __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  end if

  ! no more needed, can be deallocated
  call qrm_adealloc(front%aiptr)
  call qrm_adealloc(front%ajcn)
  call qrm_adealloc(front%aval)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)


10 continue

  if (father .ne. 0) then
     front => fdata%front_list(father)
     !$ if(par) call omp_set_lock(front%lock)
     front%status = front%status+1
     !$ if(par) call omp_unset_lock(front%lock)
  end if


  call qrm_err_act_restore(err_act)
  return

9999 continue
  ! We get to this point if some allocation failed. Deallocate
  ! everything and return
  call qrm_adealloc(first)
  if(present(work)) then
     nullify(gcolmap)
  else
     call qrm_pdealloc(gcolmap)
  end if
  call qrm_adealloc(front%cols)
  call qrm_adealloc(front%ptable)
  call qrm_adealloc(front%stair)
  call qrm_adealloc(front%rows)
  call qrm_adealloc(front%front)
  call qrm_adealloc(front%tau)
  call qrm_adealloc(front%t)
  call qrm_adealloc(front%colmap)
  call qrm_adealloc(front%rowmap)

  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_init_front
