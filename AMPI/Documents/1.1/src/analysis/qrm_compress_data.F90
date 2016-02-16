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
!> @file qrm_compress_data.F90
!! this file contains the subroutine that compresses the result of several operations done
!! during the analysis to a size that is proportional to the number of nodes in the elimination
!! tree. These data are of size ~n in input
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This routine compresses the results of a number of operations in the analysis phase.
!! Basically, the input data is of size n and the output of size adata%nnodes (which
!! is the number of nodes in the elimination tree).

!! The arrays on input are all of size n (where n is the column size of
!! the original matrix). The purpose of this subroutine is to store the
!! same information into arrays of size nsteps, where nsteps is the
!! number of nodes in the assembly tree. Only porder will remain of size n
!!
!! @param[out] adata  on output the following fields will be
!!          modified:
!!          - cp_ptr
!!          - rc
!!          - parent
!!          - nnodes
!!          - stair
!!          - child. for each node, it contains the list of its children
!!          - childptr. pointers to the list of children. chil(childptr(i),...,childptr(i+1)-1)
!!            contains all the children of node i
!!          - icperm. the inverse column permutation is built
!!          
!! @param[in] porder  contains the postorder. porder(i)=k means that the i-th column
!!                    in the computed ordering is column k in the original matrix.
!!                    Inthis postorder principal variables always come before the
!!                    corresponding subordinate variables.
!!          
!! @param[in] parent  contains the assembly tree. parent(i)=k:
!!                    - k>0: means that i is a principal variable in a node and k
!!                      is the principal variable in the father's node
!!                    - k<0: means that i is a subordinate variable in a supernode
!!                      whose principal variable is k.
!!          
!! @param[in] rc      this array contains the rowcount, i.e., rc(i)=k means that
!!                    in the R factor the rows corresponding to the node whose
!!                    principal variable is i have k nonzeroes. rc(i)=0 for all
!!                    subordinate varibales and rc(i)=-1 for all the column singletons
!!          
!! @param[in] stair   stair(i) contains the number of rows in the step related to
!!                    node i. see @link ::_qrm_rowperm_ @endlink
!!          
!! @param[in] n       the number of columns in the original matrix (i.e. the size of the
!!                    input arrays)
!!          
subroutine qrm_compress_data(adata, porder, parent, rc, stair, n)
  
  use qrm_mem_mod
  use qrm_adata_mod
  use qrm_common_mod, protect => qrm_compress_data
  implicit none
  
  type(qrm_adata_type) :: adata
  integer              :: porder(:), parent(:), rc(:), stair(:)
  integer              :: n
  
  integer              :: i, pnt, svar, f, ss
  integer, allocatable :: work(:)
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_compress_data'

  call qrm_err_act_save(err_act)

  adata%nnodes = 0
  do i=1, n
     if(rc(i) .ne. 0) adata%nnodes = adata%nnodes+1
  end do
  
  call qrm_aalloc(adata%cp_ptr, adata%nnodes+1)
  call qrm_aalloc(work, n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  
  work = 0
  adata%cp_ptr=0
  ! build pointers and inverse mapping
  ! work(i)=k means that principal variable i is in node number k
  adata%cp_ptr(1) = 1
  work(porder(1)) = 1
  pnt             = 2
  do i=2, n
     svar=porder(i)
     if (rc(svar) .ne. 0) then
        adata%cp_ptr(pnt) = i
        work(svar) = pnt
        pnt = pnt+1
     end if
  end do
  adata%cp_ptr(adata%nnodes+1) = n+1

  ! build adata%parent
  ! build adata%rc
  ! build adata%child and adata%childptr
  call qrm_aalloc(adata%parent, adata%nnodes)
  call qrm_aalloc(adata%rc, adata%nnodes)
  call qrm_aalloc(adata%stair, adata%nnodes)
  call qrm_aalloc(adata%child, adata%nnodes)
  call qrm_aalloc(adata%childptr, adata%nnodes+1)
  call qrm_aalloc(adata%icperm, n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  adata%childptr = 0
  ss = 0
  do i=1, adata%nnodes
     svar = porder(adata%cp_ptr(i))
     f    = parent(svar)
     adata%rc(i) = rc(svar)
     ! adata%nfrows(i) = stair(svar)-ss
     adata%stair(i) = stair(svar)
     ss = adata%stair(i)
     if(f .eq. 0) then
        adata%parent(i) = 0
     else     
        adata%parent(i) = work(f)
        adata%childptr(work(f)+1) = adata%childptr(work(f)+1)+1
     end if
  end do

  adata%childptr(1) = 1
  do i=2, adata%nnodes+1
     adata%childptr(i) = adata%childptr(i-1)+adata%childptr(i)
  end do


  adata%child=0
  work(1:adata%nnodes)=0
  do i=1, adata%nnodes
     f = adata%parent(i)
     if (f .ne. 0) then
        adata%child(adata%childptr(f)+work(f)) = i
        work(f) = work(f)+1
     end if
  end do

  do i=1, n
     adata%icperm(adata%cperm(i)) = i
  end do


  call qrm_adealloc(work)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine qrm_compress_data
