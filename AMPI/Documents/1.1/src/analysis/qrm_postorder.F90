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
!> @file qrm_postorder.F90
!! This file contains the routine that computes a postorder traversal of a tree
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!>   @brief This subroutine computes a postorder by traversing a tree in dfs.
!!   
!!   @param[in] parent  integer array of size n. parent(i)=j means that node j is the
!!                      father of node i in the tree
!!   
!!   @param[in] n       number of nodes in the tree
!!   
!!   @param[out] porder an integer array of size n containing the postorder
!!
!!   @param[in] weight  an optional array containing nodes weights. If present, the
!!                      children of each node will be sorted by increasing weight.
!!   
subroutine qrm_postorder(parent, n, porder, weight)

  use qrm_mem_mod
  implicit none

  integer           :: n
  integer           :: parent(:), porder(:)
  integer, optional :: weight(:)

  integer, allocatable   :: son(:), brother(:), stack(:)
  integer                :: i, father, br, head, hp, pp, t, w, next
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_analyse'

  call qrm_err_act_save(err_act)

  call qrm_aalloc(son, n)
  call qrm_aalloc(brother, n)
  call qrm_aalloc(stack, n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  son = 0

  ! build a tree that can be traversed top-to-bottom
  if(present(weight)) then
     ! use stack as a workspace
     stack = 0
     do i=1, n
        w = weight(i)
        brother(i) = stack(w)
        stack(w)   = i
     end do

     do w=n, 1, -1
        i = stack(w)
        do
           if(i .eq. 0) exit
           next   = brother(i)
           father = parent(i)
           if(father .ne. 0) then
              brother(i) = son(father)
              son(father) = i
           end if
           i = next
        end do
     end do
  else
     do i=n, 1, -1
        father = parent(i)
        if (father .ne. 0) then
           br          = son(father)
           brother(i)  = br
           son(father) = i
        end if
     end do
  end if


  head = 0
  hp   = 0
  pp   = 1
  ! the tree is processed in dfs. Starting from a root, we go down
  ! and put all the encountered nodes on a stack. When we reach the bottom,
  ! we pop the leaf from the stack, we put it in the postorder and we go down again
  ! along the branch starting from the brother of the node just popped
  do t=1, n
     if (parent(t) .ne. 0) cycle
     ! at this point t is the root of a tree
     hp        = hp+1
     stack(hp) = t
     head      = t
     do
        if(son(head) .eq. 0) then
           ! we reached the bottom
           porder(pp) = head
           pp = pp+1
           hp = hp-1
           if (parent(head) .ne. 0) then
              son(parent(head)) = brother(head)
           end if
           if (hp .eq. 0) then
              ! tree is exhausted
              exit
           end if
           head = stack(hp)
        else
           ! go down one more level
           hp = hp+1
           stack(hp) = son(head)
           head = son(head)
        end if
     end do
  end do

  call qrm_adealloc(son)
  call qrm_adealloc(brother)
  call qrm_adealloc(stack)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine qrm_postorder
