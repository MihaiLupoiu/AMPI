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
!> @file qrm_reorder_tree.F90
!! This file contains the routine that computes a reordering of the tree to reduce the search
!! space for task scheduling
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This subroutine reorders the assembly tree in order to reduce
!! the tasks search space.

!> A nice side effect of this is that memory consumption is also reduced. The code here
!! basically follows the idea in:
!!
!! J. W. H. Liu. <em>On the storage requirement in the out-of-core
!!                   multifrontal method for sparse factorization.</em>
!!                   ACM Transactions on Mathematical Software, 12:127–148, 1986.
!!
!! @param[in,out] adata a qrm_adata_type data containing a full
!!                caracterization of the assembly tree and fronts
!!                structure.
!! @todo Add better explication of the algorithm
!!
subroutine qrm_reorder_tree(adata)

  use qrm_adata_mod
  use qrm_mem_mod
  use qrm_sort_mod
  implicit none

  type(qrm_adata_type) :: adata

  integer, allocatable :: peaks(:), children_peaks(:), aux(:), &
       & roots(:), child(:), tnch(:)
  integer :: root, maxch, nroots, nleaves, i, nl, leaf, f, c
  integer :: nch, nnodes, node, peak

  call qrm_aalloc(peaks, adata%nnodes)

  nroots  = 0
  maxch   = 0
  nleaves = 0
  do node = 1, adata%nnodes
     if(adata%parent(node) .eq. 0) then
        nroots = nroots+1
     end if
     
     nch = adata%childptr(node+1)-adata%childptr(node)
     if(nch .gt. maxch) maxch = nch
     if(nch .eq. 0    ) nleaves = nleaves+1
     peaks(node) = -nch
  end do

  ! uncomment to bypass the reordertree
  ! adata%nleaves = 0
  ! call qrm_aalloc(adata%leaves, nleaves)
  ! do node = 1, adata%nnodes
     ! nch = adata%childptr(node+1)-adata%childptr(node)
     ! if(nch .eq. 0) then
        ! adata%nleaves = adata%nleaves+1
        ! adata%leaves(adata%nleaves) = node
     ! end if
  ! end do
  ! return
  ! until here



  maxch = max(maxch, nroots)
  call qrm_aalloc(children_peaks, maxch)
  call qrm_aalloc(aux, maxch+2, lbnd=0)
  call qrm_aalloc(roots, nroots)

  nnodes = adata%nnodes
  node   = 1
  nl     = nleaves
  nroots = 0

  main: do leaf=1, adata%nnodes

     nch = adata%childptr(leaf+1)-adata%childptr(leaf)
     if(nch .eq. 0) then
        node = leaf

        ! start going up the subtree
        do
           if(nnodes .eq. 0) exit main
           ! peaks(node).eq.0 only occurs if node is a leaf or
           ! if its peak can be computed
           if(peaks(node) .ne. 0) exit
           
           peak=0
           ! sort node's children (if any) and compute its peak
           nch = 0
           if(adata%childptr(node) .ne. adata%childptr(node+1)) then
              do i = adata%childptr(node), adata%childptr(node+1)-1
                 nch = nch+1
                 children_peaks(nch) = peaks(adata%child(i))
              end do
              ! sort
              call qrm_mergesort(nch, children_peaks(1:nch), aux(0:nch+1), order=-1)

              call qrm_mergeswap(nch, aux(0:nch+1), children_peaks(1:nch), &
                   & adata%child(adata%childptr(node):adata%childptr(node+1)-1))
              ! compute node's peak
              do i=1, nch
                 peak = max(peak, i-1+children_peaks(i))
              end do
           end if
           
           peaks(node) = max(peak,nch+1)
           nnodes = nnodes -1

           f = adata%parent(node)
           if (f .ne. 0) then
              peaks(f) = peaks(f)+1
              if(peaks(f) .eq. 0) then
                 node = f
              else
                 cycle main
              end if
           else
              nroots = nroots+1
              roots(nroots) = node
           end if
        end do

     end if

  end do main

  ! sort roots
  do i=1, nroots
     ! write(*,'("Root peak: ",i10)')peaks(roots(i))
     children_peaks(i) = peaks(roots(i))
  end do
  
  call qrm_mergesort(nroots, children_peaks(1:nroots), aux(0:nroots+1), order=-1)
  call qrm_mergeswap(nroots, aux(0:nroots+1), children_peaks(1:nroots), &
       & roots(1:nroots))

  call qrm_adealloc(children_peaks)
  call qrm_adealloc(aux)
  call qrm_aalloc(adata%leaves, nleaves)
  call qrm_aalloc(tnch, adata%nnodes)
  call move_alloc(peaks, child)
  child = 0

  ! compute the list of leaves to start facto. Even if the whole tree
  ! is reordered, only the leaves in l0 are returned; this may have to
  ! be improved. TODO
  
  tnch = 0
  do node=1, adata%nnodes
     do i=adata%childptr(node), adata%childptr(node+1)-1
        c = adata%child(i)
        if(adata%small(c) .eq. 0) tnch(node) = tnch(node)+1
     end do
  end do

  nleaves = 0
  roots_loop: do root=1, nroots
     node = roots(root)
     ! go down the subtree
     
     nodes_loop: do
        if(tnch(node) .eq. 0) then
           nleaves = nleaves+1
           adata%leaves(nleaves) = node
        else
        
           nch=adata%childptr(node+1) - adata%childptr(node)
           do
              if(child(node) .ge. nch) exit
                 child(node) = child(node)+1
                 c = adata%child(adata%childptr(node)+child(node)-1)
                 if(adata%small(c) .ne. 1) then
                    node = c
                    cycle nodes_loop
                 end if
              end do
           end if
        node = adata%parent(node)
        if (node .eq. 0) cycle roots_loop
        cycle nodes_loop

     end do nodes_loop
  end do roots_loop

  adata%nleaves = nleaves
  
  
  call qrm_adealloc(roots)
  call qrm_adealloc(child)
  call qrm_adealloc(tnch)
  

  return

end subroutine qrm_reorder_tree
