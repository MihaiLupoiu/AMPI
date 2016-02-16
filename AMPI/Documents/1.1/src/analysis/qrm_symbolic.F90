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
!> @file qrm_symbolic.F90
!! This files contains the routine that does the symbolic factorization.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the symbolic QR factorization of a matrix.


!> This routine completely characterizes the structure of fronts in
!! the elimination tree and does a number of other symbolic operations
!! that are essential for the subsequent numerical
!! factorization. Specifically:
!! - it computes the number of rows in each front (note that the
!!   number of columns was already computed by the @link
!!   _qrm_rowcount @endlink routine)
!! - it computes the column indices for each front
!! - it defines a layer in the tree below which, subtrees are treated
!!   sequentially
!! - it estimates the number of flops done in the numerical factorization
!! - it computes the number of nnz in R and H
!! 
!! @param[in] graph This is the adjacency graph of the matrix to be
!!                  factorized in CSC format. On exit, its adata
!!                  member will bemmodfied. This is the global data
!!                  structure holding all the information computed in
!!                  the analysis phase and needed for the numerical
!!                  factorization. On input adata\%rc, adata\%cperm,
!!                  adata\%parent, adata\%cp_ptr, adata\%nnodes,
!!                  adata\%icperm, adata\%rperm, adata\%child,
!!                  adata\%childptr must be as produced by @link
!!                  qrm_compress_data @endlink. On output, the result
!!                  will be sotred in the following fields:
!!                  - nfrows: the number of rows in each front
!!                  - fcol and fcol_ptr: all the column indices of
!!                    front i will be stored in fcol(fcol_ptr(i):fcol_ptr(i+1)-1)
!!                  - small: if small(i)=1 the it means that front i is the root
!!                    of a subtree that will be treated sequentially during the
!!                    numerical factorization.
!!
subroutine _qrm_symbolic(graph)

  use _qrm_spmat_mod
  use qrm_adata_mod
  use qrm_error_mod
  use qrm_sort_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type), target :: graph

  integer                         :: i, j, f, p, pp, ppp, root, node, roff, ne, np
  integer                         :: first, c, ib, nlz, nth, leaves, totleaves
  integer                         :: m, n, k, cyc, nb, fm, fn, fk
  real(kind(1.d0)), allocatable   :: n_weight(:), t_weight(:), lzero_w(:), proc_w(:)
  real(kind(1.d0))                :: rm, rk, rn, totflops, smallth
  integer, allocatable            :: col_map(:), mark(:), stair(:), lzero(:), aux(:)
  integer, pointer                :: porder(:), rc(:), parent(:), fcol(:), fcol_ptr(:)
  type(_qrm_spmat_type)           :: g_csr
  logical                         :: found
  type(qrm_adata_type), pointer   :: adata
  integer(kind=8)                 :: hsize, rsize
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_symbolic'

  call qrm_err_act_save(err_act)

  ! just to simplify
  adata    => graph%adata
  porder   => adata%cperm
  rc       => adata%rc
  parent   => adata%parent
  

  call qrm_aalloc(adata%fcol_ptr, adata%nnodes+1)
  call qrm_aalloc(adata%fcol, sum(rc))
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  ! just to simplify
  fcol     => adata%fcol
  fcol_ptr => adata%fcol_ptr

  
  ! first, determine the columns in each front.  

  ! col_map(j)=f means that global column j is a principal variable in
  ! front f
  call qrm_aalloc(col_map, graph%n)
  call qrm_aalloc(mark, adata%nnodes)
  __QRM_CHECK_RET(name,'qrm_aalloc2',9999)
  mark = 0

  fcol_ptr(1:2)=1

  ! extract the first iteration in order to initialize fcol_ptr as
  ! well
  f = 1
  do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
     j = porder(p)
     col_map(j) = f
  end do

  do f=2, adata%nnodes
     fcol_ptr(f+1) = fcol_ptr(f)+max(rc(f-1),0)
     do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
        j = porder(p)
        col_map(j) = f
     end do
  end do
  
#if defined(debug)
  if(p .ne. graph%n+1) then
     __QRM_PRNT_DBG('("Error in symbolic. i .ne. n ",i5,2x,i5)')p, graph%n
  end if
#endif
  
  ! on input the graph is in csc format. we also need a csr
  ! representation in order to determine the row-subtrees of A'*A
  call _qrm_spmat_convert(graph, g_csr, 'csr', .false.)
  __QRM_CHECK_RET(name,'qrm_spmat_convert',9999)

  do i=1,g_csr%nz
     g_csr%jcn(i) = adata%icperm(g_csr%jcn(i))
  end do

  call sort_mat(g_csr)

  do i=1,g_csr%nz
     g_csr%jcn(i) = adata%cperm(g_csr%jcn(i))
  end do

 
  ! for every front f, for every variable i in f, determine the
  ! coefficients j in (A'*A)(i,1:i) and 
  cyc=0
  do f=1, adata%nnodes
     do p=adata%cp_ptr(f), adata%cp_ptr(f+1)-1
        i = porder(p)

        ! the the coefficient j exists in row i of A'*A iff i and j
        ! are both present in a row of A. So, for each nnz k in column j of A,
        ! every nnz (k,j) in A defines a coefficient (i,j) in A'*A
        fcol(fcol_ptr(f+1)) = i
        fcol_ptr(f+1) = fcol_ptr(f+1)+1

        do pp = graph%jptr(i), graph%jptr(i+1)-1
           k = graph%irn(pp)
           
           ! for every nnz k in column i, go along the corresponding
           ! row and, in the tree, go up until node f (containing i)
           do ppp=g_csr%iptr(k), g_csr%iptr(k+1)-1
              j = g_csr%jcn(ppp)
              if(adata%icperm(j) .ge. adata%icperm(i)) exit

              ! (i,j) is in tril(A'*A). Go up the tree until node f, for
              ! every node met, add column i to the corresponding
              ! front
              node = col_map(j)
              do
                 ! go up the tree
                 if((mark(node) .eq. i) .or. (node .eq. f)) exit

                 fcol(fcol_ptr(node+1)) = i
                 fcol_ptr(node+1) = fcol_ptr(node+1)+1

                 mark(node) = i
                 node = parent(node)
              end do
           end do
        end do
        
     end do
  end do
  call qrm_adealloc(mark)
  call qrm_aalloc(adata%nfrows, adata%nnodes)
  __QRM_CHECK_RET(name,'qrm_aalloc2.5',9999)


  ! n_weight(i) is meant to hold the weight of front i
  ! t_weight(i) is meant to hold the weight of the subtree rooted at i
  call qrm_aalloc(n_weight, adata%nnodes) 
  call qrm_aalloc(t_weight, adata%nnodes)

  n_weight = 0.d0
  t_weight = 0.d0

  call qrm_aalloc(stair, maxval(rc)+1)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  call qrm_get(graph, 'qrm_ib', ib)
  call qrm_get(graph, 'qrm_nb', nb)

  hsize=0
  rsize=0

  ! determine structure and weight of all nodes
  do f=1, adata%nnodes

#if defined(debug)
     col_map=0
#endif

     ! under the assumption of postordered nodes, we can simply sweep
     ! the list of fronts

     ! build the col_map for front f. col_map(k)=j means that global
     ! column k is column j inside front f
     do j=1, rc(f)
        k = fcol(fcol_ptr(f)+j-1)
        col_map(k)=j
     end do

     if(f .eq. 1) then
        roff = 1
     else
        roff = adata%stair(f-1)+1
     end if
     
     stair(1:rc(f)) = 0
     ! assemble the rows from the original matrix
     do p=roff, adata%stair(f)
        ! i is a row of the original matrix to be assembled into front f
        i = adata%rperm(p) 
        ! sweep this row and determine its first coefficient
        first = col_map(g_csr%jcn(g_csr%iptr(i)))
        ! first = rc(f)+1
        ! do pp=g_csr%iptr(i), g_csr%iptr(i+1)-1
           ! j = col_map(g_csr%jcn(pp))
           ! ! TODO: can be optimized id g_csr is sorted
           ! if (j .lt. first) first=j
        ! end do
        stair(first) = stair(first)+1
     end do
     
     ! assemble the CBs from the children
     do ppp=adata%childptr(f), adata%childptr(f+1)-1
        c = adata%child(ppp)
        
        ! ne is the number of Householder vectors computed on the
        ! child c. np is the number of fully assembled pivots in c
        ne = min(rc(c), adata%nfrows(c))
        np = adata%cp_ptr(c+1)-adata%cp_ptr(c)
        
        ! count in all the rows on the CB of c
        do i=np+1, ne
           j = fcol(fcol_ptr(c)+i-1)
           first = col_map(j)
           stair(first) = stair(first)+1
        end do
     end do

     ! finalize stair
     do i=2, rc(f)
        stair(i) = stair(i)+stair(i-1)
     end do
     
     adata%nfrows(f) = stair(rc(f))

     ! At this point it is possible to determine the computational
     ! weight of node f
     ne = min(rc(f),adata%nfrows(f))
     do i=1, ne
        n_weight(f) = n_weight(f)+qrm_count_flops(max(stair(i)-i+1,0),rc(f)-i,1,'panel')
        n_weight(f) = n_weight(f)+qrm_count_flops(max(stair(i)-i+1,0),rc(f)-i,1,'update')
        hsize = hsize+max(stair(i)-i+1,0)
     end do

     ! do the same for R
     np = adata%cp_ptr(f+1)-adata%cp_ptr(f)
     rsize = rsize + np*(np+1)/2 + np*(rc(f)-np)

     t_weight(f) = t_weight(f)+n_weight(f)
     p = parent(f)
     if(p .ne. 0) t_weight(p) = t_weight(p)+t_weight(f)
     
  end do

  totflops = sum(n_weight)

  graph%gstats(qrm_e_facto_flops_) = floor(totflops,kind(graph%gstats(1)))
  graph%gstats(qrm_e_nnz_r_) = rsize
  graph%gstats(qrm_e_nnz_h_) = hsize

#if defined(debug)
  __QRM_PRNT_DBG('("Total estimated number of MFLOPS: ",i10)')floor(totflops)
#endif

  
  call _qrm_spmat_destroy(g_csr, all=.true.)
  __QRM_CHECK_RET(name,'qrm_spmat_destroy',9999)

  ! these are no more needed
  call qrm_adealloc(col_map)
  call qrm_adealloc(stair)
  call qrm_aalloc(lzero, adata%nnodes)
  call qrm_aalloc(adata%small, adata%nnodes)

  ! at this point we start going down the tree until we identify a set
  ! of nodes such that the subtrees rooted at them can be scheduled to
  ! threads with a good load balancing. Small nodes (or subtrees) will
  ! be pruned away during the descent.
  
  call qrm_aalloc(lzero_w, adata%nnodes)
  call qrm_aalloc(aux, adata%nnodes+2, lbnd=0)
  call qrm_get(graph, 'qrm_nthreads', nth)
  
  if(nth .gt. adata%nnodes) nth = adata%nnodes ! you never know

  call qrm_aalloc(proc_w, nth)


  smallth = 0.01
10 continue

  totleaves = 0
  adata%small = 0

  ! goto 20
  nlz = 0
  ! initialize the l0 layer with the root nodes
  do i=1, adata%nnodes
     if(parent(i) .eq. 0) then
        if(t_weight(i) .gt. smallth*totflops) then
           nlz = nlz+1
           lzero(nlz)   = i
           lzero_w(nlz) = t_weight(i)
        else
           adata%small(i) = 1 ! node is too smal; mark it
        end if
     end if
     if(adata%childptr(i+1) .eq. adata%childptr(i)) totleaves = totleaves+1
  end do

  leaves = 0

  ! start the loop 
  godown: do
     if(nth .eq. 1) exit ! shortcut for serial execution
     if(nlz .gt. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) exit ! exit if already too many nodes in l0
     
     proc_w = 0.d0
     ! sort the nodes in l0 in order of decreasing weight
     call qrm_mergesort(nlz, lzero_w(1:nlz), aux(0:nlz+1), order=-1)
     call qrm_mergeswap(nlz, aux(0:nlz+1), lzero(1:nlz), lzero_w(1:nlz))

     ! map subtrees to threads round-robin 
     do i=1, nlz
        ! find the least loaded proc
        p = minloc(proc_w,1)
        proc_w(p) = proc_w(p) + lzero_w(i)
     end do

     ! all the subtrees have been mapped. Evaluate load balance
     rm = minval(proc_w)/maxval(proc_w)
     ! write(*,*)'==>',lzero_w(1:nlz)
     ! write(*,*)'-->',nth,nlz,rm,leaves,totleaves
     if((rm .gt. 0.9) .and. (nlz .ge. 1*nth)) exit ! if balance is higher than 90%, we're happy

     ! if load is not balanced, replace heaviest node with its kids (if any)
     found = .false.
     findn: do
        if(leaves .eq. totleaves) exit godown

        if(leaves .eq. nlz) then
           if(nlz .ge. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) then 
              exit godown ! all the nodes in l0 are leaves. nothing to do
           else
              smallth = smallth/2.d0
              if(smallth .lt. 0.0001) then
                 exit godown
              else
                 goto 10
              end if
           end if
        end if
        n = lzero(leaves+1) ! n is the node that must be replaced

        ! appends children of n 
        do p=adata%childptr(n), adata%childptr(n+1)-1
           c = adata%child(p)
           if(t_weight(c) .gt. smallth*totflops) then
              ! this child is big enough, add it
              found = .true.
              nlz = nlz+1
              lzero  (nlz) = c
              lzero_w(nlz) = t_weight(c)
           else
              adata%small(c) = 1 ! node is too smal; mark it
           end if
        end do
        if(found) exit findn ! if at least one child was added then we redo the mapping
        leaves = leaves+1
     end do findn

     ! swap n with last element
     lzero  (leaves+1) = lzero  (nlz)
     lzero_w(leaves+1) = lzero_w(nlz)
     nlz = nlz-1

  end do godown

  ! mark all the children of nodes in l0
  do i=1, nlz
     n = lzero(i)
        do p=adata%childptr(n), adata%childptr(n+1)-1
           c = adata%child(p)
           adata%small(c) = 1
        end do
  end do

20 continue
  
  t_weight = t_weight/totflops * 100
  ! call qrm_print_tree('atree.dot',adata, t_weight)


  call qrm_adealloc(lzero)
  call qrm_adealloc(lzero_w)
  call qrm_adealloc(proc_w)
  call qrm_adealloc(aux)
  call qrm_adealloc(n_weight)
  call qrm_adealloc(t_weight)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return


contains

  subroutine sort_mat(mat)
    use qrm_sort_mod
    implicit none

    type(_qrm_spmat_type) :: mat
    integer, allocatable :: aux(:)
    integer :: i, n, k1, k2



    allocate(aux(0:mat%n+1))

    do i=1,mat%m
       k1= mat%iptr(i)
       k2= mat%iptr(i+1)-1
       n = k2-k1+1
       call qrm_mergesort(n, mat%jcn(k1:k2), aux(0:n+1))
       call qrm_mergeswap(n, aux(0:n+1), mat%jcn(k1:k2))
    end do

    deallocate(aux)


    return
  end subroutine sort_mat

end subroutine _qrm_symbolic
