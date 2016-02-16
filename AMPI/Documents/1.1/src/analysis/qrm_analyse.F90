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
!> @file qrm_analyse.F90
!! This file contains the main analysis driver
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This is the driver routine for the analysis phase.

!> This routine performa a number of symbolic operations in
!! preparation for the numerical factorization:
!! -# it computes the graph of the matrix removing duplicates
!! -# it detects the presence of singletons in the matrix
!! -# it computes a column permutation in order to reduce the fill-in
!! -# it computes the elimination tree
!! -# it postorders the elimination tree
!! -# it computes the rowcount
!! -# it performs an amalgamation of the elimination tree (with fill)
!! -# it computes a row permutation of the matrix in order to have 
!!    a global staircase structure
!! -# merges the singletons into the results of the above operations
!! -# it compresses the results of the above operations to get a size
!!    which is proportional to the number of nodes in the amalgamated tree
!! -# it doea a symbolic factorization which completely characterizes the
!!    structure of fronts etc.
!! -# it computes a tree traversal order which reduces the search space 
!!    for scheduling tasks in the numerical factorization
!!
!! @param[in,out] qrm_mat a qrm_spmat_type data which contains the input matrix.
!!                        On output qrm_mat%adata will contain the results of the 
!!                        analysis phase
!!
!! @param[in] transp      a character saying whether to do analysis on A or A'
!!
subroutine _qrm_analyse(qrm_mat, transp)

  use _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use _qrm_utils_mod
  use qrm_string_mod
  use _qrm_analysis_mod, savesym => _qrm_analyse

  implicit none

  type(_qrm_spmat_type), target   :: qrm_mat
  character, optional, intent(in) :: transp

  type(_qrm_spmat_type), target   :: graph
  integer, allocatable            :: parent(:), rc(:), srow(:), scol(:)
  integer, allocatable            :: mcperm(:), mrperm(:), stair(:), nvar(:)
  integer                         :: ncsing, nrsing
  integer, pointer                :: cperm_p(:)=>null(), rperm_p(:)=>null()
  integer, pointer                :: cperm(:)=>null(), rperm(:)=>null(), tmp(:)
  integer                         :: i, j, ii, lastrc, cnt, ierr
  real(kind(1.d0))                :: t1, t2
  logical                         :: sing
  character                       :: itransp
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_analyse'

  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Entering the analysis driver")')

  call _qrm_check_spmat(qrm_mat, qrm_analyse_)
  __QRM_CHECK_RET(name,'qrm_check_spmat',9999)
  
  call qrm_adata_destroy(qrm_mat%adata)
  __QRM_CHECK_RET(name,'qrm_adata_destroy',9999)

  ncsing=0
  nrsing=0

  nullify(cperm_p,rperm_p,cperm,rperm,tmp)

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
  end if

  call _qrm_compute_graph(qrm_mat, graph)
  __QRM_CHECK_RET(name,'qrm_compute_graph',9999)

  call qrm_get(qrm_mat, 'qrm_sing', i)
  sing = i .eq. qrm_yes_

  if(sing) then
     ! do singleton detection and prune graph
     call _qrm_detect_singletons(graph, scol, srow, mrperm, mcperm, nrsing, ncsing)
     __QRM_CHECK_RET(name,'qrm_detect_singletons',9999)

     __QRM_PRNT_DBG('("row/column singletons ",i10,x,i10)')nrsing,ncsing

  else 
     ncsing = 0
     nrsing = 0
  end if

  ! the cperm array will contain
  call qrm_aalloc(graph%adata%cperm, graph%n+1)  ! FIXME: can't rememebr why n+1, check
  call qrm_aalloc(graph%adata%rperm, graph%m)
  __QRM_CHECK_RET(name,'qrm_alloc',9999)
  cperm => graph%adata%cperm
  rperm => graph%adata%rperm

  graph%adata%ncsing = ncsing
  graph%adata%nrsing = nrsing

  ! shortcut: in case a simple permutation is enough to have a triangular
  ! matrix
  if ( (ncsing .eq. qrm_mat%n) .or. (nrsing .eq. qrm_mat%m) ) goto 9998 ! FIXME: goto 10 is not the right thing to do

  ! at this point we the graph and the singletons. Time to go for the ordering
  ! the ordering will be computed inside cperm(ncsing+1:n) and will be mapped to
  ! the columns of A only later using mcperm
  call _qrm_do_ordering(graph, cperm, qrm_mat%cperm_in) 
  __QRM_CHECK_RET(name,'qrm_do_ordering',9998)

  ! build the elimination tree
  call qrm_aalloc(parent, graph%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9998)

  call _qrm_elim_tree(graph, cperm, parent)
  __QRM_CHECK_RET(name,'qrm_elim_tree',9998)

  ! compute a postorder traversal of the tree
  call qrm_postorder(parent, graph%n, cperm)
  __QRM_CHECK_RET(name,'qrm_postorder',9998)

  call qrm_aalloc(rc, graph%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9998)

  ! do the symbolic facto
  call _qrm_rowcount(graph, parent, cperm, rc)
  __QRM_CHECK_RET(name,'qrm_rowcount',9998)

  call qrm_aalloc(nvar, graph%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9998)

  ! amalgamate the tree
  call qrm_amalg_tree(graph%n, parent, rc, cperm, nvar, qrm_mat%icntl(3), qrm_mat%rcntl(1))
  __QRM_CHECK_RET(name,'qrm_amalg_tree',9998)
  ! call qrm_prnt_array(nvar(1:qrm_mat%n), 'nvar')

  call qrm_aalloc(stair, graph%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9998)

  ! Compute the row permutation to put the matrix in staircase form
  call _qrm_rowperm(graph, cperm, rperm, nvar, stair)
  __QRM_CHECK_RET(name,'qrm_rowperm',9998)

  ! Rework all the data and compress it to nnodes size instead of n
  call qrm_compress_data(graph%adata, cperm, parent, rc, stair, graph%n)
  __QRM_CHECK_RET(name,'qrm_compress_data',9998)

  ! Do the symbolic facto: estimate fronts size, flops, parallelism etc.
  call _qrm_symbolic(graph)
  __QRM_CHECK_RET(name,'qrm_symbolic',9998)

  ! reorder the tree in order to minimize something
  call qrm_reorder_tree(graph%adata)

  if (ncsing .gt. 0) then
     ! the result of the ordering, symbolic facto and amalgamation
     ! has to be mapped back to the original matrix
     call _qrm_attach_singletons(graph, scol, srow, mrperm, mcperm, nrsing, ncsing)
     __QRM_CHECK_RET(name,'qrm_attach_singletons',9998)

  end if

  call qrm_adata_move(graph%adata, qrm_mat%adata)
  qrm_mat%gstats(1:3) = graph%gstats(1:3)

  ! compute the number of rows in each front and some global stats
  ! call qrm_prnt_array(qrm_mat%adata%cperm(1:qrm_mat%n), 'cp')
  ! call qrm_prnt_array(qrm_mat%adata%cp_ptr(1:qrm_mat%adata%nnodes+1), 'pt')
  ! call qrm_prnt_array(qrm_mat%adata%rperm(1:qrm_mat%m), 'rp')
  

9998 continue

  ! in case transp=='t' restore the row and column indices as well as m and n
  if(qrm_str_tolower(itransp) .eq. 't') then
     tmp => qrm_mat%irn
     qrm_mat%irn => qrm_mat%jcn
     qrm_mat%jcn => tmp
     i = qrm_mat%m
     qrm_mat%m = qrm_mat%n
     qrm_mat%n = i
  end if

  call qrm_adealloc(mcperm)
  call qrm_adealloc(mrperm)
  call qrm_adealloc(srow)
  call qrm_adealloc(scol)
  call qrm_adealloc(parent)
  call qrm_adealloc(rc)
  call qrm_adealloc(nvar)
  call qrm_adealloc(stair)
  call qrm_spmat_destroy(graph, all=.true.)
  __QRM_CHECK_RET(name,'qrm_dealloc',9999)


  qrm_mat%adata%ok = .true.
  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_analyse


