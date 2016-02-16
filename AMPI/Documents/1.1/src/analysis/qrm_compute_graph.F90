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
!> @file qrm_compute_graph.F90
!! This file contains the routine that computes the graph of a matrix
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief Computes the adjacency graph of a matrix

!> This subroutine computes the column graph associated to an input matrix qrm_mat 
!! COO format. The output graph has no duplicates as well as no self-edges
!!
!! @param[in] qrm_mat the input matrix
!!
!! @param[out] graph  the adjacenc graph in CSC format
!!
subroutine _qrm_compute_graph(qrm_mat, graph)

  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod
  use _qrm_analysis_mod, savesym => _qrm_compute_graph

  implicit none

  type(_qrm_spmat_type)              :: qrm_mat
  type(_qrm_spmat_type), intent(out) :: graph

  integer              :: i, j, pnt, savepnt, dups, ii
  integer, allocatable     :: dupsmap(:)
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_compute_graph'

  call qrm_err_act_save(err_act)

  call qrm_aalloc(dupsmap, qrm_mat%m)

  ! At this moment, if the matrix is centralized, we just need to 
  ! convert the structure from COO to CSC. Otherwise, we also need
  ! to gather it from the other procs.

  call _qrm_spmat_convert(qrm_mat, graph, 'csc', values=.false.)
  __QRM_CHECK_RET(name,'alloc/convert',9999)
  
  ! make a pass in order to remove duplicates
  dupsmap = 0
  pnt     = 0
  dups    = 0
  savepnt = 1
  do j=1, graph%n
     do ii=graph%jptr(j), graph%jptr(j+1)-1
        i = graph%irn(ii)
        if(dupsmap(i) .eq. j) then
           ! duplicate entry, skip
           dups = dups+1
        else
           ! flag the entry as visited
           dupsmap(i) = j
           pnt = pnt+1
           graph%irn(pnt) = i
        end if
     end do
     graph%jptr(j) = savepnt
     savepnt = pnt+1
  end do
  graph%jptr(graph%n+1)=savepnt
  graph%nz = pnt
  graph%icntl = qrm_mat%icntl
  graph%rcntl = qrm_mat%rcntl

  __QRM_PRNT_DBG('("Number of duplicates in the matrix: ",i10)')dups
  
  
  call qrm_adealloc(dupsmap)
  __QRM_CHECK_RET(name,'qrm_pdealloc',9999)
  
  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
  
end subroutine _qrm_compute_graph
