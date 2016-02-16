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
!> @file qrm_do_metis.F90
!! This file contains the routine that computes a METIS permutation of the input matrix
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> Please refer to:
!!
!! <em>A fast and high quality multilevel scheme for partitioning
!! irregular graphs</em>. George Karypis and Vipin Kumar.  International
!! Conference on Parallel Processing, pp. 113-122, 1995
!! 
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_metis(graph, cperm)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_metis
  use qrm_mem_mod
  use iso_c_binding

  implicit none

  type(_qrm_spmat_type) :: graph
  integer               :: cperm(:)

  interface
     subroutine qrm_metis(n, iptr, jcn, cperm, iperm) bind(c, name='qrm_metis')
       use iso_c_binding
       integer(c_int)  :: n
       integer(c_int)  :: iptr(*), jcn(*), cperm(*), iperm(*)
     end subroutine qrm_metis
  end interface

  integer                     :: i, idx, cnt, tmp, alen
  type(_qrm_spmat_type)       :: ata_graph
  integer, allocatable        :: iperm(:)
  ! error management
  integer                     :: err_act
  character(len=*), parameter :: name='qrm_do_metis'

  call qrm_err_act_save(err_act)

  call _qrm_ata_graph(graph, ata_graph)
  __QRM_CHECK_RET(name,'qrm_ata_graph',9999)

  call qrm_aalloc(iperm, graph%n)
  __QRM_CHECK_RET(name,'qrm_ata_graph/alloc',9999)

  call qrm_metis(graph%n, ata_graph%iptr, ata_graph%jcn, cperm, iperm)

  call qrm_adealloc(iperm)
  call _qrm_spmat_destroy(ata_graph, all=.true.)
  __QRM_CHECK_RET(name,'dealloc/destroy',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return
end subroutine _qrm_do_metis
