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
!> @file qrm_do_scotch.F90
!! This file contains the routine that computes a SCOTCH permutation of the input matrix
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
!! <em>"SCOTCH 5.1 User's guide</em>.
!! Technical report, LaBRI, September 2008.
!! F. Pellegrini.
!!
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_scotch(graph, cperm)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_scotch
  use qrm_mem_mod
  use iso_c_binding
  implicit none

#if defined(have_scotch)
  include "scotchf.h"
#endif

  type(_qrm_spmat_type) :: graph
  integer                :: cperm(:)
  
#if defined(have_scotch)
  integer                :: i, info, cblknbr
  type(_qrm_spmat_type) :: ata_graph
  integer, allocatable   :: iperm(:)
  real(kind(1.d0))       :: grafdat(scotch_graphdim), stradat(scotch_stratdim), orderdat(scotch_orderdim)
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_do_scotch'

  call qrm_err_act_save(err_act)

  call _qrm_ata_graph(graph, ata_graph)
  __QRM_CHECK_RET(name,'qrm_ata_graph',9999)
  
  info = 0
  call scotchfgraphinit(grafdat, info)
  call scotchfstratinit(stradat, info)
  if(info .ne. 0) then
     call qrm_err_push(19,name)
     goto 9999
  end if

  call scotchfgraphbuild(grafdat, 1, ata_graph%n, ata_graph%iptr(1), &
       & ata_graph%iptr(2), ata_graph%iptr, ata_graph%iptr, ata_graph%nz, &
       & ata_graph%jcn, ata_graph%jcn, info)
  if(info .ne. 0) then
     call qrm_err_push(19,name)
     goto 9999
  end if

  call scotchfgraphorder(grafdat, stradat, grafdat, cperm, cblknbr, &
       & grafdat, grafdat, info)
  if(info .ne. 0) then
     call qrm_err_push(19,name)
     goto 9999
  end if

  call scotchfgraphexit(grafdat)
  call scotchfstratexit(stradat)

  call _qrm_spmat_destroy(ata_graph, all=.true.)
  __QRM_CHECK_RET(name,'qrm_spmat_destroy',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
#endif

  return
end subroutine _qrm_do_scotch
