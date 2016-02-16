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
!> @file qrm_do_ordering.F90
!! This file contains the routine computes (through different methods) a column permutation
!! of the input matrix in order to reduce fill-in
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This routine computes (through different methods) a column permutation
!! of the input matrix in order to reduce fill-in

!> Supported methods are, currently, COLAMD, SCOTCH and METIS. The user may also
!! provide his own permutation, in which case a check on its validity is done.
!!
!! @param[in] graph    the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm   the new column order
!!
!! @param[in] cperm_in the permutation eventually provided by the user
!!
subroutine _qrm_do_ordering(graph, cperm, cperm_in)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_ordering
  use qrm_mem_mod
  use qrm_const_mod
  use qrm_common_mod
  implicit none

  type(_qrm_spmat_type) :: graph
  integer                :: cperm(:)
  integer, pointer       :: cperm_in(:)

  integer                :: iord=0
  integer                :: i
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_do_ordering'

  call qrm_err_act_save(err_act)

  call qrm_get(graph, 'qrm_ordering', iord)

  if(iord .eq. qrm_auto_) then
     iord = qrm_choose_ordering()
  end if

  
  select case(iord)
  case(qrm_natural_)
     ! perm has to be allocated and built explicitly because
     ! we are going to use an equivalent permutation
     do i=1, graph%n
        cperm(i) = i
     end do
  case(qrm_given_)
     ! in this case we just need to check that
     ! the given permutation makes sense
     if(.not. associated(cperm_in)) then
        call qrm_err_push(8, name)
        goto 9999
     else
        call qrm_check_cperm(cperm_in, graph%n)
        __QRM_CHECK_RET(name,'qrm_check_perm',9999)
        cperm(1:graph%n) = cperm_in(1:graph%n)
     end if
  case(qrm_colamd_)
     ! COLAMD
#if defined (have_colamd)
     call _qrm_do_colamd(graph, cperm)
#else
     call qrm_err_push(16, 'qrm_do_ordering',aed='colamd')
     goto 9999
#endif
  case(qrm_metis_)
     ! METIS
#if defined (have_metis)
     call _qrm_do_metis(graph, cperm)
#else
     call qrm_err_push(16, 'qrm_do_ordering',aed='metis')
     goto 9999
#endif
  case(qrm_scotch_)
     ! SCOTCH
#if defined (have_scotch)
     call _qrm_do_scotch(graph, cperm)
#else
     call qrm_err_push(16, 'qrm_do_ordering',aed='scotch')
     goto 9999
#endif
  case default
     call qrm_err_push(9, 'qrm_do_ordering',ied=(/iord, 0, 0, 0, 0/))
     goto 9999
  end select

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

contains

  function qrm_choose_ordering()
    ! Function: qrm_choose_ordering
    ! This function sets the ordering for the case where an automatic
    ! choice is requested
    !

    integer :: qrm_choose_ordering
    integer :: iord

    iord = 1

#if defined(ccolamd)    
    iord = qrm_colamd_
#endif

#if defined(metis)    
    iord = qrm_metis_
#endif

#if defined(scotch)    
    iord = qrm_scotch_
#endif

    qrm_choose_ordering = iord
    return
    
  end function qrm_choose_ordering


end subroutine _qrm_do_ordering
