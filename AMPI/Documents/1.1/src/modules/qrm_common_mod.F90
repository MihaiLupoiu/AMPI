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
!> @file qrm_common_mod.F90
!! this module contains generic interfaces for all the untyped routines.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This module contains the interfaces of all non-typed routines
module qrm_common_mod
  use qrm_const_mod


  !> @brief Generic interface for the @link ::qrm_print_nsteps_tree @endlink,
  !! @link ::qrm_print_elim_tree @endlink and @link
  !! ::qrm_print_asm_tree @endlink routines
  interface qrm_print_tree
     subroutine qrm_print_nsteps_tree(file, adata, weight)
       use qrm_adata_mod
       type(qrm_adata_type)       :: adata
       character                  :: file*(*)
       real(kind(1.d0)), optional :: weight(:)
     end subroutine qrm_print_nsteps_tree
     subroutine qrm_print_elim_tree(file, parent, n)
       character :: file*(*)
       integer   :: parent(:)
       integer   :: n
     end subroutine qrm_print_elim_tree
     subroutine qrm_print_asm_tree(file, parent, rc, n)
       character :: file*(*)
       integer   :: parent(:), rc(:)
       integer   :: n
     end subroutine qrm_print_asm_tree
  end interface

  !> @brief Generic interface for the @link ::qrm_postorder @endlink routine
  interface qrm_postorder
     subroutine qrm_postorder(parent, n, porder, weight)
       integer           :: n
       integer           :: parent(:), porder(:)
       integer, optional :: weight(:)
     end subroutine qrm_postorder
  end interface

  !> @brief Generic interface for the @link ::qrm_amalg_tree @endlink routine
  interface
     subroutine qrm_amalg_tree(n, parent, rc, porder, nvar, min_var, fill_thresh)
       integer          :: n, min_var
       integer          :: parent(:), rc(:), porder(:), nvar(:)
       real(kind(1.d0)) :: fill_thresh
     end subroutine qrm_amalg_tree
  end interface

  !> @brief Generic interface for the @link ::qrm_compress_data @endlink routine
  interface qrm_compress_data
     subroutine qrm_compress_data(adata, porder, parent, rc, stair, n)
       use qrm_mem_mod
       use qrm_adata_mod
       type(qrm_adata_type) :: adata
       integer              :: porder(:), parent(:), rc(:), stair(:)
       integer              :: n
     end subroutine qrm_compress_data
  end interface

  !> @brief Generic interface for the @link ::qrm_reorder_tree @endlink routine
  interface qrm_reorder_tree
     subroutine qrm_reorder_tree(adata)
       use qrm_adata_mod
       type(qrm_adata_type) :: adata
     end subroutine qrm_reorder_tree
  end interface

  !> @brief Generic interface for the @link ::qrm_prnt_iarray
  !> @endlink, @link ::qrm_prnt_sarray @endlink, @link
  !> ::qrm_prnt_darray @endlink, @link ::qrm_prnt_carray @endlink and
  !> @link ::qrm_prnt_zarray @endlink routines
  interface qrm_prnt_array
     subroutine qrm_prnt_iarray(a, lab, unit)
       integer   :: a(:)
       character :: lab*(*)
       integer, optional :: unit
     end subroutine qrm_prnt_iarray
     subroutine qrm_prnt_sarray(a, lab, unit)
       real(kind(1.e0))  :: a(:)
       character         :: lab*(*)
       integer, optional :: unit
     end subroutine qrm_prnt_sarray
     subroutine qrm_prnt_darray(a, lab, unit)
       real(kind(1.d0))  :: a(:)
       character         :: lab*(*)
       integer, optional :: unit
     end subroutine qrm_prnt_darray
     subroutine qrm_prnt_carray(a, lab, unit)
       complex(kind(1.e0)) :: a(:)
       character           :: lab*(*)
       integer, optional   :: unit
     end subroutine qrm_prnt_carray
     subroutine qrm_prnt_zarray(a, lab, unit)
       complex(kind(1.d0)) :: a(:)
       character           :: lab*(*)
       integer, optional   :: unit
     end subroutine qrm_prnt_zarray
  end interface

  !> @brief Generic interface for the @link ::qrm_swtime @endlink routine
  interface qrm_swtime
     function qrm_swtime() bind(c, name='qrm_swtime')
       real(kind(1.d0)) :: qrm_swtime
     end function qrm_swtime
  end interface

  !> @brief Generic interface for the @link ::qrm_uwtime @endlink routine
  interface qrm_uwtime
     function qrm_uwtime() bind(c, name='qrm_uwtime')
       real(kind(1.d0)) :: qrm_uwtime
     end function qrm_uwtime
  end interface

  !> @brief Generic interface for the @link ::qrm_msleep @endlink routine
  interface qrm_msleep
     subroutine qrm_msleep(n) bind(c, name='qrm_msleep')
       use iso_c_binding
       integer(c_int), value :: n
     end subroutine qrm_msleep
  end interface


  !> @brief Generic interface for the @link ::qrm_check_cperm @endlink routine
  interface qrm_check_cperm
     subroutine qrm_check_cperm(cperm, n)
       integer :: cperm(:)
       integer :: n
     end subroutine qrm_check_cperm
  end interface

#if defined(hwloc)
  !> @brief Generic interface for the @link ::qrm_hwloc_bind @endlink routine
  interface
     subroutine qrm_hwloc_bind(id) bind(c, name='qrm_hwloc_bind')
       use iso_c_binding
       integer(c_int), value  :: id
     end subroutine qrm_hwloc_bind
  end interface

  !> @brief Generic interface for the @link ::qrm_hwloc_topo @endlink routine
  interface
     subroutine qrm_hwloc_topo(nnodes, topo) bind(c, name='qrm_hwloc_topo')
       use iso_c_binding
       integer(c_int) :: topo(nnodes,*)
     end subroutine qrm_hwloc_topo
  end interface

  !> @brief Generic interface for the @link ::qrm_hwloc_info @endlink routine
  interface
     subroutine qrm_hwloc_info(ncores, nnodes, cnode) bind(c, name='qrm_hwloc_info')
       use iso_c_binding
       integer(c_int) :: ncores, nnodes, cnode
     end subroutine qrm_hwloc_info
  end interface

#endif

  
  !> @brief Generic interface for the
  !! @link ::qrm_count_realflops @endlink
  !! @link ::qrm_count_pureflops @endlink
  interface qrm_count_flops
     module procedure qrm_count_realflops, qrm_count_pureflops
     ! function qrm_count_realflops(m, n, k, op)
       ! real(kind(1.d0)) :: qrm_count_realflops
       ! integer :: m, k, n
       ! character :: op*(*)
     ! end function qrm_count_realflops
     ! function qrm_count_pureflops(stair, n, j, nb)
       ! implicit none
       ! real(kind(1.d0)) :: qrm_count_pureflops
       ! integer :: stair(:)
       ! integer :: n, nb, j
     ! end function qrm_count_pureflops
  end interface

  interface qrm_set
     module procedure  qrm_gseti
  end interface

  interface qrm_get
     module procedure  qrm_ggeti, qrm_ggetii
  end interface

contains

  subroutine qrm_gseti(string, ival)
    use qrm_string_mod
    use qrm_error_mod
    use qrm_mem_mod
    implicit none

    character(len=*)     :: string
    integer :: ival

    character(len=len(string)) :: istring
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='qrm_gseti'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_eunit') .eq. 1) then
       qrm_eunit = ival
    else if(index(istring,'qrm_ounit') .eq. 1) then
       qrm_ounit = ival
    else if(index(istring,'qrm_dunit') .eq. 1) then
       qrm_dunit = ival
    else if(index(istring,'qrm_max_mem') .eq. 1) then
       qrm_max_mem(0) = ival
    else if(index(istring,'qrm_tot_mem') .eq. 1) then
       qrm_tot_mem(0) = ival
    else if(index(istring,'qrm_exact_mem') .eq. 1) then
       qrm_exact_mem = ival
    else if(index(istring,'qrm_error_action') .eq. 1) then
       call qrm_err_act_set(ival)
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if


    return
  end subroutine qrm_gseti

  subroutine qrm_ggeti(string, ival)
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    character(len=*) :: string
    integer          :: ival

    integer(kind=8)  :: iival

    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='qrm_ggeti'

    call qrm_err_act_save(err_act)

    call qrm_ggetii(string, iival)
    __QRM_CHECK_RET(name,'qrm_ggetii',9999)
    ival = iival

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return
  end subroutine qrm_ggeti

  subroutine qrm_ggetii(string, iival)
    use qrm_string_mod
    use qrm_error_mod
    use qrm_mem_mod
    implicit none

    character(len=*) :: string
    integer(kind=8)  :: iival

    character(len=len(string)) :: istring
    integer                    :: ival
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='qrm_ggeti'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_error') .eq. 1) then
       call qrm_err_get(ival)
       iival = ival
    else if(index(istring,'qrm_error_action') .eq. 1) then
       iival = qrm_err_act
    else if(index(istring,'qrm_max_mem') .eq. 1) then
       iival = qrm_max_mem(0)
    else if(index(istring,'qrm_tot_mem') .eq. 1) then
       iival = qrm_tot_mem(0)
    else if(index(istring,'qrm_ounit') .eq. 1) then
       iival = qrm_ounit
    else if(index(istring,'qrm_eunit') .eq. 1) then
       iival = qrm_eunit
    else if(index(istring,'qrm_dunit') .eq. 1) then
       iival = qrm_dunit
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return
  end subroutine qrm_ggetii




!> @brief Used for counting the actual flops
function qrm_count_realflops(m, n, k, op)
  real(kind(1.d0)) :: qrm_count_realflops
  integer :: m, k, n
  character :: op*(*)
  
  real(kind(1.d0)) :: rk, rm, rn
  
  rk = real(k, kind(1.d0))
  rm = real(m, kind(1.d0))
  rn = real(n, kind(1.d0))
  
  qrm_count_realflops = 0.d0
  
  select case(op)
  case('panel')
     if(m .gt. k) then
        qrm_count_realflops = 2*rk*rk*(rm-rk/3.d0)
     else
        qrm_count_realflops = 2*rm*rm*(rk-rm/3.d0)
     end if
     ! qrm_count_realflops = qrm_count_realflops + rk*rk*(rm-rk/3.d0)
  case('update')
     qrm_count_realflops = rk*rn*(4*rm-rk)
  end select
  
  return
end function qrm_count_realflops


!> @brief Used for counting the real flops (i.e., it ignores the zeros
!> included by the blocking)
function qrm_count_pureflops(stair, n, j, nb)
  implicit none
  real(kind(1.d0)) :: qrm_count_pureflops
  integer :: stair(:)
  integer :: n, nb, j
  
  integer :: i
  
  qrm_count_pureflops=0
  do i=j, min(j+nb-1, n)
     qrm_count_pureflops = qrm_count_pureflops+(stair(i)-i+1)*(3 + 4*(n-i))
  end do
  
  return
end function qrm_count_pureflops

end module qrm_common_mod
