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
!> @file qrm_do_subtree.F90
!! This file contains the routines that take care of factorizing sequentially an entire subtree
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine does the sequential factorization of an
!! entire subtree
!!
!! @param[in,out] qrm_mat the whole problem data structure
!!
!! @param[in] fnum the id of the root of the subtree
!! 
!! @param[in,out] flops a counter for the flops which is updated
!! 
subroutine _qrm_do_subtree(qrm_mat, fnum, flops)

  use qrm_mem_mod
  use _qrm_spmat_mod
  use _qrm_fdata_mod
  use qrm_adata_mod
  use _qrm_utils_mod
  use qrm_common_mod
  use qrm_sort_mod
  use _qrm_factorization_mod, protect => _qrm_do_subtree
  implicit none
  
  type(_qrm_spmat_type), target :: qrm_mat
  integer                        :: fnum
  real(kind(1.d0))    :: flops

  type(_qrm_front_type), pointer :: front, cfront
  type(qrm_adata_type), pointer   :: adata
  type(_qrm_fdata_type), pointer :: fdata
  
  integer, allocatable    :: iwork(:)
  _qrm_data, allocatable  :: rwork(:)
  integer                 :: node, p, c, i, j, cnt

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_do_subtree'

  call qrm_err_act_save(err_act)

  ! simplify 
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  node = fnum
  front => fdata%front_list(node)
  ! write(*,*)'--->',fnum,20*front%n*front%nb
  call qrm_aalloc(rwork, adata%rc(fnum)*qrm_mat%icntl(qrm_nb_))
  
  cnt = 0
  call qrm_aalloc(iwork, qrm_mat%n)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)

  subtree: do

     front => fdata%front_list(node)
     
     if (front%status .eq. qrm_ready_) then
        ! the front is ready to be assembled and factorized

        ! assemble
        call _qrm_init_front(qrm_mat, node, .false., work=iwork)
        __QRM_CHECK_RET(name,'qrm_init_front',9999)

        ! clean
        do p = adata%childptr(node), adata%childptr(node+1)-1
           c = adata%child(p)
           cfront => fdata%front_list(c)
           call _qrm_clean_front(qrm_mat, cfront%num)
           __QRM_CHECK_RET(name,'qrm_clean_front',9999)
           cnt = cnt+1
           cfront%status = qrm_done_
        end do

        ! factorize
        call sequential_fct(front, flops)
        __QRM_CHECK_RET(name,'sequential_fct',9999)

        front%status = qrm_factorized_

        if(node .eq. fnum) exit subtree  ! we're done
        
        node = adata%parent(node)
        if(node .eq. 0) exit subtree

     else
        ! go down the subtree
        node = adata%child(adata%childptr(node+1)+front%status)
        cycle subtree
     end if


  end do subtree

  !$ call omp_set_lock(qrm_mat%fdata%lock)
  qrm_mat%fdata%done = qrm_mat%fdata%done+cnt
  !$ call omp_unset_lock(qrm_mat%fdata%lock)


  call qrm_adealloc(rwork)
  call qrm_adealloc(iwork)
  __QRM_CHECK_RET(name,'qrm_adealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

contains








!==========================================================================================
!==========================================================================================
  subroutine sequential_fct(front, flops)
    ! this routine does the sequential factorization of a front

    implicit none
    
    type(_qrm_front_type)   :: front
    real(kind(1.d0))        :: flops

    integer                 :: ib, i, j, m, n, k, thn, info

    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='sequential_fct'
    
    call qrm_err_act_save(err_act)


    call qrm_get(qrm_mat, 'qrm_ib', ib)

    call qrm_arealloc(rwork, front%n*front%nb)
    __QRM_CHECK_RET(name,'qrm_aalloc',9999)

    do i=1, front%ne, ib
       n = min(ib, front%n-i+1)
       m = max(front%stair(min(i+n-1,front%ne))- i+1,0) ! the height of the panel
       k = min(m,n)
       j = i+n

       call _xgeqrf(m, n, front%front(i,i), front%m, front%tau(i), rwork(1), &
            & front%nb**2, info)
       call _xlarft( 'f', 'c', m, k, front%front(i,i), front%m, front%tau(i), &
            & front%t(1,i), front%ib )
       
       flops = flops+qrm_count_flops(m, n, n, 'panel')
       
       n = front%n-i-n+1
       if(n.gt.0) call _xlarfb('l', 't', 'f', 'c', m, n, k, front%front(i,i), &
            & front%m, front%t(1,i), front%ib, &
            & front%front(i,j), front%m, rwork(1), front%n)
       flops = flops+qrm_count_flops(m, n, k, 'update')

       ! do j=i+n,front%n, front%nb
          ! n = min(front%nb,front%n-j+1)
          ! call _xlarfb('l', 't', 'f', 'c', m, n, k, front%front(i,i), &
               ! & front%m, front%t(1,1,(i+ib-1)/ib), front%nb, &
               ! & front%front(i,j), front%m, rwork(1), front%nb)
          ! flops = flops+qrm_count_flops(m, n, k, 'update')
       ! end do

    end do

    __QRM_CHECK_RET(name,'qrm_adealloc',9999)

    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    
    return
    
  end subroutine sequential_fct


end subroutine _qrm_do_subtree
