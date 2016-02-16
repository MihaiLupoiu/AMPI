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
!> @file qrm_fdata_mod.F90
!! This file contains the module which holds the data types used during the factorization.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains the definition of all the data related to the
!! factorization phase. 
!! 
module _qrm_fdata_mod

#if defined(_OPENMP)
  use omp_lib
#endif

  integer, parameter :: qrm_busy_       = -huge(0)
  integer, parameter :: qrm_ready_      = 0 ! this can never be changed
  integer, parameter :: qrm_activable_  = 1
  integer, parameter :: qrm_active_     = 2
  integer, parameter :: qrm_factorized_ = 3 ! the front has been completely factorized
  integer, parameter :: qrm_free_       = 4 ! the front has been freed
  integer, parameter :: qrm_done_       = huge(0)

 !> @brief This type defines a data structure containing all the data related to a front.
 type _qrm_front_type
     !> the front number
     integer                 :: num=0
     !> the number of rows in the front
     integer                 :: m=0
     !> the number of columns in the front
     integer                 :: n=0
     !> the number of pivots (fully assembled variables in the front
     integer                 :: npiv=0
     !> the list of row indices in the front
     integer, allocatable    :: rows(:)
     !> the list of column indices in the front
     integer, allocatable    :: cols(:)
     !> The pointer to the beginning of rows of coefficients from
     !> the original matrix. 
     integer, allocatable    :: aiptr(:)
     !> The columns indices of coefficients from
     !> the original matrix. 
     integer, allocatable    :: ajcn(:)
     !> The values of coefficients from original matrix. These
     !! coefficients will be assembled into
     !! the front at the moment of its activation
     _qrm_data, allocatable :: aval(:)
     !> The initial data from the problem matrix is stored in CSR format. This
     !  integer contains the number of rows in this front from the original matrix
     integer                 :: anrows
     !> this integer array of size front%n contains a mapping between the front's
     !! column indices and father front's column indices, i.e. colmap(i)=k means
     !! that column i of front goes into k-th column of its parent
     integer, allocatable    :: colmap(:)
     !> this integer array is of size equal to the number of rows in
     !! the contribution block produced by this front and contains a
     !! mapping between the rows of the CB and the front of the
     !! parent, i.e., rowmap(i)=k means that i-th row of this
     !! contribution block goes into row k of the front of the father.
     integer, allocatable    :: rowmap(:)
     !> this array of size front%n will ultimately define the staircase
     !! structure of the front. stair(i)=k means that the lowest
     !! element in column i of the front is in row k
     integer, allocatable    :: stair(:)
     !> a 2D array containing the front. This area will be dynamically
     !! allocated at front assembly time
     _qrm_data, allocatable :: front(:,:)
     !> after the front is factorized, the corresponding elemnts of R will be
     !! copied into this array
     _qrm_data, allocatable :: r(:)
     !> after the front is factorized, the corresponding elemnts of H will be
     !! copied into this array
     _qrm_data, allocatable :: h(:)
     !> in the case where the front is completely factorized but assembly
     !! has not yet started, memory can be saved by copying the contribution block
     !! in here are freeing the front matrix
     _qrm_data, allocatable :: cb(:)
     !> this array holds the tau values resulting from the factorization. Its size is
     !! equal to the number of eliminations done in the front (i.e., the number of columns in H)
     _qrm_data, allocatable :: tau(:)
     !> this array holds the T matrices resulting from each panel factorization. If
     !! ib is the inner blocking size, then T is (ib,n).
     _qrm_data, allocatable :: t(:,:)
     !> this 2-D array will be used to gather the elements of a vaector concerned by the
     !! application of Q' or Q
     _qrm_data, allocatable :: b(:,:)
     !> this integer array contains a progress table for the factorization.
     !! Each element of the table tracks the status of a block-column like this:
     !! - qrm_busy_: means that the column is busy. This is to prevent that multiple
     !!            threads work on the same column
     !! - -nc: where nc is the number of contributions coming from children
     !!        to this block column. This value is incremented each time a contribution
     !!        is assembled into the block column until...
     !! - qrm_ready_: ...the column is fully assembled and it's ready to be treated
     !! - >0: the index of the panel wrt which the column is up to date
     !! - qrm_done_: nothing left to do on this column
     integer, allocatable    :: ptable(:)
     !> this integer holds the status of the front. The status is coded like this:
     !! - qrm_busy_: means the front is being activated/cleared. This is mostly to prevent
     !!   concurrent work of the same front
     !! - - k, where k is the number of children of the front. This is the initial front
     !!   state before factorization starts. Each time one of its children is activated,
     !!   the status is incremented by 1 until...
     !! - qrm_ready_: ...it reaches the "ready for activation state"
     !! - qrm_activable_: once the front is ready it may be pushed
     !!   onto the "ready" queue in which case it becomes activable
     !! - qrm_active_: means that the front has been activated and thus we can start assemblying
     !!   values into it
     !! - qrm_done_: done. The front has been factorized
     integer                 :: status
     !> this field holds the id of the thread that owns the front. A thread becomes the owner
     !! of a front the moment it activates it
     integer                 :: owner 
     !> an integer containing the index of the block-column that has been factorized last 
     integer                 :: lastpnl=0
     !> the block-column sizes for this panel
     integer                 :: nb, ib
     !> the number of block-columns in this front
     integer                 :: nc
     !> the number of block-panel to be factorized
     integer                 :: np
     !> the number of eliminations to be performed on the front.
     !!  This basically corresponds to min(front%m,front%n)
     integer                 :: ne 
     !> the number of entries in R and H in this front
     integer(kind=8)         :: rsize, hsize
#if defined(_OPENMP)
     !> an openmp lock to prevent simultaneous access to critical data in the front
     integer(kind=omp_lock_kind) :: lock 
#endif
  end type _qrm_front_type



  !> @brief The data structure meant to store all the results of the
  !! factorization phase.
  type _qrm_fdata_type
     !> an integer containing the number of fronts assigned to the process
     integer                             :: nfronts=0
     !> just a counter to count the number of fronts done during the
     !> factorization. It is used to stop the factorization operation
     integer                             :: done=0
     !> an array of type _qrm_front_type containing the list of fronts
     !! assigned to the process
     type(_qrm_front_type), allocatable  :: front_list(:)
#if defined(_OPENMP)
     !> a lock to synchronize all dangerous accesses to a qrm_fdata_type
     integer(kind=omp_lock_kind)         :: lock 
#endif
     !> it is set to .true. if the factorization is done (with success), to .false. otherwise
     logical                             :: ok=.false.
  end type _qrm_fdata_type

  interface qrm_fdata_destroy
     module procedure _qrm_fdata_destroy
  end interface qrm_fdata_destroy



contains

  ! TODO: add copy!

  !> @brief Frees a qrm_front_type instance
  !!
  !! @param[in,out] qrm_front the data to be freed
  !!
  subroutine _qrm_front_destroy(qrm_front)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(_qrm_front_type) :: qrm_front

      ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='qrm_front_destroy'
    
    call qrm_err_act_save(err_act)
    
    call qrm_adealloc(qrm_front%aiptr)
    call qrm_adealloc(qrm_front%ajcn)
    call qrm_adealloc(qrm_front%aval)
    call qrm_adealloc(qrm_front%front)
    call qrm_adealloc(qrm_front%r)
    call qrm_adealloc(qrm_front%h)
    call qrm_adealloc(qrm_front%rows)
    call qrm_adealloc(qrm_front%cols)
    call qrm_adealloc(qrm_front%rowmap)
    call qrm_adealloc(qrm_front%colmap)
    call qrm_adealloc(qrm_front%stair)
    call qrm_adealloc(qrm_front%ptable)
    call qrm_adealloc(qrm_front%cb)
    call qrm_adealloc(qrm_front%tau)
    call qrm_adealloc(qrm_front%t)
    call qrm_adealloc(qrm_front%b)
#if defined(_OPENMP)
    call omp_destroy_lock(qrm_front%lock)
#endif
    __QRM_CHECK_RET(name,'qrm_adelloc',9999)

    qrm_front%m = 0
    qrm_front%n = 0
    
    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return
    
  end subroutine _qrm_front_destroy
  
  
  
  !> @brief Destroys a @link _qrm_fdata_type @endlink instance
  !! 
  !! @param[in,out] qrm_fdata the instace to be freed
  !! 
  subroutine _qrm_fdata_destroy(qrm_fdata)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none
    
    type(_qrm_fdata_type) :: qrm_fdata
    integer :: i
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='qrm_fdata_destroy'
    
    call qrm_err_act_save(err_act)
    
    
    if(allocated(qrm_fdata%front_list)) then
       do i=1, qrm_fdata%nfronts
          call _qrm_front_destroy(qrm_fdata%front_list(i))
       end do
       deallocate(qrm_fdata%front_list)
    end if
    __QRM_CHECK_RET(name,'qrm_front_destroy',9999)

    qrm_fdata%done    = 0
    qrm_fdata%nfronts = 0
    qrm_fdata%ok      = .false.
#if defined(_OPENMP)
    ! FIXME: how to check whether a lock is initialized or not?
    ! call omp_destroy_lock(qrm_fdata%lock)
#endif
    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return
    
  end subroutine _qrm_fdata_destroy
  


end module _qrm_fdata_mod
