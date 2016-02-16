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
!> @file qrm_task_mod.F90
!! This file contains a module that implements all the task handling facilities
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

#include "qrm_common.h"

!> @brief This module contains the definition of a task type that is used
!! for scheduling tasks during the factorization and solve, and the
!! associated methods

!> The parallelism in qr_mumps is based on a dynamic execution of
!! computational tasks. These computational tasks and their
!! dependencies are defined by a set of rules that are handled in the
!! factorization or solve routines. Once all the dependencies of a
!! task are satisfied, the task is ready for being executed. The
!! execution of all the ready tasks is managed by the methods in this
!! module. Specifically, a task-queue handle is shared among all the
!! threads in a family; each task queue groups a set of queues, one
!! for each thread in the family. Each thread can access its queue
!! through the handle of the family it belongs to and its rank.
!!
!! Usage example:
!! @code
!!  
!! type(qrm_task_queue_handle) :: h
!! type(qrm_task_type)         :: task
!! integer                     :: iam
!! logical                     :: got
!!  
!! ! The handle is always shared!!!
!! !$omp parallel shared(h)
!!
!! ! initialize the family queues
!! call qrm_init_task_queue(h)
!!
!! iam = omp_get_thread_num()
!! ! init the task using its constructor
!! task = qrm_task_type(qrm_act_, iam, 0, 0, .true.)
!!  
!! ! push the task on my own queue
!! call qrm_sched_task(h, task, 't', iam)
!!  
!! ! push the task on the queue
!! got = qrm_get_task(h, task)
!!  
!! if(got) then
!!    ! execute the task
!! end if
!!  
!! !$omp end parallel
!!
!! @endcode
!!
!! In the small example above, each thread pushes on its own queue a
!! task corresponding to the activation of the front whose ID is the
!! same as the thread's rank. Note that the task is bound, i.e., it
!! can only be executed by the thread associated to the queue where
!! the task was pushed. Then each thread pops a task and executes it.
!!
module qrm_task_mod

#if defined (_OPENMP)  
  use omp_lib
#endif
  use qrm_common_mod
  use qrm_error_mod 

  ! parameter variables follow the convention of having an underscore at the end
  ! 0=exit, 1=panel, 2=update, 3=activate, 4=assemble, 5=free, 6=clean
  integer, parameter :: qrm_task_exit_  = 0      ! exit
  integer, parameter :: qrm_task_pnl_   = 1      ! panel
  integer, parameter :: qrm_task_upd_   = 2      ! update
  integer, parameter :: qrm_task_act_   = 3      ! activate
  integer, parameter :: qrm_task_asm_   = 4      ! assemble
  integer, parameter :: qrm_task_free_  = 5      ! free
  integer, parameter :: qrm_task_cln_   = 6      ! clean
  integer, parameter :: qrm_task_app_   = 7      ! apply (for Q or Q')
  integer, parameter :: qrm_task_sol_   = 8      ! solve (for R or R')

  !> @brief This type defines a computational task
  type qrm_task_type
     !> The type of action to execute
     integer :: id
     !> The id of the related front
     integer :: front
     !> The panel number
     integer :: pnl
     !> The column number
     integer :: col
     !> Whether the task is bound to a thread or not
     logical :: bind
  end type qrm_task_type


  !> The max size of a task queue attached to a thread
  integer, parameter   :: max_tasks = 300
  
  integer, private :: qrm_task_thn, qrm_task_nth
  !$omp threadprivate(qrm_task_thn, qrm_task_nth)

  !> @brief This type defines the task queue attached to a thread
  type qrm_task_queue
     integer :: h, t, n
     type(qrm_task_type) :: q(max_tasks)
#if defined(_OPENMP)
     integer(kind=omp_lock_kind) :: lock 
#endif
  end type qrm_task_queue
  
  !> @brief This type defines the handle for the queues attached to a
  !! family of threads
  type qrm_task_queue_handle
     type(qrm_task_queue), allocatable :: queues(:)
     integer, allocatable              :: proxy_list(:,:)
     integer                           :: stolen=0, ntsk=0, ncores, nnodes, cnode
  end type qrm_task_queue_handle

contains

#ifndef NOLOC
  !> @brief Inititalizes a set of queues attached to a family of
  !! threads referenced through the handle h
  !!
  !! @param[in] h The handle that references the queues
  !! 
  subroutine qrm_init_task_queue(h)
    implicit none
    type(qrm_task_queue_handle) :: h

    integer :: i
    
#if defined (_OPENMP)
    qrm_task_nth = omp_get_num_threads()
    qrm_task_thn = omp_get_thread_num()
#else
    qrm_task_nth = 1
    qrm_task_thn = 0
#endif

#if defined (hwloc)
    !$omp master
    call qrm_hwloc_info(h%ncores, h%nnodes, h%cnode)
    !$omp end master
    !$omp barrier
    ! bind the thread but only if we are at level 1
    !$ if (omp_get_level() .eq. 1) call qrm_hwloc_bind(mod(qrm_task_thn,h%ncores))
#else
    !$omp master
    h%ncores = qrm_task_nth
    h%nnodes = 1
    h%cnode  = qrm_task_nth
    !$omp end master
#endif

    !$omp master
    call qrm_task_proximity(h)

    allocate(h%queues(0:qrm_task_nth-1))
    do i=0, qrm_task_nth-1
       h%queues(i)%h = 0
       h%queues(i)%t = 0
       h%queues(i)%n = 0
       h%queues(i)%q = qrm_task_type(0, 0, 0, 0, .false.)
       !$ call omp_init_lock(h%queues(i)%lock)
    end do
    !$omp end master
    !$omp barrier !TODO maybe remove
    return
  end subroutine qrm_init_task_queue


  !> @brief Defines the order in which queues have to be visited by
  !! each thread
  !! 
  subroutine qrm_task_proximity(h)
    use iso_c_binding
    implicit none
    type(qrm_task_queue_handle) :: h

    integer :: i, j, ii, jj, cnt, id, c, n
    integer, allocatable :: topo(:,:), core_to_node(:)

#if defined (hwloc)
    call qrm_hwloc_info(h%ncores, h%nnodes, h%cnode)
    ! write(*,'("ncores: ",i2,"    nnodes: ",i2,"    cnode: ",i2)')h%ncores,h%nnodes,h%cnode
    allocate(topo(0:h%cnode-1,0:h%nnodes-1))
    allocate(h%proxy_list(0:qrm_task_nth-1,0:qrm_task_nth-1))
    call qrm_hwloc_topo(h%nnodes, topo)
#else
    h%cnode = min(qrm_task_nth,1)
    h%ncores = qrm_task_nth
    h%nnodes = (qrm_task_nth-1)/h%cnode+1
    allocate(topo(0:h%cnode-1,0:h%nnodes-1))
    allocate(h%proxy_list(0:qrm_task_nth-1,0:qrm_task_nth-1))
    cnt = 0
    do j=0, h%nnodes-1
       do i=0, h%cnode-1
          topo(i,j)=cnt
          cnt = cnt+1
       end do
    end do
#endif

    allocate(core_to_node(0:h%nnodes*h%cnode-1))
    do j=0, h%nnodes-1
       do i=0, h%cnode-1
          core_to_node(topo(i,j)) = j
       end do
    end do

    do id=0, qrm_task_nth-1
       cnt = 0
       ! add thread id in first position
       c = 0
       do
          if(c*h%ncores +id .ge. qrm_task_nth) exit
          h%proxy_list(id,cnt) = c*h%ncores+id
          cnt = cnt+1
          c = c+1
       end do

       ! add all the other threads
       do j=0, h%nnodes-1
          n = core_to_node(mod(id,h%ncores))
          n = mod(n+j,h%nnodes)
          do i = 0, h%cnode-1
             c = 0
             do
                if(c*h%ncores +topo(i,n) .ge. qrm_task_nth) exit
                if(c*h%ncores +topo(i,n) .eq. id) exit
                h%proxy_list(id,cnt) = c*h%ncores+topo(i,n)
                cnt = cnt+1
                c = c+1
             end do
          end do
       end do
     
    end do

    deallocate(topo, core_to_node)

    return
  end subroutine qrm_task_proximity
  
  !> @brief Pushes a task on a queue
  !!
  !! @param[in,out] h the handle to the set of queues
  !!
  !! @param[in] tsk the task to be pushed
  !!
  !! @param[in] pol the policy. It can be either 'h', in which case
  !!                the task will be pushed on the head of the queue, or 't', in
  !!                which case it will be pushed on the tail of the queue
  !!
  !! @param[in] q the queue in the set where to push the task. This
  !!              can be used to improve data locality, e.g., a task can be pushed
  !!              on the queue attached to the thread which owns the data to be
  !!              accessed by the task itself
  function qrm_sched_task(h, tsk, pol, q)
    implicit none
    type(qrm_task_queue_handle) :: h
    type(qrm_task_type)         :: tsk
    character                   :: pol
    integer, optional           :: q
    logical                     :: qrm_sched_task
                                 
    integer                     :: i, iq

    if(present(q)) then
       iq = q
    else
       iq = qrm_task_thn
    end if

    !$ call omp_set_lock(h%queues(iq)%lock)
    if(h%queues(iq)%n .lt. max_tasks) then
       qrm_sched_task = .true. ! scheduling succeeded
!       !$omp critical(numtsk)
       h%ntsk = h%ntsk+1
!       !$omp end critical(numtsk)
       if(h%queues(iq)%n .eq. 0) then
          i = 1
          h%queues(iq)%n = 1
          h%queues(iq)%h = 1
       else
          if(pol .eq. 't') then
             i = mod(h%queues(iq)%h+h%queues(iq)%n-1,max_tasks)+1
          else if(pol .eq. 'h') then
             h%queues(iq)%h = h%queues(iq)%h-1
             if (h%queues(iq)%h .le. 0) h%queues(iq)%h=max_tasks
             i = h%queues(iq)%h
          end if
          h%queues(iq)%n = h%queues(iq)%n+1
       end if
       h%queues(iq)%q(i) = tsk
    else
       qrm_sched_task = .false. ! sched failed
    end if
    !$ call omp_unset_lock(h%queues(iq)%lock)

    return
  end function qrm_sched_task
  

  !> @brief Pops a task from a queue. Tasks are always popped from the
  !! head of the queue. The return value is .true. if something was
  !! found, .false. otherwise
  !!
  !! @param[in,out] h The handle to the set of queues
  !!
  !! @param[out] tsk the returned task (if any)
  !! 
  function qrm_get_task(h, tsk)
    implicit none
    type(qrm_task_queue_handle) :: h
    type(qrm_task_type)         :: tsk
    logical                     :: qrm_get_task

    integer :: i, queue

    tsk = qrm_task_type(-10, 0, 0, 0, .false.)
    qrm_get_task = .false.
    
    do i=0, qrm_task_nth-1
       queue = h%proxy_list(qrm_task_thn,i)
       ! !$ call omp_set_lock(h%queues(queue)%lock)
#if defined(_OPENMP)
       if(omp_test_lock(h%queues(queue)%lock)) then
#endif
          if(h%queues(queue)%n .gt. 0) then
             tsk   = h%queues(queue)%q(h%queues(queue)%h)
             if(tsk%bind .and. i .ne. 0) goto 10
             qrm_get_task = .true.
             h%queues(queue)%q(h%queues(queue)%h) = qrm_task_type(0, 0, 0, 0, .false.) !TODO: REMOVE!!!
             h%queues(queue)%h = mod(h%queues(queue)%h, max_tasks)+1
             h%queues(queue)%n = h%queues(queue)%n-1
             if (h%queues(queue)%n .eq. 0) h%queues(queue)%h = 0
             if(i .gt. 0) then
                h%stolen = h%stolen+1
             end if
             !$ call omp_unset_lock(h%queues(queue)%lock)
             return
          else
             tsk = qrm_task_type(-10, 0, 0, 0, .false.)
             qrm_get_task = .false.
          end if
10        continue
#if defined(_OPENMP)
          call omp_unset_lock(h%queues(queue)%lock)
       end if
#endif
    end do


    return
  end function qrm_get_task


  !> @brief Returns the number of tasks present on a set of queues
  !! referenced by a handle
  !!
  !! @param[in] h the handle
  function qrm_task_queue_card(h)
    type(qrm_task_queue_handle) :: h
    integer                     :: qrm_task_queue_card

    ! qrm_task_queue_card = h%queues(qrm_task_thn)%n
    qrm_task_queue_card = sum(h%queues(:)%n)

    return
  end function qrm_task_queue_card

  !> @brief Tells whether one, or all, queues are empty
  !!
  !! @param[in] h the handle
  !!
  !! @param[in] who (optional) if present, the function will only look
  !!                at the queue attached to thread "who" otherwise it will look at
  !!                all queues referenced by h
  !!                
  function qrm_task_queue_empty(h, who)
    type(qrm_task_queue_handle) :: h
    integer, optional           :: who
    logical                     :: qrm_task_queue_empty

    integer :: i    

    if(present(who)) then
       !$ call omp_set_lock(h%queues(who)%lock)
       qrm_task_queue_empty = h%queues(who)%n .eq. 0
       !$ call omp_unset_lock(h%queues(who)%lock)
    else
       qrm_task_queue_empty = .true.
       ! !$ do i=0, qrm_task_thn-1
       ! !$    call omp_set_lock(h%queues(i)%lock)
       ! !$ end do
          qrm_task_queue_empty = qrm_task_queue_empty .and. (maxval(h%queues(0:qrm_task_nth-1)%n) .eq. 0)
       ! !$ do i=qrm_task_thn-1, 0, -1
       ! !$    call omp_unset_lock(h%queues(i)%lock)
       ! !$ end do
    end if

    return
  end function qrm_task_queue_empty

  !> @brief Destroyes a set of queues
  !!
  !! @param[in,out] h the handle referencing the set of queues to be
  !!                  destroyed
  subroutine qrm_clean_task_queue(h)

    type(qrm_task_queue_handle) :: h

    integer :: i
    
    !$omp master
    do i=0, qrm_task_nth-1
       !$ call omp_destroy_lock(h%queues(i)%lock)
    end do
#if defined(debug)
    __QRM_PRNT_DBG('("Queues cleaned  -- ntsk:",i6,"   stolen:",i6)')h%ntsk,h%stolen
#endif

    deallocate(h%queues, h%proxy_list)
    !$omp end master
    return
  end subroutine qrm_clean_task_queue




#else
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! the following code is for a basic scheduling with no locality at
  !! all. mostly used for testing
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  !> @brief Inititalizes a set of queues attached to a family of
  !! threads referenced through the handle h
  !!
  !! @param[in] h The handle that references the queues
  !! 
  subroutine qrm_init_task_queue(h)
    implicit none
    type(qrm_task_queue_handle) :: h

    !$omp master
    allocate(h%queues(0:0))
    h%queues(0)%h = 0
    h%queues(0)%t = 0
    h%queues(0)%n = 0
    h%queues(0)%q = qrm_task_type(0, 0, 0, 0, .false.)
    !$ call omp_init_lock(h%queues(0)%lock)
    !$omp end master
    !$omp barrier !TODO maybe remove
    return
  end subroutine qrm_init_task_queue


  
  !> @brief Pushes a task on a queue
  !!
  !! @param[in,out] h the handle to the set of queues
  !!
  !! @param[in] tsk the task to be pushed
  !!
  !! @param[in] pol the policy. It can be either 'h', in which case
  !!                the task will be pushed on the head of the queue, or 't', in
  !!                which case it will be pushed on the tail of the queue
  !!
  !! @param[in] q the queue in the set where to push the task. This
  !!              can be used to improve data locality, e.g., a task can be pushed
  !!              on the queue attached to the thread which owns the data to be
  !!              accessed by the task itself
  function qrm_sched_task(h, tsk, pol, q)
    implicit none
    type(qrm_task_queue_handle) :: h
    type(qrm_task_type)         :: tsk
    character                   :: pol
    integer, optional           :: q
    logical                     :: qrm_sched_task
               
    integer :: i
                  
    !$ call omp_set_lock(h%queues(0)%lock)
    if(h%queues(0)%n .lt. max_tasks) then
       qrm_sched_task = .true. ! scheduling succeeded
       h%ntsk = h%ntsk+1
       if(h%queues(0)%n .eq. 0) then
          i = 1
          h%queues(0)%n = 1
          h%queues(0)%h = 1
       else
          if(pol .eq. 't') then
             i = mod(h%queues(0)%h+h%queues(0)%n-1,max_tasks)+1
          else if(pol .eq. 'h') then
             h%queues(0)%h = h%queues(0)%h-1
             if (h%queues(0)%h .le. 0) h%queues(0)%h=max_tasks
             i = h%queues(0)%h
          end if
          h%queues(0)%n = h%queues(0)%n+1
       end if
       h%queues(0)%q(i) = tsk
    else
       qrm_sched_task = .false. ! sched failed
    end if
    !$ call omp_unset_lock(h%queues(0)%lock)

    return
  end function qrm_sched_task
  

  !> @brief Pops a task from a queue. Tasks are always popped from the
  !! head of the queue. The return value is .true. if something was
  !! found, .false. otherwise
  !!
  !! @param[in,out] h The handle to the set of queues
  !!
  !! @param[out] tsk the returned task (if any)
  !! 
  function qrm_get_task(h, tsk)
    implicit none
    type(qrm_task_queue_handle) :: h
    type(qrm_task_type)         :: tsk
    logical                     :: qrm_get_task

    tsk = qrm_task_type(-10, 0, 0, 0, .false.)
    qrm_get_task = .false.
    
#if defined(_OPENMP)
       ! if(omp_test_lock(h%queues(0)%lock)) then
    call omp_set_lock(h%queues(0)%lock)
#endif
    if(h%queues(0)%n .gt. 0) then
       tsk   = h%queues(0)%q(h%queues(0)%h)
       qrm_get_task = .true.
       h%queues(0)%q(h%queues(0)%h) = qrm_task_type(0, 0, 0, 0, .false.) !TODO: REMOVE!!!
       h%queues(0)%h = mod(h%queues(0)%h, max_tasks)+1
       h%queues(0)%n = h%queues(0)%n-1
       if (h%queues(0)%n .eq. 0) h%queues(0)%h = 0
       goto 10
    else
       tsk = qrm_task_type(-10, 0, 0, 0, .false.)
       qrm_get_task = .false.
    end if
10  continue
#if defined(_OPENMP)
    call omp_unset_lock(h%queues(0)%lock)
       ! end if
#endif

    return
  end function qrm_get_task


  !> @brief Returns the number of tasks present on a set of queues
  !! referenced by a handle
  !!
  !! @param[in] h the handle
  function qrm_task_queue_card(h)
    type(qrm_task_queue_handle) :: h
    integer                     :: qrm_task_queue_card

    qrm_task_queue_card = h%queues(0)%n

    return
  end function qrm_task_queue_card

  !> @brief Tells whether one, or all, queues are empty
  !!
  !! @param[in] h the handle
  !!
  !! @param[in] who (optional) if present, the function will only look
  !!                at the queue attached to thread "who" otherwise it will look at
  !!                all queues referenced by h
  !!                
  function qrm_task_queue_empty(h, who)
    type(qrm_task_queue_handle) :: h
    integer, optional           :: who
    logical                     :: qrm_task_queue_empty

    qrm_task_queue_empty =  h%queues(0)%n .eq. 0

    return
  end function qrm_task_queue_empty

  !> @brief Destroyes a set of queues
  !!
  !! @param[in,out] h the handle referencing the set of queues to be
  !!                  destroyed
  subroutine qrm_clean_task_queue(h)

    type(qrm_task_queue_handle) :: h

    !$omp master
    !$ call omp_destroy_lock(h%queues(0)%lock)
#if defined(debug)
    __QRM_PRNT_DBG('("Queues cleaned  -- ntsk:",i6,"   stolen:",i6)')h%ntsk,h%stolen
#endif

    deallocate(h%queues)
    !$omp end master
    return
  end subroutine qrm_clean_task_queue


#endif


end module qrm_task_mod
