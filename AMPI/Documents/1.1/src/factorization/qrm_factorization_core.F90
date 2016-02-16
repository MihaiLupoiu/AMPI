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
!> @file qrm_factorization_core.F90
!! This file holds the core of the QR numerical factorization.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This is the main factorization routine. It performs the factorization of all
!! the fronts that have been assigned to the node (MPI task). The whole process
!! is OpenMP parallel.

!> This factorization is completely dynamic and results in an out of order execution
!! of computational and symbolic tasks. Each front is cut into block-columns. Computational
!! tasks are scheduled on a per-front or per-block-column basis according to a graph of
!! dependencies; each task is pushed on the queue of the thread who owns the related front.
!! The ownership of a front is fixed at the moment where the front is activated. If a thread
!! runs out of tasks to perform (i.e., its queue is empty) it can steal tasks from other
!! threads (see the @link qrm_task_mod @endlink module for the details).
!!
!! Six different types of tasks can be executed by a thread:
!! - panel    : the panel factorization of a column in a front
!! - update   : the update of a block-column wrt to a panel in a front
!! - activate : activation of a front. this operation consists in computing the full
!!              structure of a front and allocating all the memory areas needed for the
!!              subsequent factorization
!! - assemble : this operation consists in assemblying one block-column in the contribution
!!              block of a front into the father node
!! - clean    : cleans up a front, i.e., it deallocates all the memory areas that are no more
!!              needed
!!
!! @param[in,out] qrm_mat This is the main data structure associated
!!                to a problem. The @link _qrm_factorization_core @endlink
!!                assumes that the analysis has been done and thus that the
!!                qrm_mat%adata field contains meaningful stuff
!!

subroutine _qrm_factorization_core(qrm_mat)

  use _qrm_spmat_mod
  use qrm_mem_mod
  use qrm_queue_mod
  use qrm_task_mod
  use qrm_error_mod
  use qrm_common_mod
  use _qrm_factorization_mod, protect => _qrm_factorization_core
#if defined(trace)
  use qrm_trace_mod
#endif
  !$ use omp_lib
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat

  integer :: qrm_nth, nth, thn, info, i

#if defined(trace)
  integer :: pnl_id, upd_id, act_id, asm_id, clean_id, sched_id
#endif
  type(qrm_task_type)             :: task
  type(qrm_task_queue_handle)     :: tq_h
  _qrm_data, allocatable          :: work(:)
  type(qrm_queue)                 :: active_q, ready_q
  real(kind(1.d0))                :: flops, times(6), t1, t2
  logical                         :: got_task
  type(qrm_adata_type), pointer   :: adata
  type(_qrm_fdata_type), pointer  :: fdata

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_factorization_core'

  call qrm_err_act_save(err_act)

#if defined(trace)
  call qrm_trace_init(1)
  call qrm_trace_create_event('Panel', pnl_id)
  call qrm_trace_create_event('Update', upd_id)
  call qrm_trace_create_event('Activate', act_id)
  call qrm_trace_create_event('Assemble', asm_id)
  call qrm_trace_create_event('Clean', clean_id)
#endif

  ! simplify 
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  ! The factorization uses a double system of queues:
  ! 1) front queues: these are used to hold the fronts that are
  !    currently being treated. There is an active front queue that contains
  !    all the fronts that are currently active and a another one holding all
  !    the fronts that can be activated (i.e., whose children are all active).
  !    When scheduling tasks, a thread first tries to schedule tasks related to
  !    the active fronts. If non is available the thread will schedule the activation
  !    of the first front in the ready queue.
  ! 2) task queues: each thread has its own queue of tasks (as defined in the
  !    qrm_task_mod). Tasks are pushed on the queue of the thread that own the front.
  !    Then, in the main loop below, each thread tries to fetch a task from its own
  !    queue but in the case where it is empty, the thread can steal a task from
  !    others.

  ! initialize front queues
  call qrm_queue_init(active_q, adata%nnodes, qrm_fifo_)
  call qrm_queue_init(ready_q , adata%nnodes, qrm_lifo_)
  ! push all the leaves in the ready queue. this queue is lifo (in order
  ! to achieve as much as possible a post-order traversal) so leaves are pushed
  ! in reverse order.
  do i=adata%nleaves, 1, -1
     call qrm_queue_push(ready_q, adata%leaves(i))
     fdata%front_list(adata%leaves(i))%status = qrm_activable_
  end do


  ! initialize to safe values for the case where openmp is not used

#if defined (_OPENMP)
  call omp_set_num_threads(1)
  qrm_nth = qrm_mat%icntl(qrm_nthreads_)
#endif

  !$omp parallel &
  !$omp & num_threads(qrm_nth) &
  !$omp & private(task, nth, thn, work, info, got_task) &
  !$omp & shared(active_q, ready_q, tq_h) &
  !$omp & reduction(+:flops)

#if defined (_OPENMP)
  nth = omp_get_num_threads()
  thn = omp_get_thread_num()
#else
  nth = 1
  thn = 0
#endif

  call qrm_par_mem_init()

  flops = 0.d0

  !$omp master
  __QRM_PRNT_DBG('("Num threads: ",i3)')nth
  !$omp end master

  ! initialize the task queues
  call qrm_init_task_queue(tq_h)

  ! allocate work area for panel and update
  call qrm_aalloc(work, qrm_mat%icntl(qrm_nb_)*qrm_mat%icntl(qrm_nb_))
  __QRM_CHECK_RET(name,'qrm_alloc',9998)

  info = 0

  !$omp barrier
  taskloop: do
     ! this is the main loop. At each iteration of this loop a thread
     ! does the following things:
     ! 1) checks if errors are present on the error stack
     ! 2) pushes some tasks on the queues if queues begin to run out of entries
     ! 3) gets a task from the queues (its own queue first and then from others)
     ! 4) executes the action corresponding to the fetched task

     if(qrm_err_stack%nelem .gt. 0) goto 9998

     ! fill up the queue if there are too few tasks
     if(qrm_task_queue_card(tq_h) .lt. nth) then
        call fill_queue()
     end if
     
     got_task = qrm_get_task(tq_h, task)

     if(.not. got_task) cycle taskloop  ! queue is empty
     
     ! perform the action corresponding to the fetched task
     select case(task%id)
     case(qrm_task_exit_) 
        !$omp barrier
        ! done, exit
        exit
     case(qrm_task_pnl_)
        ! factorize a panel
        call panel   (task, thn, work, flops)
     case(qrm_task_upd_)
        ! update a block-column
        call update  (task, thn, work, flops)
     case(qrm_task_act_)
        ! activate a front
        call activate(task, thn, flops)
     case(qrm_task_asm_)
        ! assemble column into father front
        call assemble(task, thn)
     case(qrm_task_cln_)
        ! clean the front
        call clean   (task, thn)
     end select
  end do taskloop

9998 continue
  call qrm_adealloc(work)
  call qrm_clean_task_queue(tq_h)
  call qrm_par_mem_finalize()
  !$omp end parallel

#if defined(trace)
  call qrm_trace_log_dump('trace.svg')
#endif

  ! free the front queues
  call qrm_queue_free(active_q)
  call qrm_queue_free(ready_q)

  ! check if errors were generated
  if(qrm_err_stack%nelem .gt. 0) then
     call qrm_err_push(20, name)
     goto 9999
  end if

  qrm_mat%gstats(qrm_facto_flops_) = floor(flops,kind(qrm_mat%gstats(1)))

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
  subroutine fill_queue()
    implicit none

    type(_qrm_front_type), pointer :: front
    integer   :: f
    integer   :: thn
    logical   :: found

    thn=0
    !$ thn=omp_get_thread_num()
    found = .false.
    ! loop over the queue of active nodes
    f = 0
    active_fronts: do
       f = qrm_queue_next(active_q, f)
       if(f .eq. 0) exit active_fronts
       front => fdata%front_list(f)
       ! if(front%owner .ne. thn) cycle active_fronts ! front f does not belong to me
#if defined(_OPENMP)
       if(.not. omp_test_lock(front%lock)) cycle active_fronts ! front f does not belong to me
#endif
       found = found .or. front_sched_ops(front)
#if defined(_OPENMP)
       call omp_unset_lock(front%lock)
#endif
    end do active_fronts

    ! if nothing was found above, then look for fronts to activate
    ! otherwise return
    if(found) return

    ! schedule the activation of a new front
    f = 0
    ready_fronts: do
       f = qrm_queue_next(ready_q, f)
       if(f .eq. 0) exit ready_fronts
       front => fdata%front_list(f)
       found = found .or. front_sched_act(front)
       if (found) exit ready_fronts
    end do ready_fronts
    ! if nothing was found above, then check if the facto is over
    ! otherwise return
    if(found) return

    call check_facto_over()

    return
  end subroutine fill_queue






!==========================================================================================
!==========================================================================================
  function front_sched_ops(front)
    ! schedules operations for a specific (active) front
    implicit none
    logical                :: front_sched_ops
    type(_qrm_front_type) :: front

    type(qrm_task_type)             :: tsk
    type(_qrm_front_type), pointer :: father, cfront
    integer                         :: p, c, i, j, childpt, child

    ! initialize
    front_sched_ops = .false.

    if(  (front%status .le. 0) .or. &
        &(front%status .eq. qrm_done_)) &
        & return

    ! check if the front can be cleaned
    if(minval(front%ptable) .eq. qrm_done_) then
       tsk = qrm_task_type(qrm_task_cln_, front%num, 0, 0, .false.)
       if(qrm_sched_task(tq_h, tsk, 't', front%owner)) then
          ! mark the column as assigned
          front_sched_ops = .true.
          front%status = qrm_busy_
          return ! if the front can be cleaned
                 ! there's nothing more we can schedule on it
       end if
    end if

#if defined(fullasm)
    if(maxval(front%ptable) .lt. 0) return
#endif

    ! first look for panels
    if (front%lastpnl .lt. front%np) then
       ! there are still panels left to do in this front
       ! check if it's possible to reduce panel lastpnl+1
       if (front%ptable(front%lastpnl+1) .eq. front%lastpnl) then
          ! column lastpnl+1 is up-to-date wrt panel lastpnl.
          ! that means it can be reduced
          tsk = qrm_task_type(qrm_task_pnl_, front%num, front%lastpnl+1, 0,.false.)
          if(.not. qrm_sched_task(tq_h, tsk, 'h', front%owner)) then
             ! cannot schedule the task (queue is full)
             return
          else
             ! mark the column as assigned
             front_sched_ops = .true.
             front%ptable(front%lastpnl+1) = qrm_busy_
          end if
       end if
    end if



    ! now look for updates
    do p=1, front%lastpnl
       ! for every reduced panel, check if there are updates to be done
       do c=p+1, front%nc
          if(front%ptable(c) .eq. qrm_busy_) cycle ! column is busy
          ! just need to check whether c is up-to-date wrt p-1
          if(front%ptable(c) .eq. p-1) then
             tsk = qrm_task_type(qrm_task_upd_, front%num, p, c, .false.)
             if(.not. qrm_sched_task(tq_h, tsk, 't', front%owner)) then
                ! cannot schedule the task (queue is full)
                return
             else
                ! mark the column as assigned
                front_sched_ops = .true.
                front%ptable(c)  = qrm_busy_
             end if
          end if
       end do
    end do


    ! now check if there's something to assemble
    p = adata%parent(front%num)
    if (p .eq. 0) return ! front is a root node. no need to assemble
    father => fdata%front_list(p)
    ! check is father node is active
    if(father%status .eq. qrm_active_) then
       do c=front%npiv/front%nb+1, front%nc
          if (front%ptable(c) .eq. min(c, front%np)) then
             tsk = qrm_task_type(qrm_task_asm_, front%num, 0, c, .false.)
             if(.not. qrm_sched_task(tq_h, tsk, 'h', front%owner)) then
                ! cannot schedule the task (queue is full)
                exit
             else
                ! mark the column as assigned
                front_sched_ops = .true.
                front%ptable(c) = qrm_busy_
                father%status = qrm_busy_
                ! better not schedule more than one assembly at once
                ! because we may end up executing mutiple assemblies
                ! at the same time which results in memory conflicts.
                exit
             end if
          end if
       end do
    end if


    return
  end function front_sched_ops








!==========================================================================================
!==========================================================================================
  function front_sched_act(front)
    ! check if a specific front can be activated and eventually
    ! schedules its activation
    implicit none
    logical                :: front_sched_act
    type(_qrm_front_type) :: front

    type(qrm_task_type)    :: tsk

    ! initialize
    front_sched_act = .false.

    !$ if(omp_test_lock(front%lock)) then
    if(front%status .eq. qrm_activable_) then
       tsk = qrm_task_type(qrm_task_act_, front%num, 0, 0, .false.)
       if(qrm_sched_task(tq_h, tsk, 't')) then
          ! mark the column as assigned
          front_sched_act = .true.
          front%status = qrm_busy_
       end if
    end if
    !$ call omp_unset_lock(front%lock)
    !$ end if

    return
  end function front_sched_act






!==========================================================================================
!==========================================================================================
  subroutine check_facto_over()
    ! checks whether the facto is over and eventually schedules a
    ! termination task
    implicit none

    type(qrm_task_type) :: tsk
    logical             :: found

    ! all the fronts are done
    if(fdata%done .eq. fdata%nfronts) then
       tsk = qrm_task_type(qrm_task_exit_, 0, 0, 0, .false.)
       found = qrm_sched_task(tq_h, tsk, 't')
    end if

    return
  end subroutine check_facto_over









!==========================================================================================
!==========================================================================================
  subroutine panel(task, thn, work, flops)
    ! this subroutine performs the panel reduction
    ! block-column task%pnl in front task%front
    ! 
    
    use _qrm_fdata_mod
    use qrm_common_mod
    implicit none

    real(kind(1.d0))    :: t1, t2 
    type(qrm_task_type) :: task
    integer             :: thn
    _qrm_data           :: work(:)
    real(kind(1.d0))    :: flops

    type(_qrm_front_type), pointer :: front
    integer                         :: f, np, nc, nb, &
         & m, n, k, i, info, j, ib, fcol, lcol, ii, jj
    type(qrm_task_type)             :: tsk

    front => fdata%front_list(task%front)
    f     =  adata%parent(front%num)

    call qrm_get(qrm_mat, 'qrm_ib', ib)

#if defined(debug)    
    __QRM_PRNT_DBG('(i3," =--> Panel   : ",i4," -- ",i4)')thn,task%front,task%pnl
#endif

#if defined(trace)
    call qrm_trace_event_start(pnl_id, thn)
#endif

    fcol = (task%pnl-1)*front%nb+1        ! the first pivot in the panel
    lcol = min(task%pnl*front%nb,front%n) ! the last pivot in the panel

    ! flops = flops+pureflopscount(front%stair, front%n, j, front%nb)

    do i=fcol, min(task%pnl*front%nb,front%ne), ib
       n   = min(ib, lcol-i+1)               ! the width of the panel
       m   = max(front%stair(min(i+n-1,front%ne))- i+1,0) ! the height of the panel
       k   = min(m,n)
       j   = i+n
       
       call _xgeqrf(m, n, front%front(i,i), front%m, front%tau(i), work(1), &
            & ib**2, info)
       
       call _xlarft( 'f', 'c', m, k, front%front(i,i), front%m, front%tau(i), &
            & front%t(1,i), front%ib )

       flops = flops+qrm_count_flops(m, n, n, 'panel')

       n   = lcol-j+1
       if ( n .le. 0) exit 
       call _xlarfb('l', 't', 'f', 'c', m, n, k, front%front(i,i), &
            & front%m, front%t(1,i), front%ib, &
            & front%front(i,j), front%m, work(1), front%nb)

       flops = flops+qrm_count_flops(m, n, k, 'update')

    end do

#if defined(trace)
    call qrm_trace_event_stop(pnl_id, thn)
#endif

    ! once all the computations are done, update the status of the
    ! factorization

    !$ call omp_set_lock(front%lock)
    if(task%pnl .lt. front%nc) then
       if(front%ptable(task%pnl+1) .eq. task%pnl-1) then
          tsk = qrm_task_type(qrm_task_upd_, front%num, task%pnl, task%pnl+1, .false.)
          if(qrm_sched_task(tq_h, tsk, 'h', front%owner)) then
             front%ptable(task%pnl+1)  = qrm_busy_
          end if
       end if
    end if

    if (f .eq. 0) then
       ! if column does not have to be assembled, it's marked done
       front%ptable(task%pnl) = qrm_done_
    else
       if ((task%pnl .le. front%npiv/front%nb) .or. &
            & (front%npiv .eq. front%ne)) then
          ! colums corresponding to pivots do not have to be assembled
          front%ptable(task%pnl) = qrm_done_
       else
          front%ptable(task%pnl) = task%pnl
       end if
    end if
    ! check if the front is completely factorized
    if(task%pnl .eq. front%nc) then
       front%status = qrm_factorized_
    end if
    front%lastpnl=task%pnl
    !$ call omp_unset_lock(front%lock)
    
    return

  end subroutine panel







!==========================================================================================
!==========================================================================================
  subroutine update(task, thn, work, flops)
    ! this subroutine performs the update of
    ! block-column task%col wrt panel task%pnl
    ! in front task%front
    ! 
    
    use _qrm_fdata_mod
    use qrm_common_mod
    implicit none

    type(qrm_task_type) :: task
    integer             :: thn
    _qrm_data           :: work(:)

    real(kind(1.d0))    :: flops
    type(qrm_task_type) :: tsk
    real(kind(1.d0))    :: t1, t2 


    type(_qrm_front_type), pointer :: front
    integer :: f, i, j, m, n, k, ib, fcol, lcol

    front => fdata%front_list(task%front)
    f     =  adata%parent(front%num)

#if defined(debug)
    __QRM_PRNT_DBG('(i3," =--> Update  : ",i4," -- ",i4,2x,i4)')thn,task%front,task%pnl,task%col
#endif

#if defined(trace)
    call qrm_trace_event_start(upd_id, thn)
#endif

    call qrm_get(qrm_mat, 'qrm_ib', ib)

    fcol = (task%pnl-1)*front%nb+1          ! the first pivot in the panel
    lcol = min(task%pnl*front%nb, front%ne) ! the last pivot in the panel
    j    = (task%col-1)*front%nb+1          ! the first column in the update
    n    = min(front%nb, front%n-j+1)       ! the number of columns to update

    do i=fcol, lcol, ib
       k     = min(ib, lcol-i+1)               ! the width of the reflectors set
       m     = max(front%stair(min(i+k-1,front%ne))- i+1,0) ! the height of the panel
       k     = min(k,m)
    
       call _xlarfb('l', 't', 'f', 'c', m, n , k, front%front(i,i), &
            & front%m, front%t(1,i), front%ib, &
            & front%front(i,j), front%m, work(1), front%nb)

       flops = flops+qrm_count_flops(m, n, k, 'update')

    end do

    ! once all the computations are done, update the status of the
    ! factorization

#if defined(trace)
    call qrm_trace_event_stop(upd_id, thn)
#endif

    !$ call omp_set_lock(front%lock)
    if((f .eq. 0) .and. (task%pnl .eq. front%np)) then
       ! we just updated wrt the last panel and no assembly is required:
       ! mark as done
       front%ptable(task%col) = qrm_done_
    else if ((task%pnl .eq. front%np) .and. &
            & (front%npiv .eq. front%ne)) then
       front%ptable(task%col) = qrm_done_
    else 
       ! mark the column as up-to-date wrt panel task%pnl
       front%ptable(task%col) = task%pnl
    end if
    i = minval(front%ptable)
    if(i .ge. front%np) front%status = qrm_factorized_

    if((task%col .eq. (task%pnl+1)) .and. (task%col .le. front%np)) then
       tsk = qrm_task_type(qrm_task_pnl_, front%num, task%col, 0, .false.)
       if(qrm_sched_task(tq_h, tsk, 'h', front%owner)) then
          front%ptable(task%col)  = qrm_busy_
       end if
    end if
    !$ call omp_unset_lock(front%lock)
    
    return

  end subroutine update






!==========================================================================================
!==========================================================================================
  subroutine activate(task, thn, flops)
    ! this subroutine activates front task%front

    use _qrm_fdata_mod
    use _qrm_factorization_mod
    use qrm_queue_mod
    implicit none

    type(qrm_task_type) :: task
    integer             :: thn
    real(kind(1.d0))    :: flops

    type(_qrm_front_type), pointer :: front
    integer                         :: msize, i, f
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='activate'
    
    call qrm_err_act_save(err_act)
  
    front  => fdata%front_list(task%front)

#if defined(_OPENMP)
    call omp_set_lock(front%lock)
#endif

#if defined(debug)
    __QRM_PRNT_DBG('(i3," =--> Activate: ",i4)')thn,task%front
#endif

#if defined(trace)
    call qrm_trace_event_start(act_id, thn)
#endif

    call _qrm_activate_front(qrm_mat, front%num, flops)
    __QRM_CHECK_RET(name,'qrm_activate_front',9999)

    call qrm_queue_rm(ready_q, task%front)
    call qrm_queue_push(active_q, task%front)
    if((front%m .eq. 0) .or. (front%n .eq. 0)) then
       front%status = qrm_done_
       !$ call omp_set_lock(fdata%lock)
       fdata%done = fdata%done+1
       !$ call omp_unset_lock(fdata%lock)
    else
       front%status = qrm_active_
    end if
     !$ call omp_unset_lock(front%lock)

#if defined(trace)
    call qrm_trace_event_stop(act_id, thn)
#endif

    f = adata%parent(front%num)
    if(f .ne. 0) then
       front  => fdata%front_list(f)
       !$ call omp_set_lock(front%lock)
       if (front%status .eq. qrm_ready_) then 
          call qrm_queue_push(ready_q, f) 
          front%status = qrm_activable_
       end if
       !$ call omp_unset_lock(front%lock)
    end if

    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    !$ call omp_unset_lock(front%lock)
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    
    return
    
  end subroutine activate







!==========================================================================================
!==========================================================================================
  subroutine assemble(task, thn)
    ! this subroutine performs the assembly of column
    ! task%col of front task%front into its father
    ! 
    use _qrm_fdata_mod
    implicit none

    type(qrm_task_type) :: task
    integer :: thn

    type(_qrm_front_type), pointer :: front, father
    integer :: f, i, j, nb, row, col, frow, fcol, b
    integer :: nc, mc, ptr
    type(qrm_task_type) :: tsk
    integer, allocatable :: tmp(:)
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='assemble'
    
    call qrm_err_act_save(err_act)

    front  => fdata%front_list(task%front)
    f   =  adata%parent(front%num)
    father => fdata%front_list(f)

#if defined(debug)
    __QRM_PRNT_DBG('(i3," =--> Assemble: ",i4," -- ",i4)')thn,task%front,task%col
#endif

#if defined(trace)
    call qrm_trace_event_start(asm_id, thn)
#endif
    
    call qrm_aalloc(tmp, father%nc)
    __QRM_CHECK_RET(name,'qrm_aalloc',9999)
    tmp = 0

    do col = (task%col-1)*front%nb+1, min(front%n, task%col*front%nb)
       if(col .le. front%npiv) cycle
       fcol = front%colmap(col-front%npiv)
       do row = front%npiv+1, min(col, min(front%m, front%n))
          frow = front%rowmap(row - front%npiv)
          father%front(frow, fcol) = front%front(row,col)
       end do
       b = (fcol-1)/father%nb+1
       tmp(b) = tmp(b)+1
    end do

    !$ call omp_set_lock(father%lock)
    father%ptable = father%ptable+tmp
    father%status = qrm_active_
    !$ call omp_unset_lock(father%lock)

    call qrm_adealloc(tmp)

#if defined(trace)
    call qrm_trace_event_stop(asm_id, thn)
#endif
    
    !$ call omp_set_lock(front%lock)
    front%ptable(task%col) = qrm_done_
    !$ call omp_unset_lock(front%lock)

    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    
    return
    
  end subroutine assemble



!==========================================================================================
!==========================================================================================
  subroutine clean(task, thn)
    ! cleanup: deallocates all the memory areas
    ! needed for the front factorization.
    ! The corresponding R and (eventually) H parts
    ! are stored into dedicated memory areas

    use _qrm_fdata_mod
    implicit none

    type(qrm_task_type) :: task
    integer :: thn

    type(_qrm_front_type), pointer :: front, father
    integer :: f, mn, h, i, j
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='clean'
    
    call qrm_err_act_save(err_act)

    front  => fdata%front_list(task%front)

#if defined(debug)
    __QRM_PRNT_DBG('(i3," =--> Clean   : ",i4)')thn,task%front
#endif

#if defined(trace)
    call qrm_trace_event_start(clean_id, thn)
#endif

    !$ call omp_set_lock(front%lock)
    call _qrm_clean_front(qrm_mat, task%front)
    !$ call omp_unset_lock(front%lock)
    __QRM_CHECK_RET(name,'qrm_clean_front',9999)

#if defined(trace)
    call qrm_trace_event_stop(clean_id, thn)
#endif

    ! front is done, remove it from the queue of active nodes
    call qrm_queue_rm(active_q, front%num)
    
    !$ call omp_set_lock(fdata%lock)
    fdata%done = fdata%done+1
    !$ call omp_unset_lock(fdata%lock)
    
    !$ call omp_set_lock(front%lock)
    front%status = qrm_done_
    !$ call omp_unset_lock(front%lock)

    call qrm_err_act_restore(err_act)
    return
    
9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    
    return
    
  end subroutine clean





end subroutine _qrm_factorization_core



