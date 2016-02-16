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
!> @file qrm_solve_r.F90
!! This file contains a routine that solves for R against multiple vectors
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function solves for R against multiple vectors
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in]     b     a 2d array containing the RHS vectors
!!
!! @param[out]    x     a 2d array containing the solution vectors
!!
subroutine _qrm_solve_r(qrm_mat, b, x)

  use _qrm_spmat_mod
  use _qrm_rfpf_mod
  use qrm_mem_mod
  use qrm_common_mod
  use qrm_task_mod
  use qrm_queue_mod
  use _qrm_solve_mod, protect=>_qrm_solve_r
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(inout) :: b(:,:)
  _qrm_data, intent(out)   :: x(:,:)
  
  integer                         :: nth, thn, info, f, dones
  integer                         :: i, p
  type(qrm_task_type)             :: task
  type(qrm_task_queue_handle)     :: tq_h
  type(qrm_queue)                 :: ready_q
  type(_qrm_front_type), pointer  :: front
  integer, allocatable            :: status(:)
  type(qrm_adata_type), pointer   :: adata
  type(_qrm_fdata_type), pointer  :: fdata
  logical                         :: got_task
#if defined (_OPENMP)
  integer(kind=omp_lock_kind), allocatable :: locks(:)
  integer(kind=omp_lock_kind)  :: dlock
#endif

  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_solve_r'
  
  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Solving for R")')

  ! simplify 
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  ! initialize queue
  call qrm_queue_init(ready_q , adata%nnodes, qrm_lifo_)

  ! initialize to safe values for the case where openmp is not used
  ! the status and locks arrays are needed because multiple solves may
  ! be run in parallel on the same fdata
  call qrm_aalloc(status, qrm_mat%adata%nnodes)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  !$ allocate(locks(adata%nnodes))

  do f = 1, adata%nnodes
     status(f) = qrm_ready_
     !$ call omp_init_lock(locks(f))
     p = adata%parent(f)
     if(p .eq. 0 .and. (adata%rc(f).ge.0)) call qrm_queue_push(ready_q, f)
  end do

  !$ call omp_init_lock(dlock)

  if(adata%ncsing .gt. 0) then
     dones = 1
  else
     dones = 0
  end if

  !$ call omp_set_num_threads(1)

  !$ nth = qrm_mat%icntl(qrm_nthreads_)
  !$omp parallel &
  !$omp & num_threads(nth) &
  !$omp & private(got_task, task, nth, thn, info) &
  !$omp & shared(ready_q, status, locks, dlock, dones, tq_h)

#if defined (_OPENMP)
  nth = omp_get_num_threads()
  thn = omp_get_thread_num()
#else
  nth   = 1
  thn   = 0
#endif

  info  = 0

  call qrm_par_mem_init()
  call qrm_init_task_queue(tq_h)

  taskloop: do
     if(qrm_err_stack%nelem .gt. 0) goto 9998

     ! fill up the queue if there are too few tasks
     if(qrm_task_queue_card(tq_h) .lt. nth) then
        call fill_queue_r( )
     end if
     
     got_task = qrm_get_task(tq_h, task)

     if(.not. got_task) cycle taskloop  ! queue is empty

     select case(task%id)
     case(qrm_task_exit_) 
        !$omp barrier 
        ! done, exit
        exit
     case(qrm_task_sol_)
        ! apply a set of Householder vectors
        call solve_r(task, thn)
     end select
  end do taskloop

9998 continue
  call qrm_clean_task_queue(tq_h)
  call qrm_par_mem_finalize()
  !$omp end parallel


  call qrm_adealloc(status)
  !$ deallocate(locks)
  call qrm_queue_free(ready_q)

  if(qrm_err_stack%nelem .gt. 0) then
     call qrm_err_push(22, name)
     goto 9999
  end if

  ! check if there is a singleton block and eventually solve for it
  if (qrm_mat%adata%ncsing .gt. 0) then
     call _qrm_solve_sing_front(qrm_mat, b, x, 'n')
     __QRM_CHECK_RET(name,'qrm_solve_sing_front',9999)
  end if

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
  subroutine fill_queue_r( )
    implicit none

    type(_qrm_front_type), pointer :: front
    integer   :: f
    integer   :: thn
    logical   :: found
    type(qrm_task_type)             :: tsk

#if defined(_OPENMP)
    thn=omp_get_thread_num()
#else
    thn=0
#endif

    found = .false.

    f = 0
    do

       f = qrm_queue_next(ready_q, f)
       if(f .eq. 0) exit

       front => fdata%front_list(f)
       
#if defined(_OPENMP)
       if(.not. omp_test_lock(locks(f))) cycle
#endif
       if(status(f) .eq. qrm_ready_) then
          call qrm_queue_rm(ready_q, f)
          tsk = qrm_task_type(qrm_task_sol_, front%num, 0, 0, .false.)
          if(qrm_sched_task(tq_h, tsk, 'h')) then
             ! mark the column as assigned
             found = .true.
             status(f) = qrm_busy_
          end if
       end if
#if defined(_OPENMP)
       call omp_unset_lock(locks(f))
#endif
    end do

    ! if nothing was found above, then check if the facto is over
    ! otherwise return
    if(found) return

    call check_solver_over( )

    return
  end subroutine fill_queue_r




!==========================================================================================
!==========================================================================================
  subroutine check_solver_over( )
    ! checks whether the apply_qt is over and eventually schedules a
    ! termination task
    implicit none

    type(qrm_task_type) :: tsk
    logical             :: found

    ! all the fronts are done
    !$ call omp_set_lock( dlock )
    if(dones .eq. fdata%nfronts) then
       tsk = qrm_task_type(qrm_task_exit_, 0, 0, 0, .false.)
       found = qrm_sched_task(tq_h, tsk, 't')
    end if
    !$ call omp_unset_lock( dlock )

    return
  end subroutine check_solver_over






!==========================================================================================
!==========================================================================================
  subroutine solve_r(task, thn)
    implicit none
    type(qrm_task_type) :: task
    integer         :: thn

    type(_qrm_front_type), pointer :: front
    integer                         :: f, p, c, info

    front  => null()
    
    ! to make things easier
    front => qrm_mat%fdata%front_list(task%front)
    info  = 0

    call front_r(front, info)
    status(task%front) = qrm_done_
    call qrm_queue_rm(ready_q, front%num)

#if defined(_OPENMP)
    call omp_set_lock( dlock )
#endif
    dones = dones+1
#if defined(_OPENMP)
    call omp_unset_lock( dlock )
#endif

    ! sweep over the children. Small children are treated immediately,
    ! the others are pushed on the ready_q
    do p = adata%childptr(front%num), adata%childptr(front%num+1)-1
       c = adata%child(p)
       if(adata%small(c) .eq. 1) then
          call do_subtree_r(c, info)
          if(info .ne. 0) goto 9997
       else
          call qrm_queue_push(ready_q, c)
       end if
    end do

    
9997 continue
    return
  end subroutine solve_r






!==========================================================================================
!==========================================================================================
  subroutine do_subtree_r(fnum, info)
    implicit none

    integer :: fnum, info

    type(_qrm_front_type), pointer :: front
    integer :: node, c, acc, thn, p
    type(qrm_queue) :: sub_q

    info = 0
    acc = 0

#if defined(_OPENMP)
    thn = omp_get_thread_num()
#else
    thn = 0
#endif

    call qrm_queue_init(sub_q, adata%nnodes, qrm_fifo_)
    call qrm_queue_push(sub_q, fnum)

    do
       node = qrm_queue_pop(sub_q)

       if(node .eq. 0) exit

       front => qrm_mat%fdata%front_list(node)
       
       call front_r(front, info)
       status(node) = qrm_done_

       !$ call omp_set_lock( dlock )
       dones = dones+1
       !$ call omp_unset_lock( dlock )

       ! sweep over the children. Small children are treated immediately,
       ! the others are pushed on the ready_q
       do p = adata%childptr(front%num), adata%childptr(front%num+1)-1
          c = adata%child(p)
          call qrm_queue_push(sub_q, c)
          acc = acc+1
       end do
       
    end do
    
    call qrm_queue_free(sub_q)
  
    return

  end subroutine do_subtree_r






!==========================================================================================
!==========================================================================================
  subroutine front_r(front, info)

    use _qrm_utils_mod
    use _qrm_rfpf_mod
    use _qrm_spmat_mod
    use qrm_mem_mod
    use qrm_common_mod
    implicit none

    type(_qrm_front_type) :: front
    integer                :: info

    integer :: k, m, pk, jp, n, i, thn, cnt, j, pv1, pv2, r

    ! TODO: optimize this by allocating once in the subtree
    _qrm_data, allocatable :: in_x(:,:)

    ! error management
    character(len=*), parameter     :: name='front_r'

    ! shortcut
    if (min(front%m, front%n) .le.0) goto 9999
    if (front%npiv .le. 0) goto 9999

#if defined(_OPENMP)
    thn = omp_get_thread_num()
#else
    thn = 0
#endif
    ! !$ call omp_set_lock( locks(front%num) )

    ! write(*,'(i3," =--> Solve R  : ",i4)')thn, front%num
    call qrm_aalloc(in_x, front%n, size(x,2), info)
    __QRM_CHECK_RET(name,'qrm_aalloc',9999)

    in_x(1:front%npiv,:)         = b(front%rows(1:front%npiv),:)
    in_x(front%npiv+1:front%n,:) = x(front%cols(front%npiv+1:front%n),:)


#if defined (RFPF)
    ! pv1 = 1
     
    ! ! set r to be the starting row of the R factor
    ! r = ((front%npiv-1)/front%nb)*front%nb+1
    ! do 
    !    if(r .lt. 1) exit
    !    ! k is the size of the panel
    !    k   = min(front%nb,front%npiv-r+1)
    !    ! n is the col size of the rectangular part of the panel
    !    n   = front%n - r-k+1
    !    ! pv1 is the position of the first element of the panel in the
    !    ! front%r array (where the rfpf part starts)
    !    pv1 = (front%n-r+1)*(r-1)+r*(r-1)/2+1
    !    ! pv2 is the position in the front%r array where the
    !    ! rectangular part of this panel starts
    !    if(n .le. 0) then
    !       pv2 = 1
    !    else
    !       pv2 = pv1+(k*(k+1))/2
    !    end if
    !    if(n .gt. 0) call _xgemm('n', 'n', k, size(x,2), n, -_qrm_one, front%r(pv2), &
    !         &k, in_x(r+k,1), front%n, _qrm_one, in_x(r,1), front%n)
    !    call _xrpsm('l', 'u', 'n', 'n', k, size(x,2), _qrm_one, front%r(pv1), in_x(r,1), front%n)
    !    ! go backwards one panel
    !    r = r - front%nb
    ! end do
#else
    n = size(x,2)
    cnt=size(front%r)+1
    outer: do jp = front%npiv - mod(front%npiv, front%nb)+1, 1, -front%nb
       pk = min(front%nb, front%npiv-jp+1) 
       if(pk .le. 0) cycle
       ! write(*,*)'solve',jp,pk

       inner: do j = jp+pk-mod(pk,front%ib), jp, -front%ib
          m = min(front%ib, jp+pk - j)
          if(m .le. 0) cycle
          k = front%n-j-m+1
          cnt = cnt-m*(front%n-j+1)
          ! write(*,*)size(front%r), cnt, j, k+m, m
          if(k.gt.0) call _xgemm('n', 'n', m, n, k, -_qrm_one, front%r(cnt+m*m), m, &
               & in_x(j+m,1), front%n, _qrm_one, in_x(j,1), front%n)
          call _xtrsm('l', 'u', 'n', 'n', m, n, _qrm_one, front%r(cnt), &
               & m, in_x(j,1), front%n)
       end do inner
    end do outer
#endif

    x(front%cols(1:front%npiv),:) = in_x(1:front%npiv,:)

    call qrm_adealloc(in_x)
    __QRM_CHECK_RET(name,'qrm_adealloc',9999)
    ! !$ call omp_unset_lock( locks(front%num) )

9999 continue !error
    return

  end subroutine front_r

end subroutine _qrm_solve_r
