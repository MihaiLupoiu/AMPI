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
!> @file qrm_apply_q.F90
!! This file contains a routine that applies Q to a vector/matrix
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function applies Q to a vector/matrix
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in,out] b     a 1d array containing the vector to which Q will 
!!                      be applied. 
!!
subroutine _qrm_apply_q(qrm_mat, b)

  use _qrm_spmat_mod
  use _qrm_rfpf_mod
  use qrm_mem_mod
  use qrm_common_mod
  use qrm_task_mod
  use qrm_queue_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(inout)  :: b(:,:)
  

  integer                         :: nth, thn, info, f, dones, qrm_nth
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
  character(len=*), parameter     :: name='qrm_apply_q'
  
  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Applying Q")')

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

#if defined(_OPENMP)
  call omp_set_num_threads(1)
  qrm_nth = qrm_mat%icntl(qrm_nthreads_)
#endif

  !$omp parallel &
  !$omp & num_threads(qrm_nth) &
  !$omp & private(got_task, task, nth, thn, info) &
  !$omp & shared(ready_q, status, locks, dlock, dones, tq_h)

#if defined(_OPENMP)
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
        call fill_queue_q( )
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
        call apply_q(task, thn)
     end select
  end do taskloop

9998 continue
  call qrm_clean_task_queue(tq_h)
  call qrm_par_mem_finalize()
  !$omp end parallel

  call qrm_adealloc(status)
  call qrm_queue_free(ready_q)

#if defined (_OPENMP)
  deallocate(locks)
#endif

  if(qrm_err_stack%nelem .gt. 0) then
     call qrm_err_push(22, name)
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

contains



!==========================================================================================
!==========================================================================================
  subroutine fill_queue_q( )
    implicit none

    type(_qrm_front_type), pointer :: front
    integer                        :: f
    integer                        :: thn
    logical                        :: found
    type(qrm_task_type)            :: tsk

#if defined (_OPENMP)
    thn = omp_get_thread_num()
#else
    thn = 0
#endif

    found = .false.

    f = 0
    do

       f = qrm_queue_next(ready_q, f)
       if(f .eq. 0) exit

       front => fdata%front_list(f)
       
#if defined (_OPENMP)
       if(.not. omp_test_lock(locks(f))) cycle
       ! call omp_set_lock(locks(f))
#endif
       if(status(f) .eq. qrm_ready_) then
          tsk = qrm_task_type(qrm_task_sol_, front%num, 0, 0, .false.)
          if(qrm_sched_task(tq_h, tsk, 'h')) then
             ! mark the column as assigned
             found = .true.
             status(f) = qrm_busy_
          end if
          
       end if
#if defined (_OPENMP)
       call omp_unset_lock(locks(f))
#endif
    end do

    ! if nothing was found above, then check if the facto is over
    ! otherwise return
    if(found) return

    call check_applyq_over( )

    return
  end subroutine fill_queue_q




!==========================================================================================
!==========================================================================================
  subroutine check_applyq_over( )
    ! checks whether the apply_qt is over and eventually schedules a
    ! termination task
    implicit none

    type(qrm_task_type) :: tsk
    logical             :: found

    ! all the fronts are done
#if defined (_OPENMP)
    call omp_set_lock( dlock )
#endif
    if(dones .eq. fdata%nfronts) then
       tsk = qrm_task_type(qrm_task_exit_, 0, 0, 0, .false.)
       found = qrm_sched_task(tq_h, tsk, 't')
    end if
#if defined (_OPENMP)
    call omp_unset_lock( dlock )
#endif

    return
  end subroutine check_applyq_over






!==========================================================================================
!==========================================================================================
  subroutine apply_q(task, thn)
    implicit none
    type(qrm_task_type) :: task
    integer             :: thn

    type(_qrm_front_type), pointer :: front
    integer                         :: f, p, c, info

    front  => null()
    
    ! to make things easier
    front => qrm_mat%fdata%front_list(task%front)
    info  = 0

    ! write(*,'(i3," =--> Apply Q  : ",i4)')thn, front%num

    call front_q(front, info)
    status(task%front) = qrm_done_
    call qrm_queue_rm(ready_q, front%num)

    !$ call omp_set_lock( dlock )
    dones = dones+1
    !$ call omp_unset_lock( dlock )

    ! sweep over the children. Small children are treated immediately,
    ! the others are pushed on the ready_q
    do p = adata%childptr(front%num), adata%childptr(front%num+1)-1
       c = adata%child(p)
       if(adata%small(c) .eq. 1) then
          call do_subtree_q(c, info)
          if(info .ne. 0) goto 9997
       else
          call qrm_queue_push(ready_q, c)
       end if
    end do

    
9997 continue
    return
  end subroutine apply_q






!==========================================================================================
!==========================================================================================
  subroutine do_subtree_q(fnum, info)
    implicit none

    integer :: fnum, info

    type(_qrm_front_type), pointer :: front
    integer :: node, c, acc, thn, p
    type(qrm_queue) :: sub_q

    info = 0
    acc = 0

    call qrm_queue_init(sub_q, adata%nnodes, qrm_fifo_)
    call qrm_queue_push(sub_q, fnum)

    do
       node = qrm_queue_pop(sub_q)

       if(node .eq. 0) exit

       front => qrm_mat%fdata%front_list(node)
       
       call front_q(front, info)

       !$ call omp_set_lock( dlock )
       dones = dones+1
       !$ call omp_unset_lock( dlock )

       status(node) = qrm_done_

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

  end subroutine do_subtree_q




!==========================================================================================
!==========================================================================================
  subroutine front_q(front, info)
    use _qrm_utils_mod
    use _qrm_rfpf_mod
    use _qrm_spmat_mod
    use qrm_mem_mod
    use qrm_common_mod
    implicit none

    type(_qrm_front_type) :: front
    integer                :: info

    integer :: thn, i
    integer :: pv1, c, k, m, pv2, n, cnt, j, jp, pk
    _qrm_data, allocatable :: work(:,:), in_b(:,:), t(:,:)
    ! error management
    character(len=*), parameter     :: name='front_q'
    
    ! shortcut
    if (min(front%m, front%n) .le. 0) goto 9999

    !$ call omp_set_lock( locks(front%num) )

    thn = 0
    !$ thn = omp_get_thread_num()

    ! allocate everything that's needed
    n = size(b,2)
    call qrm_aalloc(in_b, front%m, n)
    call qrm_aalloc(t, front%nb, front%nb)
    call qrm_aalloc(work, n, front%nb)
    __QRM_CHECK_RET(name,'qrm_aalloc',9999)
    
    ! gather b into front%b
    in_b = b(front%rows,:)
    
#if defined (RFPF)
    ! ! do all the computations here
    ! pv1 = 1
    ! do c = 1, ((front%ne-1)/front%nb)*front%nb, front%nb
       ! k = min(front%nb, front%ne-c+1) 
       ! pv1 = pv1 + k*(k+1)/2 + (max(front%stair(c+k-1),c+k-1) -k-c+1)*k
    ! end do

    ! ! set r to be the starting col of the Q factor
    ! c = ((front%ne-1)/front%nb)*front%nb+1
    ! do 

       ! k = min(front%nb, front%ne-c+1)
       ! m = max(front%stair(c+k-1)-c-k+1,0)

       ! if(m .le. 0) then
          ! pv2 = 1
       ! else
          ! pv2 = pv1+(k*(k+1))/2
       ! end if

       ! call _xrprft('f', 'c', m, k, front%h(pv1), front%h(pv2), m, front%tau(c), t(1,1), front%nb)
       
       ! call _xrprfb('l', 'n', 'f', 'c', m+k, n, k, front%h(pv1), front%h(pv2), m, &
            ! & t(1,1), front%nb, in_b(c,1), front%m, work(1,1), n)

       ! if(c .le. 1) exit
       ! pv1 = pv1-(max(front%stair(c-1),c-1)-c+front%nb+1)*front%nb + front%nb*(front%nb-1)/2
       ! c = c-front%nb

    ! end do

#else

    cnt=size(front%h)+1
    outer: do jp = front%ne - mod(front%ne, front%nb)+1, 1, -front%nb
       pk = min(front%nb, front%ne-jp+1) 
       if(pk .le. 0) cycle
       ! write(*,*)'solve',jp,pk

       inner: do j = jp+pk-mod(pk,front%ib), jp, -front%ib
          k = min(front%ib, jp+pk - j)
          if(k .le. 0) cycle
          m = max(front%stair(j+k-1),j+k-1) - j+1
          cnt = cnt-m*k

          call _xlarfb('l', 'n', 'f', 'c', m, n, k, front%h(cnt), m, &
               & front%h(cnt), m, in_b(j,1), front%m, work(1,1), n)
       end do inner
    end do outer
#endif
    
    ! scatter front%b into b
    b(front%rows,:) = in_b

    ! cleanup
    call qrm_adealloc(in_b)
    call qrm_adealloc(t)
    call qrm_adealloc(work)
    __QRM_CHECK_RET(name,'qrm_adelloc',9999)
    !$ call omp_unset_lock( locks(front%num) )

9999 continue
    return 
    
  end subroutine front_q
  


end subroutine _qrm_apply_q
