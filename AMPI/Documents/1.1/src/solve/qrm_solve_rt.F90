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
!> @file qrm_solve_rt.F90
!! This file contains a routine that solves for R' against multiple vectors
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This function solves for R' against multiple vectors
!! 
!! @param[in] qrm_mat   the main qrm data structure after factorization.
!!
!! @param[in]     b     a 2d array containing the RHS vectors
!!
!! @param[out]    x     a 2d array containing the solution vectors
!!
subroutine _qrm_solve_rt(qrm_mat, b, x)

  use _qrm_spmat_mod
  use _qrm_rfpf_mod
  use qrm_mem_mod
  use qrm_common_mod
  use qrm_task_mod
  use qrm_queue_mod
  use _qrm_solve_mod, protect => _qrm_solve_rt
  use qrm_error_mod
  implicit none

  type(_qrm_spmat_type), target :: qrm_mat
  _qrm_data, intent(inout)  :: b(:,:)
  _qrm_data, intent(out)    :: x(:,:)

  integer                         :: qrm_nth, nth, thn, info, f, dones
  integer                         :: i, c
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
  character(len=*), parameter     :: name='qrm_aply_rt'
  
  call qrm_err_act_save(err_act)

  __QRM_PRNT_DBG('("Solving for R^T")')

  ! simplify 
  adata => qrm_mat%adata
  fdata => qrm_mat%fdata

  ! check if there is a singleton block and eventually solve for it
  if (qrm_mat%adata%ncsing .gt. 0) then
     call _qrm_solve_sing_front(qrm_mat, b, x, 't')
     __QRM_CHECK_RET(name,'qrm_solve_sing_front',9999)
  end if

  ! initialize queue
  call qrm_queue_init(ready_q , adata%nnodes, qrm_lifo_)

  ! initialize to safe values for the case where openmp is not used
  ! the status and locks arrays are needed because multiple solves may
  ! be run in parallel on the same fdata
  call qrm_aalloc(status, qrm_mat%adata%nnodes)
  __QRM_CHECK_RET(name,'qrm_aalloc',9999)
  !$ allocate(locks(adata%nnodes))

  do f = 1, adata%nnodes
     status(f) = 0
     !$ call omp_init_lock(locks(f))
     do i=adata%childptr(f), adata%childptr(f+1)-1
        c = adata%child(i)
        if(adata%small(c) .eq. 0) status(f) = status(f)-1
     end do
  end do

  ! push all the leaves in the ready queue
  do i=adata%nleaves, 1, -1
     call qrm_queue_push(ready_q, adata%leaves(i))
     status(adata%leaves(i)) = qrm_ready_
  end do
  !$ call omp_init_lock(dlock)

  if(adata%ncsing .gt. 0) then
     dones = 1
  else
     dones = 0
  end if

#if defined (_OPENMP)
  call omp_set_num_threads(1)
  qrm_nth = qrm_mat%icntl(qrm_nthreads_)
#endif

  !$omp parallel &
  !$omp & num_threads(qrm_nth) &
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
        call fill_queue_rt( )
     end if
     
     got_task = qrm_get_task(tq_h, task)
     if(.not. got_task) cycle taskloop  ! queue is empty

     select case(task%id)
     case(qrm_task_exit_) 
        !$omp barrier 
        ! done, exit
        exit
     case(qrm_task_app_)
        ! apply a set of Householder vectors
        call solve_rt(task, thn)
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
     call qrm_err_push(21, name)
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
  subroutine fill_queue_rt( )
    implicit none

    type(_qrm_front_type), pointer :: front
    integer   :: f
    integer   :: thn
    logical   :: found
    type(qrm_task_type)             :: tsk

#if defined (_OPENMP)
    thn=omp_get_thread_num()
#else
    thn=0
#endif

    found = .false.

    f = 0
    do
       ! write(*,*)thn,' queue next'
       f = qrm_queue_next(ready_q, f)
       if(f .eq. 0) exit

       front => fdata%front_list(f)

       !$ if(.not. omp_test_lock(locks(f))) cycle
       if(status(f) .eq. qrm_ready_) then
          tsk = qrm_task_type(qrm_task_app_, front%num, 0, 0, .false.)
          if(qrm_sched_task(tq_h, tsk, 'h')) then
             ! mark the column as assigned
             found = .true.
             status(f) = qrm_busy_
          end if
          
       end if
       !$ call omp_unset_lock(locks(f))
    end do

    ! if nothing was found above, then check if the facto is over
    ! otherwise return
    if(found) return

    call check_solvert_over( )

    return
  end subroutine fill_queue_rt




!==========================================================================================
!==========================================================================================
  subroutine check_solvert_over( )
    ! checks whether the apply_qt is over and eventually schedules a
    ! termination task
    implicit none

    type(qrm_task_type) :: tsk
    logical             :: found

    ! all the fronts are done
    if(dones .eq. fdata%nfronts) then
       tsk = qrm_task_type(qrm_task_exit_, 0, 0, 0, .false.)
       found = qrm_sched_task(tq_h, tsk, 't')
    end if

    return
  end subroutine check_solvert_over



  
!==========================================================================================
!==========================================================================================
  subroutine solve_rt(task, thn)
    implicit none
    type(qrm_task_type) :: task
    integer         :: thn

    type(_qrm_front_type), pointer :: front
    integer                         :: f, p, c, info

    front  => null()
    
    ! to make things easier
    front => qrm_mat%fdata%front_list(task%front)
    f     =  qrm_mat%adata%parent(task%front)
    info  = 0
      
    ! sweep over the children to check whether there are small subtrees
    ! that have to be treated
    do p = adata%childptr(front%num), adata%childptr(front%num+1)-1
       c = adata%child(p)
       if(adata%small(c) .eq. 1) call do_subtree_rt(c, info)
       if(info .ne. 0) goto 9997
    end do

    call front_rt(front, info)
    status(task%front) = qrm_done_

    !$ call omp_set_lock( dlock )
    dones = dones+1
    !$ call omp_unset_lock( dlock )

    if(f .ne. 0) then
       !$ call omp_set_lock( locks(f) )
       status(f) = status(f)+1
       if(status(f) .eq. qrm_ready_) call qrm_queue_push(ready_q, f)
       !$ call omp_unset_lock( locks(f) )
    end if
    call qrm_queue_rm(ready_q, task%front)

9997 continue
    return
  end subroutine solve_rt

  

!==========================================================================================
!==========================================================================================
  subroutine do_subtree_rt(fnum, info)
    implicit none

    integer :: fnum, info

    type(_qrm_front_type), pointer :: front
    integer :: node, f
    
    info = 0
    
    node = fnum
    
    subtree: do
     
       front => fdata%front_list(node)
     
       if (status(node) .eq. qrm_ready_) then
          ! the front is ready to be assembled and factorized
          
          ! factorize
          call front_rt(front, info)
          if(info .ne. 0) goto 9998

          !$ call omp_set_lock( dlock )
          dones = dones+1
          !$ call omp_unset_lock( dlock )

          f = adata%parent(node)
          if(f .ne. 0) then
             !$ call omp_set_lock( locks(f) )
             status(f) = status(f)+1
             !$ call omp_unset_lock( locks(f) )
          end if

          status(node) = qrm_done_
          
          if(node .eq. fnum) exit subtree  ! we're done
          
          node = adata%parent(node)
          if(node .eq. 0) exit subtree
          
       else
          ! go down the subtree
          node = adata%child(adata%childptr(node+1)+status(node))
          cycle subtree
       end if
       
       
    end do subtree
    
9998 continue
    return
  end subroutine do_subtree_rt




!==========================================================================================
!==========================================================================================
  subroutine front_rt(front, info)
    use _qrm_utils_mod
    use _qrm_rfpf_mod
    use _qrm_spmat_mod
    use qrm_mem_mod
    use qrm_common_mod
    implicit none

    type(_qrm_front_type) :: front
    integer                :: info


    integer :: thn

#if defined (sprec) || defined (dprec)
    character, parameter :: tr='t', notr='n'
#elif defined (cprec) || defined (zprec)
    character, parameter :: tr='c', notr='n'
#endif

    integer :: pv1, c, k, m, pv2, n, r, i, cnt, jp, pk, j
    _qrm_data, allocatable :: in_b(:,:)
    ! error management
    character(len=*), parameter     :: name='solve_rt'

#if defined (_OPENMP)    
    thn=omp_get_thread_num()
#else
    thn = 0
#endif
    
    ! shortcut
    if (min(front%m, front%n) .le. 0) goto 10
    if (front%npiv .le. 0) goto 10

    !$ call omp_set_lock( locks(front%num) )
    ! write(*,'(i3," =--> Solve R''  : ",i4)')thn, front%num

    ! allocate everything that's needed
    call qrm_aalloc(in_b, front%n, size(x,2))
    __QRM_CHECK_RET(name,'qrm_aalloc',9999)

    ! do computations
    in_b(1:front%npiv,:)         = b(front%cols(1:front%npiv),:)
    in_b(front%npiv+1:front%n,:) = _qrm_zero

#if defined (RFPF)
    ! ! pv1 = 1
    ! r=1
    ! do 
    !    if (r .gt. front%npiv) then
    !       exit
    !    else
    !       k = min(front%nb, front%npiv-r+1)
    !    end if

    !    n = front%n - r-k+1
    !    if(n .le. 0) then
    !       pv2 = 1
    !    else
    !       pv2 = pv1+(k*(k+1))/2
    !    end if

    !    call _xrpsm('l', 'u', tr, 'n', k, size(b,2), _qrm_one, front%r(pv1), in_b(r,1), front%n)

    !    !FIXME: why _qrm_one in the BETA parameter???
    !    if(n .gt. 0) call _xgemm(tr, notr, n, size(b,2), k, -_qrm_one, front%r(pv2), &
    !         &k, in_b(r,1), front%n, _qrm_one, in_b(r+k,1), front%n)
       
    !    pv1 = pv1+(k*(k+1))/2 + k*n
    !    r = r+k
    ! end do
#else
    n = size(in_b,2)
    cnt=1
    outer: do jp = 1, front%npiv, front%nb
       pk = min(front%nb, front%npiv-jp+1) 
       if(pk .le. 0) exit outer
       
       inner: do j = jp, jp+pk-1, front%ib
          m = min(front%ib, jp+pk - j)
          if(m .le. 0) exit inner
          k = front%n-j-m+1

          ! write(*,*)thn !,front%num,'++>',cnt, j, k+m, m
          call _xtrsm('l', 'u', tr, 'n', m, n, _qrm_one, front%r(cnt), m, &
               & in_b(j,1), front%n)
          if(k.gt.0) call _xgemm(tr, notr, k, n, m, -_qrm_one, front%r(cnt+m*m), m, &
               & in_b(j,1), front%n, _qrm_one, in_b(j+m,1), front%n)
          cnt = cnt+m*(front%n-j+1)
       end do inner
    end do outer
#endif  


    !FIXME check if this critical region can be replaced with a lock of some sort
    !$omp critical (sc_gt)
    x(front%rows(1:front%npiv),:) = in_b(1:front%npiv,:)
    b(front%cols(front%npiv+1:front%n),:) = b(front%cols(front%npiv+1:front%n),:)+in_b(front%npiv+1:front%n,:)
    !$omp end critical (sc_gt)

    call qrm_adealloc(in_b)
    __QRM_CHECK_RET(name,'qrm_adealloc',9999)

    !$ call omp_unset_lock( locks(front%num) )

10  continue

9999 continue ! error
    return
  end subroutine front_rt


end subroutine _qrm_solve_rt
