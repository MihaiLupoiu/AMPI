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
!> @file qrm_queue_mod.F90
!! This file contains the module that handles the front queues
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains all the facilities for front queues.

!> It is actually meant to be more general than just for
!! fronts. However each queue can contain only a list of elements with
!! unique indices<maxelems, e.g. it can contain this sequence of
!! elements:
!! @verbatim {5 2 8 7 3}, maxelems=8 @endverbatim
!! but not this one
!! @verbatim {5 2 8 5 3}, maxelems=8 @endverbatim
!! nor this one
!! @verbatim {5 2 8 7 3}, maxelems=7 @endverbatim
!!
module qrm_queue_mod

#if defined(_OPENMP)  
   use omp_lib
#endif

   !> parameter to define the policy of the queue: FIFO
  integer, parameter :: qrm_fifo_=0

   !> parameter to define the policy of the queue: LIFO
  integer, parameter :: qrm_lifo_=1 

  !> @brief A data type meant to to define a queue
  type qrm_queue
     integer, allocatable :: elems(:) !> contains the element of the queue
     integer :: nelems                !> number of elements present in the queue
     integer :: h, t                  !> the head/tail of the queue
     integer :: maxelems              !> the maximum number of values allowed in the queue
     integer :: pol                   !> the policy of the queue
#if defined (_OPENMP)
     integer(omp_lock_kind) :: lock !> a lock to prevent simultaneous access to the queue
#endif
  end type qrm_queue

  interface qrm_queue_init
     module procedure qrm_queue_init
  end interface


contains

  !> @brief Initializes a queue
  !!
  !! @param[in,out] q  the queue to be initialized
  !!
  !! @param[in] nelems the max number of elements the queu can
  !!                   hold. This actually corresponds to the range of
  !!                   IDs of the elements that can be pushed on the
  !!                   queue
  !!
  !! @param[in] pol    the policy (either qrm_fifo_ or qrm_lifo_)
  !! 
  subroutine qrm_queue_init(q, nelems, pol)
    use qrm_mem_mod
    implicit none

    type(qrm_queue) :: q
    integer :: nelems, pol

    call qrm_aalloc(q%elems, nelems)
    q%elems    = 0
    q%h        = 0
    q%t        = 0
    q%nelems   = 0
    q%maxelems = nelems
    q%pol      = pol
    !$ call omp_init_lock(q%lock)
    return

  end subroutine qrm_queue_init


  !> @brief Frees a queue
  !!
  !! @param[in,out] q  the queue to be freed
  !!
  subroutine qrm_queue_free(q)
    use qrm_mem_mod
    implicit none

    type(qrm_queue) :: q

    call qrm_adealloc(q%elems)

    q%nelems = 0
    q%h      = 0
    q%t      = 0
    !$ call omp_destroy_lock(q%lock)

    return
  end subroutine qrm_queue_free



  !> @brief Pushes an element on a queue
  !!
  !! @param[in,out] q  the queue where to push
  !!
  !! @param[in]  elem  the element to be pushed
  !! 
  subroutine qrm_queue_push(q, elem)
    implicit none

    type(qrm_queue) :: q
    integer         :: elem
    
    !$ call omp_set_lock(q%lock)
    if (q%nelems .eq. 0) then
       q%h = elem
       q%t = elem
       q%elems(elem) = -1
    else if (q%nelems .eq. q%maxelems) then
       write(*,'("Cannot push anymore", i4,x,i4)')elem, q%nelems
       goto 10
    else
       if(q%pol .eq. qrm_fifo_) then
          q%elems(q%t)  = elem
          q%elems(elem) = -1
          q%t           = elem
       else if(q%pol .eq. qrm_lifo_) then
          q%elems(elem) = q%h
          q%h           = elem
       end if
    end if

    q%nelems = q%nelems+1
    ! write(*,'("Pushed ",i4,x,i4,x,i4)')elem,q%nelems,q%h,q%t
10 continue
    !$ call omp_unset_lock(q%lock)

    return

  end subroutine qrm_queue_push




  !> @brief Prints the content of a queue
  !!
  !! @param[in,out] q  the queue to be printed
  !!
  subroutine qrm_queue_prnt(q)
    implicit none

    type(qrm_queue) :: q

    integer i, elem

    !$ call omp_set_lock(q%lock)
    elem = q%h
    do i=1, q%nelems
       write(*,'(i3,"->")',advance='no')elem
       elem = q%elems(elem)
    end do
    write(*,'(" ")')
    !$ call omp_unset_lock(q%lock)
    return
  end subroutine qrm_queue_prnt


  !> @brief Pops an element from a queue
  !!
  !! @param[in,out] q  the queue where to pop from
  !!
  function qrm_queue_pop(q)
    implicit none

    type(qrm_queue) :: q
    integer         :: qrm_queue_pop

    !$ call omp_set_lock(q%lock)
    if (q%nelems .eq. 0) then
       qrm_queue_pop = 0
    else
       qrm_queue_pop = q%h
       q%h           = q%elems(q%h)
       q%elems(qrm_queue_pop)  = 0
       q%nelems      = q%nelems-1
    end if
    !$ call omp_unset_lock(q%lock)

    return

  end function qrm_queue_pop


  !> @brief Removes (without returning it) an element from a queue
  !!
  !! @param[in,out] q the queue where to rmove from
  !!
  !! @param[in]     n the element to be removed
  !! 
  subroutine qrm_queue_rm(q, n)
    implicit none

    type(qrm_queue) :: q
    integer :: n

    integer :: tmp, i, s


    !$ call omp_set_lock(q%lock)
    if (q%nelems .eq. 0) then
       goto 10
    end if

    if(n .eq. q%h) then
       q%h = q%elems(n)
       q%nelems = q%nelems-1
    else
       tmp = q%h
       do i =1, q%nelems-1
          if (q%elems(tmp) .eq. n) then
             q%elems(tmp) = q%elems(n)
             q%nelems = q%nelems-1
             if(n .eq. q%t) q%t = tmp
             exit
          end if
          tmp = q%elems(tmp)
       end do
       
    end if
    q%elems(n) = 0

10 continue
    !$ call omp_unset_lock(q%lock)

    return

  end subroutine qrm_queue_rm


  !> @brief Returns the element that follows n in the queue q. Very
  !! useful for sweeping through a queue. Example:
  !! @code
  !!  n = 0
  !!  do
  !!     n = qrm_queue_next(q,n)
  !!     if (n .eq. 0) exit ! we went though the whole queue
  !! 
  !!     ! do something on n
  !! 
  !!  end do
  !! @endcode 
  !!
  !! @param[in,out] q the queue where to get the element from
  !!
  !! @param[in]     n the element whose next has to be returned
  !! 
  function qrm_queue_next(q, n)
    implicit none

    type(qrm_queue) :: q
    integer :: qrm_queue_next, n

    integer :: i

#if defined(_OPENMP)
    call omp_set_lock(q%lock)
#endif
    if (q%nelems .eq. 0) then
       qrm_queue_next = 0
    else if (n .eq. q%t) then
       qrm_queue_next = 0
    else

       if(n.eq.0) then
          qrm_queue_next = q%h
       else
          qrm_queue_next = q%elems(n)
       end if
    end if
#if defined(_OPENMP)
    call omp_unset_lock(q%lock)
#endif
    return
  end function qrm_queue_next

end module qrm_queue_mod
