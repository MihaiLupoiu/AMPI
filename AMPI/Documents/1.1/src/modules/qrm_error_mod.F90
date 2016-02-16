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
!> @file qrm_error_mod.F90
!! This file contains the module that implements the error management
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains all the error management routines and 
!! data. 

!> The error management is based on a stack of error messages.
!! Every time an error is detected, the related error code and eventually 
!! a message are pushed onto the errors stack. No error will be raised until one 
!! or all the processor perform a check on the error stack. At this time
!! if there is something on the error stack, two possible actions can be done
!! -# if @link qrm_err_act @endlink=qrm_return_, control is
!!    returned to the calling routine
!! -# if @link qrm_err_act @endlink=qrm_abort_, the execution
!!    is aborted and the content of the stack is dumped on @link qrm_eunit @endlink
!!
!! Having a stack of errors allows to descend the chain of function
!! calls until the leaf level where the error was generated.
!!
!! Basically, upon entrance every routine suba has to save the @link
!! qrm_err_act @endlink value and set it to return so that if a
!! subroutine subb called by suba will return control to suba; at this
!! point suba restores @link qrm_err_act @endlink and performs the
!! corresponding action. Better with an example:
!! @code
!!    subroutine suba()
!!      ! error management
!!      integer                         :: err_act
!!      character(len=*), parameter     :: name='suba'
!!      
!!      call qrm_err_act_save(err_act)
!!     
!!     
!!      call subb()
!!      __QRM_CHECK_RET(name,'subb',9999)
!!      
!!     
!!      call qrm_err_act_restore(err_act)
!!      return
!!      
!!    9999 continue ! error management
!!      call qrm_err_act_restore(err_act)
!!      if(err_act .eq. qrm_abort_) then
!!         call qrm_err_check()
!!      end if
!!     
!!      return
!!     
!!    end subroutine suba
!! @endcode
!!
!! In this example __QRM_CHECK_RET(name,'subb',9999) checks whether an
!! error is present on the error stack and eventually pushed an "error
!! on return" on the stack and then goes to label 9999. This macro is
!! defined in the @link qrm_common.h @endlink fine.
!!

module qrm_error_mod

  !> @brief This is the basic type for error message
  type qrm_err_type

     !> the error code
     integer                     :: code
     !> the name of the subroutine where the error was detected
     character(len=30)           :: sub=''
     !> an array that contains optional integer data for the
     !! error message
     integer, dimension(5)       :: i_err_data=0  
     !> an optional string for the error message
     character(len=40)           :: a_err_data=''      
     !> pointer to the message below in the stack
     type(qrm_err_type), pointer :: next=>null()
  end type qrm_err_type

  !> @brief This type is to represent the errors stack
  type qrm_err_stack_type
     !> Pointer to the top of the stack
     type(qrm_err_type), pointer :: top=>null()
     !> The number of error messages on the stack
     integer                     :: nelem=0
  end type qrm_err_stack_type

  !> @brief The errors stack
  type(qrm_err_stack_type), save  :: qrm_err_stack ! the error messages stack
  !> @brief Possible actions to be performed upon detection of an error
  integer, parameter              :: qrm_abort_=0, qrm_return_=1
  !> @brief Default action
  integer                         :: qrm_err_act=qrm_abort_

contains
  !> @brief Saves a copy of the @link qrm_err_act @endlink variable
  !! @param[out] err_act The variable where to save the value of 
  !!                     @link qrm_err_act @endlink 
  subroutine qrm_err_act_save(err_act)
    integer :: err_act
    err_act = qrm_err_act
    qrm_err_act = qrm_return_
    return
  end subroutine qrm_err_act_save


  !> @brief Sets the default error action
  !! @param[in] err_act The new default error action
  subroutine qrm_err_act_set(err_act)
    integer :: err_act
    if((err_act .ne. qrm_abort_) .and. &
         & (err_act .ne. qrm_return_)) then
       call qrm_err_push(26, 'qrm_err_act_save',ied=(/err_act,0,0,0,0/))
       if(qrm_err_act .eq. qrm_abort_) call qrm_err_check()
    else
       qrm_err_act = err_act
    end if

    return
  end subroutine qrm_err_act_set



  !> @brief Restores the value of the @link qrm_err_act @endlink variable
  !! @param[in] err_act The value to be restored in @link qrm_err_act @endlink
  subroutine qrm_err_act_restore(err_act)
    integer :: err_act
    qrm_err_act = err_act
    return
  end subroutine qrm_err_act_restore

  !> @brief This subroutine pushes an error on top of the stack
  !!
  !! @param[in] code the error code
  !! @param[in] sub  (optional) the name of the subroutine
  !! @param[in] ied  (optional) an array of integers of size 5 
  !!                 containing optional data for the error message
  !! @param[in] aed  (optional) an array of integers of size 5 
  !!                 containing optional data for the error message
  !!
  subroutine qrm_err_push(code, sub, ied, aed)
    
    implicit none

    integer                         :: code
    character(len=*), optional      :: sub
    integer, dimension(5), optional :: ied
    character(len=*), optional      :: aed

    type(qrm_err_type), pointer     :: new_err

    !$omp critical (error)
    allocate(new_err)

    new_err%code = code
    if(present(sub)) new_err%sub        = sub
    if(present(ied)) new_err%i_err_data = ied
    if(present(aed)) new_err%a_err_data = aed

    new_err%next  => qrm_err_stack%top
    qrm_err_stack%top => new_err
    qrm_err_stack%nelem = qrm_err_stack%nelem+1
    !$omp end critical (error)

    return

  end subroutine qrm_err_push


  !> @brief Pushes an error on the stack and the flushes the stack itself.
  !! Basically does err_push and err_check at the same time
  !!
  !! @param[in] code the error code
  !! @param[in] sub  (optional) the name of the subroutine
  !! @param[in] ied  (optional) an array of integers of size 5 
  !!                 containing optional data for the error message
  !! @param[in] aed  (optional) an array of integers of size 5 
  !!                 containing optional data for the error message
  !!
  subroutine qrm_err_raise(code, sub, ied, aed)
    implicit none

    integer                         :: code
    character(len=*), optional      :: sub
    integer, dimension(5), optional :: ied
    character(len=*), optional      :: aed

    call qrm_err_push(code, sub, ied, aed)
    call qrm_err_check()

    return
    
  end subroutine qrm_err_raise
  

  !> @brief This subroutine return the code of the first error on the stack
  !! or zero if the stack is empty
  subroutine qrm_err_get(info)

    implicit none

    integer :: info
    type(qrm_err_type), pointer :: curr
    integer :: i
    ! local check
    
    info = 0
    curr => qrm_err_stack%top
    
    do i=1, qrm_err_stack%nelem-1
       curr => curr%next
    end do

    if(associated(curr)) info = curr%code
       
    return
  end subroutine qrm_err_get


  !> @brief This subroutine checks the errors stack. If something is found
  !! all the entries in the stack are popped and an abort is executed.
  !!
  subroutine qrm_err_check( )
    use qrm_const_mod
    implicit none

    ! local check
    if(qrm_err_stack%nelem .gt. 0) then
       ! errors were previously detected
       __QRM_PRNT_ERR('(" ")')
       __QRM_PRNT_ERR('("==================== Error detected in qr_mumps ====================")')
       call qrm_flush_err_stack(.true.)
       __QRM_PRNT_ERR('("====================================================================")')
       stop
    end if
    
    return

  end subroutine qrm_err_check


  !> @brief This subroutine flushes the errors stack optionally printing
  !! all the messages on the eunit output unit.
  !!
  subroutine qrm_flush_err_stack(prnt)
    implicit none 
    logical, optional :: prnt

    type(qrm_err_type), pointer :: curr, tmp
    logical :: iprnt

    if(present(prnt)) then
       iprnt = prnt
    else
       iprnt = .true.
    end if

    curr => qrm_err_stack%top

    do while (associated(curr))
       if(iprnt) call qrm_process_msg(curr)
       tmp => curr
       curr => curr%next
       deallocate(tmp)
       qrm_err_stack%nelem = qrm_err_stack%nelem-1
    end do

    qrm_err_stack%top => null()

    return

  end subroutine qrm_flush_err_stack


  !> @brief This routine prints out a message on the error unit
  !!
  !! @param[in] msg a qrm_err_type data containing info on the error
  !!                message to be printed
  !!
  subroutine qrm_process_msg(msg)
    use qrm_const_mod
    implicit none

    type(qrm_err_type), pointer :: msg

    __QRM_PRNT_ERR('("Error in subroutine ",a30 " :")')msg%sub
    select case(msg%code)
    case(0)
       __QRM_PRNT_ERR('("Generic error")')
    case(1)
       __QRM_PRNT_ERR('("Sparse matrix format ",a3," is not (yet) supported.")')msg%a_err_data(1:3)
    case(2)
       __QRM_PRNT_ERR('("Symmetric matrices are not supported.")')
    case(3)
       __QRM_PRNT_ERR('("qrm_spmat%cntl is not associated/valid.")')
    case(4)
       __QRM_PRNT_ERR('("Trying to allocate an already allocated array.")')
    case(5)
       __QRM_PRNT_ERR('("Memory allocation problem. Size required: ",i30)')msg%i_err_data(1)
    case(6) 
       __QRM_PRNT_ERR('("Trying to deallocate an unallocated array.")')
    case(7)
       __QRM_PRNT_ERR('("Memory deallocation problem. ",i30)')msg%i_err_data(1)
    case(8)
       __QRM_PRNT_ERR('("Input column permutation not provided/valid")')
    case(9)
       __QRM_PRNT_ERR('("Requested ordering method unknown: ",i3)')msg%i_err_data(1)
    case(10)
       __QRM_PRNT_ERR('("Insufficient size for array: ",a20)')msg%a_err_data(1:20)
    case(11)
       __QRM_PRNT_ERR('("Error in lapack routine: ",i3)')msg%i_err_data(1)
    case(12)
       __QRM_PRNT_ERR('("No more memory available")')
    case(13)
       __QRM_PRNT_ERR('("The analysis must be done before the factorization")')
    case(14)
       __QRM_PRNT_ERR('("The factorization must be done before the solve")')
    case(15)
       __QRM_PRNT_ERR('("This type of norm is not implemented.")')
    case(16)
       __QRM_PRNT_ERR('("Requested ordering method not available: ",a20)')msg%a_err_data
    case(17)
       __QRM_PRNT_ERR('("Error from call to subroutine ",a30)')msg%a_err_data
    case(18)
       __QRM_PRNT_ERR('("COLAMD error ")')
    case(19)
       __QRM_PRNT_ERR('("SCOTCH error ")')
    case(20)
       __QRM_PRNT_ERR('("Factorization error ")')
    case(21)
       __QRM_PRNT_ERR('("Apply error ")')
    case(22)
       __QRM_PRNT_ERR('("Solve error ")')
    case(23)
       __QRM_PRNT_ERR('("Incorrect argument to qrm_set/qrm_get ",a30)')msg%a_err_data
    case(25)
       __QRM_PRNT_ERR('("Problem opening file ",a30)')msg%a_err_data
    case(26)
       __QRM_PRNT_ERR('("Unknown error action ",i10)')msg%i_err_data(1)
    case(27)
       __QRM_PRNT_ERR('("Incompatible values in qrm_spmat%icntl ",i2,2x,i2)')msg%i_err_data(1:2)
    case(28)
       __QRM_PRNT_ERR('("Incorrect value for qrm_nb_/qrm_ib_ ",i2,2x,i2)')msg%i_err_data(1)
    case(29)
       __QRM_PRNT_ERR('("Incorrect value for qrm_spmat%m/n/nz ",i6,2x,i6,2x,i16)')msg%i_err_data(1:3)
    case(30)
       __QRM_PRNT_ERR('("qrm_apply cannot be called if the H matrix is discarded.")')
    case default
       __QRM_PRNT_ERR('("Unknown error code",i4)')msg%code
    end select

    __QRM_PRNT_ERR('(" ")')
    return

  end subroutine qrm_process_msg



end module qrm_error_mod
