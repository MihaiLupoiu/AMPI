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
!> @file qrm_coverage.F90
!! This file contains coverage tests for the solve operation
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

subroutine _qrm_test_slv(c, m)
  use _qrm_coverage_mod, protect => _qrm_test_slv
  implicit none

  integer :: c, m

  integer, parameter :: ncases=15
  logical :: cases(ncases)

  if(c .eq. -1) then
     cases = .true.
  else if (c .le. ncases) then
     cases = .false.
     cases(c) = .true.
  end if

  if(cases(1 )) call case1 (m)

  return


contains

  !> @brief Test the multiple rhs parallelism
  subroutine case1(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:,:), b(:,:), r(:,:)
    _qrm_real                      :: err
    integer                        :: i, j, nb(4)

    nb=(/1, 2, 4, 5/)

    if((m.ne.-1) .and. &
         & (m.lt.11) .and. &
         & (m.gt.size(matrices))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. &
            & (m.ne.i)) cycle
       nullify(x, r, b)
       
       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       call qrm_palloc(x,qrm_spmat%n,20)
       call qrm_palloc(b,qrm_spmat%m,20)
       call qrm_palloc(r,qrm_spmat%m,20)
       b = _qrm_one
       ! set the number of threads
       call qrm_set(qrm_spmat,'qrm_nthreads',2)

       do j=1, 4

          call qrm_set(qrm_spmat,'qrm_rhsnthreads',1)
          call qrm_set(qrm_spmat,'qrm_rhsnb',nb(j))

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)
       
          ! print message
          call _qrm_prnt_testmesg(6, 1, j, i, (info .eq. 0) .and. (err .lt. eps))
       
          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)
       end do

       call qrm_pdealloc(x)
       call qrm_pdealloc(b)
       call qrm_pdealloc(r)
    end do
    
    return

  end subroutine case1
  

end subroutine _qrm_test_slv
