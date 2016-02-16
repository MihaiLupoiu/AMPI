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
!! This file contains coverage tests for the C interface
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

subroutine _qrm_test_c(c, m)
  use _qrm_coverage_mod, protect => _qrm_test_c
  use _qrm_c_interface
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

  !> @brief Test C problem solve
  subroutine case1(m)
    use _qrm_mod
    use _qrm_c_interface
    implicit none

    interface
       function _qrm_test_solve_c(qrm_spmat_c, b_c, x_c, r_c) bind(c)
         use iso_c_binding
         use _qrm_c_interface
         _qrm_real_fc            :: _qrm_test_solve_c
         type(_qrm_spmat_type_c) :: qrm_spmat_c
         type(c_ptr), value      :: b_c, x_c, r_c
       end function _qrm_test_solve_c
    end interface
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr)                    :: b_c, x_c, r_c


    if((m.ne.-1) .and. ((m.lt.11).and.(m.gt.size(matrices)))) then
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end if

    ! loop over all the file-matrices
    do i=11, size(matrices)
       if((m.ne.-1) .and. (i.ne.m)) cycle

       ! get the data
       qrm_spmat => _qrm_get_test_mat(i)
       b         => _qrm_get_test_b(i)
       x         => _qrm_get_test_x(i)
       r         => _qrm_get_test_r(i)
          
       call _qrm_spmat_init_c(qrm_spmat_c)

       qrm_spmat_c%m   = qrm_spmat%m
       qrm_spmat_c%n   = qrm_spmat%n
       qrm_spmat_c%nz  = qrm_spmat%nz
       ! FIXME: dangerous, verify
       qrm_spmat_c%irn = c_loc(qrm_spmat%irn(1))
       qrm_spmat_c%jcn = c_loc(qrm_spmat%jcn(1))
       qrm_spmat_c%val = c_loc(qrm_spmat%val(1))
       b_c             = c_loc(b(1))
       x_c             = c_loc(x(1))
       r_c             = c_loc(r(1))
          
       ! solve and get back the error
       err = _qrm_test_solve_c(qrm_spmat_c, b_c, x_c, r_c)
       
       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)
       
       ! print message
       call _qrm_prnt_testmesg(3, 1, 1, i, (info .eq. 0) .and. (err .lt. eps))
          
       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)

    end do

    return

  end subroutine case1
  



end subroutine _qrm_test_c
