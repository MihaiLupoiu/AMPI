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
!! This file contains coverage tests for the orderings
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

subroutine _qrm_test_ord(c, m)
  use _qrm_coverage_mod, protect => _qrm_test_ord
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
  if(cases(2 )) call case2 (m)
  if(cases(3 )) call case3 (m)
  if(cases(4 )) call case4 (m)
  if(cases(5 )) call case5 (m)



  return


contains

  !> @brief Test natural ordering
  subroutine case1(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i

    select case (m)
    case (-1)
       ! loop over all the file-matrices
       do i=11, size(matrices)
          ! get the data
          qrm_spmat => _qrm_get_test_mat(i)
          b         => _qrm_get_test_b(i)
          x         => _qrm_get_test_x(i)
          r         => _qrm_get_test_r(i)
          
          ! set the ordering
          call qrm_set(qrm_spmat,'qrm_ordering',qrm_natural_)

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)
          
          ! print message
          call _qrm_prnt_testmesg(2, 1, 1, i, (info .eq. 0) .and. (err .lt. eps))

          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)
       end do
       
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       ! get the data
       qrm_spmat => _qrm_get_test_mat(m)
       b         => _qrm_get_test_b(m)
       x         => _qrm_get_test_x(m)
       r         => _qrm_get_test_r(m)

       ! set the ordering
       call qrm_set(qrm_spmat,'qrm_ordering',qrm_natural_)

       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)

       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)
       ! print message
       call _qrm_prnt_testmesg(2, 1, 1, m, (info .eq. 0) .and. (err .lt. eps))
       write(*,*)err, eps
       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)
    case default
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end select

    return

  end subroutine case1
  


  !> @brief Test given ordering
  subroutine case2(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i, s

    select case (m)
    case (-1)
       ! loop over all the file-matrices
       do i=11, size(matrices)
          ! get the data
          qrm_spmat => _qrm_get_test_mat(i)
          b         => _qrm_get_test_b(i)
          x         => _qrm_get_test_x(i)
          r         => _qrm_get_test_r(i)

          ! just a simple column permutation
          s = min(qrm_spmat%n,qrm_spmat%m)
          call qrm_palloc(qrm_spmat%cperm_in, s)
          qrm_spmat%cperm_in = (/(i,i=1,s)/)
          qrm_spmat%cperm_in(1) = s
          qrm_spmat%cperm_in(s) = 1

          ! set the ordering
          call qrm_set(qrm_spmat,'qrm_ordering',qrm_given_)

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)
          
          ! print message
          call _qrm_prnt_testmesg(2, 2, 1, i, (info .eq. 0) .and. (err .lt. eps))

          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)

          ! deallocate cperm_in
          call qrm_pdealloc(qrm_spmat%cperm_in)
       end do
       
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       ! get the data
       qrm_spmat => _qrm_get_test_mat(m)
       b         => _qrm_get_test_b(m)
       x         => _qrm_get_test_x(m)
       r         => _qrm_get_test_r(m)

       ! just a simple column permutation
       s = min(qrm_spmat%n,qrm_spmat%m)
       call qrm_palloc(qrm_spmat%cperm_in, s)
       qrm_spmat%cperm_in = (/(i,i=1,s)/)
       qrm_spmat%cperm_in(1) = s
       qrm_spmat%cperm_in(s) = 1

       ! set the ordering
       call qrm_set(qrm_spmat,'qrm_ordering',qrm_given_)

       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)

       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)

       ! print message
       call _qrm_prnt_testmesg(2, 2, 1, m, (info .eq. 0) .and. (err .lt. eps))

       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)

       ! deallocate cperm_in
       call qrm_pdealloc(qrm_spmat%cperm_in)
    case default
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end select

    return

  end subroutine case2
  


  !> @brief Test colamd ordering
  subroutine case3(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i

    select case (m)
    case (-1)
       ! loop over all the file-matrices
       do i=11, size(matrices)
          ! get the data
          qrm_spmat => _qrm_get_test_mat(i)
          b         => _qrm_get_test_b(i)
          x         => _qrm_get_test_x(i)
          r         => _qrm_get_test_r(i)
          
          ! set the ordering
          call qrm_set(qrm_spmat,'qrm_ordering',qrm_colamd_)

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)

          ! print message
#if defined(have_colamd)
          call _qrm_prnt_testmesg(2, 3, 1, i, (info .eq. 0) .and. (err .lt. eps))
#else
          call _qrm_prnt_testmesg(2, 3, 1, i, (info .eq. 16))
#endif
          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)
       end do
       
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       ! get the data
       qrm_spmat => _qrm_get_test_mat(m)
       b         => _qrm_get_test_b(m)
       x         => _qrm_get_test_x(m)
       r         => _qrm_get_test_r(m)

       ! set the ordering
       call qrm_set(qrm_spmat,'qrm_ordering',qrm_colamd_)

       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)

       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)


       ! print message
#if defined(have_colamd)
       call _qrm_prnt_testmesg(2, 3, 1, m, (info .eq. 0) .and. (err .lt. eps))
#else
       call _qrm_prnt_testmesg(2, 3, 1, m, (info .eq. 16))
#endif

       ! call _qrm_spmat_destroy(qrm_spmat)
       call qrm_adata_destroy(qrm_spmat%adata)
       call _qrm_fdata_destroy(qrm_spmat%fdata)
       ! deallocate(matrices(m)%qrm_spmat)

       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)
    case default
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end select

    return

  end subroutine case3
  


  !> @brief Test metis ordering
  subroutine case4(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i

    select case (m)
    case (-1)
       ! loop over all the file-matrices
       do i=11, size(matrices)
          ! get the data
          qrm_spmat => _qrm_get_test_mat(i)
          b         => _qrm_get_test_b(i)
          x         => _qrm_get_test_x(i)
          r         => _qrm_get_test_r(i)
          
          ! set the ordering
          call qrm_set(qrm_spmat,'qrm_ordering',qrm_metis_)

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)
          
          ! print message
#if defined(have_metis)
          call _qrm_prnt_testmesg(2, 4, 1, i, (info .eq. 0) .and. (err .lt. eps))
#else
          call _qrm_prnt_testmesg(2, 4, 1, i, (info .eq. 16))
#endif
          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)
       end do
       
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       ! get the data
       qrm_spmat => _qrm_get_test_mat(m)
       b         => _qrm_get_test_b(m)
       x         => _qrm_get_test_x(m)
       r         => _qrm_get_test_r(m)

       ! set the ordering
       call qrm_set(qrm_spmat,'qrm_ordering',qrm_metis_)

       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)

       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)

       ! print message
#if defined(have_metis)
          call _qrm_prnt_testmesg(2, 4, 1, m, (info .eq. 0) .and. (err .lt. eps))
#else
          call _qrm_prnt_testmesg(2, 4, 1, m, (info .eq. 16))
#endif

       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)
    case default
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end select

    return

  end subroutine case4
  



  !> @brief Test scotch oerdering
  subroutine case5(m)
    use _qrm_mod
    implicit none
    
    integer :: m

    type(_qrm_spmat_type), pointer :: qrm_spmat
    integer                        :: info
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: err
    integer                        :: i

    select case (m)
    case (-1)
       ! loop over all the file-matrices
       do i=11, size(matrices)
          ! get the data
          qrm_spmat => _qrm_get_test_mat(i)
          b         => _qrm_get_test_b(i)
          x         => _qrm_get_test_x(i)
          r         => _qrm_get_test_r(i)
          
          ! set the ordering
          call qrm_set(qrm_spmat,'qrm_ordering',qrm_scotch_)

          ! solve and get back the error
          err = _qrm_test_solve(qrm_spmat, b, x, r)
          
          ! check if anything went wrong
          call qrm_err_get(info)
          call qrm_flush_err_stack(.false.)
          
          ! print message
#if defined(have_scotch)
          call _qrm_prnt_testmesg(2, 5, 1, i, (info .eq. 0) .and. (err .lt. eps))
#else
          call _qrm_prnt_testmesg(2, 5, 1, i, (info .eq. 16))
#endif
          ! put cntl back in its initial state
          call _qrm_cntl_init(qrm_spmat)
       end do
       
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       ! get the data
       qrm_spmat => _qrm_get_test_mat(m)
       b         => _qrm_get_test_b(m)
       x         => _qrm_get_test_x(m)
       r         => _qrm_get_test_r(m)

       ! set the ordering
       call qrm_set(qrm_spmat,'qrm_ordering',qrm_scotch_)

       ! solve and get back the error
       err = _qrm_test_solve(qrm_spmat, b, x, r)

       ! check if anything went wrong
       call qrm_err_get(info)
       call qrm_flush_err_stack(.false.)

          ! print message
#if defined(have_scotch)
          call _qrm_prnt_testmesg(2, 5, 1, m, (info .eq. 0) .and. (err .lt. eps))
#else
          call _qrm_prnt_testmesg(2, 5, 1, m, (info .eq. 16))
#endif

       ! put cntl back in its initial state
       call _qrm_cntl_init(qrm_spmat)
    case default
       write(*,'("Matrix ",i2," is not available for this test")')m
       return
    end select

    return

  end subroutine case5

end subroutine _qrm_test_ord
