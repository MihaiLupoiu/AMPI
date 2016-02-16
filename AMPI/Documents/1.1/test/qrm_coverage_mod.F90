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
!> @file qrm_coverage_mod.F90
!! Various tools for the coverage test
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

module _qrm_coverage_mod
  use _qrm_spmat_mod

  interface 
     subroutine _qrm_test_err(c)
       integer :: c
     end subroutine _qrm_test_err
  end interface

  interface 
     subroutine _qrm_test_ord(c, m)
       integer :: c, m
     end subroutine _qrm_test_ord
  end interface

  interface 
     subroutine _qrm_test_c(c, m)
       integer :: c, m
     end subroutine _qrm_test_c
  end interface

  interface 
     subroutine _qrm_test_sing(c, m)
       integer :: c, m
     end subroutine _qrm_test_sing
  end interface

  interface 
     subroutine _qrm_test_facto(c, m)
       integer :: c, m
     end subroutine _qrm_test_facto
  end interface

  interface 
     subroutine _qrm_test_slv(c, m)
       integer :: c, m
     end subroutine _qrm_test_slv
  end interface

  interface _qrm_test_solve
     module procedure _qrm_test_solve1d, _qrm_test_solve2d
  end interface _qrm_test_solve

  type tmat
     character(len=50)                  :: mfile
     type(_qrm_spmat_type), allocatable :: qrm_spmat
     _qrm_data, allocatable             :: x(:), b(:), r(:)
  end type tmat

  type(tmat), target, allocatable :: matrices(:)

#if defined(dprec) || defined(zprec)
  _qrm_real, parameter :: eps=10d-9
#elif defined(sprec) || defined(cprec)
  _qrm_real, parameter :: eps=10e-3
#endif


contains

  subroutine _qrm_prnt_testmesg(t, c, s, m, res)
    integer :: t, m, c, s
    logical :: res

    if(res) then
       write(*,'("TEST ",i2,"     CASE ",i2,"     SUBCASE ",i2,"     MATRIX ",i2," : OK")')t,c,s,m
    else
       write(*,'("TEST ",i2,"     CASE ",i2,"     SUBCASE ",i2,"     MATRIX ",i2," : FAILED")')t,c,s,m
    end if

    return
  end subroutine _qrm_prnt_testmesg
  


  subroutine read_matfile()
    use _qrm_spmat_mod
    implicit none

    integer :: info, nmats, i

    open(4,file='matfile.txt', status='OLD', action='READ', iostat=info)
    read(4,*)nmats
    
    allocate(matrices(nmats+10))

    do i=1, nmats
       read(4,*)matrices(10+i)%mfile
    end do

    return

  end subroutine read_matfile
  

  function _qrm_get_test_mat(m)
    use _qrm_utils_mod
    implicit none
    integer :: m
    type(_qrm_spmat_type), pointer :: _qrm_get_test_mat
    integer :: iseed(4)

    iseed = (/1,1,1,1/)
    select case (m)
    case(1)
       if(.not. allocated(matrices(1)%qrm_spmat)) then
          allocate(matrices(1)%qrm_spmat)
          call _qrm_spmat_init(matrices(1)%qrm_spmat)
          call _qrm_spmat_alloc(matrices(1)%qrm_spmat, 5, 4, 4, 'coo')
          matrices(1)%qrm_spmat%irn=(/1, 2, 3, 4, 4/)
          matrices(1)%qrm_spmat%jcn=(/1, 2, 3, 4, 3/)
          call _xlarnv(2,iseed,matrices(1)%qrm_spmat%nz, matrices(1)%qrm_spmat%val)
          ! call random_number(matrices(1)%qrm_spmat%val)
       end if
       _qrm_get_test_mat => matrices(1)%qrm_spmat
    case(11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if
       if(.not. allocated(matrices(m)%qrm_spmat)) then
          allocate(matrices(m)%qrm_spmat)
          call _qrm_spmat_init(matrices(m)%qrm_spmat)
          call qrm_readmat(matrices(m)%mfile, matrices(m)%qrm_spmat, .true.)
       end if
       _qrm_get_test_mat => matrices(m)%qrm_spmat
    case default
       write(*,'("Matrix ",i2," does not exist")')m
       return
    end select

    return

  end function _qrm_get_test_mat
  

  function _qrm_get_test_b(m)
    use _qrm_utils_mod
    implicit none
    integer :: m
    _qrm_data, pointer :: _qrm_get_test_b(:)
    type(_qrm_spmat_type), pointer :: qrm_spmat

    select case (m)
    case(1,11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       qrm_spmat => _qrm_get_test_mat(m)
       if(.not. allocated(matrices(m)%b)) then
          call qrm_aalloc(matrices(m)%b, qrm_spmat%m)
          ! call random_number(matrices(m)%b)
          matrices(m)%b = _qrm_one
       end if
       _qrm_get_test_b => matrices(m)%b
    case default
    end select

    return

  end function _qrm_get_test_b
  
  

  function _qrm_get_test_x(m)
    use _qrm_utils_mod
    implicit none
    integer :: m
    _qrm_data, pointer :: _qrm_get_test_x(:)
    type(_qrm_spmat_type), pointer :: qrm_spmat

    select case (m)
    case(1,11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       qrm_spmat => _qrm_get_test_mat(m)
       if(.not. allocated(matrices(m)%x)) then
          call qrm_aalloc(matrices(m)%x, qrm_spmat%n)
          matrices(m)%x = _qrm_zero
       end if
       _qrm_get_test_x => matrices(m)%x
    case default
    end select

    return

  end function _qrm_get_test_x



  function _qrm_get_test_r(m)
    use _qrm_utils_mod
    implicit none
    integer :: m
    _qrm_data, pointer :: _qrm_get_test_r(:)
    type(_qrm_spmat_type), pointer :: qrm_spmat

    select case (m)
    case(1,11:)
       if(m .gt. size(matrices)) then
          write(*,'("Matrix ",i2," does not exist")')m
          return
       end if

       qrm_spmat => _qrm_get_test_mat(m)
       if(.not. allocated(matrices(m)%r)) then
          call qrm_aalloc(matrices(m)%r, qrm_spmat%m)
          matrices(m)%r = _qrm_zero
       end if
       _qrm_get_test_r => matrices(m)%r
    case default
    end select

    return

  end function _qrm_get_test_r


  function _qrm_test_solve1d(qrm_spmat, b, x, r)
    use _qrm_mod
    implicit none 

    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, pointer             :: x(:), b(:), r(:)
    _qrm_real                      :: _qrm_test_solve1d

    integer                        :: info

    r = b
    if(qrm_spmat%m .ge. qrm_spmat%n) then
       call qrm_least_squares(qrm_spmat, r, x)
    else
       call qrm_min_norm(qrm_spmat, r, x)
    end if

    call qrm_err_get(info)
    if(info .ne. 0) return
    
    _qrm_test_solve1d  = _qrm_rzero

    r = b
    ! FIXME: check this 
    call qrm_residual_norm(qrm_spmat, r, x, _qrm_test_solve1d)
    if(_qrm_test_solve1d .gt. eps) then
       call qrm_residual_orth(qrm_spmat, r, _qrm_test_solve1d)
    end if
    
    
    return

  end function _qrm_test_solve1d

  function _qrm_test_solve2d(qrm_spmat, b, x, r)
    use _qrm_mod
    implicit none 

    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, pointer             :: x(:,:), b(:,:), r(:,:)
    _qrm_real                      :: _qrm_test_solve2d
    
    _qrm_real, allocatable         :: norms(:)
    integer                        :: info, i

    r = b
    if(qrm_spmat%m .ge. qrm_spmat%n) then
       call qrm_least_squares(qrm_spmat, r, x)
    else
       call qrm_min_norm(qrm_spmat, r, x)
    end if

    ! check if anything went wrong
    call qrm_err_get(info)
    if(info .ne. 0) return
    
    _qrm_test_solve2d  = _qrm_rzero

    call qrm_aalloc(norms, size(x,2))

    r = b

    ! FIXME: check this 
    ! do i=1, size(x,2)
       call qrm_residual_norm(qrm_spmat, r(:,:), x(:,:), norms)
    ! end do
    _qrm_test_solve2d = maxval(norms)

    if(_qrm_test_solve2d .gt. eps) then
       ! do i=1, size(x,2)
          call qrm_residual_orth(qrm_spmat, r(:,:), norms)
       ! end do
       _qrm_test_solve2d = maxval(norms)
    end if

    call qrm_adealloc(norms)
    return

  end function _qrm_test_solve2d


end module _qrm_coverage_mod
