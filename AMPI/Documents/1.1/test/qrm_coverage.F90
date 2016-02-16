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
!! This file contains a coverage test program to check as many features and lines of
!! code as possible.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


program _qrm_coverage
  use _qrm_mod
  use _qrm_coverage_mod
  implicit none

  integer :: t, c, m
  integer :: nargs, i
  integer, parameter :: ntests=6
  logical :: tests(ntests)


  call get_args(t,c,m)
  call read_matfile()

  call qrm_set('qrm_dunit',-1)
  call qrm_set('qrm_eunit',-1)
  call qrm_set('qrm_ounit',-1)

  write(*,'("Matrices used for the test")')
  do i=11, size(matrices)
     write(*,'(i2," -- ",a30)')i,matrices(i)%mfile
  end do

  if(t .eq. -1) then
     tests = .true.
  else
     tests = .false.
     tests(t) = .true.
  end if

  
  call qrm_err_act_set(qrm_return_)
  ! call qrm_err_act_set(qrm_abort_)

  if(tests(1)) call _qrm_test_err(c)
  if(tests(2)) call _qrm_test_ord(c, m)
  if(tests(3)) call _qrm_test_c(c, m)
  if(tests(4)) call _qrm_test_sing(c, m)
  if(tests(5)) call _qrm_test_facto(c, m)
  if(tests(6)) call _qrm_test_slv(c, m)

  stop

contains
  
  subroutine get_args(t, c, m)
    
    integer :: t, c, m
    
    character(len=50) :: str
    integer :: idx, len, i

    t = -1
    m = -1
    c = -1
    
    nargs = command_argument_count()
    if(nargs .gt. 0) then
       call get_command_argument(1,value=str,length=len)
       idx = index(str,'-h')
       if(idx .eq. 1) then
          write(*,'("============= _qrm_coverage usage =============")')
          stop
       end if
    end if

    i = 1
    do
       if(i .gt. nargs) exit

       call get_command_argument(i,value=str,length=len)
       select case(str(1:2))
       case('-t')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)t
       case('-m')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)m
       case('-c')
          i = i+1
          if(i .gt. nargs) exit
          call get_command_argument(i,value=str,length=len)
          read(str,*)c
       case default
          write(*,'("Unrecognized option (try with -h)")')
       end select
       i = i+1
    end do

    return

  end subroutine get_args
  


end program _qrm_coverage
