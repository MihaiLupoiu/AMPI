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
!> @file qrm_trace_mod.F90
!! This file contains a module for visualizing the execution profile of a parallel code
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains all the facilities for visualizing the
!! execution profile of a parallel code
!! 
module qrm_trace_mod
  use qrm_common_mod
  use qrm_error_mod

  public qrm_trace_init, qrm_trace_create_event, &
       & qrm_trace_event_start, qrm_trace_event_stop,&
       & qrm_trace_log_dump

  interface qrm_trace_init
     module procedure qrm_trace_init
  end interface

  interface qrm_trace_create_event
     module procedure qrm_trace_create_event
  end interface

  interface qrm_trace_event_start
     module procedure qrm_trace_event_start
  end interface

  interface qrm_trace_event_stop
     module procedure qrm_trace_event_stop
  end interface

  interface qrm_trace_log_dump
     module procedure qrm_trace_log_dump
  end interface

  private

  type event_type
     integer :: id, thn
     real(kind(1.d0)) :: start, stop
  end type event_type

  real(kind(1.d0)), save :: timezero, start, stop
  integer, parameter :: maxevents=15000, maxtypes=20, maxth=32
  logical, save :: pendings(0:maxth-1)
  real(kind(1.d0)), save :: starts(0:maxth-1), stops(0:maxth-1), ttimes(1:maxtypes)

  type(event_type), allocatable, save :: events(:,:)
  character(len=20), save :: labels(maxtypes)
  integer, save :: nevtype, nodeid
  integer, save :: nevents(0:maxth-1)
  character(len=7) :: colors(maxtypes)

contains

  subroutine qrm_trace_init(node)
    integer :: node

    nodeid = node
    pendings = .false.
    nevtype  = 0
    nevents  = 0
    ttimes   = 0
    allocate(events(0:maxth-1,maxevents))
    colors(1:7) = (/'#d38d5f', '#ffdd55', '#8dd35f', '#80b3ff', '#e580ff', '#ac9d93', '#bcd35f'/)
    timezero = qrm_swtime()
    return

  end subroutine qrm_trace_init
  

  
  subroutine qrm_trace_create_event(label, id)
    character :: label*(*)
    integer   :: id

    nevtype = nevtype+1
    id      = nevtype
    labels(id) = label

    return
    
  end subroutine qrm_trace_create_event
  
  
  subroutine qrm_trace_event_start(id, thn)
    integer :: id, thn
    if(pendings(thn)) then
       __QRM_PRNT_ERR('("Tracing error!!! events nesting not supported")')
       return
    end if
    pendings(thn) = .true.
    starts(thn) = qrm_swtime()
    return
  end subroutine qrm_trace_event_start
  
  subroutine qrm_trace_event_stop(id, thn)
    integer :: id, thn
    
    stops(thn) = qrm_swtime()
    nevents(thn) = nevents(thn)+1
    events(thn, min(nevents(thn),maxevents)) = event_type(id, thn, starts(thn)-timezero, stops(thn)-timezero)
    ttimes(id) = ttimes(id)+stops(thn)-starts(thn)
    pendings(thn) = .false.

    return
  end subroutine qrm_trace_event_stop
  

  subroutine qrm_trace_log_dump(ofile)

    character :: ofile*(*)

    integer :: i, j
    real(kind(1.d0)) :: tottime
    real(kind(1.d0)), parameter :: h=20.d0, scale=10000

    open(4, file=ofile, action='write')
    
    write(4,'("<svg>")')
    do i=0, maxth-1
       do j=1, min(nevents(i),maxevents)
          write(4,'("<rect style=""fill:",a7,";stroke:#000000;stroke-width:0;fill-opacity:1""&
               & height=""",f5.1,""" width=""",f10.2""" y=""",f10.2,""" x=""",f10.2,""" />")')&
               & colors(events(i,j)%id), h, &
               & (events(i,j)%stop-events(i,j)%start)*scale, &
!                & h*(i-1), events(i,j)%start*scale
               & real(h)*real(maxth-i), events(i,j)%start*scale
       end do
    end do
    
    tottime = sum(ttimes)

    do i=1, nevtype
       write(4,'("<text x=""0"" y=""",f10.2,""" font-size=""",i3,""" fill=""",a7,""">",a20," -- ",f4.1,"%</text>")')&
            & real(h)*real(maxth+i+2),floor(h),colors(i),labels(i),(ttimes(i)/tottime)*100
    end do

    write(4,'("</svg>")')
    close(4)


    deallocate(events)

    return
    
  end subroutine qrm_trace_log_dump
  
end module qrm_trace_mod
