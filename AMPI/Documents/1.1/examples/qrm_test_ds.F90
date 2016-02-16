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
!!##############################################################################################
!> @file testboth.F90
!! FIXME: add comments
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!!##############################################################################################


program qrm_test_ds

  use sqrm_mod
  use dqrm_mod
  use qrm_string_mod

  implicit none

  type(sqrm_spmat_type)         :: sqrm_mat
  type(dqrm_spmat_type)         :: dqrm_mat
  integer                       :: ierr, nargs, i, nrhs
  character                     :: matfile*30, transp
  real(kind(1.e0)), allocatable :: sb(:), sx(:), sr(:)
  real(kind(1.d0)), allocatable :: db(:), dx(:), dr(:)
  real(kind(1.e0))  :: srnrm, sonrm
  real(kind(1.d0))  :: drnrm, donrm

  write(*,'(30("="))')
  call qrm_spmat_init(dqrm_mat)
  call qrm_spmat_init(sqrm_mat)

  call qrm_set(sqrm_mat,'qrm_ordering',qrm_natural_)
  call qrm_set(dqrm_mat,'qrm_ordering',qrm_natural_)

  ! allocate arrays for the input matrix
  call qrm_palloc(sqrm_mat%irn, 13)
  call qrm_palloc(sqrm_mat%jcn, 13)
  call qrm_palloc(sqrm_mat%val, 13)

  ! allocate arrays for the input matrix
  call qrm_palloc(dqrm_mat%irn, 13)
  call qrm_palloc(dqrm_mat%jcn, 13)
  call qrm_palloc(dqrm_mat%val, 13)

  ! initialize the input matrix
  sqrm_mat%jcn = (/1,1,1,2,2,3,3,3,3,4,4,5,5/)
  sqrm_mat%irn = (/2,3,6,1,6,2,4,5,7,2,3,2,4/)
  sqrm_mat%val = (/0.7,0.6,0.4,0.1,0.1,0.3,0.6,0.7,0.2,0.5,0.2,0.1,0.6/)
  sqrm_mat%m   = 7
  sqrm_mat%n   = 5
  sqrm_mat%nz  = 13

  ! initialize the input matrix
  dqrm_mat%jcn = (/1,1,1,2,2,3,3,3,3,4,4,5,5/)
  dqrm_mat%irn = (/2,3,6,1,6,2,4,5,7,2,3,2,4/)
  dqrm_mat%val = (/0.7,0.6,0.4,0.1,0.1,0.3,0.6,0.7,0.2,0.5,0.2,0.1,0.6/)
  dqrm_mat%m   = 7
  dqrm_mat%n   = 5
  dqrm_mat%nz  = 13

  call qrm_aalloc(sb, sqrm_mat%m)
  call qrm_aalloc(sr, sqrm_mat%m)
  call qrm_aalloc(sx, sqrm_mat%n)

  call qrm_aalloc(db, dqrm_mat%m)
  call qrm_aalloc(dr, dqrm_mat%m)
  call qrm_aalloc(dx, dqrm_mat%n)

  sb = 1.e0
  db = 1.d0

  sr = sb
  dr = db

  write(*,'("Data initialized in both arithmetics")')
  write(*,'("Solving...")')
  
  call qrm_least_squares(sqrm_mat, sb, sx)
  call qrm_least_squares(dqrm_mat, db, dx)
  

  write(*,'("Done")')
  write(*,'("Error analysis for sp")')
  call qrm_residual_norm(sqrm_mat, sr, sx, srnrm)
  call qrm_residual_orth(sqrm_mat, sr, sonrm)   
  write(*,'("||r||/||A||    = ",e10.2)')srnrm
  write(*,'("||A^tr||/||r|| = ",e10.2)')sonrm

  write(*,'("Error analysis for dp")')
  call qrm_residual_norm(dqrm_mat, dr, dx, drnrm)
  call qrm_residual_orth(dqrm_mat, dr, donrm)   
  write(*,'("||r||/||A||    = ",e10.2)')drnrm
  write(*,'("||A^tr||/||r|| = ",e10.2)')donrm
  

  call qrm_spmat_destroy(sqrm_mat, all=.true.)
  call qrm_spmat_destroy(dqrm_mat, all=.true.)


  stop

end program qrm_test_ds
