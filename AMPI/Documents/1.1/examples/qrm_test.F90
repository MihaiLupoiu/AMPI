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
!> @file test_qrm.F90
!!
!! An example program which reads parameters from a file, a matrix from a MM file and
!! solves the problem. Run this program as
!!
!! ./_qrm_test < input.txt
!!
!! (where _ is d, s, c or z) using the provided input.txt file.
!!
!! $Date: 2012-03-20 14:23:11 +0100 (mar., 20 mars 2012) $
!! $Author: abuttari $
!! $Version 0.0.1$
!! $Revision: 346 $
!!
!!##############################################################################################


program _qrm_test

  use _qrm_mod
  use _qrm_methods_mod
  use qrm_string_mod
  implicit none

  type(_qrm_spmat_type)       :: qrm_mat

  integer                      :: ierr, nargs, i, nrhs, ounit
  character                    :: matfile*30='', transp 
  _qrm_data, allocatable, target :: b(:), x(:), r(:)
  integer, pointer :: tmp(:)
  real(kind(1.d0)) :: t1, ta, tf, ts
  _qrm_real :: rnrm, onrm

  ! initialize the control data structure. 
  call qrm_spmat_init(qrm_mat)

  ! read the input file
  call qrm_read_parms(qrm_mat, matfile)
  call qrm_get('qrm_ounit', ounit)

  if(ounit.gt.0) write(ounit,'(30("="))')

  ! read the matrix. This subroutine has an overloaded interface
  ! and thus the type/precision of the data read will match that
  ! of qrm_mat
  call qrm_readmat(matfile, qrm_mat, .true.)

  if(qrm_mat%m .ge. qrm_mat%n) then
     transp='n'
  else
     if(ounit.gt.0) write(ounit,'("Transpose")')
     transp='t'
  end if

  if(ounit.gt.0) write(ounit,'("Starting Analysis")')
  t1 = qrm_swtime()
  call qrm_analyse(qrm_mat, transp)
  ta = qrm_swtime()-t1
  if(ounit.gt.0) write(ounit,'("  Estimated nonzeroes in R           : ",i20)')qrm_mat%gstats(qrm_e_nnz_r_)
  if(ounit.gt.0) write(ounit,'("  Estimated nonzeroes in H           : ",i20)')qrm_mat%gstats(qrm_e_nnz_h_)
  if(ounit.gt.0) write(ounit,'("  Estimated total flops at facto     : ",i20)')qrm_mat%gstats(qrm_e_facto_flops_)

  if(ounit.gt.0) write(ounit,'("Starting Factorization")')
  t1 = qrm_swtime()
  call qrm_factorize(qrm_mat, transp)
  tf = qrm_swtime()-t1

  if(qrm_mat%icntl(5) .gt. 0) then
     if(ounit.gt.0) write(ounit,'("Starting Solve")')
     call qrm_aalloc(b, qrm_mat%m)
     call qrm_aalloc(r, qrm_mat%m)
     call qrm_aalloc(x, qrm_mat%n)
     
     b = _qrm_one
     r = b

     t1 = qrm_swtime()
     if(transp .eq. 'n') then
        call qrm_apply(qrm_mat, 't', b)
        call qrm_solve(qrm_mat, 'n', b, x)
     else if(transp .eq. 't') then
        call qrm_solve(qrm_mat, 't', b, x)
        call qrm_apply(qrm_mat, 'n', x)
     end if
     ts = qrm_swtime()-t1
     
     ! compute the residual
     call qrm_residual_norm(qrm_mat, r, x, rnrm)
     call qrm_residual_orth(qrm_mat, r, onrm)   
     if(ounit.gt.0) write(ounit,'("||r||/||A||    = ",e10.2)')rnrm
     if(ounit.gt.0) write(ounit,'("||A^tr||/||r|| = ",e10.2)')onrm

     call qrm_adealloc(b)
     call qrm_adealloc(r)
     call qrm_adealloc(x)
  end if

10 continue

  call qrm_spmat_destroy(qrm_mat, all=.true.)
  if(ounit.gt.0) write(ounit,'("Done.")')
  if(ounit.gt.0) write(ounit,'(" ")')
  if(ounit.gt.0) write(ounit,'(" ")')
 
  if(ounit.gt.0) write(ounit,'("  Time to do the analysis  : ",es10.3)')ta
  if(ounit.gt.0) write(ounit,'("  Time to do the facto     : ",es10.3)')tf
  if(ounit.gt.0) write(ounit,'("  Time to compute solution : ",es10.3)')ts
  if(ounit.gt.0) write(ounit,'("  Nonzeroes in R           : ",i20)')qrm_mat%gstats(qrm_nnz_r_)
  if(ounit.gt.0) write(ounit,'("  Nonzeroes in H           : ",i20)')qrm_mat%gstats(qrm_nnz_h_)
  if(ounit.gt.0) write(ounit,'("  Total flops at facto     : ",i20)')qrm_mat%gstats(qrm_facto_flops_)
  if(ounit.gt.0) write(ounit,'("  Total unallocated memory : ",i20)')qrm_tot_mem
  if(ounit.gt.0) write(ounit,'("  Memory peak              : ",f9.3," MB")') &
       &real(qrm_max_mem,kind(1.d0))/1024.d0/1024.d0

  stop

contains

  subroutine qrm_read_parms(qrm_mat, matfile)
    ! Function: read_parms
    !
    ! *Input*:
    !
    ! *Output*:
    !
    use _qrm_spmat_mod
    implicit none

    type(_qrm_spmat_type)  :: qrm_mat
    character(len=*)  :: matfile

    character :: line*100, key*30, str*90
    integer :: ival
    real(kind(1.d0)) :: rval


    do
       read(*,'(a)')line
       read(line,*)key
       select case(key)
       case('end')
          exit
       case ('matfile')
          read(line,*)key, str
          str = adjustl(str)
          matfile = str(1:len_trim(str))
       case ('qrm_ordering','qrm_sing','qrm_keeph','qrm_nb',&
            & 'qrm_ib','qrm_nthreads','qrm_rhsnb','qrm_rhsnthreads')
          read(line,*)key,ival
          call qrm_set(qrm_mat, key, ival)
       case ('qrm_ounit','qrm_eunit','qrm_error_action')
          read(line,*)key,ival
          call qrm_set(key, ival)
       case default
          if(ounit.gt.0) write(ounit,'("Unknown parameter ",a20)')key
       end select

    end do

    return
    
  end subroutine qrm_read_parms
  

end program _qrm_test







