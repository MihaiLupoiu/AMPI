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
!> @file qrm_readmat.F90
!! Contains a routine which can read a matrix from a Matrix-Market file
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine reads a Matrix Market matrix from a file
!! and stores it on the host processor
!!
!! @param[in] matfile   a string containing the name of the matrix file
!!
!! @param[in] fakec     an optional logical argument which controls how a complex
!!                      matrix is built when reading a file with a real one.
!!                      - .true. : the imaginary part of the matrix coefficients is set
!!                                 to be equal to the real part. This is mostly used for
!!                                 testing purposes due to the lack or complex matrices
!!                      - .false. : the imaginary part is set to zero. This is the default
!!                                  when fakec is not present
!!
!! *Output*:
!! qrm_mat - a _qrm_spmat_type data where the matrix will be stored
!!
subroutine _qrm_readmat(matfile, qrm_mat, fakec)

  use _qrm_spmat_mod
  use qrm_error_mod

  implicit none

  character(len=*), intent(in)         :: matfile
  type(_qrm_spmat_type), intent(inout) :: qrm_mat
  logical, optional                    :: fakec

  character(len=20)   :: rep, field, symm, typ, fmt
  integer             :: m, n, nnz, ival, i, myid, info, nprocs
  logical             :: values, ifakec
  _qrm_real           :: rnds(2), re, im
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_read_mat'

  call qrm_err_act_save(err_act)

  __QRM_PRNT_MSG('("Reading Matrix: ",a20)')matfile

  rep   = ''
  field = ''
  symm  = ''
  typ   = ''
  fmt   = ''

  if(present(fakec)) then
     ifakec = fakec
  else
     ifakec = .false.
  end if

  open(4,file=matfile, status='OLD', action='READ', iostat=info)
  if(info .gt. 0) then
     call qrm_err_push(25, name, aed=matfile)
     goto 9999
  end if

  read(4,*)rep,typ,fmt,field,symm

  read(4,*)rep
  do
     if(rep(1:1) .ne. '%') exit
     read(4,*)rep
  end do

  backspace(4)

  read(4,*)m,n,nnz

  values = field .ne. 'pattern'

  qrm_mat%m  = m
  qrm_mat%n  = n
  qrm_mat%nz = nnz

  call qrm_palloc( qrm_mat%irn, qrm_mat%nz )
  call qrm_palloc( qrm_mat%jcn, qrm_mat%nz )
  call qrm_palloc( qrm_mat%val, qrm_mat%nz )
  __QRM_CHECK_RET(name,'qrm_palloc',9999)
  if(values) then
     do i=1, nnz
#if defined (cprec) || defined(zprec)     
        if(field .eq. 'complex') then
           read(4,*)qrm_mat%irn(i), qrm_mat%jcn(i), re, im
           qrm_mat%val(i) = cmplx(re,im,kind(_qrm_one))
        else if((field.eq.'real') .or. (field.eq.'integer')) then
           read(4,*)qrm_mat%irn(i), qrm_mat%jcn(i), re
           if(ifakec) then
              qrm_mat%val(i)  = cmplx(re,re,kind(_qrm_one))
           else
              qrm_mat%val(i) = cmplx(re,_qrm_rzero,kind(_qrm_one))
           end if
        end if
#elif defined (sprec) || defined(dprec)     
        if(field .eq. 'complex') then
           read(4,*)qrm_mat%irn(i), qrm_mat%jcn(i), qrm_mat%val(i), im
        else if((field.eq.'real') .or. (field.eq.'integer')) then
           read(4,*)qrm_mat%irn(i), qrm_mat%jcn(i), qrm_mat%val(i)
        end if
#endif
     end do
  else
     do i=1, nnz
        read(4,*)qrm_mat%irn(i),qrm_mat%jcn(i)
     end do
#if defined (cprec) || defined(zprec)     
     qrm_mat%val = _qrm_one
#elif defined (sprec) || defined(dprec)     
     qrm_mat%val = _qrm_one
#endif
  end if

  close(4)
  __QRM_PRNT_MSG('("Matrix read.")')

  qrm_mat%fmt= 'coo'

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if
  return

end subroutine _qrm_readmat


