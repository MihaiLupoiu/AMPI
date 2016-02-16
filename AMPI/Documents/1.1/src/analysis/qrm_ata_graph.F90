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
!> @file qrm_ata_graph.F90
!! This file contains the routine that computes the graph of A^T * A
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This subroutine computes the fill reducing ordering using METIS
!! 
!! @param[in] g_csc  the graph (in CSC format) associated to A
!! 
!! @param[out] ata_graph  the graph of A'*A in CSR format
!! 
subroutine _qrm_ata_graph(g_csc, ata_graph)
  use _qrm_spmat_mod
  use qrm_error_mod
  use qrm_mem_mod

  implicit none

  type(_qrm_spmat_type), intent(in)  :: g_csc
  type(_qrm_spmat_type), intent(out) :: ata_graph

  type(_qrm_spmat_type) :: g_csr
  integer :: i, j, row1, col1, row2, col2
  integer, allocatable :: mark(:)
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_ata_graph'

  call qrm_err_act_save(err_act)

  call _qrm_spmat_convert(g_csc, g_csr, 'csr', .false.)

  call qrm_palloc(ata_graph%iptr,g_csc%n+2)
  ata_graph%iptr = 0
  ata_graph%iptr(1:2) = 1

  call qrm_aalloc(mark, g_csc%n)
  __QRM_CHECK_RET(name,'convert/alloc',9999)

  mark = 0

  do col1=1, g_csc%n
     ! loop over all the columns of A
     do i=g_csc%jptr(col1), g_csc%jptr(col1+1)-1
        ! for each nnz in the column, we go through the corresponding
        ! row and count all the i,j pairs
        row1 = g_csc%irn(i)

        do j=g_csr%iptr(row1), g_csr%iptr(row1+1)-1
           col2 = g_csr%jcn(j)
           ! now, the element (col1,col2) is present in A'*A. Check if
           ! it was already found, otherwise add it.
           ! skip the diagonal
           if(col1 .eq. col2) cycle
           if(mark(col2) .lt. col1) then
              mark(col2) = col1
              ata_graph%iptr(col1+2) = ata_graph%iptr(col1+2)+1
           end if
        end do
     end do
  end do

  do i=3, g_csc%n+2
     ata_graph%iptr(i) = ata_graph%iptr(i)+ata_graph%iptr(i-1)
  end do

  ata_graph%nz = ata_graph%iptr(g_csc%n+2)
  call qrm_palloc(ata_graph%jcn,ata_graph%nz)
  __QRM_CHECK_RET(name,'qrm_palloc',9999)

  ! second pass to fill up
  mark = 0
  do col1=1, g_csc%n
     ! loop over all the columns of A
     do i=g_csc%jptr(col1), g_csc%jptr(col1+1)-1
        ! for each nnz in the column, we go through the corresponding
        ! row and count all the i,j pairs
        row1 = g_csc%irn(i)

        do j=g_csr%iptr(row1), g_csr%iptr(row1+1)-1
           col2 = g_csr%jcn(j)
           ! now, the element (col1,col2) is present in A'*A. Check if
           ! it was already found, otherwise add it.
           ! skip the diagonal
           if(col1 .eq. col2) cycle
           if(mark(col2) .lt. col1) then
              mark(col2) = col1
              ata_graph%jcn(ata_graph%iptr(col1+1)) = col2
              ata_graph%iptr(col1+1) = ata_graph%iptr(col1+1)+1
           end if
        end do
     end do
  end do

  ata_graph%n = g_csc%n
  ata_graph%m = g_csc%n

  call _qrm_spmat_destroy(g_csr, all=.true.)
  call qrm_adealloc(mark)
  __QRM_CHECK_RET(name,'destroy/dealloc',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_ata_graph
