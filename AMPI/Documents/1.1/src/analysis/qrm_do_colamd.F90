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
!> @file qrm_do_colamd.F90
!! This file contains the routine that computes a COLAMD permutation of the input matrix
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This subroutine computes the fill reducing ordering using COLAMD

!> Please refer to:
!!
!! <em>"A column approximate minimum degree ordering algorithm"</em>,
!! T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng, ACM Transactions on
!! Mathematical Software, vol 30, no. 3, Sept. 2004, pp. 353-376.
!! 
!! for the details of the reordering method.
!!
!! @param[in] graph  the graph associated to the matrix to be ordered.
!!
!! @param[out] cperm an integer array containing the new column order
!!
subroutine _qrm_do_colamd(graph, cperm)

  use _qrm_spmat_mod
  use qrm_error_mod
  use _qrm_analysis_mod, savesym => _qrm_do_ordering, savesym2 => _qrm_do_colamd
  use qrm_mem_mod
  use iso_c_binding

  implicit none

  type(_qrm_spmat_type) :: graph
  integer, target        :: cperm(:)

  interface
     subroutine qrm_colamd(n_row, n_col, Alen, A, p, err) bind(c, name='qrm_colamd')
       use iso_c_binding
       integer(c_int), value :: n_row, n_col, Alen
       integer(c_int)        :: A(*), p(*), err
     end subroutine qrm_colamd
  end interface

  interface
     subroutine qrm_colamd_recommended(Alen, nnz, n_row, n_col) bind(c, name='qrm_colamd_recommended')
       use iso_c_binding
       integer(c_int), value :: nnz, n_row, n_col
       integer(c_int)        :: Alen
     end subroutine qrm_colamd_recommended
  end interface

  type(_qrm_spmat_type) :: gcopy
  integer                :: i, idx, cnt, tmp, alen, err
  ! error management
  integer                         :: err_act
  character(len=*), parameter     :: name='qrm_do_colamd'

  call qrm_err_act_save(err_act)


  ! at this point we have to make a copy of the graph.
  ! this is a huge waste of mem but we don't have a choice
  ! since ccolamd destroys the graph which, instead, we want to
  ! save for successive computations
  
  ! compute the memory required by ccolamd (a lot) and allocate
  call qrm_colamd_recommended(alen, graph%nz, graph%m, graph%n)
  call qrm_palloc(gcopy%irn, alen)
  __QRM_CHECK_RET(name,'qrm_palloc',9999)
  gcopy%jptr => cperm

  ! make the copy
  call _qrm_spmat_copy(graph, gcopy, values=.false.)
  __QRM_CHECK_RET(name,'qrm_spmat_copy',9999)

  ! ccolamd wants 0 based indices (argh!)
  gcopy%irn(1:gcopy%nz) = gcopy%irn(1:gcopy%nz)-1
  gcopy%jptr(1:gcopy%n+1) = gcopy%jptr(1:gcopy%n+1)-1

  ! call ccolamd
  call qrm_colamd(gcopy%m, gcopy%n, alen, gcopy%irn, gcopy%jptr, err)
  if(err .eq. 0) then
     call qrm_err_push(18,name)
     goto 9999
  end if

  nullify(gcopy%jptr)
  ! ccolamd return a 0-based permutation (re-argh!)
  cperm = cperm+1

  call _qrm_spmat_destroy(gcopy, all=.true.)
  __QRM_CHECK_RET(name,'qrm_spmat_destroy',9999)

  call qrm_err_act_restore(err_act)
  return

9999 continue ! error management
  call qrm_err_act_restore(err_act)
  if(err_act .eq. qrm_abort_) then
     call qrm_err_check()
  end if

  return

end subroutine _qrm_do_colamd
