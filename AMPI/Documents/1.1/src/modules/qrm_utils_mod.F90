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
!> @file qrm_utils_mod.F90
!! This file contains a module with generic interfaces for a number of auxiliary tools
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains generic interfaces for a number of auxiliary tools
module _qrm_utils_mod
  use iso_c_binding

  !> @brief Generic interface for the @link ::_qrm_readmat @endlink
  !> routine
  interface qrm_readmat
     subroutine _qrm_readmat(matfile, qrm_mat, fakec)
       use _qrm_spmat_mod
       character, intent(in)                :: matfile*30
       type(_qrm_spmat_type), intent(inout) :: qrm_mat
       logical, optional                    :: fakec
     end subroutine _qrm_readmat
  end interface

  !> @brief Generic interface for the @link ::_qrm_matmul2d @endlink
  !> and @link ::_qrm_matmul1d @endlink routines
  interface _qrm_matmul
     module procedure _qrm_matmul2dw, _qrm_matmul1dw
  end interface _qrm_matmul

  !> @brief Generic interface for the @link ::_qrm_matmul2d @endlink
  !> and @link ::_qrm_matmul1d @endlink routines
  interface qrm_matmul
     subroutine _qrm_matmul2d(qrm_mat, transp, alpha, x, beta, y)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(out)        :: y(:,:)
       _qrm_data, intent(in)         :: x(:,:)
       _qrm_data, intent(in)         :: alpha, beta
       character(len=*)              :: transp
     end subroutine _qrm_matmul2d
     subroutine _qrm_matmul1d(qrm_mat, transp, alpha, x, beta, y)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)         :: qrm_mat
       _qrm_data, intent(out)        :: y(:)
       _qrm_data, intent(in)         :: x(:)
       _qrm_data, intent(in)         :: alpha, beta
       character(len=*)              :: transp
     end subroutine _qrm_matmul1d
  end interface


  !> @brief Generic interface for the @link ::_qrm_vecnrm2d @endlink
  !> and @link ::_qrm_vecnrm1d @endlink routines
  interface _qrm_vecnrm
     module procedure _qrm_vecnrm2dw, _qrm_vecnrm1dw
  end interface _qrm_vecnrm

  !> @brief Generic interface for the @link ::_qrm_vecnrm2d @endlink
  !> and @link ::_qrm_vecnrm1d @endlink routines
  interface qrm_vecnrm
     subroutine _qrm_vecnrm2d(vec, n, ntype, nrm)
       _qrm_data, intent(in)  :: vec(:,:)
       _qrm_real              :: nrm(:)
       integer, intent(in)    :: n
       character              :: ntype
     end subroutine _qrm_vecnrm2d
     subroutine _qrm_vecnrm1d(vec, n, ntype, nrm)
       _qrm_data, intent(in)  :: vec(:)
       _qrm_real              :: nrm
       integer, intent(in)    :: n
       character              :: ntype
     end subroutine _qrm_vecnrm1d
  end interface

  !> @brief Generic interface for the @link ::_qrm_remap_pnt @endlink
  !> routine
  interface qrm_remap_pnt
     subroutine _qrm_remap_pnt(arr1d, pnt2d, n)
       integer :: n
       _qrm_data, target  :: arr1d(1:n)
       _qrm_data, pointer :: pnt2d(:,:)
     end subroutine _qrm_remap_pnt
  end interface


  
  interface qrm_matnrm
     subroutine _qrm_matnrm(qrm_mat, ntype, nrm)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: qrm_mat
       _qrm_real                         :: nrm
       character                         :: ntype
     end subroutine _qrm_matnrm
  end interface qrm_matnrm



contains

  subroutine _qrm_matmul2dw(qrm_mat, transp, alpha, x, beta, y)
    use _qrm_spmat_mod
    type(_qrm_spmat_type), target :: qrm_mat
    _qrm_data, intent(out)        :: y(:,:)
    _qrm_data, intent(in)         :: x(:,:)
    _qrm_data, intent(in)         :: alpha, beta
    character(len=*)              :: transp
    call _qrm_matmul2d(qrm_mat, transp, alpha, x, beta, y)
    return
  end subroutine _qrm_matmul2dw

  subroutine _qrm_matmul1dw(qrm_mat, transp, alpha, x, beta, y)
    use _qrm_spmat_mod
    type(_qrm_spmat_type)         :: qrm_mat
    _qrm_data, intent(out)        :: y(:)
    _qrm_data, intent(in)         :: x(:)
    _qrm_data, intent(in)         :: alpha, beta
    character(len=*)              :: transp
    call _qrm_matmul1d(qrm_mat, transp, alpha, x, beta, y)
    return
  end subroutine _qrm_matmul1dw


  subroutine _qrm_vecnrm2dw(vec, n, ntype, nrm)
    _qrm_data, intent(in)  :: vec(:,:)
    _qrm_real              :: nrm(:)
    integer, intent(in)    :: n
    character              :: ntype
    call _qrm_vecnrm2d(vec, n, ntype, nrm)
    return
  end subroutine _qrm_vecnrm2dw

  subroutine _qrm_vecnrm1dw(vec, n, ntype, nrm)
    _qrm_data, intent(in)  :: vec(:)
    _qrm_real              :: nrm
    integer, intent(in)    :: n
    character              :: ntype
    call _qrm_vecnrm1d(vec, n, ntype, nrm)
    return
  end subroutine _qrm_vecnrm1dw

end module _qrm_utils_mod
