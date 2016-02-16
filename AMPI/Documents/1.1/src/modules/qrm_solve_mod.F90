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
!> @file qrm_solve_mod.F90
!! This file contains a module with generic interfaces for all the solve typed routines
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains all the interfaces for the typed routines
!! in the solve phase
module _qrm_solve_mod

  !> @brief Generic interface for the @link ::_qrm_apply_q @endlink routine
  interface qrm_apply_q
     subroutine _qrm_apply_q(qrm_mat, b)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)  :: b(:,:)
     end subroutine _qrm_apply_q
  end interface

  !> @brief Generic interface for the @link ::_qrm_apply_qt @endlink routine
  interface qrm_apply_qt
     subroutine _qrm_apply_qt(qrm_mat, b)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       _qrm_data :: b(:,:)
     end subroutine _qrm_apply_qt
  end interface


  !> @brief Generic interface for the
  !! @link ::_qrm_apply @endlink and
  !! @link ::_qrm_apply1d @endlink routines
  !!
  interface _qrm_apply
     module procedure _qrm_apply2dw, _qrm_apply1dw
  end interface _qrm_apply


  ! !> @brief Generic interface for the
  ! !! @link ::_qrm_apply @endlink and
  ! !! @link ::_qrm_apply1d @endlink routines
  ! !!
  interface qrm_apply
     subroutine _qrm_apply2d(qrm_mat, transp, b)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: qrm_mat
       _qrm_data, intent(inout)          :: b(:,:)
       character(len=*), intent(in)       :: transp
     end subroutine _qrm_apply2d
     subroutine _qrm_apply1d(qrm_mat, transp, b)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: qrm_mat
       _qrm_data, intent(inout)          :: b(:)
       character(len=*), intent(in)       :: transp
     end subroutine _qrm_apply1d
  end interface qrm_apply

  !> @brief Generic interface for the @link ::_qrm_solve_r @endlink routine
  interface qrm_solve_r
     subroutine _qrm_solve_r(qrm_mat, b, x)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)  :: b(:,:)
       _qrm_data, intent(out)    :: x(:,:)
     end subroutine _qrm_solve_r
  end interface

  !> @brief Generic interface for the @link ::_qrm_solve_rt @endlink routine
  interface qrm_solve_rt
     subroutine _qrm_solve_rt(qrm_mat, b, x)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)  :: b(:,:)
       _qrm_data, intent(out)    :: x(:,:)
     end subroutine _qrm_solve_rt
  end interface

  !> @brief Generic interface for the
  !! @link ::_qrm_solve @endlink and
  !! @link ::_qrm_solve1d @endlink routines
  !!
  interface _qrm_solve
     module procedure _qrm_solve1dw, _qrm_solve2dw
  end interface

  !> @brief Generic interface for the
  !! @link ::_qrm_solve @endlink and
  !! @link ::_qrm_solve1d @endlink routines
  !!
  interface qrm_solve
     subroutine _qrm_solve2d(qrm_mat, transp, b, x)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)  :: b(:,:)
       _qrm_data, intent(out)    :: x(:,:)
       character(len=*) :: transp
     end subroutine _qrm_solve2d
     subroutine _qrm_solve1d(qrm_mat, transp, b, x)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)  :: b(:)
       _qrm_data, intent(out)    :: x(:)
       character(len=*) :: transp
     end subroutine _qrm_solve1d
  end interface

  !> @brief Generic interface for the @link ::_qrm_solve_sing_front @endlink routine
  interface qrm_solve_sing_front
     subroutine _qrm_solve_sing_front(qrm_mat, b, x, trans)
       use _qrm_spmat_mod
       use _qrm_fdata_mod
       type(_qrm_spmat_type), target :: qrm_mat
       _qrm_data, intent(inout)      :: b(:,:)
       _qrm_data, intent(inout)      :: x(:,:)
       character                      :: trans
     end subroutine _qrm_solve_sing_front
  end interface

contains

  subroutine _qrm_apply2dw(qrm_mat, transp, b)
    use _qrm_spmat_mod
    type(_qrm_spmat_type), intent(in) :: qrm_mat
    _qrm_data, intent(inout)          :: b(:,:)
    character(len=*), intent(in)      :: transp
    call _qrm_apply2d(qrm_mat, transp, b)
    return
  end subroutine _qrm_apply2dw

  subroutine _qrm_apply1dw(qrm_mat, transp, b)
    use _qrm_spmat_mod
    type(_qrm_spmat_type), intent(in) :: qrm_mat
    _qrm_data, intent(inout)          :: b(:)
    character(len=*), intent(in)      :: transp
    call _qrm_apply1d(qrm_mat, transp, b)
    return
  end subroutine _qrm_apply1dw

  subroutine _qrm_solve2dw(qrm_mat, transp, b, x)
    use _qrm_spmat_mod
    type(_qrm_spmat_type), target :: qrm_mat
    _qrm_data, intent(inout)  :: b(:,:)
    _qrm_data, intent(out)    :: x(:,:)
    character(len=*) :: transp
    call  _qrm_solve2d(qrm_mat, transp, b, x)
    return
  end subroutine _qrm_solve2dw

  subroutine _qrm_solve1dw(qrm_mat, transp, b, x)
    use _qrm_spmat_mod
    type(_qrm_spmat_type), target :: qrm_mat
    _qrm_data, intent(inout)  :: b(:)
    _qrm_data, intent(out)    :: x(:)
    character(len=*) :: transp
    call _qrm_solve1d(qrm_mat, transp, b, x)
    return
  end subroutine _qrm_solve1dw


end module _qrm_solve_mod
