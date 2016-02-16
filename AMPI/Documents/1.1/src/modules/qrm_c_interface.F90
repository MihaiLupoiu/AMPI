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
!> @file qrm_cintface.F90
!! This file contains the C interface for qr_mumps.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> @brief This module contains the definition of the qr_mumps C interface.
module _qrm_c_interface
  use iso_c_binding
  use _qrm_spmat_mod

  !> @brief This is the main qr_mumps data type which is meant to hold
  !! all the information related to a matrix. It is equivalent to the
  !! @link _qrm_spmat_mod::_qrm_spmat_type @endlink Fortran type.
  type, bind(c) :: _qrm_spmat_type_c
     !> This array contains the list of row indices of the nonzeroes
     !! in the matrix.
     type(c_ptr)      :: irn
     !> This array contains the list of column indices of the nonzeroes
     !! in the matrix.
     type(c_ptr)      :: jcn
     !> This array contains the list of values of the nonzeroes
     !! in the matrix.
     type(c_ptr)      :: val
     !> The number of rows in the matrix.
     integer(c_int)   :: m
     !> The number of columns in the matrix.
     integer(c_int)   :: n
     !> The number of nonzeroes in the matrix.
     integer(c_int)   :: nz
     !> A pointer to an array containing a column permutation provided
     !! by the user
     type(c_ptr)      :: cperm_in
     !> The integer control parameters
     integer(c_int)   :: icntl(20)
     !> The real control parameters
     real(c_double)   :: rcntl(10)
     !> The stats resulting from various operations
     integer(c_long)  :: gstats(10)
     !> The handle to the @link _qrm_spmat_mod::_qrm_spmat_type
     !! @endlink instance which will be used internally by qr_mumps
     integer(c_int)   :: h
     type(c_ptr)      :: mat_ptr
  end type _qrm_spmat_type_c
  
  type inst 
     type(_qrm_spmat_type), pointer :: fmat => null()
  end type inst
  
  integer, parameter :: max_inst=10
  type(inst), save :: spmat_array(max_inst)
  integer, save :: ninst=0

contains

  !> @brief C equivalent of the @link _qrm_spmat_mod::_qrm_spmat_init @endlink
  !> routine
  subroutine _qrm_spmat_init_c(qrm_spmat_c) bind(c)
    use _qrm_spmat_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    type(_qrm_spmat_type), pointer :: fmat
    integer :: i
    
    allocate(fmat)
    call _qrm_spmat_init(fmat)

    qrm_spmat_c%icntl  = fmat%icntl
    qrm_spmat_c%rcntl  = fmat%rcntl
    qrm_spmat_c%gstats = fmat%gstats

    qrm_spmat_c%mat_ptr = c_loc(fmat)
    
    nullify(fmat)

    return

  end subroutine _qrm_spmat_init_c


  !> @brief C equivalent of the @link _qrm_spmat_mod::_qrm_spmat_destroy @endlink
  !> routine
  subroutine _qrm_spmat_destroy_c(qrm_spmat_c) bind(c)
    use _qrm_spmat_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)
 
    call _qrm_spmat_destroy(fmat, all=.false.)

    qrm_spmat_c%icntl  = fmat%icntl
    qrm_spmat_c%rcntl  = fmat%rcntl
    qrm_spmat_c%gstats = fmat%gstats

    qrm_spmat_c%mat_ptr = c_null_ptr
    return

  end subroutine _qrm_spmat_destroy_c

  
  !> @brief C equivalent of the @link ::_qrm_analyse @endlink
  !> routine
  subroutine _qrm_analyse_c(qrm_spmat_c, transp) bind(c)
    use _qrm_analysis_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    character(kind=c_char), value :: transp

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%cperm_in, fmat%cperm_in,(/qrm_spmat_c%n/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_analyse(fmat, transp)

    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_analyse_c

  !> @brief C equivalent of the @link ::_qrm_factorize @endlink
  !> routine
  subroutine _qrm_factorize_c(qrm_spmat_c, transp) bind(c)
    use _qrm_factorization_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    character(kind=c_char), value :: transp

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_factorize(fmat, transp)

    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_factorize_c


  !> @brief C equivalent of the @link ::_qrm_solve @endlink
  !> routine
  subroutine _qrm_solve_c(qrm_spmat_c, transp, b, x, nrhs) bind(c)
    use _qrm_solve_mod
    implicit none
    type(_qrm_spmat_type_c)      :: qrm_spmat_c
    character(kind=c_char), value :: transp
    type(c_ptr), value            :: b, x
    integer(c_int), value         :: nrhs

    _qrm_data, pointer           :: ib(:,:), ix(:,:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(b, ib,(/qrm_spmat_c%m, nrhs/))
    call c_f_pointer(x, ix,(/qrm_spmat_c%n, nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_solve(fmat, transp, ib, ix)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_solve_c


  !> @brief C equivalent of the @link ::_qrm_apply @endlink
  !> routine
  subroutine _qrm_apply_c(qrm_spmat_c, transp, b, nrhs) bind(c)
    use _qrm_solve_mod
    implicit none
    type(_qrm_spmat_type_c)      :: qrm_spmat_c
    character(kind=c_char), value :: transp
    type(c_ptr), value            :: b
    integer(c_int), value         :: nrhs

    _qrm_data, pointer           :: ib(:,:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    if(transp.eq.'t') then
       call c_f_pointer(b, ib,(/qrm_spmat_c%m, nrhs/))
    else if(transp.eq.'n') then
       call c_f_pointer(b, ib,(/qrm_spmat_c%n, nrhs/))
    end if

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_apply(fmat, transp, ib)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_apply_c



  !> @brief C equivalent of the @link ::_qrm_matmul @endlink
  !> routine
  subroutine _qrm_matmul_c(qrm_spmat_c, transp, alpha, x, beta, y, nrhs) bind(c)
    use _qrm_utils_mod
    implicit none
    type(_qrm_spmat_type_c)      :: qrm_spmat_c
    character(kind=c_char), value :: transp
    _qrm_data_fc, value          :: alpha, beta
    type(c_ptr), value            :: x, y
    integer(c_int), value         :: nrhs

    _qrm_data, pointer           :: ix(:,:), iy(:,:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(x, ix,(/qrm_spmat_c%n, nrhs/))
    call c_f_pointer(y, iy,(/qrm_spmat_c%m, nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_matmul(fmat, transp, alpha, ix, beta, iy)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_matmul_c
  

  !> @brief C equivalent of the @link ::_qrm_matnrm @endlink
  !> routine
  subroutine _qrm_matnrm_c(qrm_spmat_c, ntype, nrm) bind(c)
    use _qrm_utils_mod
    implicit none
    type(_qrm_spmat_type_c)       :: qrm_spmat_c
    character(kind=c_char), value :: ntype
    _qrm_real_fc                  :: nrm

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call _qrm_matnrm(fmat, ntype, nrm)

    return
  end subroutine _qrm_matnrm_c

  !> @brief C equivalent of the @link ::_qrm_vecnrm @endlink
  !> routine
  subroutine _qrm_vecnrm_c(x, n, nrhs, ntype, nrm) bind(c)
    use _qrm_utils_mod
    implicit none
    type(c_ptr), value            :: x
    integer(c_int), value         :: nrhs
    integer(c_int), value         :: n
    character(kind=c_char), value :: ntype
    type(c_ptr), value            :: nrm
    _qrm_data, pointer            :: ix(:,:)
    _qrm_real, pointer            :: inrm(:)

    call c_f_pointer(x, ix,(/n,nrhs/))
    call c_f_pointer(nrm, inrm,(/nrhs/))

    call _qrm_vecnrm(ix, n, ntype, inrm)

    return
  end subroutine _qrm_vecnrm_c


  !> @brief C equivalent of the @link ::_qrm_least_squares @endlink
  !> routine
  subroutine _qrm_least_squares_c(qrm_spmat_c, b, x, nrhs) bind(c)
    use _qrm_spmat_mod
    use _qrm_methods_mod
    implicit none
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr), value             :: b, x
    integer(c_int), value          :: nrhs
 
    _qrm_data, pointer             :: ib(:,:), ix(:,:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(b, ib,(/qrm_spmat_c%m, nrhs/))
    call c_f_pointer(x, ix,(/qrm_spmat_c%n, nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call qrm_least_squares(fmat, ib, ix)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_least_squares_c

  !> @brief C equivalent of the @link ::_qrm_min_norm @endlink
  !> routine
  subroutine _qrm_min_norm_c(qrm_spmat_c, b, x, nrhs) bind(c)
    use _qrm_spmat_mod
    use _qrm_methods_mod
    implicit none
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr), value             :: b, x
    integer(c_int), value          :: nrhs

    _qrm_data, pointer             :: ib(:,:), ix(:,:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(b, ib,(/qrm_spmat_c%m, nrhs/))
    call c_f_pointer(x, ix,(/qrm_spmat_c%n, nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call qrm_min_norm(fmat, ib, ix)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_min_norm_c


  !> @brief C equivalent of the @link ::_qrm_residual_norm @endlink
  !> routine
  subroutine _qrm_residual_norm_c(qrm_spmat_c, b, x, nrhs, nrm) bind(c)
    use _qrm_spmat_mod
    use _qrm_methods_mod
    implicit none
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    integer(c_int), value          :: nrhs
    type(c_ptr), value             :: b, x
    type(c_ptr), value             :: nrm

    _qrm_data, pointer             :: ib(:,:), ix(:,:)
    _qrm_real, pointer             :: inrm(:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(b, ib,(/qrm_spmat_c%m,nrhs/))
    call c_f_pointer(x, ix,(/qrm_spmat_c%n,nrhs/))
    call c_f_pointer(nrm, inrm,(/nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call qrm_residual_norm(fmat, ib, ix, inrm)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_residual_norm_c

  !> @brief C equivalent of the @link ::_qrm_residual_orth @endlink
  !> routine
  subroutine _qrm_residual_orth_c(qrm_spmat_c, r, nrhs, nrm) bind(c)
    use _qrm_spmat_mod
    use _qrm_methods_mod
    implicit none
    type(_qrm_spmat_type_c)        :: qrm_spmat_c
    type(c_ptr), value             :: r
    integer(c_int), value          :: nrhs
    type(c_ptr), value             :: nrm

    _qrm_data, pointer             :: ir(:,:)
    _qrm_real, pointer             :: inrm(:)

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)

    fmat%m  = qrm_spmat_c%m   
    fmat%n  = qrm_spmat_c%n   
    fmat%nz = qrm_spmat_c%nz  

    call c_f_pointer(qrm_spmat_c%irn, fmat%irn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%jcn, fmat%jcn,(/qrm_spmat_c%nz/))
    call c_f_pointer(qrm_spmat_c%val, fmat%val,(/qrm_spmat_c%nz/))

    call c_f_pointer(r, ir,(/qrm_spmat_c%m,nrhs/))
    call c_f_pointer(nrm, inrm,(/nrhs/))

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 

    call qrm_residual_orth(fmat, ir, inrm)
    
    qrm_spmat_c%gstats = fmat%gstats

    return
  end subroutine _qrm_residual_orth_c




  !> @brief C equivalent of the @link _qrm_spmat_mod::_qrm_pseti @endlink
  !> routine (only for global, integer type)
  subroutine _qrm_pseti_c(qrm_spmat_c, string, val) bind(c)
    use _qrm_spmat_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    character(kind=c_char)  :: string(40)
    integer(c_int), value   :: val
    
    character(len=40) :: a

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)
    
    write(a,'(40a)')string
    
    call _qrm_pseti(fmat, a, val)

    qrm_spmat_c%icntl = fmat%icntl  
    qrm_spmat_c%rcntl = fmat%rcntl 

    return
    
  end subroutine _qrm_pseti_c


  !> @brief C equivalent of the @link _qrm_spmat_mod::_qrm_pgeti @endlink
  !> routine (only for global, integer type)
  subroutine _qrm_pgeti_c(qrm_spmat_c, string, val) bind(c)
    use _qrm_spmat_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    character(kind=c_char)  :: string(40)
    integer(c_int)          :: val
    
    character(len=40) :: a

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)
    
    write(a,'(40a)')string

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 
    
    call _qrm_pgeti(fmat, a, val)

    return
    
  end subroutine _qrm_pgeti_c

  !> @brief C equivalent of the @link _qrm_spmat_mod::_qrm_pgetii @endlink
  !> routine (only for global, integer type)
  subroutine _qrm_pgetii_c(qrm_spmat_c, string, val) bind(c)
    use _qrm_spmat_mod
    implicit none
    type(_qrm_spmat_type_c) :: qrm_spmat_c
    character(kind=c_char)  :: string(40)
    integer(c_long_long)    :: val
    
    character(len=40) :: a

    type(_qrm_spmat_type), pointer :: fmat
    call c_f_pointer(qrm_spmat_c%mat_ptr, fmat)
    
    write(a,'(40a)')string

    fmat%icntl = qrm_spmat_c%icntl 
    fmat%rcntl = qrm_spmat_c%rcntl 
    
    call _qrm_pgetii(fmat, a, val)

    return
    
  end subroutine _qrm_pgetii_c


end module _qrm_c_interface

