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

!> @brief This module contains the definition of the qr_mumps C interface common
!! to all precisions/types.
module qrm_c_comm_interface
  use iso_c_binding
  use qrm_error_mod
  use qrm_common_mod

contains

  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_gseti @endlink
  !> routine (only for global, integer type)
  subroutine qrm_gseti_c(string, val) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_int), value  :: val
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_gseti(a, val)

    return
    
  end subroutine qrm_gseti_c

  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_ggeti @endlink
  !> routine (only for global, integer type)
  subroutine qrm_ggeti_c(string, val) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_int)         :: val
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_ggeti(a, val)

    return
    
  end subroutine qrm_ggeti_c

  !> @brief C equivalent of the @link _qrm_common_mod::_qrm_ggetii @endlink
  !> routine (only for global, integer type)
  subroutine qrm_ggetii_c(string, val) bind(c)
    character(kind=c_char) :: string(40)
    integer(c_long_long)   :: val
    
    character(len=40) :: a
    
    write(a,'(40a)')string
    
    call qrm_ggetii(a, val)

    return
    
  end subroutine qrm_ggetii_c

  !> @brief C equivalent of the @link _qrm_error_mod::_qrm_error_check @endlink
  !> routine
  subroutine qrm_err_check_c() bind(c)
    use qrm_error_mod

    call qrm_err_check()

    return
  end subroutine qrm_err_check_c




end module qrm_c_comm_interface

