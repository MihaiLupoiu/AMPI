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
!> @file qrm_sdata_mod.F90
!! This file contains a module that defines the data type for storing
!! the results of the solve phase (not used at the moment)
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains the definition of all the data related to the
!! solution phase. 
!! 
module _qrm_sdata_mod


  ! Type: _qrm_sdata_type
  !
  ! This type defines a data structure containing all the data related to the solve phase.
  !
  ! Fields:
  ! rhs  - a pointer pointing to the location where the rhs is stored
  ! x    - a pointer pointing to the location where the solution must be stored
  !
  type _qrm_sdata_type
     _qrm_data, pointer :: rhs(:), x(:)

  end type _qrm_sdata_type


contains


  subroutine _qrm_sdata_destroy(sdata)
    ! Function: qrm_sdata_destroy
    ! This subroutine clean sup an sdata data structure.
    ! 
    ! *Input/Output*:
    ! sdata - the data structure to be cleaned
    !
  
    implicit none
    type(_qrm_sdata_type) :: sdata

    nullify(sdata%rhs)
    nullify(sdata%x)
  
    return
  
  end subroutine _qrm_sdata_destroy
  

end module _qrm_sdata_mod
