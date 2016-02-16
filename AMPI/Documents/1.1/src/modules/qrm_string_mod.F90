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
!> @file qrm_string_mod.F90
!! This file contains a module that implements string handling tools
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains various string handling routines
module qrm_string_mod


  character(len=*), parameter   :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter   :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface qrm_str_tolower
     module procedure qrm_str_tolower
  end interface

  interface qrm_str_toupper
     module procedure qrm_str_toupper
  end interface

contains 

  function  qrm_str_tolower(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: qrm_str_tolower
    integer  :: i,k
    
    do i=1,len(string)
       k = index(ucase, string(i:i))
       if (k /= 0) then 
          qrm_str_tolower(i:i) = lcase(k:k)
       else          
          qrm_str_tolower(i:i) = string(i:i)
       end if
    enddo
   
    return
    
  end function qrm_str_tolower
 
  function  qrm_str_toupper(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: qrm_str_toupper
    integer  :: i,k
    
    do i=1,len(string)
       k = index(lcase, string(i:i))
       if (k /= 0) then 
          qrm_str_toupper(i:i) = ucase(k:k)
       else          
          qrm_str_toupper(i:i) = string(i:i)
       end if
    enddo
   
    return
    
  end function qrm_str_toupper
 
end module qrm_string_mod
