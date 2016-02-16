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
!> @file qrm_factorization_mod.F90
!! This file contains the @link _qrm_factorization_mod @endlink with  the generic interfaces
!! for all the factorization related routines.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains all the generic interfaces for the typed routines
!! in the factorization phase
module _qrm_factorization_mod

  !> @brief Generic interface for the @link ::_qrm_factorization_init
  !> @endlink routine
  interface qrm_factorization_init
     subroutine _qrm_factorization_init(qrm_mat)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
     end subroutine _qrm_factorization_init
  end interface

  !> @brief Generic interface for the @link ::_qrm_activate_front
  !> @endlink routine
  interface qrm_activate_front
     subroutine _qrm_activate_front(qrm_mat, fnum, flops)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       integer                :: fnum
       real(kind(1.d0))       :: flops
     end subroutine _qrm_activate_front
  end interface

  !> @brief Generic interface for the @link ::_qrm_factorize
  !> @endlink routine
  interface qrm_factorize
     subroutine _qrm_factorize(qrm_mat, transp)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       character, optional, intent(in) :: transp
     end subroutine _qrm_factorize
  end interface
  
  !> @brief Generic interface for the @link ::_qrm_factorization_core
  !> @endlink routine
  interface qrm_factorization_core
     subroutine _qrm_factorization_core(qrm_mat)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
     end subroutine _qrm_factorization_core
  end interface
  
  !> @brief Generic interface for the @link ::_qrm_init_front
  !> @endlink routine
  interface qrm_init_front
     subroutine _qrm_init_front(qrm_mat, fnum, par, work)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       integer                        :: fnum
       logical                        :: par
       integer, optional              :: work(:)
     end subroutine _qrm_init_front
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_subtree
  !> @endlink routine
  interface qrm_do_subtree
     subroutine _qrm_do_subtree(qrm_mat, fnum, flops)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       integer                        :: fnum
       real(kind(1.d0))               :: flops
     end subroutine _qrm_do_subtree
  end interface

  !> @brief Generic interface for the @link ::_qrm_clean_front
  !> @endlink routine
  interface qrm_clean_front
     subroutine _qrm_clean_front(qrm_mat, fnum)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), target :: qrm_mat
       integer                        :: fnum
     end subroutine _qrm_clean_front
  end interface

  !> @brief Generic interface for the @link ::_qrm_store_h
  !> @endlink routine
  interface qrm_store_h
     subroutine _qrm_store_h(front)
       use _qrm_fdata_mod
       type(_qrm_front_type) :: front
     end subroutine _qrm_store_h
  end interface

  !> @brief Generic interface for the @link ::_qrm_store_r
  !> @endlink routine
  interface qrm_store_r
     subroutine _qrm_store_r(front)
       use _qrm_fdata_mod
       type(_qrm_front_type) :: front
     end subroutine _qrm_store_r
  end interface

end module _qrm_factorization_mod
