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
!> @file qrm_coverage.F90
!! This file contains coverage tests for the error handling
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################

subroutine _qrm_test_err(c)
  use _qrm_coverage_mod, protect => _qrm_test_err
  implicit none

  integer :: c

  integer, parameter :: ncases=15
  logical :: cases(ncases)

  if(c .eq. -1) then
     cases = .true.
  else if (c .le. ncases) then
     cases = .false.
     cases(c) = .true.
  end if

  if(cases(1 )) call case1 ()
  if(cases(2 )) call case2 ()
  if(cases(3 )) call case3 ()
  if(cases(4 )) call case4 ()
  if(cases(5 )) call case5 ()
  if(cases(6 )) call case6 ()
  if(cases(7 )) call case7 ()
  if(cases(8 )) call case8 ()
  if(cases(9 )) call case9 ()  
  if(cases(10)) call case10()
  if(cases(11)) call case11()
  if(cases(12)) call case12()
  if(cases(13)) call case13()
  if(cases(14)) call case14()


  return


contains

  !> @brief Check for error 1: unsupported matrix format
  subroutine case1()
    use _qrm_spmat_mod
    implicit none
    
    type(_qrm_spmat_type) :: qrm_spmat
    integer               :: info

    call _qrm_spmat_alloc(qrm_spmat, 1, 1, 1, 'xyz')
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 1, 1, -1, info .eq. 1)
    call qrm_flush_err_stack(.false.)
    call _qrm_spmat_destroy(qrm_spmat)
    return

  end subroutine case1
  

  !> @brief Check for error 4: allocating an already allocated
  subroutine case2()
    use qrm_mem_mod
    implicit none
    
    _qrm_data, allocatable :: a(:)
    integer :: info

    info = 0
    call qrm_aalloc(a, 10)
    call qrm_aalloc(a, 10)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 2, 1, -1, info .eq. 4)
    call qrm_adealloc(a)
    call qrm_flush_err_stack(.false.)

    return

  end subroutine case2
  


  !> @brief Check for error 8: invalid or unprovided colperm
  subroutine case3()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info, i

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_set(qrm_spmat,'qrm_ordering',qrm_given_)
    ! subcase one, cperm_in not provided
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 3, 1, 1, info .eq. 8)
    call qrm_flush_err_stack(.false.)



    ! subcase two, cperm_in is not valid
    info = 0
    call qrm_palloc(qrm_spmat%cperm_in, qrm_spmat%n)
    qrm_spmat%cperm_in = (/(i,i=0,qrm_spmat%n-1)/)
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 3, 2, 1, info .eq. 8)
    call qrm_flush_err_stack(.false.)

    call qrm_pdealloc(qrm_spmat%cperm_in)
    nullify(qrm_spmat)

    return

  end subroutine case3
  

  !> @brief Check for error 9: unknown ordering method
  subroutine case4()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_set(qrm_spmat,'qrm_ordering',10)

    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 4, 1, 1, info .eq. 9)
    call qrm_flush_err_stack(.false.)

    nullify(qrm_spmat)
    return
  end subroutine case4



  !> @brief Check for error 13: facto before analysis
  subroutine case5()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)

    info = 0
    call qrm_factorize(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 5, 1, 1, info .eq. 13)
    call qrm_flush_err_stack(.false.)

    nullify(qrm_spmat)
    return
  end subroutine case5


  !> @brief Check for error 14: solve before facto
  subroutine case6()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, allocatable :: b(:), x(:)
    integer   :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_aalloc(b, qrm_spmat%m)
    call qrm_aalloc(x, qrm_spmat%n)

    info = 0
    call qrm_apply(qrm_spmat, 'n', b)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 6, 1, 1, info .eq. 14)
    call qrm_flush_err_stack(.false.)

    info = 0
    call qrm_solve(qrm_spmat, 'n', b, x)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 6, 2, 1, info .eq. 14)
    call qrm_flush_err_stack(.false.)

    nullify(qrm_spmat)
    call qrm_adealloc(b)
    call qrm_adealloc(x)
    return
  end subroutine case6


  !> @brief Check for error 15: norm not implemented
  subroutine case7()
    use _qrm_mod
    implicit none
    
    _qrm_data, allocatable :: b(:), x(:)
    integer   :: info
    _qrm_real :: nrm

    call qrm_aalloc(x, 10)
    info = 0
    call _qrm_vecnrm(x, 10, 'a', nrm)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 7, 1, -1, info .eq. 15)
    call qrm_flush_err_stack(.false.)
    call qrm_adealloc(x)

    return
  end subroutine case7



  !> @brief Check for error 16: unknown ordering method
  subroutine case8()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_colamd_)
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
#if defined(have_colamd)
    call _qrm_prnt_testmesg(1, 8, 1, 1, info .eq. 0)
#else
    call _qrm_prnt_testmesg(1, 8, 1, 1, info .eq. 16)
#endif
    call qrm_flush_err_stack(.false.)

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_metis_)
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
#if defined(have_metis)
    call _qrm_prnt_testmesg(1, 8, 2, 1, info .eq. 0)
#else
    call _qrm_prnt_testmesg(1, 8, 2, 1, info .eq. 16)
#endif
    call qrm_flush_err_stack(.false.)

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_scotch_)
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
#if defined(have_metis)
    call _qrm_prnt_testmesg(1, 8, 3, 1, info .eq. 0)
#else
    call _qrm_prnt_testmesg(1, 8, 3, 1, info .eq. 16)
#endif
    call qrm_flush_err_stack(.false.)

    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    return
  end subroutine case8


  !> @brief Check for error 23: wrong argument to set/get
  subroutine case9()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)
    info = 0
    call qrm_set(qrm_spmat,'qrm_pippo',qrm_colamd_)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 9, 1, 1, info .eq. 23)
    call qrm_flush_err_stack(.false.)

    nullify(qrm_spmat)
    return
  end subroutine case9


  !> @brief Check for error 26: unknown error action
  subroutine case10()
    use _qrm_mod
    implicit none
    
    integer :: info

    info = 0
    call qrm_err_act_set(10)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 10, 1, -1, info .eq. 26)
    call qrm_flush_err_stack(.false.)

    return
  end subroutine case10



  !> @brief Check for error 27: Incompatible icntl
  subroutine case11()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)

    call qrm_set(qrm_spmat,'qrm_ordering',qrm_given_)
    call qrm_set(qrm_spmat,'qrm_sing',qrm_yes_)
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 10, 1, 1, info .eq. 27)
    call qrm_flush_err_stack(.false.)
    call qrm_adata_destroy(qrm_spmat%adata)


    call qrm_set(qrm_spmat,'qrm_ordering',qrm_natural_)
    call qrm_set(qrm_spmat,'qrm_sing',qrm_no_)
    qrm_spmat%icntl(qrm_nb_) = 10
    qrm_spmat%icntl(qrm_ib_) = 15
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_factorize(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 11, 2, 1, info .eq. 27)
    call qrm_flush_err_stack(.false.)
    
    call _qrm_fdata_destroy(qrm_spmat%fdata)
    call qrm_adata_destroy(qrm_spmat%adata)
    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    return
  end subroutine case11



  !> @brief Check for error 28: incorrect ib/nb
  subroutine case12()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info

    qrm_spmat => _qrm_get_test_mat(1)

    call qrm_set(qrm_spmat,'qrm_nb',-10)
    call qrm_set(qrm_spmat,'qrm_ib',-15)
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 12, 1, 1, info .eq. 28)
    call qrm_flush_err_stack(.false.)
    
    call _qrm_fdata_destroy(qrm_spmat%fdata)
    call qrm_adata_destroy(qrm_spmat%adata)
    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    return
  end subroutine case12



  !> @brief Check for error 29: wrong m/n/nz
  subroutine case13()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat

    integer :: info, m, n, nz

    qrm_spmat => _qrm_get_test_mat(1)

    m  = qrm_spmat%m
    n  = qrm_spmat%n
    nz = qrm_spmat%nz

    qrm_spmat%m  = -1

    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 13, 1, 1, info .eq. 29)
    call qrm_flush_err_stack(.false.)
    call qrm_adata_destroy(qrm_spmat%adata)

    qrm_spmat%m  = m
    qrm_spmat%n  = n
    qrm_spmat%nz = m*n+1
    info = 0
    call qrm_analyse(qrm_spmat)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 13, 2, 1, info .eq. 29)
    call qrm_flush_err_stack(.false.)
    call qrm_adata_destroy(qrm_spmat%adata)

    qrm_spmat%m  = m
    qrm_spmat%n  = n
    qrm_spmat%nz = nz

    nullify(qrm_spmat)
    return
  end subroutine case13



  !> @brief Check for error 30: apply when H is discarded
  subroutine case14()
    use _qrm_mod
    implicit none
    
    type(_qrm_spmat_type), pointer :: qrm_spmat
    _qrm_data, allocatable :: b(:)
    integer   :: info

    qrm_spmat => _qrm_get_test_mat(1)
    call qrm_aalloc(b, qrm_spmat%m)

    info = 0
    call qrm_set(qrm_spmat, 'qrm_keeph', qrm_no_)
    call qrm_analyse(qrm_spmat)
    call qrm_factorize(qrm_spmat)
    call qrm_apply(qrm_spmat, 't', b)
    call qrm_err_get(info)
    call _qrm_prnt_testmesg(1, 14, 1, 1, info .eq. 30)
    call qrm_flush_err_stack(.false.)

    call qrm_cntl_init(qrm_spmat)
    nullify(qrm_spmat)
    call qrm_adealloc(b)
    return
  end subroutine case14


end subroutine _qrm_test_err
