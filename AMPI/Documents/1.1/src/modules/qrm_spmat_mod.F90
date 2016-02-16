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
!> @file qrm_spmat_mod.F90
!! This file contains the module that implements the main qr_mumps data structure
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1 $
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This module contains the definition of the basic sparse matrix type
!! and of the associated methods
module _qrm_spmat_mod
  use qrm_common_mod
  use qrm_error_mod
  use qrm_mem_mod
  use qrm_adata_mod
  use _qrm_fdata_mod
  use qrm_error_mod

  !> @brief Generif interface for the @link ::_qrm_spmat_alloc @endlink routine
  interface qrm_spmat_alloc
     module procedure _qrm_spmat_alloc
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_init @endlink routine
  interface qrm_spmat_init
     module procedure _qrm_spmat_init
  end interface

  !> @brief Generif interface for the @link ::_qrm_cntl_init @endlink routine
  interface qrm_cntl_init
     module procedure _qrm_cntl_init
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_convert @endlink routine
  interface qrm_spmat_convert
     module procedure _qrm_spmat_convert
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_copy @endlink routine
  interface qrm_spmat_copy
     module procedure _qrm_spmat_copy
  end interface

  !> @brief Generif interface for the @link ::_qrm_spmat_destroy @endlink routine
  interface qrm_spmat_destroy
     module procedure _qrm_spmat_destroy
  end interface


  !> @brief Generif interface for the
  !! @link ::_qrm_pseti @endlink,
  !! @link ::_qrm_psetr @endlink and
  ! ! !! @link ::_qrm_psetl @endlink routines
  !!
  interface qrm_set
     module procedure _qrm_pseti, _qrm_psetr!, _qrm_psetl
  end interface

  !> @brief Generif interface for the
  !! @link ::_qrm_pgeti @endlink,
  !! @link ::_qrm_pgetr @endlink and
  ! !! @link ::_qrm_pgetl @endlink routines
  !!
  interface qrm_get
     module procedure _qrm_pgeti, _qrm_pgetii, _qrm_pgetr!, _qrm_pgetl
  end interface



  !> @brief This type defines the data structure used to store a matrix. 

  !> A matrix can be represented either in COO, CSR or CSC
  !! format. Following the qr_mumps convention any array
  !! visible/usable from the users interface will be a pointer;
  !! otherwise allocatbles will be used because they normally provide
  !! better performance.
  !!
  type _qrm_spmat_type
     !> an integer array containing control parameters
     !! expressed as an integer value. Undocumented entries are
     !! assigned as nicntl, nicntl-1 etc.
     !! 
     !! Meaning:
     !! - icntl(qrm_ordering_=1)     : the ordering method to be used:
     !!             - qrm_auto_=0    : automatic choice
     !!             - qrm_natural_=1 : natural order
     !!             - qrm_given_=2   : given order
     !!             - qrm_colamd_=3  : COLAMD
     !!             - qrm_metis_=4   : METIS
     !!             - qrm_scotch_=5  : SCOTCH
     !! - icntl(qrm_sing_=2)         : singleton detection switch:
     !!             - qrm_no_=0  : no singleton detection
     !!             - qrm_yes_=1 : singleton detection
     !! - icntl(qrm_minamalg_=3)     : minimum node size for amalgamation. Node i can be amalgamated
     !!               to its father iff they both have less than icntl(3) variables
     !!               and the fill in introduced is less than rcntl(1).
     !! - icntl(qrm_nb_=4)           : the block size to be used during factorization
     !! - inctl(qrm_keeph_=5)        : whether H (the Householder vectors) has to be stored or not 
     !!             - qrm_no_=0  : the Householder vectors are discarded
     !!             - qrm_yes_=1 : the Householder vectors are kept
     !! - inctl(qrm_ib_=6)           : the internal blocking size to reduce flops
     !! - icntl(qrm_rhsnb_=7)        : the blocking parameter to handle multiple right-hand-sides
     !! - icntl(qrm_nthreads_=8)     : the number of threads to use in facto
     !! - icntl(qrm_rhsnthreads_=9)  : the number of threads to use for the outer parallel loop on
     !!               blocks of rhss
     !! - icntl(qrm_nlz_=nicntl)     : the minimum number of subtrees in L0 is this times the
     !!                                number of threads
     !! - icntl(qrm_cnode_=nicntl-1) : the number of cores per NUMA node (or cache, whatever)
     
     integer                           :: icntl(20)=0
     !> an double precision array containing control parameters
     !! expressed as a real value. Undocumented entries are
     !! assigned as nrcntl, nrcntl-1 etc.
     !!
     !! Meaning:
     !! - rcntl(qrm_amalgth_=1)      : fill-in threshold for amalgamation. Node i can be amalgamated
     !!                                to its father iff they both have less than icntl(3) variables
     !!                                and the fill in introduced is less than rcntl(1).
     !! - rcntl(qrm_rweight_=nrcntl) : subtrees with relative weight below this threashold will
     !!                                be put in L0
     real(kind(1.d0))                  :: rcntl(10)=0.d0
     !> an array containing local global stats. some of its content will only
     !! be relevant on the master node. Meaning:
     !! - gstats(1) = total number of flops
     !! - gstats(2) = total number of nonzeroes in R (estimated)
     !! - gstats(3) = total number of nonzeroes in H (estimated)
     !! - gstats(4) = maximum front size (MxN)
     !! - gstats(5) = number of column singletons
     integer(kind=8)                   :: gstats(10)=0
     !> Pointer to the beginning of each row in CSR format
     integer, pointer,    dimension(:) :: iptr => null()
     !> Pointer to the beginning of each column in CSC format 
     integer, pointer,    dimension(:) :: jptr => null()
     !> Row indices
     integer, pointer,    dimension(:) :: irn => null()
     !> Column indices
     integer, pointer,    dimension(:) :: jcn => null()
     !> Numerical values
     _qrm_data, pointer, dimension(:)  :: val => null()
     !> Number of rows
     integer                           :: m=0
     !> Number of columns
     integer                           :: n=0
     !> Number of nonzero elements
     integer                           :: nz=0
     !> A pointer to an array containing a column permutation provided
     !! by the user
     integer, pointer,    dimension(:) :: cperm_in => null()
     !> a @link qrm_adata_mod::qrm_adata_type @endlink data which is meant to
     !> contain all the data related to the analysis phase         
     type(qrm_adata_type)              :: adata
     !> a @link _qrm_fdata_mod::_qrm_fdata_type @endlink data which is meant to
     !> contain all the data related to the factorization phase         
     type(_qrm_fdata_type)             :: fdata
     !> Storage format; can be either 'COO', 'CSR' or 'CSC'
     character(len=3)                  :: fmt='coo'
  end type _qrm_spmat_type

  
contains

  !> @brief This subroutine allocates memory for a sparse matrix.
  !!
  !! @param[in,out] qrm_spmat A @link _qrm_spmat_mod::_qrm_spmat_type @endlink data
  !!           structure. The memory for storing the matrix is
  !!           allocated according to the storage format. Also
  !!           qrm_spmat%nz, qrm_spmat%m and qrm_spmat%n are set to
  !!           nz, m and n respectively.
  !!           These are the sizes of the arrays in output
  !!           * coo: irn(nz), jcn(nz), val(nz)
  !!           * csr: iptr(m+1), jcn(nz), val(nz)
  !!           * csc: irn(nz), jptr(n+1), val(nz)
  !!
  !! @param[in] nz  The number of nonzeroes contained in the matrix
  !!
  !! @param[in] m   The number of rows in the matrix
  !!
  !! @param[in] n   The number of columns in the matrix
  !!
  !! @param[in] fmt The matrix storage format. Can be either "coo" or "csr"
  !!           or "csc"
  !!
  subroutine _qrm_spmat_alloc(qrm_spmat, nz, m, n, fmt)
    use qrm_mem_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type), intent(inout) :: qrm_spmat
    integer, intent(in)                  :: nz, m, n
    character, intent(in)                :: fmt*(*)

    logical :: lsamen ! LAPACK subroutine to test strings
    ! error management
    integer                              :: err_act
    character(len=*), parameter          :: name='_qrm_spmat_alloc'

    call qrm_err_act_save(err_act)

#if defined(debug)
    __QRM_PRNT_DBG('("Allocating Matrix")')
#endif

    if(fmt .eq. 'coo') then
       call qrm_palloc(qrm_spmat%irn, nz)
       call qrm_palloc(qrm_spmat%jcn, nz)
       call qrm_palloc(qrm_spmat%val, nz)
       __QRM_CHECK_RET(name,'qrm_palloc',9999)
    else if(fmt .eq. 'csr') then
       call qrm_palloc(qrm_spmat%iptr, m+1)
       call qrm_palloc(qrm_spmat%jcn , nz)
       call qrm_palloc(qrm_spmat%val , nz)
       __QRM_CHECK_RET(name,'qrm_palloc',9999)
    else if(fmt .eq. 'csc') then
       call qrm_palloc(qrm_spmat%irn , nz)
       call qrm_palloc(qrm_spmat%jptr, n+1)
       call qrm_palloc(qrm_spmat%val , nz)
       __QRM_CHECK_RET(name,'qrm_palloc',9999)
    else
       call qrm_err_push(1, '_qrm_spmat_convert',aed=fmt)
       goto 9999
    end if

    qrm_spmat%nz = nz
    qrm_spmat%m  = m
    qrm_spmat%n  = n

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_spmat_alloc


  !> @brief This subroutine initializes a qrm_spmat_type instance setting
  !! default values into the control parameters
  !!
  !! @param[in,out] qrm_spmat The matrix to be initialized
  !! 
  subroutine _qrm_spmat_init(qrm_spmat)

    implicit none
    type(_qrm_spmat_type), intent(inout) :: qrm_spmat

    character(LEN=10) :: str
    integer :: ierr

    call _qrm_cntl_init(qrm_spmat)


    nullify(qrm_spmat%iptr, qrm_spmat%jptr, qrm_spmat%irn, qrm_spmat%jcn, &
         & qrm_spmat%val, qrm_spmat%cperm_in)

    return

  end subroutine _qrm_spmat_init


  !> @brief This subroutine initializes a qrm_spmat_type instance setting
  !! default values into the control parameters
  !!
  !! @param[in,out] qrm_spmat The matrix to be initialized
  !! 
  subroutine _qrm_cntl_init(qrm_spmat)

    use qrm_common_mod
    implicit none
    type(_qrm_spmat_type), intent(inout) :: qrm_spmat

    character(LEN=10) :: str
    integer :: ierr

    ! set default values for icntl and rcntl
    qrm_spmat%icntl(qrm_ordering_)      = qrm_auto_
    qrm_spmat%icntl(qrm_minamalg_)      = 4
    qrm_spmat%icntl(qrm_nb_)            = 120
    qrm_spmat%icntl(qrm_keeph_)         = qrm_yes_
    qrm_spmat%icntl(qrm_ib_)            = 120
    qrm_spmat%icntl(qrm_rhsnb_)         = -1
    qrm_spmat%icntl(qrm_rhsnthreads_)   = 1
    qrm_spmat%icntl(qrm_nlz_)           = 8
    qrm_spmat%icntl(qrm_cnode_)         = 1
    qrm_spmat%icntl(qrm_sing_)          = qrm_no_
    
    qrm_spmat%rcntl(qrm_amalgth_)       = 0.05
    qrm_spmat%rcntl(qrm_rweight_)       = 0.001
    qrm_spmat%fmt = 'coo'

    call get_environment_variable(name="QRM_NUM_THREADS",value=str, status=ierr)
    if(ierr .eq. 1) then
       qrm_spmat%icntl(qrm_nthreads_) = 1
    else
       read(str,*)qrm_spmat%icntl(qrm_nthreads_)
    end if

    return

  end subroutine _qrm_cntl_init





  !> This subroutine converts an input matrix into a different
  !! storage format. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in] fmt         the format of the output matrix
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_spmat_convert(in_mat, out_mat, fmt, values)
    implicit none

    type(_qrm_spmat_type), intent(in) :: in_mat
    type(_qrm_spmat_type)             :: out_mat
    character, intent(in)              :: fmt*(*)
    logical, optional                  :: values
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_spmat_convert'

    call qrm_err_act_save(err_act)

    select case(in_mat%fmt)
    case('csc')
       select case(fmt)
       case('csr')
          call _qrm_csc_to_csr(in_mat, out_mat, values)
       case default
          call qrm_err_push(1, name,aed=in_mat%fmt)
          goto 9999
       end select
    case('coo')
       select case(fmt)
       case('csc')
          call _qrm_coo_to_csc(in_mat, out_mat, values)
       case default
          call qrm_err_push(1, name,aed=in_mat%fmt)
          goto 9999
       end select
    case default
       call qrm_err_push(1, name,aed=in_mat%fmt)
       goto 9999
    end select

    out_mat%icntl    = in_mat%icntl
    out_mat%rcntl    = in_mat%rcntl

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_spmat_convert


  !> This subroutine converts a COO matrix into a CSC
  !! matrix. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_coo_to_csc(in_mat, out_mat, values)

    implicit none

    type(_qrm_spmat_type), intent(in) :: in_mat
    type(_qrm_spmat_type)             :: out_mat
    logical, optional                  :: values

    integer, allocatable :: work(:)
    logical :: ivalues, ob
    integer :: i, j, idx, k, m, n
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_coo_to_csc'

    call qrm_err_act_save(err_act)

    if(present(values)) then
       ivalues = values
    else
       ivalues = .true.
    end if

    call qrm_aalloc(work, in_mat%n+1)

    call qrm_prealloc(out_mat%jptr, in_mat%n+1)
    call qrm_prealloc(out_mat%irn , in_mat%nz)
    if(ivalues) call qrm_prealloc(out_mat%val , in_mat%nz)
    __QRM_CHECK_RET(name,'qrm_alloc',9999)

    work=0
    ob = .false.

    m = in_mat%m
    n = in_mat%n

    ! first loop to calculate # of nz per column
    do k=1, in_mat%nz
       j = in_mat%jcn(k)
       i = in_mat%irn(k)
       if((j.gt.0) .and. (j.le. n) .and. (i.gt.0) .and. (i.le. m) ) then
          work(j) = work(j)+1
       else
          ! out of bounds coefficients. ignore and print a warning at the end
          ob = .true.
       end if
    end do

    if(ob) then
       __QRM_PRNT_DBG('("** Out-of-bounds coefficients present **")')
    end if
    
    ! loop to convert the counts into ptrs
    out_mat%jptr(1) = 1
    do j=2, n+1
       out_mat%jptr(j) = out_mat%jptr(j-1)+work(j-1)
    end do


    ! last loop to put things into place
    work=0
    ! instead of putting an "if" inside the loop
    ! I put it here to gain some speed
    if(ivalues) then
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%jptr(j)+work(j)
          out_mat%irn(idx) = i
          out_mat%val(idx) = in_mat%val(k)
          work(j) = work(j)+1
       end do
    else
       do k=1, in_mat%nz
          j = in_mat%jcn(k)
          i = in_mat%irn(k)
          if((j.le.0) .or. (j.gt. n) .or. (i.le.0) .or. (i.gt. m) ) cycle
          idx = out_mat%jptr(j)+work(j)
          out_mat%irn(idx) = i
          work(j) = work(j)+1
       end do
    end if

    call qrm_adealloc(work)
    __QRM_CHECK_RET(name,'qrm_adelloc',9999)

    out_mat%m   = in_mat%m
    out_mat%n   = in_mat%n
    out_mat%nz  = in_mat%nz
    out_mat%fmt = 'csc'

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_coo_to_csc

  !> This subroutine converts a CSC matrix into a CSR
  !! matrix. Optionally the values may be ignored
  !! (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_csc_to_csr(in_mat, out_mat, values)

    type(_qrm_spmat_type), intent(in) :: in_mat
    type(_qrm_spmat_type)             :: out_mat
    logical, optional                  :: values

    integer, allocatable :: work(:)

    logical :: ivalues, ob
    integer :: i, j, idx, ii, m, n
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_csc_to_csr'

    call qrm_err_act_save(err_act)

    if(present(values)) then
       ivalues=values
    else
       ivalues = .true.
    end if

    ob = .false.

    m = in_mat%m
    n = in_mat%n

    call qrm_aalloc(work, m+1)

    call qrm_prealloc(out_mat%iptr, m+1)
    call qrm_prealloc(out_mat%jcn , in_mat%nz)
    if(ivalues) call qrm_prealloc(out_mat%val , in_mat%nz)
    __QRM_CHECK_RET(name,'qrm_alloc',9999)

    work=0
    ! first loop to calculate # of nz per row
    do j = 1, n
       do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
          i = in_mat%irn(ii)
          if((i.gt.0) .and. (i.le.m)) then
             work(i) = work(i)+1
          else
             ob = .true.
          end if
       end do
    end do

    if(ob) then
       __QRM_PRNT_DBG('("** Out-of-bounds coefficients present **")')
    end if

    ! loop to convert the counts into ptrs
    out_mat%iptr(1) = 1
    do j=2, m+1
       out_mat%iptr(j) = out_mat%iptr(j-1)+work(j-1)
    end do


    ! last loop to put things into place
    work=0
    ! instead of putting an "if" inside the loop
    ! I put it here to gain some speed
    if(ivalues) then
       do j = 1, n
          do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
             i = in_mat%irn(ii)
             if((i.le.0) .or. (i.gt.m)) cycle
             idx = out_mat%iptr(i)+work(i)
             out_mat%jcn(idx) = j
             out_mat%val(idx) = in_mat%val(ii)
             work(i) = work(i)+1
          end do
       end do
    else
       do j = 1, n
          do ii= in_mat%jptr(j), in_mat%jptr(j+1)-1
             i = in_mat%irn(ii)
             if((i.le.0) .or. (i.gt.m)) cycle
             idx = out_mat%iptr(i)+work(i)
             out_mat%jcn(idx) = j
             work(i) = work(i)+1
          end do
       end do
    end if

    call qrm_adealloc(work)
    __QRM_CHECK_RET(name,'qrm_adelloc',9999)

    out_mat%m   = in_mat%m
    out_mat%n   = in_mat%n
    out_mat%nz  = in_mat%nz
    out_mat%fmt = 'csr'

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_csc_to_csr






  !> This subroutine makes a copy of a matrix. Optionally the values
  !! may be ignored (this comes handy during the analysis)
  !!
  !! @param[in] in_mat      the input matrix
  !!
  !! @param[in,out] out_mat the output matrix in fmt format
  !!
  !! @param[in]  values      (optional) if values=.true. the output matrix will include
  !!                        numerical values, otherwise only the structure
  !!
  subroutine _qrm_spmat_copy(in_mat, out_mat, values)

    type(_qrm_spmat_type), intent(in) :: in_mat
    type(_qrm_spmat_type)             :: out_mat
    logical, optional                  :: values

    logical :: ivalues=.true.
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_spmat_copy'

    call qrm_err_act_save(err_act)

    ! TODO complete with other types

    if(present(values)) ivalues=values

    select case(in_mat%fmt)
    case('csc')
       call qrm_prealloc(out_mat%jptr, in_mat%n+1)
       call qrm_prealloc(out_mat%irn,  in_mat%nz)
       __QRM_CHECK_RET(name,'qrm_prelloc',9999)

       do i=1, in_mat%n+1
          out_mat%jptr(i) = in_mat%jptr(i)
       end do
       do i=1, in_mat%nz
          out_mat%irn(i)   = in_mat%irn(i)
       end do
       if(ivalues) then
          call qrm_prealloc(out_mat%val,  in_mat%nz)
          __QRM_CHECK_RET(name,'qrm_prealloc',9999)
          out_mat%val = in_mat%val
       end if
    case('coo')
       call qrm_prealloc(out_mat%jcn, in_mat%nz)
       call qrm_prealloc(out_mat%irn, in_mat%nz)
       __QRM_CHECK_RET(name,'qrm_prealloc',9999)
       do i=1, in_mat%nz
          out_mat%jcn(i) = in_mat%jcn(i)
          out_mat%irn(i) = in_mat%irn(i)
       end do
       if(ivalues) then
          call qrm_prealloc(out_mat%val,  in_mat%nz)
          __QRM_CHECK_RET(name,'qrm_realloc',9999)
          out_mat%val = in_mat%val
       end if
    case default
       call qrm_err_push(1, name,aed=in_mat%fmt)
       goto 9999
    end select

    out_mat%n        = in_mat%n
    out_mat%m        = in_mat%m
    out_mat%nz       = in_mat%nz
    out_mat%fmt      = in_mat%fmt
    out_mat%icntl    = in_mat%icntl
    out_mat%rcntl    = in_mat%rcntl

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_spmat_copy

  !> @brief This subroutine destroyes a qrm_spmat instance
  !!
  !! @param[in,out] qrm_spmat the matrix to be destroyed
  !!
  !! @param[in] all whether to deallocate all the memory or not
  !!
  subroutine _qrm_spmat_destroy(qrm_spmat, all)

    use qrm_mem_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    logical, optional      :: all

    logical                         :: iall
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_spmat_destroy'

    call qrm_err_act_save(err_act)

    if(present(all)) then
       iall = all
    else
       iall = .false.
    end if

    if(iall) then
       call qrm_pdealloc(qrm_spmat%irn)
       call qrm_pdealloc(qrm_spmat%jcn)
       call qrm_pdealloc(qrm_spmat%iptr)
       call qrm_pdealloc(qrm_spmat%jptr)
       call qrm_pdealloc(qrm_spmat%val)
       call qrm_pdealloc(qrm_spmat%cperm_in)
       __QRM_CHECK_RET(name,'qrm_pdealloc',9999)
    end if

    qrm_spmat%n       = 0
    qrm_spmat%m       = 0
    qrm_spmat%nz      = 0
    qrm_spmat%fmt     = ''

    call qrm_adata_destroy(qrm_spmat%adata)
    __QRM_CHECK_RET(name,name,9999)
    call _qrm_fdata_destroy(qrm_spmat%fdata)
    __QRM_CHECK_RET(name,name,9999)

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_spmat_destroy



  ! The following subroutine set or get control parameters from the
  ! cntl or rcntl control arrays. All the set and get routines are
  ! gathered under the same, overloaded interface, respectively


  !> @brief This subroutine is meant to set the integer control parameters
  !!
  !! @param[in,out] qrm_spmat The qrm_spmat instance concerned by the setting
  !!
  !! @param[in] string a string describing the parameter to be
  !!            set. Accepted values are:
  !!            - "qrm_ordering" : to set a method for the fill-reducing
  !!                               column permutation. Accepted values are:
  !!                               - qrm_auto_=0    for automatic choice
  !!                               - qrm_natural_=1 for natural ordering
  !!                               - qrm_given_=2   for given ordering (through
  !!                                                the qrm_spmat%cperm_in pointer)
  !!                               - qrm_colamd_=3  for COLAMD
  !!                               - qrm_metis_=4   for METIS
  !!                               - qrm_scotch_=5  for SCOTCH
  !!            - "qrm_minamalg" : fronts whose size is smaller than this will be
  !!                               systematically amalgamated to their father
  !!            - "qrm_nb" : the block-size that defines the granularity of parallel tasks
  !!            - "qrm_ib" : the block-size for computations
  !!            - "qrm_rhsnb" : the block-size for grouping RHSs
  !!            - "qrm_nthreads" : the number of threads to be used in the factorization
  !!                               (this can also be controlled through the
  !!                               QRM_NUM_THREADS environment variable)
  !!            - "qrm_rhsnthreads" : the number of threads for the solve phase
  !!            - "qrm_keeph" : whether to store or not the Householder vectors. Accepted
  !!                            values are qrm_yes_ and qrm_no_
  !!            - "qrm_sing" : whether or not to detect the presence of singletons. Accepted
  !!                           values are qrm_yes_ and qrm_no_
  !!            - "qrm_nlz"  : the number of subtrees in L0 will be at least this times
  !!                           the number of threads
  !!            - "qrm_cnode" : the number of cores per node
  !!
  !! @param[in] ival Any of the accepted values described above
  !!
  subroutine _qrm_pseti(qrm_spmat, string, ival)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    character(len=*)     :: string
    integer              :: ival

    character(len=len(string)) :: istring
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_pseti'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)
    if(index(istring,'qrm_ordering') .eq. 1) then
       qrm_spmat%icntl(qrm_ordering_) = ival
    else if (index(istring,'qrm_minamalg') .eq. 1) then
       qrm_spmat%icntl(qrm_minamalg_) = ival
    else if (index(istring,'qrm_nb') .eq. 1) then
       qrm_spmat%icntl(qrm_nb_) = ival
       if (qrm_spmat%icntl(qrm_ib_).gt.qrm_spmat%icntl(qrm_nb_)) then
          __QRM_PRNT_MSG('("Warning: qrm_ib is being set equal to qrm_nb")')
          qrm_spmat%icntl(qrm_ib_) = qrm_spmat%icntl(qrm_nb_)
       end if
    else if (index(istring,'qrm_ib') .eq. 1) then
       qrm_spmat%icntl(qrm_ib_) = ival
       if (qrm_spmat%icntl(qrm_ib_).gt.qrm_spmat%icntl(qrm_nb_)) then
          __QRM_PRNT_MSG('("Warning: qrm_nb is being set equal to qrm_ib")')
          qrm_spmat%icntl(qrm_nb_) = qrm_spmat%icntl(qrm_ib_)
       end if
    else if (index(istring,'qrm_rhsnb') .eq. 1) then
       qrm_spmat%icntl(qrm_rhsnb_) = ival
    else if (index(istring,'qrm_nthreads') .eq. 1) then
       qrm_spmat%icntl(qrm_nthreads_) = ival
    else if (index(istring,'qrm_rhsnthreads') .eq. 1) then
       qrm_spmat%icntl(qrm_rhsnthreads_) = ival
    else if (index(istring,'qrm_keeph') .eq. 1) then
       if(ival .eq. qrm_yes_) then
          qrm_spmat%icntl(qrm_keeph_) = ival
       else
          qrm_spmat%icntl(qrm_keeph_) = qrm_no_
       end if
    else if (index(istring,'qrm_sing') .eq. 1) then
       if(ival .eq. qrm_yes_) then
          qrm_spmat%icntl(qrm_sing_) = ival
       else
          qrm_spmat%icntl(qrm_sing_) = qrm_no_
       end if
    else if (index(istring,'qrm_nlz') .eq. 1) then
       qrm_spmat%icntl(qrm_nlz_) = ival
    else if (index(istring,'qrm_cnode') .eq. 1) then
       qrm_spmat%icntl(qrm_cnode_) = ival
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return

  end subroutine _qrm_pseti


  !> @brief This subroutine is meant to set the real control parameters
  !!
  !! @param[in,out] qrm_spmat The qrm_spmat instance concerned by the setting
  !!
  !! @param[in] string a string describing the parameter to be
  !!            set. Accepted values are:
  !!            - "qrm_amalgth" : the threshold that controls the amalgamation.
  !!                              A higher threshold means more fill-in but also
  !!                              more BLAS-3. Any real value is accepted
  !!
  !! @param[in] rval Any of the accepted values described above
  !!
  subroutine _qrm_psetr(qrm_spmat, string, rval)

    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    character(len=*)     :: string
    real(kind(1.d0))     :: rval

    character(len=len(string)) :: istring
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_psetr'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_amalgth') .eq. 1) then
       qrm_spmat%rcntl(qrm_amalgth_) = rval
    else if(index(istring,'qrm_rweight') .eq. 1) then
       qrm_spmat%rcntl(qrm_rweight_) = rval
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return

  end subroutine _qrm_psetr



  ! subroutine _qrm_psetl(qrm_spmat, string, lval)

    ! use qrm_string_mod
    ! use qrm_error_mod
    ! implicit none

    ! type(_qrm_spmat_type) :: qrm_spmat
    ! character(len=*)     :: string
    ! logical              :: lval

    ! character(len=len(string)) :: istring
    ! ! error management
    ! integer                         :: err_act
    ! character(len=*), parameter     :: name='_qrm_psetl'

    ! call qrm_err_act_save(err_act)

    ! istring = qrm_str_tolower(string)

    ! if(index(istring,'qrm_keeph') .eq. qrm_yes_) then
       ! if(lval) then
          ! qrm_spmat%icntl(qrm_keeph_) = qrm_yes_
       ! else
          ! qrm_spmat%icntl(qrm_keeph_) = qrm_no_
       ! end if
    ! else if(index(istring,'qrm_sing') .eq. 1) then
       ! if(lval) then
          ! qrm_spmat%icntl(qrm_sing_) = 1
       ! else
          ! qrm_spmat%icntl(qrm_sing_) = 0
       ! end if
    ! else
       ! call qrm_err_push(23, name, aed=string)
       ! goto 9999
    ! end if

    ! call qrm_err_act_restore(err_act)
    ! return

! 9999 continue ! error management
    ! call qrm_err_act_restore(err_act)
    ! if(err_act .eq. qrm_abort_) then
       ! call qrm_err_check()
    ! end if

    ! return

  ! end subroutine _qrm_psetl_


  !> @brief Gets the values of an integer control parameter.
  !!        This is the dual of the @link ::_qrm_pseti @endlink
  !!        routine; the parameters and accepted values are the same.
  subroutine _qrm_pgeti(qrm_spmat, string, ival)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    character(len=*)    :: string
    integer             :: ival

    character(len=len(string)) :: istring
    integer(kind=8) :: iival
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_pgeti'

    call qrm_err_act_save(err_act)

    call _qrm_pgetii(qrm_spmat, string, iival)
    __QRM_CHECK_RET(name,'qrm_pgetii',9999)

    ival = iival

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return

  end subroutine _qrm_pgeti

  !> @brief Gets the values of an integer control parameter.
  !!        This is the dual of the @link ::_qrm_pseti @endlink
  !!        routine; the parameters and accepted values are the same.
  subroutine _qrm_pgetii(qrm_spmat, string, ival)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    character(len=*  )    :: string
    integer(kind=8)       :: ival

    character(len=len(string)) :: istring
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_pgetii'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_ordering') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_ordering_)
    else if (index(istring,'qrm_minamalg') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_minamalg_)
    else if (index(istring,'qrm_nb') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_nb_)
    else if (index(istring,'qrm_ib') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_ib_)
    else if (index(istring,'qrm_rhsnb') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_rhsnb_)
    else if (index(istring,'qrm_nthreads') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_nthreads_)
    else if (index(istring,'qrm_rhsnthreads') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_rhsnthreads_)
    else if (index(istring,'qrm_keeph') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_keeph_)
    else if (index(istring,'qrm_sing') .eq. 1) then
       ival = qrm_spmat%icntl(qrm_sing_)
    else if (index(istring,'qrm_e_nnz_r') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_nnz_r_)
    else if (index(istring,'qrm_e_nnz_h') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_nnz_h_)
    else if (index(istring,'qrm_e_facto_flops') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_e_facto_flops_)
    else if (index(istring,'qrm_nnz_r') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_nnz_r_)
    else if (index(istring,'qrm_nnz_h') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_nnz_h_)
    else if (index(istring,'qrm_facto_flops') .eq. 1) then
       ival = qrm_spmat%gstats(qrm_facto_flops_)
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return

  end subroutine _qrm_pgetii



  !> @brief Gets the values of a real control parameter.
  !!        This is the dual of the @link ::_qrm_psetr @endlink
  !!        routine; the parameters and accepted values are the same.
  !!
  subroutine _qrm_pgetr(qrm_spmat, string, rval)

    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    character(len=*)     :: string
    real(kind(1.d0))     :: rval

    character(len=len(string)) :: istring
    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_pgetr'

    call qrm_err_act_save(err_act)

    istring = qrm_str_tolower(string)

    if(index(istring,'qrm_amalgth') .eq. 1) then
       rval = qrm_spmat%rcntl(qrm_amalgth_)
    else
       call qrm_err_push(23, name, aed=string)
       goto 9999
    end if

    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if

    return

  end subroutine _qrm_pgetr


  ! subroutine _qrm_pgetl(qrm_spmat, string, lval)

    ! use qrm_string_mod
    ! use qrm_error_mod
    ! implicit none

    ! type(_qrm_spmat_type) :: qrm_spmat
    ! character(len=*)     :: string
    ! logical              :: lval

    ! character(len=len(string)) :: istring
    ! ! error management
    ! integer                         :: err_act
    ! character(len=*), parameter     :: name='_qrm_pgetl'

    ! call qrm_err_act_save(err_act)

    ! istring = qrm_str_tolower(string)

    ! if(index(istring,'qrm_keeph') .eq. 1) then
       ! if (qrm_spmat%icntl(qrm_keeph_) .eq. qrm_yes_) then
          ! lval = .true.
       ! else
          ! lval = .false.
       ! end if
    ! else if(index(istring,'qrm_sing') .eq. 1) then
       ! if(qrm_spmat%icntl(qrm_sing_) .eq. 1) then
          ! lval = .true.
       ! else
          ! lval = .false.
       ! end if
    ! else
       ! call qrm_err_push(23, name, aed=string)
       ! goto 9999
    ! end if

    ! call qrm_err_act_restore(err_act)
    ! return

! 9999 continue ! error management
    ! call qrm_err_act_restore(err_act)
    ! if(err_act .eq. qrm_abort_) then
       ! call qrm_err_check()
    ! end if

    ! return

  ! end subroutine _qrm_pgetl


  !> @brief Check the compatibility and correctness of icntl and rcntl
  !! parameters
  subroutine _qrm_check_spmat(qrm_spmat, op)
    use qrm_common_mod
    use qrm_string_mod
    use qrm_error_mod
    implicit none

    type(_qrm_spmat_type) :: qrm_spmat
    integer, optional     :: op

    integer :: iop

    ! error management
    integer                         :: err_act
    character(len=*), parameter     :: name='_qrm_check_spmat'

    call qrm_err_act_save(err_act)

    if(present(op)) then
       iop = op
    else
       iop = qrm_allop_
    end if
    
    if((qrm_spmat%m .lt. 0) .or. (qrm_spmat%n .lt. 0) .or. &
         & (qrm_spmat%nz .lt. 0) .or. &
         & (qrm_spmat%nz .gt. (int(qrm_spmat%n,kind=8)*int(qrm_spmat%m,kind=8)))) then
       call qrm_err_push(29, name,ied=(/qrm_spmat%m,qrm_spmat%n,qrm_spmat%nz,0,0/))
       goto 9999
    end if


    if((iop.eq.qrm_allop_) .or. (iop.eq.qrm_analyse_)) then
       
       ! all the potential cases of conflict with the orderings
       select case(qrm_spmat%icntl(qrm_ordering_))
       case(:-1,6:)
          call qrm_err_push(9, name,ied=(/qrm_spmat%icntl(qrm_ordering_),0,0,0,0/))
          goto 9999
       case (qrm_given_)
          if(qrm_spmat%icntl(qrm_sing_) .eq. qrm_yes_) then
             call qrm_err_push(27, name,ied=(/qrm_ordering_,qrm_sing_,0,0,0/))
             goto 9999
          end if
       end select
       
       ! all the potential cases of conflict with the orderings
       select case(qrm_spmat%icntl(qrm_nb_))
       case(:-1)
          call qrm_err_push(28, name, ied=(/qrm_spmat%icntl(qrm_nb_),0,0,0,0/))
          goto 9999
       case default
          if(qrm_spmat%icntl(qrm_nb_) .lt. qrm_spmat%icntl(qrm_ib_)) then
             call qrm_err_push(27, name,ied=(/qrm_nb_,qrm_ib_,0,0,0/))
             goto 9999
          end if
       end select
       
       select case(qrm_spmat%icntl(qrm_ib_))
       case(:-1)
          call qrm_err_push(28, name, ied=(/qrm_spmat%icntl(qrm_ib_),0,0,0,0/))
          goto 9999
       end select
    end if
    
    call qrm_err_act_restore(err_act)
    return

9999 continue ! error management
    call qrm_err_act_restore(err_act)
    if(err_act .eq. qrm_abort_) then
       call qrm_err_check()
    end if
    return

  end subroutine _qrm_check_spmat
  








end module _qrm_spmat_mod
