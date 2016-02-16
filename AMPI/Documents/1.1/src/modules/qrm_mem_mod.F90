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
!> @file qrm_mem_mod.F90
!! This file contains the module that implements all the memory management.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module implements the memory handling
!> routines. Pretty mucch allocations and deallocations.

module qrm_mem_mod

  use qrm_error_mod
  use qrm_const_mod
  implicit none

  !> @brief Generic interface for the
  !! @link qrm_palloc_i @endlink, @link qrm_palloc_2i @endlink, 
  !! @link qrm_palloc_s @endlink, @link qrm_palloc_2s @endlink,
  !! @link qrm_palloc_d @endlink, @link qrm_palloc_2d @endlink,
  !! @link qrm_palloc_c @endlink, @link qrm_palloc_2c @endlink,
  !! @link qrm_palloc_z @endlink, @link qrm_palloc_2z @endlink,
  !! routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they allocate a
  !! pointer of size n. An error is generated or returned in the case
  !! where the pointer is already associated or if the allocation
  !! failed.
  interface qrm_palloc
     module procedure qrm_palloc_i, qrm_palloc_2i 
     module procedure qrm_palloc_s, qrm_palloc_2s
     module procedure qrm_palloc_d, qrm_palloc_2d
     module procedure qrm_palloc_c, qrm_palloc_2c
     module procedure qrm_palloc_z, qrm_palloc_2z
  end interface

  !> @brief Generic interface for the
  !! @link qrm_aalloc_i @endlink, @link qrm_aalloc_2i @endlink, 
  !! @link qrm_aalloc_s @endlink, @link qrm_aalloc_2s @endlink, @link qrm_aalloc_3s @endlink,
  !! @link qrm_aalloc_d @endlink, @link qrm_aalloc_2d @endlink, @link qrm_aalloc_3d @endlink,
  !! @link qrm_aalloc_c @endlink, @link qrm_aalloc_2c @endlink, @link qrm_aalloc_3c @endlink,
  !! @link qrm_aalloc_z @endlink, @link qrm_aalloc_2z @endlink, @link qrm_aalloc_3z @endlink,
  !! routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e. they allocate an
  !! alocatable of size n.  An error is generated or returned in the
  !! case where the allocatable is already allocated or if the
  !! allocation fails.
  interface qrm_aalloc
     module procedure qrm_aalloc_i, qrm_aalloc_2i
     module procedure qrm_aalloc_s, qrm_aalloc_2s, qrm_aalloc_3s
     module procedure qrm_aalloc_d, qrm_aalloc_2d, qrm_aalloc_3d
     module procedure qrm_aalloc_c, qrm_aalloc_2c, qrm_aalloc_3c
     module procedure qrm_aalloc_z, qrm_aalloc_2z, qrm_aalloc_3z
  end interface

  !> @brief Generic interface for the
  !! @link qrm_pdealloc_i @endlink, @link qrm_pdealloc_2i @endlink, 
  !! @link qrm_pdealloc_s @endlink, @link qrm_pdealloc_2s @endlink,
  !! @link qrm_pdealloc_d @endlink, @link qrm_pdealloc_2d @endlink,
  !! @link qrm_pdealloc_c @endlink, @link qrm_pdealloc_2c @endlink,
  !! @link qrm_pdealloc_z @endlink, @link qrm_pdealloc_2z @endlink,
  !! routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they deallocate
  !! a pointer. An error is generated/returned in the case where the
  !! deallocation fails
  interface qrm_pdealloc
     module procedure qrm_pdealloc_i, qrm_pdealloc_2i
     module procedure qrm_pdealloc_s, qrm_pdealloc_2s
     module procedure qrm_pdealloc_d, qrm_pdealloc_2d
     module procedure qrm_pdealloc_c, qrm_pdealloc_2c
     module procedure qrm_pdealloc_z, qrm_pdealloc_2z
  end interface

  !> @brief Generic interface for the
  !! @link qrm_adealloc_i @endlink, @link qrm_adealloc_2i @endlink, 
  !! @link qrm_adealloc_s @endlink, @link qrm_adealloc_2s @endlink, @link qrm_adealloc_3s @endlink,
  !! @link qrm_adealloc_d @endlink, @link qrm_adealloc_2d @endlink, @link qrm_adealloc_3d @endlink,
  !! @link qrm_adealloc_c @endlink, @link qrm_adealloc_2c @endlink, @link qrm_adealloc_3c @endlink,
  !! @link qrm_adealloc_z @endlink, @link qrm_adealloc_2z @endlink, @link qrm_adealloc_3z @endlink,
  !! routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they deallocate
  !! an alocatable. An error is generated/returned in the case where the
  !! deallocation fails
  interface qrm_adealloc
     module procedure qrm_adealloc_i, qrm_adealloc_2i
     module procedure qrm_adealloc_s, qrm_adealloc_2s, qrm_adealloc_3s
     module procedure qrm_adealloc_d, qrm_adealloc_2d, qrm_adealloc_3d
     module procedure qrm_adealloc_c, qrm_adealloc_2c, qrm_adealloc_3c
     module procedure qrm_adealloc_z, qrm_adealloc_2z, qrm_adealloc_3z
  end interface

  !> @brief Generic interface for the 
  !! @link qrm_prealloc_i @endlink
  !! @link qrm_prealloc_s @endlink
  !! @link qrm_prealloc_d @endlink
  !! @link qrm_prealloc_c @endlink
  !! @link qrm_prealloc_z @endlink, routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they reallocate
  !! a pointer. An error is generated/returned in the case where the
  !! reallocation fails
  interface qrm_prealloc
     module procedure qrm_prealloc_i
     module procedure qrm_prealloc_s
     module procedure qrm_prealloc_d
     module procedure qrm_prealloc_c
     module procedure qrm_prealloc_z
  end interface

  !> @brief Generic interface for the 
  !! @link qrm_arealloc_i @endlink
  !! @link qrm_arealloc_s @endlink
  !! @link qrm_arealloc_d @endlink
  !! @link qrm_arealloc_c @endlink
  !! @link qrm_arealloc_z @endlink, routines

  !> All the routines under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they reallocate
  !! an allocatable. An error is generated/returned in the case where the
  !! reallocation fails
  interface qrm_arealloc
     module procedure qrm_arealloc_i
     module procedure qrm_arealloc_s
     module procedure qrm_arealloc_d
     module procedure qrm_arealloc_c
     module procedure qrm_arealloc_z
  end interface

  !> @brief Generic interface for the 
  !! @link qrm_asize_i @endlink, 
  !! @link qrm_asize_s @endlink, @link qrm_asize_2s @endlink, @link qrm_asize_3s @endlink,
  !! @link qrm_asize_d @endlink, @link qrm_asize_2d @endlink, @link qrm_asize_3d @endlink,
  !! @link qrm_asize_c @endlink, @link qrm_asize_2c @endlink, @link qrm_asize_3c @endlink,
  !! @link qrm_asize_z @endlink, @link qrm_asize_2z @endlink, @link qrm_asize_3z @endlink routines

  !> All the functions under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they return the
  !! size of an allocatble. The returned size is 0 if the allocatable
  !! is not allocated
  interface qrm_asize
     module procedure qrm_asize_i
     module procedure qrm_asize_s, qrm_asize_2s, qrm_asize_3s
     module procedure qrm_asize_d, qrm_asize_2d, qrm_asize_3d
     module procedure qrm_asize_c, qrm_asize_2c, qrm_asize_3c
     module procedure qrm_asize_z, qrm_asize_2z, qrm_asize_3z
  end interface

  !> @brief Generic interface for the 
  !! @link qrm_psize_i @endlink
  !! @link qrm_psize_s @endlink
  !! @link qrm_psize_d @endlink
  !! @link qrm_psize_c @endlink
  !! @link qrm_psize_z @endlink, routines

  !> All the functions under this generic interface have the same list
  !! of arguments and perform the same actions, i.e.  they return the
  !! size of a pointer. The returned size is 0 if the pointer
  !! is not associated
  interface qrm_psize
     module procedure qrm_psize_i
     module procedure qrm_psize_s
     module procedure qrm_psize_d
     module procedure qrm_psize_c
     module procedure qrm_psize_z
  end interface

  integer :: qrm_mem_thn=0, qrm_mem_nth=1
  !$omp threadprivate(qrm_mem_thn)

  !>  a counter to keep track of the currently allocated memory, per thread
  integer(kind=8) :: qrm_tot_mem(0:qrm_maxthreads-1)=0

  !>  a counter to keep track of the peak memory, per thread
  integer(kind=8) :: qrm_max_mem(0:qrm_maxthreads-1)=0

  !>  scratchpad variable to store the memory peak on entry to a
  !!  parallel region
  integer(kind=8) :: qrm_seq_peak=0

  integer :: qrm_exact_mem = qrm_no_

  integer(kind=8), parameter :: qrm_sizeof_i_=4
  integer(kind=8), parameter :: qrm_sizeof_s_=4
  integer(kind=8), parameter :: qrm_sizeof_d_=8
  integer(kind=8), parameter :: qrm_sizeof_c_=8
  integer(kind=8), parameter :: qrm_sizeof_z_=16


#if defined(memlim)
  integer(kind=8)       :: qrm_mem_lim=500000000
#endif

contains

  !> @brief This routine has to be called at the beginning of a
  !! parallel section. Afterwards, each thread will update its own
  !! tot_mem and max_mem variables. This is done to avoid
  !! synchronizations on the update of statistics on memory
  !! consumption
  subroutine qrm_par_mem_init()
    !$ use omp_lib
    if(qrm_exact_mem .eq. qrm_yes_) then
       ! do nothing
    else
       !$ qrm_mem_thn = omp_get_thread_num()
       !$omp master
       !$ qrm_mem_nth = omp_get_num_threads()
       qrm_seq_peak = qrm_max_mem(0)
       qrm_tot_mem(1:qrm_mem_nth-1)=0
       qrm_max_mem(0:qrm_mem_nth-1)=0
       !$omp end master
       !$omp barrier
    end if
    return
  end subroutine qrm_par_mem_init
    
  subroutine qrm_par_mem_finalize()
    if(qrm_exact_mem .eq. qrm_yes_) then
       ! do nothing
    else
       !$omp barrier
       !$omp master
       qrm_tot_mem(0) = sum(qrm_tot_mem(0:qrm_mem_nth-1))
       qrm_max_mem(0) = max(sum(qrm_max_mem(0:qrm_mem_nth-1)), &
            & qrm_seq_peak)
       qrm_seq_peak=0
       qrm_tot_mem(1:qrm_mem_nth-1)=0
       qrm_max_mem(1:qrm_mem_nth-1)=0
       !$omp end master
    end if
    return
  end subroutine qrm_par_mem_finalize

  !> @brief updates memory statistics
  !!
  !! @param[in] n the amount of memory to be added (can be negative
  !!            for deallocations)
  subroutine qrm_mem_upd(n)
    integer(kind=8) :: n
    if(qrm_exact_mem .eq. qrm_yes_) then
       !$omp critical(mem)
       qrm_tot_mem(0) = qrm_tot_mem(0)+n
       if(qrm_tot_mem(0) .gt. qrm_max_mem(0)) &
            & qrm_max_mem(0) = qrm_tot_mem(0)
       !$omp end critical(mem)
    else
       qrm_tot_mem(qrm_mem_thn) = qrm_tot_mem(qrm_mem_thn)+n
       if(qrm_tot_mem(qrm_mem_thn) .gt. qrm_max_mem(qrm_mem_thn)) &
            & qrm_max_mem(qrm_mem_thn) = qrm_tot_mem(qrm_mem_thn)
    end if
    return
  end subroutine qrm_mem_upd



  !> @param[in] n the size of the pointer
  !!
  !! @param[in,out] a the pointer to be allocated.
  !! 
  !! @param[out] info  (optional). if this argument is present, then the callee wants
  !!        to handle the error on his side 
  !! 
  subroutine qrm_palloc_d(a, n, info)

    real(kind(1.d0)), pointer, dimension(:) :: a
    integer, intent(in)                     :: n
    integer, optional                       :: info

    integer :: err, disp

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_d')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_d_ .gt. disp ) then
          err = 1
       else
          allocate(a(n), stat=err)
       end if
#else
       allocate(a(n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_d',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_d_)
       end if
    end if

    return

  end subroutine qrm_palloc_d


  !> @param[in] n the size of the pointer
  !!
  !! @param[in,out] a the pointer to be allocated.
  !! 
  !! @param[out] info  (optional). if this argument is present, then the callee wants
  !!        to handle the error on his side 
  !! 
  subroutine qrm_palloc_s(a, n, info)

    real(kind(1.e0)), pointer, dimension(:) :: a
    integer, intent(in)                     :: n
    integer, optional                       :: info

    integer :: err, disp

    if(n .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_s')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1)
       !$omp end critical(mem)
       if( n*qrm_sizeof_s_ .gt. disp ) then
          err = 1
       else
          allocate(a(n), stat=err)
       end if
#else
       allocate(a(n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_s',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_s_)
       end if

    end if

    return
    
  end subroutine qrm_palloc_s


  !> @param[in] n the size of the pointer
  !!
  !! @param[in,out] a the pointer to be allocated.
  !! 
  !! @param[out] info  (optional). if this argument is present, then the callee wants
  !!        to handle the error on his side 
  !! 
  subroutine qrm_palloc_i(a, n, info)

    integer, pointer, dimension(:) :: a
    integer, intent(in)            :: n
    integer, optional              :: info

    integer :: err, disp

    if(n .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_i')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1)
       !$omp end critical(mem)
       if( n*qrm_sizeof_i_ .gt. disp ) then
          err = 1
       else
          allocate(a(n), stat=err)
       end if
#else
       allocate(a(n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_i',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_i_)
       end if

    end if

    return

  end subroutine qrm_palloc_i



  !> @param[in,out] a the pointer to be allocated
  !!
  !! @param[in] m     the rank-1 size of the pointer
  !! 
  !! @param[in] n     the rank-2 size of the pointer
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_palloc_2d(a, m, n, info)

    real(kind(1.d0)), pointer, dimension(:,:) :: a
    integer, intent(in)                       :: m, n
    integer, optional                         :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_2d')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_d_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_2d',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_d_)
       end if

    end if

    return

  end subroutine qrm_palloc_2d


  !> @param[in,out] a the pointer to be allocated
  !!
  !! @param[in] m     the rank-1 size of the pointer
  !! 
  !! @param[in] n     the rank-2 size of the pointer
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_palloc_2s(a, m, n, info)

    real(kind(1.e0)), pointer, dimension(:,:) :: a
    integer, intent(in)                       :: m, n
    integer, optional                         :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_2s')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_s_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_2s',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_s_)
       end if

    end if

    return

  end subroutine qrm_palloc_2s


  !> @param[in,out] a the pointer to be allocated
  !!
  !! @param[in] m     the rank-1 size of the pointer
  !! 
  !! @param[in] n     the rank-2 size of the pointer
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_palloc_2i(a, m, n, info)

    integer, pointer, dimension(:,:) :: a
    integer, intent(in)              :: m, n
    integer, optional                :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2i')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_i_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2i',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_i_)
       end if

    end if

    return

  end subroutine qrm_palloc_2i

  !> @param[in,out] a the pointer to be allocated
  !!
  !! @param[in] m     the rank-1 size of the pointer
  !! 
  !! @param[in] n     the rank-2 size of the pointer
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_palloc_2z(a, m, n, info)

    complex(kind(1.d0)), pointer, dimension(:,:) :: a
    integer, intent(in)                          :: m, n
    integer, optional                            :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_2z')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_z_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2z',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_z_)
       end if

    end if

    return

  end subroutine qrm_palloc_2z

  !> @param[in,out] a the pointer to be allocated
  !!
  !! @param[in] m     the rank-1 size of the pointer
  !! 
  !! @param[in] n     the rank-2 size of the pointer
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_palloc_2c(a, m, n, info)

    complex(kind(1.e0)), pointer, dimension(:,:) :: a
    integer, intent(in)                          :: m, n
    integer, optional                            :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_2c')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_c_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2c',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_c_)
       end if

    end if

    return

  end subroutine qrm_palloc_2c





  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] n     the size of the allocatable
  !! 
  !! @param[in] lbnd  optional integer specifying the lower bound
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_d(a, n, lbnd, info)

    real(kind(1.d0)), allocatable, dimension(:) :: a
    integer, intent(in)                         :: n
    integer, optional                           :: lbnd
    integer, optional                           :: info

    integer :: err, ilbnd, disp

    if(n .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_d')
    else
       if(present(lbnd)) then
          ilbnd = lbnd
       else
          ilbnd = 1
       end if

#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_d_ .gt. disp ) then
          err = 1
       else
          allocate(a(ilbnd: ilbnd+n-1), stat=err)
       end if
#else
       allocate(a(ilbnd: ilbnd+n-1), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_d',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_d_)
       end if

    end if

    return

  end subroutine qrm_aalloc_d


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] n     the size of the allocatable
  !! 
  !! @param[in] lbnd  optional integer specifying the lower bound
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_s(a, n, lbnd, info)

    real(kind(1.e0)), allocatable, dimension(:) :: a
    integer, intent(in)                         :: n
    integer, optional                           :: lbnd
    integer, optional                           :: info

    integer :: err, ilbnd, disp

    if(n .lt. 0) return
    err = 0
    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_s')
    else
       if(present(lbnd)) then
          ilbnd = lbnd
       else
          ilbnd = 1
       end if

#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_s_ .gt. disp ) then
          err = 1
       else
          allocate(a(ilbnd: ilbnd+n-1), stat=err)
       end if
#else
       allocate(a(ilbnd: ilbnd+n-1), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_s',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_s_)
       end if

    end if

    return

  end subroutine qrm_aalloc_s


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] n     the size of the allocatable
  !! 
  !! @param[in] lbnd  optional integer specifying the lower bound
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_i(a, n, lbnd, info)

    integer, allocatable, dimension(:) :: a
    integer, intent(in)                :: n
    integer, optional                  :: lbnd
    integer, optional                  :: info

    integer :: err, ilbnd, disp

    if(n .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_i')
    else
       if(present(lbnd)) then
          ilbnd = lbnd
       else
          ilbnd = 1
       end if

#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_i_ .gt. disp ) then
          err = 1
       else
          allocate(a(ilbnd: ilbnd+n-1), stat=err)
       end if
#else
       allocate(a(ilbnd: ilbnd+n-1), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_i',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_i_)
       end if

    end if

    return

  end subroutine qrm_aalloc_i


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_2d(a, m, n, info)

    real(kind(1.d0)), allocatable, dimension(:,:) :: a
    integer, intent(in)                           :: m, n
    integer, optional                             :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2d')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_d_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2d',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_d_)
       end if

    end if

    return

  end subroutine qrm_aalloc_2d


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_2s(a, m, n, info)

    real(kind(1.e0)), allocatable, dimension(:,:) :: a
    integer, intent(in)                           :: m, n
    integer, optional                             :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2s')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_s_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2s',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_s_)
       end if

    end if

    return

  end subroutine qrm_aalloc_2s


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_2i(a, m, n, info)

    integer, allocatable, dimension(:,:) :: a
    integer, intent(in)                  :: m, n
    integer, optional                    :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2i')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_i_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2i',ied=(/m*n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*qrm_sizeof_i_)
       end if

    end if

    return

  end subroutine qrm_aalloc_2i



  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[in] k     the rank-3 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_3d(a, m, n, k, info)

    real(kind(1.d0)), allocatable, dimension(:,:,:) :: a
    integer, intent(in)                             :: m, n, k
    integer, optional                               :: info

    integer :: err, disp

    if(min(min(m,n),k) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_3d')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*k*qrm_sizeof_d_ .gt. disp ) then
          err = 1
       else
       allocate(a(m,n,k), stat=err)
       end if
#else
       allocate(a(m,n,k), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_3d',ied=(/m*n*k,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_d_)
       end if

    end if

    return

  end subroutine qrm_aalloc_3d


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[in] k     the rank-3 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_3s(a, m, n, k, info)

    real(kind(1.e0)), allocatable, dimension(:,:,:) :: a
    integer, intent(in)                             :: m, n, k
    integer, optional                               :: info

    integer :: err, disp

    if(min(min(m,n),k) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_3s')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*k*qrm_sizeof_s_ .gt. disp ) then
          err = 1
       else
       allocate(a(m,n,k), stat=err)
       end if
#else
       allocate(a(m,n,k), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_3s',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(m*int(n,8)*int(k,8)*qrm_sizeof_s_)
       end if

    end if

    return

  end subroutine qrm_aalloc_3s



  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_d(a)

    real(kind(1.d0)), pointer, dimension(:) :: a

    integer :: err=0, n
    
    if(associated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_d',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_d_)
    end if

    return

  end subroutine qrm_pdealloc_d

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_s(a)

    real(kind(1.e0)), pointer, dimension(:) :: a

    integer :: err=0, n

    if(associated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else 
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_s',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_s_)
    end if
 
    return
    
  end subroutine qrm_pdealloc_s

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_2d(a)

    real(kind(1.d0)), pointer, dimension(:,:) :: a

    integer :: err=0, n, m

    if(associated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_2d',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_d_)
    end if

    return
    
  end subroutine qrm_pdealloc_2d

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_2s(a)

    real(kind(1.e0)), pointer, dimension(:,:) :: a

    integer :: err=0, n, m

    if(associated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_2s',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_s_)
    end if

    return
    
  end subroutine qrm_pdealloc_2s

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_i(a)

    integer, pointer, dimension(:) :: a

    integer :: err=0, n

    if(associated(a)) then
       n=size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_i',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_i_)
    end if

    return

  end subroutine qrm_pdealloc_i


  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_2i(a)

    integer, pointer, dimension(:,:) :: a

    integer :: err=0, n, m

    if(associated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_2i',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_i_)
    end if

    return
    
  end subroutine qrm_pdealloc_2i



  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_d(a)

    real(kind(1.d0)), allocatable, dimension(:) :: a

    integer :: err=0, n


    if(allocated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_d',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_d_)
    end if

    return

  end subroutine qrm_adealloc_d

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_s(a)

    real(kind(1.e0)), allocatable, dimension(:) :: a

    integer :: err=0, n

    if(allocated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_s',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_s_)
    end if

    return
    
  end subroutine qrm_adealloc_s

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_i(a)

    integer, allocatable, dimension(:) :: a

    integer :: err=0, n

    if(allocated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_i', ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_i_)
    end if

    return

  end subroutine qrm_adealloc_i

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_2d(a)

    real(kind(1.d0)), allocatable, dimension(:,:) :: a

    integer :: err=0, n, m

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_2d',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_d_)
    end if

    return
    
  end subroutine qrm_adealloc_2d

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_2s(a)

    real(kind(1.e0)), allocatable, dimension(:,:) :: a

    integer :: err=0, n, m

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_2s',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_s_)
    end if

    return
    
  end subroutine qrm_adealloc_2s


  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_3d(a)

    real(kind(1.d0)), allocatable, dimension(:,:,:) :: a

    integer :: err=0, n, m, k

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       k = size(a,3)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_3d',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_d_)
    end if

    return
    
  end subroutine qrm_adealloc_3d

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_3s(a)

    real(kind(1.e0)), allocatable, dimension(:,:,:) :: a

    integer :: err=0, n, m, k

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       k = size(a,3)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_3s',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_s_)
    end if

    return
    
  end subroutine qrm_adealloc_3s



  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_2i(a)

    integer, allocatable, dimension(:,:) :: a

    integer :: err=0, n, m

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_2i',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_i_)
    end if

    return
    
  end subroutine qrm_adealloc_2i


  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_prealloc_d(a, n, force, copy)
  
    real(kind(1.d0)), pointer, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    real(kind(1.d0)), pointer, dimension(:) :: tmp=>null()

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not associated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. associated(a))
    
    if(associated(a)) then
       ! a is associated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(icopy) then
             ! we must save a copy
             tmp => a
             nullify(a)
          else
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_prealloc_d',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_d_)
             end if
          end if
          ! reallocate a
          allocate(a(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_palloc_d',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(int(n,8)*qrm_sizeof_d_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(a(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_palloc_d',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(int(n,8)*qrm_sizeof_d_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          a(i) = tmp(i)
       end do
       deallocate(tmp, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_prealloc_d',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_d_)
       end if
    end if
    
    return

  end subroutine qrm_prealloc_d
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_prealloc_s(a, n, force, copy)
  
    real(kind(1.e0)), pointer, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    real(kind(1.e0)), pointer, dimension(:) :: tmp=>null()

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not associated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. associated(a))
    
    if(associated(a)) then
       ! a is associated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(icopy) then
             ! we must save a copy
             tmp => a
             nullify(a)
          else
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_prealloc_s',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_s_)
             end if
          end if
          ! reallocate a
          allocate(a(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_prealloc_s',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_s_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(a(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_prealloc_s',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_s_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          a(i) = tmp(i)
       end do
       deallocate(tmp, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_prealloc_s',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_s_)
       end if
    end if
    
    return

  end subroutine qrm_prealloc_s
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_prealloc_i(a, n, force, copy)
  
    integer, pointer, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    integer, pointer, dimension(:) :: tmp=>null()

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not associated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. associated(a))
    
    if(associated(a)) then
       ! a is associated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(icopy) then
             ! we must save a copy
             tmp => a
             nullify(a)
          else
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_prealloc_i',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_i_)
             end if
          end if
          ! reallocate a
          allocate(a(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_prealloc_i',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_i_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(a(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_prealloc_i',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_i_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          a(i) = tmp(i)
       end do
       deallocate(tmp, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_prealloc_i',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_i_)
       end if
    end if
    
    return

  end subroutine qrm_prealloc_i
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_arealloc_d(a, n, force, copy)
  
    real(kind(1.d0)), allocatable, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    real(kind(1.d0)), allocatable, dimension(:) :: tmp

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not allocated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. allocated(a))
    
    if(allocated(a)) then
       ! a is allocated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(.not. icopy) then
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_arealloc_d',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_d_)
             end if
          end if
          ! allocate tmp
          allocate(tmp(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_arealloc_d',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_d_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(tmp(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_arealloc_d',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_d_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          tmp(i) = a(i)
       end do
       deallocate(a, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_arealloc_d',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_d_)
       end if
    end if
    
    call move_alloc(from=tmp, to=a)

    return

  end subroutine qrm_arealloc_d
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_arealloc_s(a, n, force, copy)
  
    real(kind(1.e0)), allocatable, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    real(kind(1.e0)), allocatable, dimension(:) :: tmp

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not allocated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. allocated(a))
    
    if(allocated(a)) then
       ! a is allocated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(.not. icopy) then
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_arealloc_s',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_s_)
             end if
          end if
          ! allocate tmp
          allocate(tmp(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_arealloc_s',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_s_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(tmp(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_arealloc_s',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_s_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          tmp(i) = a(i)
       end do
       deallocate(a, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_arealloc_s',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_s_)
       end if
    end if
    
    call move_alloc(from=tmp, to=a)

    return

  end subroutine qrm_arealloc_s
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_arealloc_i(a, n, force, copy)
  
    integer, allocatable, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    integer, allocatable, dimension(:) :: tmp

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not allocated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. allocated(a))
    
    if(allocated(a)) then
       ! a is allocated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(.not. icopy) then
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_arealloc_i',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-asize*qrm_sizeof_i_)
             end if
          end if
          ! allocate tmp
          allocate(tmp(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_arealloc_i',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_i_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(tmp(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_arealloc_i',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_i_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          tmp(i) = a(i)
       end do
       deallocate(a, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_arealloc_i',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_i_)
       end if
    end if
    
    call move_alloc(from=tmp, to=a)

    return

  end subroutine qrm_arealloc_i
  



  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_asize_i(a)
  
    implicit none
    
    integer :: qrm_asize_i
    integer, allocatable :: a(:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_i = size(a)
    else
       qrm_asize_i = 0
    end if
  
    return
  
  end function qrm_asize_i
  

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_asize_s(a)

    implicit none
    
    integer :: qrm_asize_s
    real(kind(1.e0)), allocatable :: a(:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_s = size(a)
    else
       qrm_asize_s = 0
    end if
  
    return
  
  end function qrm_asize_s

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_asize_d(a)
  
    implicit none
    
    integer :: qrm_asize_d
    real(kind(1.d0)), allocatable :: a(:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_d = size(a)
    else
       qrm_asize_d = 0
    end if
  
    return
  
  end function qrm_asize_d


  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_2s(a)
  
    implicit none
    
    integer :: qrm_asize_2s
    real(kind(1.e0)), allocatable :: a(:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_2s = size(a)
    else
       qrm_asize_2s = 0
    end if
  
    return
  
  end function qrm_asize_2s

  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_2d(a)
  
    implicit none
    
    integer :: qrm_asize_2d
    real(kind(1.d0)), allocatable :: a(:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_2d = size(a)
    else
       qrm_asize_2d = 0
    end if
  
    return
  
  end function qrm_asize_2d



  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_3s(a)
  
    implicit none
    
    integer :: qrm_asize_3s
    real(kind(1.e0)), allocatable :: a(:,:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_3s = size(a)
    else
       qrm_asize_3s = 0
    end if
  
    return
  
  end function qrm_asize_3s


  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_3d(a)
  
    implicit none
    
    integer :: qrm_asize_3d
    real(kind(1.d0)), allocatable :: a(:,:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_3d = size(a)
    else
       qrm_asize_3d = 0
    end if
  
    return
  
  end function qrm_asize_3d


  !> @param[in] n the size of the pointer
  !!
  !! @param[in,out] a the pointer to be allocated.
  !! 
  !! @param[out] info  (optional). if this argument is present, then the callee wants
  !!        to handle the error on his side 
  !! 
  subroutine qrm_palloc_z(a, n, info)

    complex(kind(1.d0)), pointer, dimension(:) :: a
    integer, intent(in)                     :: n
    integer, optional                       :: info

    integer :: err, disp

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_z')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_z_ .gt. disp ) then
          err = 1
       else
          allocate(a(n), stat=err)
       end if
#else
       allocate(a(n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_z',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
       end if
    end if

    return

  end subroutine qrm_palloc_z


  !> @param[in] n the size of the pointer
  !!
  !! @param[in,out] a the pointer to be allocated.
  !! 
  !! @param[out] info  (optional). if this argument is present, then the callee wants
  !!        to handle the error on his side 
  !! 
  subroutine qrm_palloc_c(a, n, info)

    complex(kind(1.e0)), pointer, dimension(:) :: a
    integer, intent(in)                     :: n
    integer, optional                       :: info

    integer :: err, disp

    if(n .lt. 0) return

    if(associated(a)) then
       call qrm_err_push(4,sub='qrm_palloc_c')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_c_ .gt. disp ) then
          err = 1
       else
          allocate(a(n), stat=err)
       end if
#else
       allocate(a(n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_palloc_c',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
       end if

    end if

    return
    
  end subroutine qrm_palloc_c


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] n     the size of the allocatable
  !! 
  !! @param[in] lbnd  optional integer specifying the lower bound
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_z(a, n, lbnd, info)

    complex(kind(1.d0)), allocatable, dimension(:) :: a
    integer, intent(in)                         :: n
    integer, optional                           :: lbnd
    integer, optional                           :: info

    integer :: err, ilbnd, disp

    if(n .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_z')
    else
       if(present(lbnd)) then
          ilbnd = lbnd
       else
          ilbnd = 1
       end if

#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_z_ .gt. disp ) then
          err = 1
       else
          allocate(a(ilbnd: ilbnd+n-1), stat=err)
       end if
#else
       allocate(a(ilbnd: ilbnd+n-1), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_z',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
       end if

    end if

    return

  end subroutine qrm_aalloc_z

  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] n     the size of the allocatable
  !! 
  !! @param[in] lbnd  optional integer specifying the lower bound
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_c(a, n, lbnd, info)

    complex(kind(1.e0)), allocatable, dimension(:) :: a
    integer, intent(in)                         :: n
    integer, optional                           :: lbnd
    integer, optional                           :: info

    integer :: err, ilbnd, disp

    if(n .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_c')
    else
       if(present(lbnd)) then
          ilbnd = lbnd
       else
          ilbnd = 1
       end if

#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*qrm_sizeof_c_ .gt. disp ) then
          err = 1
       else
          allocate(a(ilbnd: ilbnd+n-1), stat=err)
       end if
#else
       allocate(a(ilbnd: ilbnd+n-1), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_c',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
       end if

    end if

    return

  end subroutine qrm_aalloc_c


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_2z(a, m, n, info)

    complex(kind(1.d0)), allocatable, dimension(:,:) :: a
    integer, intent(in)                           :: m, n
    integer, optional                             :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2z')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_z_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2z',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(m,8)*int(n,8)*qrm_sizeof_z_)
       end if

    end if

    return

  end subroutine qrm_aalloc_2z

  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_2c(a, m, n, info)

    complex(kind(1.e0)), allocatable, dimension(:,:) :: a
    integer, intent(in)                           :: m, n
    integer, optional                             :: info

    integer :: err, disp

    if(min(m,n) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_2c')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*qrm_sizeof_c_ .gt. disp ) then
          err = 1
       else
          allocate(a(m,n), stat=err)
       end if
#else
       allocate(a(m,n), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_2c',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(m,8)*int(n,8)*qrm_sizeof_c_)
       end if

    end if

    return

  end subroutine qrm_aalloc_2c




  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[in] k     the rank-3 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_3z(a, m, n, k, info)

    complex(kind(1.d0)), allocatable, dimension(:,:,:) :: a
    integer, intent(in)                             :: m, n, k
    integer, optional                               :: info

    integer :: err, disp

    if(min(min(m,n),k) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_3z')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*k*qrm_sizeof_z_ .gt. disp ) then
          err = 1
       else
       allocate(a(m,n,k), stat=err)
       end if
#else
       allocate(a(m,n,k), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_3z',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_z_)
       end if

    end if

    return

  end subroutine qrm_aalloc_3z


  !> @param[in,out] a the allocatable to be allocated
  !!
  !! @param[in] m     the rank-1 size of the allocatable
  !! 
  !! @param[in] n     the rank-2 size of the allocatable
  !! 
  !! @param[in] k     the rank-3 size of the allocatable
  !! 
  !! @param[out] info (optional). if this argument is present, then the callee wants
  !!                  to handle the error on his side 
  !! 
  subroutine qrm_aalloc_3c(a, m, n, k, info)

    complex(kind(1.e0)), allocatable, dimension(:,:,:) :: a
    integer, intent(in)                             :: m, n, k
    integer, optional                               :: info

    integer :: err, disp

    if(min(min(m,n),k) .lt. 0) return

    if(allocated(a)) then
       call qrm_err_push(4,sub='qrm_aalloc_3c')
    else
#if defined(memlim)
       !$omp critical(mem)
       disp = qrm_mem_lim - sum(qrm_tot_mem(0:qrm_mem_nth-1))
       !$omp end critical(mem)
       if( n*m*k*qrm_sizeof_c_ .gt. disp ) then
          err = 1
       else
       allocate(a(m,n,k), stat=err)
       end if
#else
       allocate(a(m,n,k), stat=err)
#endif
       if(err .ne. 0) then
          if(present(info)) then
             info = err
          else
             call qrm_err_push(5,sub='qrm_aalloc_3c',ied=(/n,0,0,0,0/))
          end if
       else
          call qrm_mem_upd(+int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_c_)
       end if

    end if

    return

  end subroutine qrm_aalloc_3c



!!!!!!!!!!!!!!!!
!! Deallocation
!!!!!!!!!!!!!!!!


  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_z(a)

    complex(kind(1.d0)), pointer, dimension(:) :: a

    integer :: err=0, n
    
    if(associated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_z',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_z_)
    end if

    return

  end subroutine qrm_pdealloc_z

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_c(a)

    complex(kind(1.e0)), pointer, dimension(:) :: a

    integer :: err=0, n

    if(associated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else 
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_c',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_c_)
    end if
 
    return
    
  end subroutine qrm_pdealloc_c


  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_2z(a)

    complex(kind(1.d0)), pointer, dimension(:,:) :: a

    integer :: err=0, n, m

    if(associated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_2z',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_z_)
    end if

    return
    
  end subroutine qrm_pdealloc_2z

  !> @param[in,out] a the pointer to be deallocated. 
  !! 
  subroutine qrm_pdealloc_2c(a)

    complex(kind(1.e0)), pointer, dimension(:,:) :: a

    integer :: err=0, n, m

    if(associated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_pdealloc_2c',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_c_)
    end if

    return
    
  end subroutine qrm_pdealloc_2c


  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_z(a)

    complex(kind(1.d0)), allocatable, dimension(:) :: a

    integer :: err=0, n


    if(allocated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_z',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_z_)
    end if

    return

  end subroutine qrm_adealloc_z

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_c(a)

    complex(kind(1.e0)), allocatable, dimension(:) :: a

    integer :: err=0, n

    if(allocated(a)) then
       n = size(a)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_c',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(n,8)*qrm_sizeof_c_)
    end if

    return
    
  end subroutine qrm_adealloc_c

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_2z(a)

    complex(kind(1.d0)), allocatable, dimension(:,:) :: a

    integer :: err=0, n, m

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_2z',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_z_)
    end if

    return
    
  end subroutine qrm_adealloc_2z

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_2c(a)

    complex(kind(1.e0)), allocatable, dimension(:,:) :: a

    integer :: err=0, n, m

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_2c',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*qrm_sizeof_c_)
    end if

    return
    
  end subroutine qrm_adealloc_2c


  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_3z(a)

    complex(kind(1.d0)), allocatable, dimension(:,:,:) :: a

    integer :: err=0, n, m, k

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       k = size(a,3)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_3z',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_z_)
    end if

    return
    
  end subroutine qrm_adealloc_3z

  !> @param[in,out] a the allocatable to be deallocated. 
  !! 
  subroutine qrm_adealloc_3c(a)

    complex(kind(1.e0)), allocatable, dimension(:,:,:) :: a

    integer :: err=0, n, m, k

    if(allocated(a)) then
       m = size(a,1)
       n = size(a,2)
       k = size(a,3)
       deallocate(a, stat=err)
    else
       return
    end if
    if(err .ne. 0) then
       call qrm_err_push(7,sub='qrm_adealloc_3c',ied=(/err,0,0,0,0/))
    else
       call qrm_mem_upd(-int(m,8)*int(n,8)*int(k,8)*qrm_sizeof_c_)
    end if

    return
    
  end subroutine qrm_adealloc_3c



!!!!!!!!!!!!!!!
!! Reallocation
!!!!!!!!!!!!!!!
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_prealloc_z(a, n, force, copy)
  
    complex(kind(1.d0)), pointer, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    complex(kind(1.d0)), pointer, dimension(:) :: tmp=>null()

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not associated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. associated(a))
    
    if(associated(a)) then
       ! a is associated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(icopy) then
             ! we must save a copy
             tmp => a
             nullify(a)
          else
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_prealloc_z',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-int(n,8)*qrm_sizeof_z_)
             end if
          end if
          ! reallocate a
          allocate(a(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_palloc_z',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(a(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_palloc_z',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          a(i) = tmp(i)
       end do
       deallocate(tmp, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_prealloc_z',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_z_)
       end if
    end if
    
    return

  end subroutine qrm_prealloc_z
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_prealloc_c(a, n, force, copy)
  
    complex(kind(1.e0)), pointer, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    complex(kind(1.e0)), pointer, dimension(:) :: tmp=>null()

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not associated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. associated(a))
    
    if(associated(a)) then
       ! a is associated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(icopy) then
             ! we must save a copy
             tmp => a
             nullify(a)
          else
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_prealloc_c',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-int(n,8)*qrm_sizeof_c_)
             end if
          end if
          ! reallocate a
          allocate(a(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_prealloc_c',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(a(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_prealloc_c',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          a(i) = tmp(i)
       end do
       deallocate(tmp, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_prealloc_c',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_c_)
       end if
    end if
    
    return

  end subroutine qrm_prealloc_c
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_arealloc_z(a, n, force, copy)
  
    complex(kind(1.d0)), allocatable, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    complex(kind(1.d0)), allocatable, dimension(:) :: tmp

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not allocated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. allocated(a))
    
    if(allocated(a)) then
       ! a is allocated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(.not. icopy) then
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_arealloc_z',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-int(n,8)*qrm_sizeof_z_)
             end if
          end if
          ! allocate tmp
          allocate(tmp(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_arealloc_z',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(tmp(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_arealloc_z',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_z_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          tmp(i) = a(i)
       end do
       deallocate(a, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_arealloc_z',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_z_)
       end if
    end if
    
    call move_alloc(from=tmp, to=a)

    return

  end subroutine qrm_arealloc_z
  

  !> @param[in,out] a the array that has to be reallocated
  !! @param[in] n     the minimum size of a on exit
  !! @param[in] force if force=.true. a is always reallocated
  !! @param[in] copy  if copy=.true. and force=.false. the content of the old area is
  !!                  copied into the new
  !!
  subroutine qrm_arealloc_c(a, n, force, copy)
  
    complex(kind(1.e0)), allocatable, dimension(:) :: a
    integer                                 :: n
    logical, optional                       :: force, copy

    integer :: err=0, asize, i
    logical :: iforce, icopy
    complex(kind(1.e0)), allocatable, dimension(:) :: tmp

    iforce=.false.
    if(present(force)) iforce=force

    ! if iforce=.true. we don't make any copies. Also, if a is not allocated
    ! then there's nothing to copy.
    icopy = .false.
    if(present(copy)) icopy=(copy .and. (.not. iforce) .and. allocated(a))
    
    if(allocated(a)) then
       ! a is allocated
       asize = size(a)
       if(iforce .or. (size(a) .lt. n)) then
          ! we need to reallocate
          if(.not. icopy) then
             deallocate(a, stat=err)
             if(err .ne. 0) then
                call qrm_err_push(7, sub='qrm_arealloc_c',ied=(/err,0,0,0,0/))
             else
                call qrm_mem_upd(-int(n,8)*qrm_sizeof_c_)
             end if
          end if
          ! allocate tmp
          allocate(tmp(n), stat=err)
          if(err .ne. 0) then
             call qrm_err_push(5,sub='qrm_arealloc_c',ied=(/n,0,0,0,0/))
          else
             call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
          end if
       else
          return
       end if
    else
       ! allocate a
       allocate(tmp(n), stat=err)
       if(err .ne. 0) then
          call qrm_err_push(5,sub='qrm_arealloc_c',ied=(/n,0,0,0,0/))
       else
          call qrm_mem_upd(+int(n,8)*qrm_sizeof_c_)
       end if
    end if
       
    ! check if copy is to be done
    if(icopy) then
       do i=1, asize
          tmp(i) = a(i)
       end do
       deallocate(a, stat=err)
       if(err .ne. 0) then
          call qrm_err_push(7, sub='qrm_arealloc_c',ied=(/err,0,0,0,0/))
       else
          call qrm_mem_upd(-asize*qrm_sizeof_c_)
       end if
    end if
    
    call move_alloc(from=tmp, to=a)

    return

  end subroutine qrm_arealloc_c
  

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_asize_c(a)
  
    implicit none
    
    integer :: qrm_asize_c
    complex(kind(1.e0)), allocatable :: a(:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_c = size(a)
    else
       qrm_asize_c = 0
    end if
  
    return
  
  end function qrm_asize_c

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_asize_z(a)
  
    implicit none
    
    integer :: qrm_asize_z
    complex(kind(1.d0)), allocatable :: a(:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_z = size(a)
    else
       qrm_asize_z = 0
    end if
  
    return
  
  end function qrm_asize_z

 
  !> @param[in]  a            the arrays whose size must be computed
  !!
 function qrm_asize_2c(a)
  
    implicit none
    
    integer :: qrm_asize_2c
    complex(kind(1.e0)), allocatable :: a(:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_2c = size(a)
    else
       qrm_asize_2c = 0
    end if
  
    return
  
  end function qrm_asize_2c

  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_2z(a)
  
    implicit none
    
    integer :: qrm_asize_2z
    complex(kind(1.d0)), allocatable :: a(:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_2z = size(a)
    else
       qrm_asize_2z = 0
    end if
  
    return
  
  end function qrm_asize_2z



  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_3c(a)
  
    implicit none
    
    integer :: qrm_asize_3c
    complex(kind(1.e0)), allocatable :: a(:,:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_3c = size(a)
    else
       qrm_asize_3c = 0
    end if
  
    return
  
  end function qrm_asize_3c


  !> @param[in]  a            the arrays whose size must be computed
  !!
  function qrm_asize_3z(a)
  
    implicit none
    
    integer :: qrm_asize_3z
    complex(kind(1.d0)), allocatable :: a(:,:,:)

    if(allocated(a)) then
       ! a is allocated    
       qrm_asize_3z = size(a)
    else
       qrm_asize_3z = 0
    end if
  
    return
  
  end function qrm_asize_3z

! ===========================================================================================



  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_psize_i(a)
  
    implicit none
    
    integer :: qrm_psize_i
    integer, pointer :: a(:)

    if(associated(a)) then
       ! a is allocated    
       qrm_psize_i = size(a)
    else
       qrm_psize_i = 0
    end if
  
    return
  
  end function qrm_psize_i
  

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_psize_s(a)
  
    implicit none
    
    integer :: qrm_psize_s
    real(kind(1.e0)), pointer :: a(:)

    if(associated(a)) then
       ! a is allocated    
       qrm_psize_s = size(a)
    else
       qrm_psize_s = 0
    end if
  
    return
  
  end function qrm_psize_s



  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_psize_c(a)
  
    implicit none
    
    integer :: qrm_psize_c
    complex(kind(1.e0)), pointer :: a(:)

    if(associated(a)) then
       ! a is allocated    
       qrm_psize_c = size(a)
    else
       qrm_psize_c = 0
    end if
  
    return
    
  end function qrm_psize_c
  
  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_psize_z(a)
  
    implicit none
    
    integer :: qrm_psize_z
    complex(kind(1.d0)), pointer :: a(:)

    if(associated(a)) then
       ! a is allocated    
       qrm_psize_z = size(a)
    else
       qrm_psize_z = 0
    end if
  
    return
  
  end function qrm_psize_z

  !> @param[in]  a           the arrays whose size must be computed
  !!
  function qrm_psize_d(a)

    implicit none
    
    integer :: qrm_psize_d
    real(kind(1.d0)), pointer :: a(:)

    if(associated(a)) then
       ! a is allocated    
       qrm_psize_d = size(a)
    else
       qrm_psize_d = 0
    end if
  
    return
    
  end function qrm_psize_d
  

  subroutine qrm_get_mem_stats(totmem, maxmem)
    ! Function: qrm_get_mem_stats
    ! This subroutine returns stats on memory.
    ! 
    ! *Output*:
    ! totmem - the amount of memory currently allocated
    ! maxmem - the peak memory allocated
    
    integer(kind=8) :: totmem, maxmem

    totmem = qrm_tot_mem(0)
    maxmem = qrm_max_mem(0)
  
    return
  
  end subroutine qrm_get_mem_stats
  


end module qrm_mem_mod
