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
!> @file qrm_adata_mod.F90
!! This file contains the module that holds the main data type for the analysis phase.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"


!> @brief This module contains the definition of the analysis data type.



module qrm_adata_mod

!> @brief The main data type for the analysis phase.

!> This data type is meant to store all the results of the anaysis
!! phase (permutations, assembly tree etc.). Moreover it is
!! not typed which means that it will be exactly the same for dp/sp
!! real or complex matrices. However, because of its untyped nature,
!! this data type cannot be included into qrm_analysis_mod which, instead
!! is typed due to the fact that many routines there take as input
!! argument a qrm_spmat_type.

  type qrm_adata_type
     !> An integer array of size n (column dimension of the problem matrix)
     !! containing the column permutation resulting from fill-in
     !! minimization orderings (COLAMD etc). It has to be a pointer to
     !! accomodate the case where the user wants to provide his own
     !! permutation.
     integer, allocatable            :: cperm(:)
     !> Inverse column permutation. icperm(i)=k means that column i is eliminated
     !! at k-th step
     integer, allocatable            :: icperm(:)
     !> This is an integer array containing a row permutation computed in order
     !! to have the original matrix in a "stair" format. See @link _qrm_rowperm_ @endlink
     integer, allocatable            :: rperm(:)
     !> This integer array holds pointers to the beginning of nodes in the
     !! column permutation. The variables contained in node i are
     !! cperm( cp_ptr(i) : cp_ptr(i+1)-1 ). It is of size nnodes+1
     integer, allocatable            :: cp_ptr(:)
     !> rc(i) contains the rowcount for the rows in node number i. This array
     !! is of size nnodes
     integer, allocatable            :: rc(:)
     !> An integer array of size nnodes containing the assembly tree.
     !! parent(i)=k means that node k is the father of node i
     integer, allocatable            :: parent(:)
     !> For each node, it contains the list of its children
     integer, allocatable            :: child(:)
     !> Pointers to the list of children. child(childptr(i),...,childptr(i+1)-1)
     !! contains all the children of node i
     integer, allocatable            :: childptr(:)
     !> nfrows(i)=k means that the frontal matrix related to node i has k rows
     !! (the number of columns is rc(i))
     integer, allocatable            :: nfrows(:)
     !> defines the staircase structure of the original matrix. Assuming A(rperm,cperm)
     !! @verbatim
     !!            |x       |
     !!            |x       |     in this case stair would be:
     !!            | x      |     stair=(/ 2, 5, 7, 7, 9 /)
     !!            | x      |
     !!            | x      |
     !!            |  xx    |
     !!            |  xx    |
     !!            |    x   |
     !!            |    x   |
     !! @endverbatim         
     integer, allocatable            :: stair(:)
     !> An array of size nnodes which flags (with 1 whereas everything is =0) the roots
     !! of subtrees that are treated sequentially in the numerical factorization
     integer, allocatable            :: small(:) 
     !> Contains the list of columns indices for each front
     integer, allocatable            :: fcol(:) 
     !> Contains pointers to the list of column indices fcol. Specifically, the list of
     !! column indices for front i is fcol(fcol_ptr(i):fcol_ptr(i+1)-1)
     integer, allocatable            :: fcol_ptr(:)
     !> A list of leaf nodes where to start the numerical factorization from. A leaf is
     !! defined as a node which only has small children
     integer, allocatable            :: leaves(:) 
     !> The number of leaves present in the tree
     integer                         :: nleaves=0
     !> The number of nodes in the elimination tree, i.e., the number of frontal
     !! matrices.
     integer                         :: nnodes=0
     !> The number of column singletons found
     integer                         :: ncsing=0 
     !> The number of row singletons found
     integer                         :: nrsing=0 
     !> it is set to .true. if the analysis is done (with success), to .false. otherwise
     logical                         :: ok=.false.
  end type qrm_adata_type

  contains


    !> @brief make s a copy of an @link qrm_adata_type @endlink instance
    !! @param[in]  adata_in  the instance to be copied
    !! @param[out] adata_out the produced copy
    subroutine qrm_adata_copy(adata_in, adata_out)

      use qrm_mem_mod
      use qrm_error_mod
      implicit none

      type(qrm_adata_type), intent(in)  :: adata_in
      type(qrm_adata_type), intent(out) :: adata_out

      ! error management
      integer                         :: err_act
      character(len=*), parameter     :: name='qrm_adata_copy'
      
      call qrm_err_act_save(err_act)
      
      call qrm_aalloc(adata_out%cperm,    qrm_asize(adata_in%cperm))
      call qrm_aalloc(adata_out%icperm,   qrm_asize(adata_in%icperm))
      call qrm_aalloc(adata_out%rperm,    qrm_asize(adata_in%rperm))
      call qrm_aalloc(adata_out%cp_ptr,   qrm_asize(adata_in%cp_ptr))
      call qrm_aalloc(adata_out%rc,       qrm_asize(adata_in%rc))
      call qrm_aalloc(adata_out%parent,   qrm_asize(adata_in%parent))
      call qrm_aalloc(adata_out%child,    qrm_asize(adata_in%child))
      call qrm_aalloc(adata_out%childptr, qrm_asize(adata_in%childptr))
      call qrm_aalloc(adata_out%nfrows,   qrm_asize(adata_in%nfrows))
      call qrm_aalloc(adata_out%stair,    qrm_asize(adata_in%stair))
      call qrm_aalloc(adata_out%leaves,   qrm_asize(adata_in%leaves))
      call qrm_aalloc(adata_out%fcol,     qrm_asize(adata_in%fcol))
      call qrm_aalloc(adata_out%fcol_ptr, qrm_asize(adata_in%fcol_ptr))
      call qrm_aalloc(adata_out%small,    qrm_asize(adata_in%small))
      __QRM_CHECK_RET(name,'qrm_alloc',9999)

      if(allocated(adata_in%cperm))adata_out%cperm       = adata_in%cperm 
      if(allocated(adata_in%icperm))adata_out%icperm     = adata_in%icperm
      if(allocated(adata_in%rperm))adata_out%rperm       = adata_in%rperm
      if(allocated(adata_in%cp_ptr))adata_out%cp_ptr     = adata_in%cp_ptr
      if(allocated(adata_in%rc))adata_out%rc             = adata_in%rc
      if(allocated(adata_in%parent))adata_out%parent     = adata_in%parent
      if(allocated(adata_in%child))adata_out%child       = adata_in%child
      if(allocated(adata_in%childptr))adata_out%childptr = adata_in%childptr
      if(allocated(adata_in%nfrows))adata_out%nfrows     = adata_in%nfrows
      if(allocated(adata_in%stair))adata_out%stair       = adata_in%stair
      if(allocated(adata_in%leaves))adata_out%leaves     = adata_in%leaves
      if(allocated(adata_in%fcol))adata_out%fcol         = adata_in%fcol
      if(allocated(adata_in%fcol_ptr))adata_out%fcol_ptr = adata_in%fcol_ptr
      if(allocated(adata_in%small))adata_out%small       = adata_in%small
      adata_out%nnodes     = adata_in%nnodes
      adata_out%ncsing     = adata_in%ncsing
      adata_out%nrsing     = adata_in%nrsing
      adata_out%nleaves    = adata_in%nleaves
      adata_out%ok         = adata_in%ok

      call qrm_err_act_restore(err_act)
      return
      
9999  continue ! error management
      call qrm_err_act_restore(err_act)
      if(err_act .eq. qrm_abort_) then
         call qrm_err_check()
      end if
      return
      
    end subroutine qrm_adata_copy


    !> @brief make s a move of an @link qrm_adata_type @endlink instance
    !! @param[in]  adata_in  the instance to be moved
    !! @param[out] adata_out the produced copy
    subroutine qrm_adata_move(adata_in, adata_out)

      use qrm_mem_mod
      use qrm_error_mod
      implicit none

      type(qrm_adata_type), intent(inout)  :: adata_in
      type(qrm_adata_type), intent(inout) :: adata_out

      ! error management
      integer                         :: err_act
      character(len=*), parameter     :: name='qrm_adata_move'
      
      call qrm_err_act_save(err_act)
      
      call move_alloc(adata_in%cperm   ,  adata_out%cperm     ) 
      call move_alloc(adata_in%icperm  ,  adata_out%icperm    ) 
      call move_alloc(adata_in%rperm   ,  adata_out%rperm     ) 
      call move_alloc(adata_in%cp_ptr  ,  adata_out%cp_ptr    ) 
      call move_alloc(adata_in%rc      ,  adata_out%rc        ) 
      call move_alloc(adata_in%parent  ,  adata_out%parent    ) 
      call move_alloc(adata_in%child   ,  adata_out%child     ) 
      call move_alloc(adata_in%childptr,  adata_out%childptr  ) 
      call move_alloc(adata_in%nfrows  ,  adata_out%nfrows    ) 
      call move_alloc(adata_in%stair   ,  adata_out%stair     ) 
      call move_alloc(adata_in%leaves  ,  adata_out%leaves    ) 
      call move_alloc(adata_in%fcol    ,  adata_out%fcol      ) 
      call move_alloc(adata_in%fcol_ptr,  adata_out%fcol_ptr  ) 
      call move_alloc(adata_in%small   ,  adata_out%small     ) 

      adata_out%nnodes     = adata_in%nnodes
      adata_out%ncsing     = adata_in%ncsing
      adata_out%nrsing     = adata_in%nrsing
      adata_out%nleaves    = adata_in%nleaves
      adata_out%ok         = adata_in%ok

      call qrm_err_act_restore(err_act)
      return
      
9999  continue ! error management
      call qrm_err_act_restore(err_act)
      if(err_act .eq. qrm_abort_) then
         call qrm_err_check()
      end if
      return
      
    end subroutine qrm_adata_move


    !> @brief Frees an @link qrm_adata_type @endlink instance
    !! @param[in,out] adata the instance to be freed
    subroutine qrm_adata_destroy(adata)

      use qrm_mem_mod
      use qrm_error_mod
      implicit none

      type(qrm_adata_type) :: adata

      ! error management
      integer                         :: err_act
      character(len=*), parameter     :: name='qrm_adata_destroy'
      
      call qrm_err_act_save(err_act)

      call qrm_adealloc(adata%cperm)
      call qrm_adealloc(adata%icperm)
      call qrm_adealloc(adata%rperm)
      call qrm_adealloc(adata%cp_ptr)
      call qrm_adealloc(adata%rc)
      call qrm_adealloc(adata%parent)
      call qrm_adealloc(adata%nfrows)
      call qrm_adealloc(adata%stair)
      call qrm_adealloc(adata%small)
      call qrm_adealloc(adata%childptr)
      call qrm_adealloc(adata%child)
      call qrm_adealloc(adata%leaves)
      call qrm_adealloc(adata%fcol)
      call qrm_adealloc(adata%fcol_ptr)
      __QRM_CHECK_RET(name,'qrm_dealloc',9999)
  
      adata%nleaves = 0
      adata%nnodes  = 0
      adata%ncsing  = 0
      adata%nrsing  = 0

      adata%ok = .false.
      call qrm_err_act_restore(err_act)
      return
      
9999  continue ! error management
      call qrm_err_act_restore(err_act)
      if(err_act .eq. qrm_abort_) then
         call qrm_err_check()
      end if
      return
      
    end subroutine qrm_adata_destroy

  
end module qrm_adata_mod
