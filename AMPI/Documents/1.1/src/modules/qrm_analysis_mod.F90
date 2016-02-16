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
!> @file qrm_analysis_mod.F90
!! This file contains a module with all the generic interfaces for the typed analysis routines.
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains the generic interfaces for all the analysis routines.
module _qrm_analysis_mod

  !> @brief Generic interface for the @link ::_qrm_analyse @endlink routine
  interface qrm_analyse
     subroutine _qrm_analyse(qrm_mat, transp)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)          :: qrm_mat
       character, optional, intent(in)  :: transp
     end subroutine _qrm_analyse
  end interface

  !> @brief Generic interface for the @link ::_qrm_compute_graph @endlink routine
  interface qrm_compute_graph
     subroutine _qrm_compute_graph(qrm_mat, graph, work)
       use _qrm_spmat_mod
       type(_qrm_spmat_type)              :: qrm_mat
       type(_qrm_spmat_type), intent(out) :: graph
       integer, target, optional            :: work(:)
     end subroutine _qrm_compute_graph
  end interface

  !> @brief Generic interface for the @link ::_qrm_detect_singletons @endlink routine
  interface qrm_detect_singletons
     subroutine _qrm_detect_singletons(graph, scol, srow, mrperm, mcperm, nrsing, ncsing)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: graph
       integer, allocatable                :: scol(:), srow(:)
       integer, allocatable                :: mcperm(:), mrperm(:)
       integer                             :: nrsing, ncsing
     end subroutine _qrm_detect_singletons
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_colamd @endlink routine
  interface qrm_do_colamd
     subroutine _qrm_do_colamd(graph, cperm)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: graph
       integer                 :: cperm(:)
     end subroutine _qrm_do_colamd
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_metis @endlink routine
  interface qrm_do_metis
     subroutine _qrm_do_metis(graph, cperm)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: graph
       integer                 :: cperm(:)
     end subroutine _qrm_do_metis
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_scotch @endlink routine
  interface qrm_do_scotch
     subroutine _qrm_do_scotch(graph, cperm)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: graph
       integer                 :: cperm(:)
     end subroutine _qrm_do_scotch
  end interface

  !> @brief Generic interface for the @link ::_qrm_ata_graph @endlink routine
  interface qrm_ata_graph
     subroutine _qrm_ata_graph(g_csc, ata_graph)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in)  :: g_csc
       type(_qrm_spmat_type), intent(out) :: ata_graph
     end subroutine _qrm_ata_graph
  end interface

  !> @brief Generic interface for the @link ::_qrm_do_ordering @endlink routine
  interface qrm_do_ordering
     subroutine _qrm_do_ordering(graph, cperm, cperm_in)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: graph
       integer                 :: cperm(:)
       integer, pointer        :: cperm_in(:)
     end subroutine _qrm_do_ordering
  end interface

  !> @brief Generic interface for the @link ::_qrm_elim_tree @endlink routine
  interface qrm_elim_tree
     subroutine _qrm_elim_tree(graph, cperm, parent)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(in) :: graph
       integer, intent(in)                 :: cperm(:)
       integer, allocatable                :: parent(:)
     end subroutine _qrm_elim_tree
  end interface

  !> @brief Generic interface for the @link ::_qrm_rowcount @endlink routine
  interface qrm_rowcount
     subroutine _qrm_rowcount(graph, parent, porder, rc)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: graph
       integer                 :: parent(:), porder(:), rc(:)
     end subroutine _qrm_rowcount
  end interface

  !> @brief Generic interface for the @link ::_qrm_symbolic @endlink routine
  interface qrm_symbolic
     subroutine _qrm_symbolic(graph)
       use _qrm_spmat_mod
       type(_qrm_spmat_type), intent(inout) :: graph
     end subroutine _qrm_symbolic
  end interface
  
  !> @brief Generic interface for the @link ::_qrm_rowperm @endlink routine
  interface qrm_rowperm
     subroutine _qrm_rowperm(qrm_mat, cperm, rperm, nvar, stair)
       use _qrm_spmat_mod
       type(_qrm_spmat_type) :: qrm_mat
       integer                 :: cperm(:), rperm(:), nvar(:), stair(:)
     end subroutine _qrm_rowperm
  end interface

  !> @brief Generic interface for the @link ::_qrm_attach_singletons @endlink routine
  interface qrm_attach_singletons
     subroutine _qrm_attach_singletons(qrm_mat, scol, srow, mrperm, mcperm, nrsing, ncsing)
       use _qrm_spmat_mod
       type(_qrm_spmat_type):: qrm_mat
       integer                :: scol(:), srow(:), mcperm(:), mrperm(:)
       integer                :: nrsing, ncsing
     end subroutine _qrm_attach_singletons
  end interface

end module _qrm_analysis_mod
   
