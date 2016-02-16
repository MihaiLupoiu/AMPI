module qrm_const_mod
  ! The maximum number of tasks allowed
  integer, parameter :: qrm_maxthreads=256

  !> @brief The unit for output messages
  integer, save :: qrm_ounit=6
  !> @brief The unit for error messages
  integer, save :: qrm_eunit=0
  !> @brief The unit for debug messages
  integer, save :: qrm_dunit=-1


  !> @brief The number of entries in icntl
  integer, parameter :: nicntl=20

 !> @brief The number of entries in rcntl
  integer, parameter :: nrcntl=10

  !> @brief constants to set control parameters
  integer, parameter :: qrm_auto_           = 0
  integer, parameter :: qrm_ordering_       = 1
  integer, parameter :: qrm_sing_           = 2
  integer, parameter :: qrm_minamalg_       = 3
  integer, parameter :: qrm_nb_             = 4
  integer, parameter :: qrm_keeph_          = 5
  integer, parameter :: qrm_ib_             = 6
  integer, parameter :: qrm_rhsnb_          = 7
  integer, parameter :: qrm_nthreads_       = 8
  integer, parameter :: qrm_rhsnthreads_    = 9
  integer, parameter :: qrm_nlz_            = nicntl
  integer, parameter :: qrm_cnode_          = nicntl-1
                                            
  integer, parameter :: qrm_amalgth_        = 1
  integer, parameter :: qrm_rweight_        = nrcntl
                                            
                                            
  integer, parameter :: qrm_natural_        = 1
  integer, parameter :: qrm_given_          = 2
  integer, parameter :: qrm_colamd_         = 3
  integer, parameter :: qrm_metis_          = 4
  integer, parameter :: qrm_scotch_         = 5

  integer, parameter :: qrm_e_facto_flops_  = 1
  integer, parameter :: qrm_e_nnz_r_        = 2
  integer, parameter :: qrm_e_nnz_h_        = 3
  integer, parameter :: qrm_facto_flops_    = 4
  integer, parameter :: qrm_nnz_r_          = 5
  integer, parameter :: qrm_nnz_h_          = 6

  integer, parameter :: qrm_no_             = 0
  integer, parameter :: qrm_yes_            = 1
                                           
                                           
  integer, parameter :: qrm_allop_          = 0
  integer, parameter :: qrm_analyse_        = 1
  integer, parameter :: qrm_factorize_      = 2
  integer, parameter :: qrm_solve_          = 3
  integer, parameter :: qrm_apply_          = 4


end module qrm_const_mod
