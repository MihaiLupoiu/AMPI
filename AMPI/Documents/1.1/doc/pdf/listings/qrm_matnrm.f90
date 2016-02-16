interface qrm_matnrm

   subroutine sqrm_matnrm(qrm_mat, ntype, nrm)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: nrm
     character             :: ntype
   end subroutine sqrm_matnrm

end interface qrm_matnrm
