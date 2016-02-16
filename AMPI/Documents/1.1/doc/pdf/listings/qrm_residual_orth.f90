interface qrm_residual_orth

   subroutine sqrm_residual_orth1d(qrm_mat, r, nrm)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: r(:)
     real                  :: nrm
   end subroutine sqrm_residual_orth1d

   subroutine sqrm_residual_orth2d(qrm_mat, r, nrm)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: r(:,:)
     real                  :: nrm
   end subroutine sqrm_residual_orth2d

end interface qrm_residual_orth
