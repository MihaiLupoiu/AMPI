interface qrm_min_norm
   subroutine sqrm_min_norm1d(qrm_mat, b, x)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: x(:)
     real                  :: b(:)
   end subroutine sqrm_min_norm1d

   subroutine sqrm_min_norm2d(qrm_mat, b, x)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: x(:,:)
     real                  :: b(:,:)
   end subroutine sqrm_min_norm2d

end interface qrm_min_norm
