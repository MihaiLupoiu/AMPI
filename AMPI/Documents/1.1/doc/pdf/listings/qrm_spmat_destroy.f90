interface qrm_spmat_destroy

   subroutine sqrm_spmat_destroy(qrm_mat, all)
     type(sqrm_spmat_type) :: qrm_mat
     logical, optional     :: all
   end subroutine sqrm_spmat_destroy

end interface qrm_spmat_destroy
