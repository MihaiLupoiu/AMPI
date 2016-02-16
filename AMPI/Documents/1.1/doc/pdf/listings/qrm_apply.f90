interface qrm_apply
   
   subroutine sqrm_apply1d(qrm_mat, transp, b)
     type(sqrm_spmat_type) :: qrm_mat
     character             :: transp
     real                  :: b(:)
   end subroutine sqrm_apply1d
   
   subroutine sqrm_apply2d(qrm_mat, transp, b)
     type(sqrm_spmat_type) :: qrm_mat
     character             :: transp
     real                  :: b(:,:)
   end subroutine sqrm_apply2d
   
end interface qrm_apply
