interface qrm_solve

   subroutine sqrm_solve1d(qrm_mat, transp, b, x)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:)
     real                  :: x(:)
     character             :: transp
   end subroutine sqrm_solve1d

   subroutine sqrm_solve2d(qrm_mat, transp, b, x)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:,:)
     real                  :: x(:,:)
     character             :: transp
   end subroutine sqrm_solve2d

end interface qrm_solve
