interface qrm_factorize

   subroutine sqrm_factorize(qrm_mat, transp)
     type(sqrm_spmat_type):: qrm_mat
     character, optional  :: transp
   end subroutine sqrm_factorize

end interface qrm_factorize
