interface qrm_analyse

   subroutine sqrm_analyse(qrm_mat, transp)
     type(sqrm_spmat_type):: qrm_mat
     character, optional  :: transp
   end subroutine sqrm_analyse

end interface qrm_analyse
