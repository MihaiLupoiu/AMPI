interface qrm_get

   subroutine sqrm_pgeti(qrm_mat, string, ival)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer               :: ival
   end subroutine sqrm_pgeti
   
   subroutine sqrm_pgetii(qrm_mat, string, ival)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer(kind=8)       :: ival
   end subroutine sqrm_pgetii
   
   subroutine qrm_ggeti(string, ival)
     character(len=*)      :: string
     integer               :: ival
   end subroutine qrm_ggeti

   subroutine qrm_ggetii(string, ival)
     character(len=*)      :: string
     integer(kind=8)       :: ival
   end subroutine qrm_ggetii

end interface qrm_get
