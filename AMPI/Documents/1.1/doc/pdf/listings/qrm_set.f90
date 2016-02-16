interface qrm_set

   subroutine sqrm_pseti(qrm_mat, string, ival)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer               :: ival
   end subroutine sqrm_pseti

   subroutine qrm_gseti(string, ival)
     character(len=*)      :: string
     integer               :: ival
   end subroutine qrm_gseti

end interface qrm_set
