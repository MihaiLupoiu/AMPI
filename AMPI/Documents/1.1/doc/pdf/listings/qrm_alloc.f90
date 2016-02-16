interface qrm_aalloc

   subroutine qrm_aalloc_s(a, m)
     real(kind(1.e0)), allocatable  :: a(:)
     integer                        :: m
   end subroutine qrm_aalloc_s
   
   subroutine qrm_aalloc_2s(a, m, n)
     real(kind(1.e0)), allocatable  :: a(:,:)
     integer                        :: m, n
   end subroutine qrm_aalloc_2s
   
   subroutine qrm_adealloc_s(a)
     real(kind(1.e0)), allocatable  :: a(:)
   end subroutine qrm_adealloc_s
   
   subroutine qrm_adealloc_2s(a)
     real(kind(1.e0)), allocatable  :: a(:,:)
   end subroutine qrm_adealloc_2s
   
end interface qrm_aalloc
