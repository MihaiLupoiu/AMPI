
  type(sqrm_spmat_type) :: qrm_mat
  real, allocatable     :: b(:), x(:)

  ! initialize the control data structure. 
  call qrm_spmat_init(qrm_mat)
  ...
  ! allocate arrays for the input matrix
  call qrm_palloc(qrm_mat%irn, nz)
  call qrm_palloc(qrm_mat%jcn, nz)
  call qrm_palloc(qrm_mat%val, nz)
  call qrm_aalloc(b, m)
  call qrm_aalloc(x, n)

  ! initialize the data
  ...

  ! solve the problem 
  call qrm_least_squares(qrm_mat, b, x)
  ...
