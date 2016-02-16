!! ##############################################################################################
!!
!! Copyright 2012 CNRS, INPT
!!  
!! This file is part of qr_mumps.
!!  
!! qr_mumps is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as 
!! published by the Free Software Foundation, either version 3 of 
!! the License, or (at your option) any later version.
!!  
!! qr_mumps is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!  
!! You can find a copy of the GNU Lesser General Public License
!! in the qr_mumps/doc directory.
!!
!! ##############################################################################################


!! ##############################################################################################
!> @file qrm_rfpf_mod.F90
!! This file contains a module that implements operations on RFPF matrices
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


#include "qrm_common.h"

!> This module contains an implementation of some operations
!! on triangular/trapezoidal matrices stored in Rectangular Full Packed Format.

!> For more details on this format please refer to:
!!
!! Fred G. Gustavson, Jerzy Wasniewski, and Jack J. Dongarra.
!! <em>"Rectangular Full Packed Format for Cholesky’s Algorithm: Factorization, Solution and Inversion."</em>
!!  LAPACK Working Note 199, April 2008.
!!
!! According to this format, a triangular matrix can be stored in a 2D array like described below.
!!
!! The RFPF format is obtained by "cutting" the upper-leftmost corner of a triangular matrix,
!! transposing it and the appending it on the right or at the bottom for lower or upper
!! triangular matrices, respectively.
!! 
!! A lower triangular matrix with odd dimension:
!!
!! @verbatim
!! 01
!! 02 03             
!! 04 05 06           is stored like    04 05 06 01 02
!! 07 08 09 10                          07 08 09 10 03
!! 11 12 13 14 15                       11 12 13 14 15
!! @endverbatim
!!
!! A lower triangular matrix with even dimension:
!!
!! @verbatim
!! 01
!! 02 03
!! 04 05 06
!! 07 08 09 10        is stored like    07 08 09 10 01 02 04
!! 11 12 13 14 15                       11 12 13 14 15 03 05
!! 16 17 18 19 20 21                    16 17 18 19 20 21 06
!! @endverbatim
!!
!! An upper triangular matrix with odd dimension
!!
!! @verbatim
!! 01 02 03 04 05                       03 04 05
!!    06 07 08 09                       07 08 09
!!       10 11 12     is stored like    10 11 12
!!          13 14                       01 13 14
!!             15                       02 06 15
!! @endverbatim
!!
!! An upper triangular matrix with even dimension
!!
!! @verbatim
!! 01 02 03 04 05 06                    04 05 06
!!    07 08 09 10 11                    09 10 11
!!       12 13 14 15                    13 14 15
!!          16 17 18  is stored like    16 17 18
!!             19 20                    01 19 20
!!                21                    02 07 21
!! @endverbatim
!!                                      03 08 12
!!
!! In general a matrix in RFPF can be split into three components
!! (e.g., in the last even, upper triangular case):
!!
!! @verbatim
!!         04 05 06
!!         09 10 11
!!         13 14 15
!!
!!           16 17 18
!!     01       19 20
!!     02 07       21
!!     03 08 12
!! @endverbatim
!!
!! and thus an operation like a TRMM (triangular matrix times a dense matrix)
!! can be translated into three sub-operations: TRMM, GEMM and TRMM.
!!
!! A general, trapezoidal matrix can be split into a rectangular block, stored in
!! conventional column-major format and a triangular one stored in RFPF.
!!
!! Routines in this module will have an old-style f77 interface a la LAPACK.
!! This is risky because the fortran standard doesn't say a word about storing
!! arrays/allocatables in contiguous memory areas although compilers always do it.
!! However, in such a case all the BLAS and LAPACK routines won't work so we're
!! ruined anyway.
!!
module _qrm_rfpf_mod

  use qrm_error_mod

#if defined (sprec) || defined (dprec)
  character, parameter :: tran='t', notran='n'
#elif defined (cprec) || defined (zprec)
  character, parameter :: tran='c', notran='n'
#endif


contains

  subroutine _qrm_to_rfpf(uplo, unit, n, a, lda, b)

    use qrm_string_mod
    implicit none

    character :: uplo, unit
    integer   :: n, lda
    _qrm_data :: a(lda,*)
    _qrm_data :: b(*)
    integer :: mb, nb, row, col, i, j
    logical :: lo, ud, up
    
    lo = qrm_str_tolower(uplo) .eq. 'l'
    up = qrm_str_tolower(uplo) .eq. 'u'
    ud = qrm_str_tolower(unit) .eq. 'u'

    if(lo) then
       mb = (n+1)/2
       nb = n+(n/2-mb+1)

       do j=1, nb
          row = j-nb+n/2 
          do i=1, j-nb+n/2
             col = i
             if(row .eq. col .and. ud) then
                b((j-1)*mb+i) = _qrm_one
             else
                b((j-1)*mb+i) = _conjg(a(row,col))
             end if
          end do
          col = j
          do i= max(j-nb+n/2+1,1), mb
             row = i+n/2
             if(row .eq. col .and. ud) then
                b((j-1)*mb+i) = _qrm_one
             else
                b((j-1)*mb+i) = a(row,col)
             end if
          end do
       end do

    else if (up) then
       nb  = (n+1)/2
       mb  = n+(n/2-nb+1)

       do j=1, nb
          col = j+n/2
          do i=1, mb-(n/2-j+1)
             row = i
             if(row .eq. col .and. ud) then
                b((j-1)*mb+i) = _qrm_one
             else
                b((j-1)*mb+i) = a(row,col)
             end if
          end do
          row = j
          do i= max(mb-(n/2-j+1)+1,1), mb
             col = i - (mb-n/2)
             if(row .eq. col .and. ud) then
                b((j-1)*mb+i) = _qrm_one
             else
                b((j-1)*mb+i) = _conjg(a(row,col))
             end if
          end do
       end do

    end if

    return

  end subroutine _qrm_to_rfpf



  subroutine _xrpmv(uplo, trans, diag, n, a, x)

    implicit none
    character :: uplo, trans, diag
    integer   :: n
    _qrm_data :: a(*), x(*)

    integer :: ma, na
    character :: nuplo, ntrans

    if(uplo .eq. 'l') then
       ma = (n+1)/2
       na = n+(n/2-ma+1)

       if((trans .eq. 't') .or. (trans .eq. 'c')) then

          call _xtrmv('u', notran, diag, n/2, a((na-n/2)*ma+1), ma, x(1), 1)
          call _xgemv(tran, ma, na-n/2-1, _qrm_one, a(1), ma, &
               & x(n/2+1), 1, _qrm_one, x(1), 1)
          call _xtrmv('l', tran, diag, n-n/2, a((na-n/2-1)*ma+1), ma, x(n/2+1), 1)

       else

          call _xtrmv('l', notran, diag, n-n/2, a((na-n/2-1)*ma+1), ma, x(n/2+1), 1)
          call _xgemv(notran, ma, na-n/2-1, _qrm_one, a(1), ma, &
               & x(1), 1, _qrm_one, x(n/2+1), 1)
          call _xtrmv('u', tran, diag, n/2, a((na-n/2)*ma+1), ma, x(1), 1)

       end if

    else if (uplo .eq. 'u') then
       na  = (n+1)/2
       ma  = n+(n/2-na+1)

       if((trans .eq. 't') .or. (trans .eq. 'c')) then

          call _xtrmv('u', tran, diag, na, a((ma-n/2-1)+1), ma, x(n/2+1), 1)
          call _xgemv(tran, ma-n/2-1, na, _qrm_one, a(1), ma, &
               & x(1), 1, _qrm_one, x(n/2+1), 1)
          call _xtrmv('l', notran, diag, n/2, a((ma-n/2)+1), ma, x(1), 1)

       else

          call _xtrmv('l', tran, diag, n/2, a((ma-n/2)+1), ma, x(1), 1)
          call _xgemv(notran, ma-n/2-1, na, _qrm_one, a(1), ma, &
               & x(n/2+1), 1, _qrm_one, x(1), 1)
          call _xtrmv('u', notran, diag, na, a((ma-n/2-1)+1), ma, x(n/2+1), 1)

       end if

    end if

    return

  end subroutine _xrpmv





  subroutine _xrpsv(uplo, trans, diag, n, a, x)

    implicit none
    character :: uplo, trans, diag
    integer   :: n
    _qrm_data :: a(*), x(*)

    integer :: ma, na
    character :: nuplo, ntrans

    if(uplo .eq. 'l') then
       ma = (n+1)/2
       na = n+(n/2-ma+1)

       if((trans .eq. 't') .or. (trans .eq. 'c')) then

          call _xtrsv('l', tran, diag, n-n/2, a((na-n/2-1)*ma+1), ma, x(n/2+1), 1)
          call _xgemv(tran, ma, na-n/2-1, -_qrm_one, a(1), ma, &
               & x(n/2+1), 1, _qrm_one, x(1), 1)
          call _xtrsv('u', notran, diag, n/2, a((na-n/2)*ma+1), ma, x(1), 1)

       else

          call _xtrsv('u', tran, diag, n/2, a((na-n/2)*ma+1), ma, x(1), 1)
          call _xgemv(notran, ma, na-n/2-1, -_qrm_one, a(1), ma, &
               & x(1), 1, _qrm_one, x(n/2+1), 1)
          call _xtrsv('l', notran, diag, n-n/2, a((na-n/2-1)*ma+1), ma, x(n/2+1), 1)

       end if

    else if (uplo .eq. 'u') then
       na  = (n+1)/2
       ma  = n+(n/2-na+1)

       if((trans .eq. tran) .or. (trans .eq. 'c')) then

          call _xtrsv('l', notran, diag, n/2, a((ma-n/2)+1), ma, x(1), 1)
          call _xgemv(tran, ma-n/2-1, na, -_qrm_one, a(1), ma, &
               & x(1), 1, _qrm_one, x(n/2+1), 1)
          call _xtrsv('u', tran, diag, na, a((ma-n/2-1)+1), ma, x(n/2+1), 1)

       else

          call _xtrsv('u', notran, diag, na, a((ma-n/2-1)+1), ma, x(n/2+1), 1)
          call _xgemv(notran, ma-n/2-1, na, -_qrm_one, a(1), ma, &
               & x(n/2+1), 1, _qrm_one, x(1), 1)
          call _xtrsv('l', tran, diag, n/2, a((ma-n/2)+1), ma, x(1), 1)

       end if

    end if

    return

  end subroutine _xrpsv







  subroutine _xrprft(direct, storev, m, k, v1, v2, ldv2, tau, t, ldt)
    ! This routine is different from LAPACK's one because the V matrix
    ! already has ones on the diagonal.
    use qrm_const_mod
    implicit none
    integer :: m, k, ldv2, ldt
    character :: direct, storev
    _qrm_data :: v1(*), v2(ldv2,*), tau(*), t(ldt,*)
    _qrm_data, parameter :: zero=_qrm_zero, one=_qrm_one
    _qrm_data :: vii, mult
    logical, external :: lsame
    integer :: mv1, nv1, i, j, off1, off2, im, in, ii, jj
    character :: tr

    if( k.eq.0 ) return

    mv1 = (k+1)/2
    nv1 = k+(k/2-mv1+1)

    if( lsame( direct, 'f' ) ) then
       do  i = 1, k
          if( tau( i ).eq. zero ) then
             ! 
             ! h(i)  =  i
             ! 
             do  j = 1, i
                t( j, i ) = zero
             end do
          else
             if(i .le. k/2) then
                mult = _qrm_one
             else
                mult = _qrm_zero
             end if

             if( lsame( storev, 'c' ) ) then
                off1 = (k/2+i)*mv1+1
                off2 = (k/2+i)*mv1+i
                im   = i-1
                in   = k/2-i+1


                ! FIXME: check this (seems to work). last good revision 293 (non tread safe)
#if defined (cprec) || defined (zprec)
                if(min(im,in) .gt. 0) then
                   ! using the lower part ot T as a scratch memory
                   call _xcopy( in, v1(off2), mv1, t(ldt-in+1,1), 1 )
                   call _xlacgv( in, t(ldt-in+1,1), 1)
                   call _xgemv(notran, im, in, -tau(i), &
                        & v1(off1), mv1, t(ldt-in+1,1), 1, zero, t(1,i),1)
                end if
#else
                if(min(im,in) .gt. 0) call _xgemv(notran, im, in, -tau(i), &
                     & v1(off1), mv1, v1(off2), mv1, zero, t(1,i),1)
#endif

                off1 = max(i-k/2,1)
                off2 = max((i-1)*mv1 + max(i-k/2,1),1)
                im   = min(mv1,mv1-(i-k/2)+1)
                in   = i-1


                if(min(im,in) .gt. 0) call _xgemv(tran, im, in, -tau(i), &
                     & v1(off1), mv1, v1(off2), 1, mult, t(1,i), 1)

                if(min(m,i-1) .gt. 0) call _xgemv(tran, m, i-1, -tau(i), &
                     & v2(1,1), ldv2, v2(1,i), 1, one, t(1,i), 1)

             else 
                __QRM_PRNT_ERR('("_RPRFT: method not supported")')
                return
             end if


             call _xtrmv( 'upper', notran, 'non-unit', i-1, t, ldt, t( 1, i ), 1 )

             t( i, i ) = tau( i )

             
          end if
       end do
    else 
       __QRM_PRNT_ERR('("_RPRFT: method not supported")')
       return
    end if
    return

  end subroutine _xrprft



  subroutine _xrprfb(side, trans, direct, storev, m, n, k, v1, v2, ldv2, t, ldt, c, ldc, work, ldwork)
    use qrm_const_mod
    implicit none

    character :: side, trans, direct, storev
    integer   :: m, n, k, ldv2, ldt, ldc, ldwork, ii
    _qrm_data :: v1(*), v2(ldv2,*), t(ldt,*), c(ldc,*), work(ldwork, *)

    integer :: i, j
    _qrm_data, parameter :: one=_qrm_one
    logical, external :: lsame
    character :: transt

    if( m.le.0 .or. n.le.0 )  return

    if( lsame( trans, 'n' ) ) then
#if defined (sprec) || defined (dprec)
       transt = 't'
#elif defined (cprec) || defined (zprec)
       transt = 'c'
#endif
    else
       transt = 'n'
    end if

    if( lsame( storev, 'c' ) ) then
       if( lsame( direct, 'f' ) ) then
          !           let  v =  ( v1 )    (first k rows)
          !                     ( v2 )
          !           where  v1  is unit lower triangular.
          if( lsame( side, 'l' ) ) then
             !              form  h * c  or  h' * c  where  c = ( c1 )
             !                                                  ( c2 )
             !
             !              w := c' * v  =  (c1'*v1 + c2'*v2)  (stored in work)
             !              w := c1'
             do j = 1, k
                call _xcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
#if defined(cprec) || defined(zprec)
                call _xlacgv(n, work(1,j), 1)
#endif
             end do
             !              w := w * v1
             call _xrpmm( 'right', 'lower', notran, 'unit', n, k, one, v1, work, ldwork )
             if( m.gt.k ) then
                !                 w := w + c2'*v2
                call _xgemm( tran, notran, n, k, m-k, one, c( k+1, 1 ), ldc, &
                     & v2( 1, 1 ), ldv2, one, work, ldwork )
             end if
             !              w := w * t'  or  w * t
             !FIXME not thread-safe (for some reason)
             call _xtrmm( 'right', 'upper', transt, 'non-unit', n, k, one, t, ldt, work, ldwork )
             !              c := c - v * w'
             if( m.gt.k ) then
                !                 c2 := c2 - v2 * w'
                call _xgemm( notran, tran, m-k, n, k, -one, v2( 1, 1 ), ldv2, work, &
                     & ldwork, one, c( k+1, 1 ), ldc )
             end if
             !              w := w * v1'
             call _xrpmm( 'right', 'lower', tran, 'unit', n, k, one, v1, work, ldwork )

             !              c1 := c1 - w'
             do  j = 1, k
                do  i = 1, n
                   c( j, i ) = c( j, i ) - _conjg(work( i, j ))
                end do
             end do

          else if( lsame( side, 'r' ) ) then
             __QRM_PRNT_ERR('("_RPRFB: method not supported")')
          end if
       else
          __QRM_PRNT_ERR('("_RPRFB: method not supported")')
       end if

    else
       __QRM_PRNT_ERR('("_RPRFB: method not supported")')
    end if

    return
  end subroutine _xrprfb





  subroutine _xrpmm(side, uplo, trans, diag, m, n, alpha, a, x, ldx)

    implicit none
    character :: uplo, trans, diag, side
    integer   :: m, n, ldx
    _qrm_data :: a(*), x(ldx,*), alpha

    integer :: ma, na, p1, p2, p3, off
    logical, external :: lsame

    if (lsame(side, 'l')) then
       if(uplo .eq. 'l') then
          ma = (m+1)/2
          na = m+(m/2-ma+1)
          p1 = 1
          p2 = (m/2)*ma+1
          p3 = (m/2+1)*ma+1
          off = m/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrmm('l', 'u', notran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(tran, notran, m/2, n, ma, alpha, a(p1), ma, x(off,1), ldx, _qrm_one, x(1,1), ldx)
             call _xtrmm('l', 'l', tran, diag, ma, n, alpha, a(p2), ma, x(off,1), ldx)

          else

             call _xtrmm('l', 'l', notran, diag, ma, n, alpha, a(p2), ma, x(off,1), ldx)
             call _xgemm(notran, notran, ma, n, m/2, alpha, a(p1), ma, x(1,1), ldx, _qrm_one, x(off,1), ldx)
             call _xtrmm('l', 'u', tran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)

          end if

       else if (uplo .eq. 'u') then
          na  = (m+1)/2
          ma  = m+(m/2-na+1)
          p1  = 1
          p2  = m/2+1
          p3  = m/2+2
          off = m/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrmm('l', 'u', tran, diag, na, n, alpha, a(p2), ma, x(off,1), ldx)
             call _xgemm(tran, notran, na, n, m/2, alpha, a(p1), ma, x(1,1), ldx, _qrm_one, x(off,1), ldx)
             call _xtrmm('l', 'l', notran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)

          else

             call _xtrmm('l', 'l', tran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, notran, m/2, n, na, alpha, a(p1), ma, x(off,1), ldx, _qrm_one, x(1,1), ldx)
             call _xtrmm('l', 'u', notran, diag, na, n, alpha, a(p2), ma, x(off, 1), ldx)

          end if

       end if

    else

       if(uplo .eq. 'l') then
          ma = (n+1)/2
          na = n+(n/2-ma+1)
          p1 = 1
          p2 = (n/2)*ma+1
          p3 = (n/2+1)*ma+1
          off = n/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrmm('r', 'l', tran, diag, m, ma, alpha, a(p2), ma, x(1,off), ldx)
             call _xgemm(notran, tran, m, ma, n/2, alpha, x(1,1), ldx, a(p1), ma, _qrm_one, x(1,off), ldx)
             call _xtrmm('r', 'u', notran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)

          else
             call _xtrmm('r', 'u', tran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, notran, m, n/2, ma, alpha, x(1,off), ldx, a(p1), ma, _qrm_one, x(1,1), ldx)
             call _xtrmm('r', 'l', notran, diag, m, ma, alpha, a(p2), ma, x(1,off), ldx)

          end if

       else if (uplo .eq. 'u') then
          na  = (n+1)/2
          ma  = n+(n/2-na+1)
          p1  = 1
          p2  = n/2+1
          p3  = n/2+2
          off = n/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrmm('r', 'l', notran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, tran, m, n/2, na, alpha, x(1,off), ldx, a(p1), ma, _qrm_one, x(1,1), ldx)
             call _xtrmm('r', 'u', tran, diag, m, na, alpha, a(p2), ma, x(1,off), ldx)

          else

             call _xtrmm('r', 'u', notran, diag, m, na, alpha, a(p2), ma, x(1, off), ldx)
             call _xgemm(notran, notran, m, na, n/2, alpha, x(1,1), ldx, a(p1), ma, _qrm_one, x(1,off), ldx)
             call _xtrmm('r', 'l', tran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)

          end if

       end if

    end if



    return

  end subroutine _xrpmm



  subroutine _xrpsm(side, uplo, trans, diag, m, n, alpha, a, x, ldx)

    implicit none
    character :: uplo, trans, diag, side
    integer   :: m, n, ldx
    _qrm_data :: a(*), x(ldx,*), alpha

    integer :: ma, na, p1, p2, p3, off, ii
    logical, external :: lsame

    if (lsame(side, 'l')) then
       if(uplo .eq. 'l') then
          ma = (m+1)/2
          na = m+(m/2-ma+1)
          p1 = 1
          p2 = (m/2)*ma+1
          p3 = (m/2+1)*ma+1
          off = m/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrsm('l', 'l', tran, diag, ma, n, alpha, a(p2), ma, x(off,1), ldx)
             call _xgemm(tran, notran, m/2, n, ma, -_qrm_one, a(p1), ma, x(off,1), ldx, alpha, x(1,1), ldx)
             call _xtrsm('l', 'u', notran, diag, m/2, n, _qrm_one, a(p3), ma, x(1,1), ldx)

          else

             call _xtrsm('l', 'u', tran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, notran, ma, n, m/2, -_qrm_one, a(p1), ma, x(1,1), ldx, alpha, x(off,1), ldx)
             call _xtrsm('l', 'l', notran, diag, ma, n, _qrm_one, a(p2), ma, x(off,1), ldx)

          end if

       else if (uplo .eq. 'u') then
          na  = (m+1)/2
          ma  = m+(m/2-na+1)
          p1  = 1
          p2  = m/2+1
          p3  = m/2+2
          off = m/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrsm('l', 'l', notran, diag, m/2, n, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(tran, notran, na, n, m/2, -_qrm_one, a(p1), ma, x(1,1), ldx, alpha, x(off,1), ldx)
             call _xtrsm('l', 'u', tran, diag, na, n, _qrm_one, a(p2), ma, x(off,1), ldx)

          else

             call _xtrsm('l', 'u', notran, diag, na, n, alpha, a(p2), ma, x(off, 1), ldx)
             call _xgemm(notran, notran, m/2, n, na, -_qrm_one, a(p1), ma, x(off,1), ldx, alpha, x(1,1), ldx)
             call _xtrsm('l', 'l', tran, diag, m/2, n, _qrm_one, a(p3), ma, x(1,1), ldx)

          end if

       end if

    else

       if(uplo .eq. 'l') then
          ma = (n+1)/2
          na = n+(n/2-ma+1)
          p1 = 1
          p2 = (n/2)*ma+1
          p3 = (n/2+1)*ma+1
          off = n/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrsm('r', 'u', notran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, tran, m, ma, n/2, -_qrm_one, x(1,1), ldx, a(p1), ma, alpha, x(1,off), ldx)
             call _xtrsm('r', 'l', tran, diag, m, ma, _qrm_one, a(p2), ma, x(1,off), ldx)

          else

             call _xtrsm('r', 'l', notran, diag, m, ma, alpha, a(p2), ma, x(1,off), ldx)
             call _xgemm(notran, notran, m, n/2, ma, -_qrm_one, x(1,off), ldx, a(p1), ma, alpha, x(1,1), ldx)
             call _xtrsm('r', 'u', tran, diag, m, n/2, _qrm_one, a(p3), ma, x(1,1), ldx)

          end if

       else if (uplo .eq. 'u') then
          na  = (n+1)/2
          ma  = n+(n/2-na+1)
          p1  = 1
          p2  = n/2+1
          p3  = n/2+2
          off = n/2+1

          if((trans .eq. 't') .or. (trans .eq. 'c')) then

             call _xtrsm('r', 'u', tran, diag, m, na, alpha, a(p2), ma, x(1,off), ldx)
             call _xgemm(notran, tran, m, n/2, na, -_qrm_one, x(1,off), ldx, a(p1), ma, alpha, x(1,1), ldx)
             call _xtrsm('r', 'l', notran, diag, m, n/2, _qrm_one, a(p3), ma, x(1,1), ldx)

          else

             call _xtrsm('r', 'l', tran, diag, m, n/2, alpha, a(p3), ma, x(1,1), ldx)
             call _xgemm(notran, notran, m, na, n/2, -_qrm_one, x(1,1), ldx, a(p1), ma, alpha, x(1,off), ldx)
             call _xtrsm('r', 'u', notran, diag, m, na, _qrm_one, a(p2), ma, x(1, off), ldx)

          end if

       end if

    end if

    return

  end subroutine _xrpsm


  subroutine _qrm_print_rfpf(uplo, n, b)

    character :: uplo
    integer   :: n
    _qrm_data :: b(*)

    integer :: i, j, mb, nb

    if(uplo .eq. 'l') then
       mb = (n+1)/2
       nb = n+(n/2-mb+1)
    else if (uplo .eq. 'u') then
       nb  = (n+1)/2
       mb  = n+(n/2-nb+1)
    end if


    do i=1, mb
       do j=1, nb
          write(*,'(f6.1, x)', advance='no')b((j-1)*mb + i)
       end do
       write(*,'(" ")')
    end do

  end subroutine _qrm_print_rfpf



end module _qrm_rfpf_mod
