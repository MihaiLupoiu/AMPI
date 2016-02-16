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
!> @file qrm_sort_mod.F90
!! This file contains a module with sorting facilities
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This module contains routines for sorting
!!

!> Mostly plain implementations from:
!!
!!      D. E. Knuth <em>"The Art of Computer Programming,"</em>
!!      vol.3: Sorting and Searching, Addison-Wesley, 1973

module qrm_sort_mod

  interface qrm_mergesort
     module procedure  qrm_mergesorti, qrm_mergesortd
  end interface


  interface qrm_mergeswap
     module procedure qrm_mergeswapii, qrm_mergeswapi, qrm_mergeswapd
     module procedure qrm_mergeswapis, qrm_mergeswapid, qrm_mergeswapic, qrm_mergeswapiz 
  end interface

contains

  subroutine qrm_mergeswapii(n, l, a1, a2)
    implicit none
    integer   :: n
    integer,intent(inout)   :: l(0:n+1), a1(n), a2(n)
    integer   :: i, lp, iswap

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       iswap    = a2(lp)
       a2(lp)   = a2(i)
       a2(i)    = iswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapii

  subroutine qrm_mergeswapi(n, l, a)

    integer   :: i, lp, iswap, n
    integer   :: l(0:), a(:)

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap  = a(lp)
       a(lp)  = a(i)
       a(i)   = iswap
       iswap  = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapi


  subroutine qrm_mergeswapis(n, l, a1, a2)

    integer          :: i, lp, iswap, n
    integer          :: l(0:n+1), a1(n)
    real(kind(1.e0)) :: a2(n)

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       iswap    = a2(lp)
       a2(lp)   = a2(i)
       a2(i)    = iswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapis

  subroutine qrm_mergeswapid(n, l, a1, a2)

    integer     :: i, lp, iswap, n
    integer     :: l(0:n+1), a1(n)
    real(kind(1.d0)) :: a2(n), dswap

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       dswap    = a2(lp)
       a2(lp)   = a2(i)
       a2(i)    = dswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapid

  subroutine qrm_mergeswapd(n, l, a1)

    integer     :: i, lp, iswap, n
    integer     :: l(0:n+1)
    real(kind(1.d0)) :: a1(n)

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapd

  subroutine qrm_mergeswapic(n, l, a1, a2)

    integer     :: i, lp, iswap, n
    integer     :: l(0:n+1), a1(n)
    complex(kind(1.e0)) :: a2(n)

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       iswap    = a2(lp)
       a2(lp)   = a2(i)
       a2(i)    = iswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapic

  subroutine qrm_mergeswapiz(n, l, a1, a2)

    integer     :: i, lp, iswap, n
    integer     :: l(0:n+1), a1(n)
    complex(kind(1.d0)) :: a2(n)

    lp = l(0)
    i  = 1
    do 
       if ((lp .eq. 0).or.(i>n)) exit
       do 
          if (lp >= i) exit
          lp = l(lp)
       end do
       iswap    = a1(lp)
       a1(lp)   = a1(i)
       a1(i)    = iswap
       iswap    = a2(lp)
       a2(lp)   = a2(i)
       a2(i)    = iswap
       iswap    = l(lp)
       l(lp) = l(i)
       l(i)  = lp
       lp = iswap 
       i  = i + 1
    enddo

    return

  end subroutine qrm_mergeswapiz



  subroutine qrm_mergesorti(n, k, l, order)

    implicit none
    integer                 :: n
    integer, intent(inout)  :: k(n), l(0:n+1)
    integer, optional       :: order

    integer    :: p, q, s, t, i, iord

    iord = 1
    if(present(order)) iord = order

    if((iord .ne. 1) .and. (iord .ne. -1)) then
       write(*,'("Wrong input in mergesort")')
       return
    end if

100 continue
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
       if (iord*k(p) .le. iord*k(p+1)) then
          l(p) = p + 1
       else
          l(t) = - (p+1)
          t = p
       end if
    end do
    l(t) = 0
    l(n) = 0

    if (l(n+1) .eq. 0) then
       return 
    else
       l(n+1) = iabs(l(n+1))
    end if

200 continue
    s = 0
    t = n+1
    p = l(s)
    q = l(t)

    if(q .eq. 0) return

300 continue
    if(iord*k(p) .gt. iord*k(q)) goto 600 

400 continue
    l(s) = sign(p,l(s))
    s = p
    p = l(p)
    if (p .gt. 0) goto 300

500 continue
    l(s) = q
    s = t
    do
       t = q
       q = l(q)
       if (q .le. 0) exit
    end do
    goto 800

600 continue
    l(s) = sign(q, l(s))
    s = q
    q = l(q)
    if (q .gt. 0) goto 300

700 continue
    l(s) = p
    s = t
    do
       t = p
       p = l(p)
       if (p .le. 0) exit
    end do

800 continue
    p = -p
    q = -q
    if(q.eq.0) then
       l(s) = sign(p, l(s))
       l(t) = 0
       goto 200
    end if

    goto 300

    return

  end subroutine qrm_mergesorti



  subroutine qrm_mergesortd(n, k, l, order)

    !      Plain implementation of the merge-sort algorithm
    !      as described in:

    !      D. E. Knuth "The Art of Computer Programming,"
    !      vol.3: Sorting and Searching, Addison-Wesley, 1973
    implicit none
    integer                :: n
    integer, intent(inout) :: l(0:n+1)
    real(kind(1.d0))       :: k(n)
    integer, optional      :: order

    integer    :: p, q, s, t, i, iord

    iord = 1
    if(present(order)) iord = order

    if((iord .ne. 1) .and. (iord .ne. -1)) then
       write(*,'("Wrong input in mergesort")')
       return
    end if

100 continue
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
       if (iord*k(p) .le. iord*k(p+1)) then
          l(p) = p + 1
       else
          l(t) = - (p+1)
          t = p
       end if
    end do
    l(t) = 0
    l(n) = 0

    if (l(n+1) .eq. 0) then
       return 
    else
       l(n+1) = iabs(l(n+1))
    end if

200 continue
    s = 0
    t = n+1
    p = l(s)
    q = l(t)

    if(q .eq. 0) return

300 continue
    if(iord*k(p) .gt. iord*k(q)) goto 600 

400 continue
    l(s) = sign(p,l(s))
    s = p
    p = l(p)
    if (p .gt. 0) goto 300

500 continue
    l(s) = q
    s = t
    do
       t = q
       q = l(q)
       if (q .le. 0) exit
    end do
    goto 800

600 continue
    l(s) = sign(q, l(s))
    s = q
    q = l(q)
    if (q .gt. 0) goto 300

700 continue
    l(s) = p
    s = t
    do
       t = p
       p = l(p)
       if (p .le. 0) exit
    end do

800 continue
    p = -p
    q = -q
    if(q.eq.0) then
       l(s) = sign(p, l(s))
       l(t) = 0
       goto 200
    end if

    goto 300

    return

  end subroutine qrm_mergesortd

end module qrm_sort_mod
