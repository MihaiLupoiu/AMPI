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
!> @file qrm_print_tree.F90
!! This file contains various routine for printing trees
!!
!! $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
!! $Author: abuttari $
!! $Version: 1.1$
!! $Revision: 1980 $
!!
!! ##############################################################################################


!> @brief This subroutine prints on a file the elimination tree
!! described by a parent array. The tree is written in "dot"
!! format (see graphviz)
!!
!! @param[in] file   the file where to print the tree
!!
!! @param[in] parent an integer array. parent(i)=j means that node/variable
!!                   j is the father of node/variable i in the tree
!!
!! @param[in] n      an integer containing the number of nodes in the tree
subroutine qrm_print_elim_tree(file, parent, n)

  character :: file*(*)
  integer   :: parent(:)
  integer   :: n
  
  integer   :: i

  open(4, file=file, action='write')
  
  write(4,'("graph G {")')
  write(4,'("node [color=black,")')
  write(4,'("fillcolor=white,")')
  write(4,'("shape=circle,")')
  write(4,'("style=filled")')
  write(4,'("];")')
  do i=1, n
     if(parent(i) .ne. 0) then
        write(4,'(i6," -- ",i6)')parent(i),i
     else
        write(4,'(i6)')i
     end if
  end do

  write(4,'("}")')

  close(4)
  return

end subroutine qrm_print_elim_tree




!> @brief This subroutine prints on a file the assembly tree
!! described by a parent and a postorder arrays. The tree is written in "dot"
!! format (see graphviz)
!!
!! @param[in] file   the file where to print the tree
!!
!! @param[in] parent an integer array containing the elimination tree in input
!!                   and the assembly tree on output. The meaning of parent on output is:
!!                   - parent(i) = j>0: i is the principal variable of a node and j if the principal
!!                                      variable of its father node
!!                   - parent(i) = j=0: i is a principal variable of a root node
!!                   - parent(i) = j<0: i is a subordinate variable inside a node whose
!!                             principal variable is j.
!!
!! @param[in] rc     an array containing
!!
!! @param[in] n      an integer containing the number of nodes in the tree
!!
subroutine qrm_print_asm_tree(file, parent, rc, n)

  use qrm_mem_mod

  character :: file*(*)
  integer   :: parent(:), rc(:)
  integer   :: n
  
  integer   :: i, maxnv
  integer, allocatable :: nvar(:)
  real(kind(1.d0)) :: s

  call qrm_aalloc(nvar, n)
  nvar=1
  maxnv = 0
  do i=1, n
     if(parent(i) .lt. 0) then
        nvar(-parent(i)) = nvar(-parent(i))+1
        if(nvar(-parent(i)) .gt. maxnv) maxnv=nvar(-parent(i))
     end if
  end do


  open(4, file=file, action='write')
  
  write(4,'("graph G {")')
  write(4,'("node [color=black,")')
  write(4,'("fillcolor=white,")')
  write(4,'("shape=circle,")')
  write(4,'("style=filled")')
  write(4,'("];")')

  do i=1, n
     if(parent(i) .ge. 0) then
        s = 5.d0/(real(maxnv,kind(1.d0))/real(nvar(i),kind(1.d0)) )
        write(4,'("node",i6.6,"[label="" node:",i6,"\n rc:",i6,'//&
             &'"\n nv:",i6,""", shape=circle, height=",i6.6,",'// &
             &'fontsize=",i6.6,"];")')i,i,rc(i),nvar(i),nvar(i),20*nvar(i)
     end if
  end do
  
  do i=1, n
     if(parent(i) .gt. 0) then
        write(4,'("node",i6.6," -- node",i6.6)')parent(i),i
     else if (parent(i) .eq. 0) then
        write(4,'("node",i6.6)')i
     end if
  end do

  write(4,'("}")')

  close(4)

  call qrm_adealloc(nvar)

  return

end subroutine qrm_print_asm_tree


!> @brief prints an assembly tree in compressed format
!! 
!! @param[in] adata  a qrm_adata_type data containing the tree
!! @param[in] file   the file where to print the tree
!! @param[in] weight optional: weight(i) is the weight of the subtree rooted at i
!!
subroutine qrm_print_nsteps_tree(file, adata, weight)

  use qrm_adata_mod
  use qrm_mem_mod
  implicit none

  type(qrm_adata_type) :: adata
  character :: file*(*)
  real(kind(1.d0)), optional :: weight(:)

  integer :: i, p, n
  integer, allocatable :: tmp(:)


  open(4, file=file, action='write')
  
  write(4,'("graph G {")')
  write(4,'("node [color=black,")')
  write(4,'("fillcolor=white,")')
  write(4,'("shape=circle,")')
  write(4,'("style=filled")')
  write(4,'("];")')

  if(present(weight)) then
     ! if(adata%nnodes .gt. 1000) then
     if(.false.) then
        call qrm_aalloc(tmp, adata%nnodes+1, lbnd=0)
        tmp = 0

        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              write(4,'("node",i6.6,"[fillcolor=red, label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"\n w:",f4.1,"% ""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i)), weight(i)
              tmp(i) = 1
              p = i
              do
                 p = adata%parent(p)
                 if((p .gt. 0) .and. (tmp(p) .eq.0)) then
                    write(4,'("node",i6.6,"[label="" node:",i6,"\n m:",i6,"\n n:",'//&
                         &'i6,"\n np:",i6,"\n pv:",i6,"\n w:",f4.1,"% ""];")')p,p,adata%nfrows(p),adata%rc(p),&
                         &adata%cp_ptr(p+1)-adata%cp_ptr(p),adata%cperm(adata%cp_ptr(p)), weight(p)
                    tmp(p) = 1
                 else
                    exit
                 end if
              end do
           end if
        end do

        tmp = 0

        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              n = i
              do
                 p = adata%parent(n)
                 if (tmp(n) .eq. 0) then
                    if(p .gt. 0) then
                       write(4,'("node",i6.6," -> node",i6.6)')p,n
                       tmp(n) = 1
                       n = p
                    else if (p .eq. 0) then
                       write(4,'("node",i6.6)')n
                       tmp(n) = 1
                       exit
                    end if
                 else
                    exit
                 end if
              end do
           end if
        end do

        call qrm_adealloc(tmp)

     else
        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              write(4,'("node",i6.6,"[fillcolor=red, label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"\n w:",f4.1,"% ""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i)),weight(i)
           else
              write(4,'("node",i6.6,"[label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"\n w:",f4.1,"% ""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i)),weight(i)
           end if
        end do

        do i=1, adata%nnodes
           if(adata%parent(i) .gt. 0) then
              write(4,'("node",i6.6," -> node",i6.6)')adata%parent(i),i
           else if (adata%parent(i) .eq. 0) then
              write(4,'("node",i6.6)')i
           end if
        end do
     end if
  else
     if(adata%nnodes .gt. 10) then
        call qrm_aalloc(tmp, adata%nnodes+1, lbnd=0)
        tmp = 0

        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              write(4,'("node",i6.6,"[fillcolor=red, label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i))
              tmp(i) = 1
              p = i
              do
                 p = adata%parent(p)
                 if((p .gt. 0) .and. (tmp(p) .eq.0)) then
                    write(4,'("node",i6.6,"[label="" node:",i6,"\n m:",i6,"\n n:",'//&
                         &'i6,"\n np:",i6,"\n pv:",i6,"""];")')p,p,adata%nfrows(p),adata%rc(p),&
                         &adata%cp_ptr(p+1)-adata%cp_ptr(p),adata%cperm(adata%cp_ptr(p))
                    tmp(p) = 1
                 else
                    exit
                 end if
              end do
           end if
        end do

        tmp = 0

        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              n = i
              do
                 p = adata%parent(n)
                 if (tmp(n) .eq. 0) then
                    if(p .gt. 0) then
                       write(4,'("node",i6.6," -- node",i6.6)')p,n
                       tmp(n) = 1
                       n = p
                    else if (p .eq. 0) then
                       write(4,'("node",i6.6)')n
                       tmp(n) = 1
                       exit
                    end if
                 else
                    exit
                 end if
              end do
           end if
        end do

        call qrm_adealloc(tmp)

     else
        do i=1, adata%nnodes
           if(adata%small(i) .eq. 1) then
              write(4,'("node",i6.6,"[fillcolor=red, label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i))
           else
              write(4,'("node",i6.6,"[label="" node:",i6,"\n m:",i6,"\n n:",'//&
                   &'i6,"\n np:",i6,"\n pv:",i6,"""];")')i,i,adata%nfrows(i),adata%rc(i),&
                   &adata%cp_ptr(i+1)-adata%cp_ptr(i),adata%cperm(adata%cp_ptr(i))
           end if
        end do

        do i=1, adata%nnodes
           if(adata%parent(i) .gt. 0) then
              write(4,'("node",i6.6," -- node",i6.6)')adata%parent(i),i
           else if (adata%parent(i) .eq. 0) then
              write(4,'("node",i6.6)')i
           end if
        end do
     end if
  end if

  write(4,'("}")')

  close(4)


  return

end subroutine qrm_print_nsteps_tree




