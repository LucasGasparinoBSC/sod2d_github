program unitt_hex27_test

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Tests how well the element type HEX27 is implemented through     !
   ! geometrical tests over the reference element.                    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use elem_hex          ! Module with all basic element definitions
   use quadrature_rules  ! Module with all quadrature rules (xgp, wgp)
   use jacobian_oper     ! Module with all jacobian operators

   implicit none

   integer(4)           :: ndime, nnode, ngaus
   integer(4)           :: idime, inode, igaus
   integer(4)           :: jdime, jnode, jgaus
   real(8)              :: s, t, z
   real(8)              :: detJe
   real(8)              :: sumNa, sumsumdNa
   real(8)              :: aux01
   real(8), allocatable :: elcod(:,:)
   real(8), allocatable :: N(:), dN(:,:)
   real(8), allocatable :: Je(:,:), He(:,:)
   logical              :: flag ! Test fail indicator

   !
   ! Initialize check
   !

   flag = .true. ! Will flip to false if a test fails

   !
   ! Element geometry (ref. elem.)
   !

   ndime = 3                    ! Element dimensions
   nnode = 27                   ! Nodes on element

   allocate(elcod(ndime,nnode)) ! Nodal coordinates array

   !                  xi      eta    zeta
   elcod(1:3,1)  = [-1.0d0, -1.0d0, -1.0d0]
   elcod(1:3,2)  = [-1.0d0, -1.0d0,  1.0d0]
   elcod(1:3,3)  = [ 1.0d0, -1.0d0,  1.0d0]
   elcod(1:3,4)  = [ 1.0d0, -1.0d0, -1.0d0]
   elcod(1:3,5)  = [-1.0d0,  1.0d0, -1.0d0]
   elcod(1:3,6)  = [-1.0d0,  1.0d0,  1.0d0]
   elcod(1:3,7)  = [ 1.0d0,  1.0d0,  1.0d0]
   elcod(1:3,8)  = [ 1.0d0,  1.0d0, -1.0d0]
   elcod(1:3,9)  = [-1.0d0, -1.0d0,  0.0d0]
   elcod(1:3,10) = [ 0.0d0, -1.0d0,  1.0d0]
   elcod(1:3,11) = [ 1.0d0, -1.0d0,  0.0d0]
   elcod(1:3,12) = [ 0.0d0, -1.0d0, -1.0d0]
   elcod(1:3,13) = [-1.0d0,  1.0d0,  0.0d0]
   elcod(1:3,14) = [ 0.0d0,  1.0d0,  1.0d0]
   elcod(1:3,15) = [ 1.0d0,  1.0d0,  0.0d0]
   elcod(1:3,16) = [ 0.0d0,  1.0d0, -1.0d0]
   elcod(1:3,17) = [-1.0d0,  0.0d0, -1.0d0]
   elcod(1:3,18) = [-1.0d0,  0.0d0,  1.0d0]
   elcod(1:3,19) = [ 1.0d0,  0.0d0,  1.0d0]
   elcod(1:3,20) = [ 1.0d0,  0.0d0, -1.0d0]
   elcod(1:3,21) = [ 0.0d0,  0.0d0, -1.0d0]
   elcod(1:3,22) = [ 0.0d0,  0.0d0,  1.0d0]
   elcod(1:3,23) = [-1.0d0,  0.0d0,  0.0d0]
   elcod(1:3,24) = [ 1.0d0,  0.0d0,  0.0d0]
   elcod(1:3,25) = [ 0.0d0, -1.0d0,  0.0d0]
   elcod(1:3,26) = [ 0.0d0,  1.0d0,  0.0d0]
   elcod(1:3,27) = [ 0.0d0,  0.0d0,  0.0d0]

   !
   ! Allocations
   !

   allocate(N(nnode))
   allocate(dN(ndime,nnode))

   !
   ! Test Na(xb) = 1 iff a==b
   !

   ! Na loop
   do inode = 1,nnode
      ! xb loop
      do jnode = 1,nnode
         s = elcod(1,jnode) ! isopar. x1
         t = elcod(2,jnode) ! isopar. x2
         z = elcod(3,jnode) ! isopar. x3
         call hex27(s,t,z,N,dN)
         ! Test
         if (jnode==inode .and. N(inode) .ne. 1.0d0) then ! Na(xa) = 1
            write(*,*) inode, jnode, N(inode)
            flag = .false.
            exit
         else if (jnode .ne. inode .and. abs(N(inode)) .ne. 0.0d0) then ! Na(xb) = 0
            write(*,*) inode, jnode, N(inode)
            flag = .false.
            exit
         end if
      end do
      ! Stop test if any errrors happen
      if (flag .eqv. .false.) then
         exit
      end if
   end do

   !
   ! Test sum(Na(xgp)) = 1 and sum(sum(dNa(xgp))) = 0
   !

   sumNa = 0.0d0          ! Initialize sum to 0
   sumsumdNa = 0.0d0      ! Initialize sum to 0
   s = 0.00d0             ! X1 gp of order 1
   t = 0.00d0             ! X2 gp of order 1
   z = 0.00d0             ! X3 gp of order 1
   call hex27(s,t,z,N,dN) ! Fill N(xgp) and dN(xgp)
   sumNa = sum(N)
   sumsumdNa = sum(dN)

   ! Test
   if (sumNa .ne. 1.0d0) then
      write(*,*) sumNa
      flag = .false.
   end if
   if (abs(sumsumdNa) .ne. 0.0d0) then
      write(*,*) sumsumdNa
      flag = .false.
   end if

   !
   ! Test Je = I and det(Je) = 1
   !

   allocate(Je(ndime,ndime)) ! Jacobian
   allocate(He(ndime,ndime)) ! Inverse Jacobian

   ! Test with 1gp
   s = 0.00d0                                           ! X1 gp of order 1
   t = 0.00d0                                           ! X2 gp of order 1
   z = 0.00d0                                           ! X3 gp of order 1
   call hex27(s,t,z,N,dN)                               ! Fill N and dN at Gauss point
   call elem_jacobian(ndime,nnode,elcod,dN,Je,detJe,He) ! Compute Je, detJe and inv(Je)

end program unitt_hex27_test
