module elem_hex

   contains  

      subroutine hex08(xi,eta,zeta,N,dN) ! HEX08 element     

         implicit none
         
         real(8), intent(in)  :: xi, eta, zeta
         real(8), intent(out) :: N(8), dN(3,8)

         ! Regarding coordinates for the isopar space:
         ! xi is left to right (-1,1)
         ! eta is bottom to top (-1,1)
         ! zeta is back to front (-1,1)
         ! These are like the Gmsh (x,y,z) system.
         ! These will apply to all orders of HEX

         N(1) = (1.0d0-xi)*(1.0d0-eta)*(1.0d0-zeta)
         N(2) = (1.0d0-xi)*(1.0d0-eta)*(1.0d0+zeta)
         N(3) = (1.0d0+xi)*(1.0d0-eta)*(1.0d0+zeta)
         N(4) = (1.0d0+xi)*(1.0d0-eta)*(1.0d0-zeta)
         N(5) = (1.0d0-xi)*(1.0d0+eta)*(1.0d0-zeta)
         N(6) = (1.0d0-xi)*(1.0d0+eta)*(1.0d0+zeta)
         N(7) = (1.0d0+xi)*(1.0d0+eta)*(1.0d0+zeta)
         N(8) = (1.0d0+xi)*(1.0d0+eta)*(1.0d0-zeta)
         N = 0.125d0*N
         
         dN(1,1) = -(1.0d0-eta)*(1.0d0-zeta)
         dN(2,1) = -(1.0d0-xi )*(1.0d0-zeta)
         dN(3,1) = -(1.0d0-xi )*(1.0d0-eta )

         dN(1,2) = -(1.0d0-eta)*(1.0d0+zeta)
         dN(2,2) = -(1.0d0-xi )*(1.0d0+zeta)
         dN(3,2) =  (1.0d0-xi )*(1.0d0-eta )

         dN(1,3) =  (1.0d0-eta)*(1.0d0+zeta)
         dN(2,3) = -(1.0d0+xi )*(1.0d0+zeta)
         dN(3,3) =  (1.0d0+xi )*(1.0d0-eta )

         dN(1,4) =  (1.0d0-eta)*(1.0d0-zeta)
         dN(2,4) = -(1.0d0+xi )*(1.0d0-zeta)
         dN(3,4) = -(1.0d0+xi )*(1.0d0-eta )

         dN(1,5) = -(1.0d0+eta)*(1.0d0-zeta)
         dN(2,5) =  (1.0d0-xi )*(1.0d0-zeta)
         dN(3,5) = -(1.0d0-xi )*(1.0d0+eta )

         dN(1,6) = -(1.0d0+eta)*(1.0d0+zeta)
         dN(2,6) =  (1.0d0-xi )*(1.0d0+zeta)
         dN(3,6) =  (1.0d0-xi )*(1.0d0+eta )

         dN(1,7) =  (1.0d0+eta)*(1.0d0+zeta)
         dN(2,7) =  (1.0d0+xi )*(1.0d0+zeta)
         dN(3,7) =  (1.0d0+xi )*(1.0d0+eta )

         dN(1,8) =  (1.0d0+eta)*(1.0d0-zeta)
         dN(2,8) =  (1.0d0+xi )*(1.0d0-zeta)
         dN(3,8) = -(1.0d0+xi )*(1.0d0+eta )
         dN = 0.125d0*dN

      end subroutine hex08
      
      subroutine hex27(xi,eta,zeta,N,dN) ! HEX27 element

         implicit none
         
         real(8), intent(in)  :: xi, eta, zeta
         real(8), intent(out) :: N(27), dN(3,27)
         real(8)              :: n1xi, n2xi, n3xi
         real(8)              :: n1eta, n2eta, n3eta
         real(8)              :: n1zeta, n2zeta, n3zeta
         real(8)              :: dn1xi, dn2xi, dn3xi
         real(8)              :: dn1eta, dn2eta, dn3eta
         real(8)              :: dn1zeta, dn2zeta, dn3zeta

         !
         ! Triquadratic brick
         !        

         ! 1D quadratic functions for xi, eta, zeta
         n1xi   = 0.5d0*xi*(xi-1.0d0)
         n1eta  = 0.5d0*eta*(eta-1.0d0)
         n1zeta = 0.5d0*zeta*(zeta-1.0d0)
         n2xi   = 1.0d0-xi*xi
         n2eta  = 1.0d0-eta*eta
         n2zeta = 1.0d0-zeta*zeta
         n3xi   = 0.5d0*xi*(xi+1.0d0)
         n3eta  = 0.5d0*eta*(eta+1.0d0)
         n3zeta = 0.5d0*zeta*(zeta+1.0d0)

         ! 1D quadratic function derivatives for xi, eta, zeta
         dn1xi   = 0.5d0*(2.0d0*xi-1.0d0)
         dn1eta  = 0.5d0*(2.0d0*eta-1.0d0)
         dn1zeta = 0.5d0*(2.0d0*zeta-1.0d0)
         dn2xi   = -2.0d0*xi
         dn2eta  = -2.0d0*eta
         dn2zeta = -2.0d0*zeta
         dn3xi   = 0.5d0*(2.0d0*xi+1.0d0)
         dn3eta  = 0.5d0*(2.0d0*eta+1.0d0)
         dn3zeta = 0.5d0*(2.0d0*zeta+1.0d0)

         !
         ! Shape functions
         !

         ! Bottom, corner
         N( 1) = n1xi*n1eta*n1zeta
         N( 2) = n1xi*n1eta*n3zeta
         N( 3) = n3xi*n1eta*n3zeta
         N( 4) = n3xi*n1eta*n1zeta

         ! Top, corner
         N( 5) = n1xi*n3eta*n1zeta
         N( 6) = n1xi*n3eta*n3zeta
         N( 7) = n3xi*n3eta*n3zeta
         N( 8) = n3xi*n3eta*n1zeta

         ! Bottom, edges
         N( 9) = n1xi*n1eta*n2zeta
         N(10) = n2xi*n1eta*n3zeta
         N(11) = n3xi*n1eta*n2zeta
         N(12) = n2xi*n1eta*n1zeta

         ! Top, edges
         N(13) = n1xi*n3eta*n2zeta
         N(14) = n2xi*n3eta*n3zeta
         N(15) = n3xi*n3eta*n2zeta
         N(16) = n2xi*n3eta*n1zeta

         ! Middle, edges
         N(17) = n1xi*n2eta*n1zeta
         N(18) = n1xi*n2eta*n3zeta
         N(19) = n3xi*n2eta*n3zeta
         N(20) = n3xi*n2eta*n1zeta

         ! Faces
         N(21) = n2xi*n2eta*n1zeta
         N(22) = n2xi*n2eta*n3zeta
         N(23) = n1xi*n2eta*n2zeta
         N(24) = n3xi*n2eta*n2zeta
         N(25) = n2xi*n1eta*n2zeta
         N(26) = n2xi*n3eta*n2zeta

         ! Volume
         N(27) = n2xi*n2eta*n2zeta

         !
         ! Shape function derivatives
         !

         ! Direction xi

         ! Bottom, corner
         dN(1, 1) = dn1xi*n1eta*n1zeta
         dN(1, 2) = dn1xi*n1eta*n3zeta
         dN(1, 3) = dn3xi*n1eta*n3zeta
         dN(1, 4) = dn3xi*n1eta*n1zeta

         ! Top, corner
         dN(1, 5) = dn1xi*n3eta*n1zeta
         dN(1, 6) = dn1xi*n3eta*n3zeta
         dN(1, 7) = dn3xi*n3eta*n3zeta
         dN(1, 8) = dn3xi*n3eta*n1zeta

         ! Bottom, edges
         dN(1, 9) = dn1xi*n1eta*n2zeta
         dN(1,10) = dn2xi*n1eta*n3zeta
         dN(1,11) = dn3xi*n1eta*n2zeta
         dN(1,12) = dn2xi*n1eta*n1zeta

         ! Top, edges
         dN(1,13) = dn1xi*n3eta*n2zeta
         dN(1,14) = dn2xi*n3eta*n3zeta
         dN(1,15) = dn3xi*n3eta*n2zeta
         dN(1,16) = dn2xi*n3eta*n1zeta

         ! Middle, edges
         dN(1,17) = dn1xi*n2eta*n1zeta
         dN(1,18) = dn1xi*n2eta*n3zeta
         dN(1,19) = dn3xi*n2eta*n3zeta
         dN(1,20) = dn3xi*n2eta*n1zeta

         ! Faces
         dN(1,21) = dn2xi*n2eta*n1zeta
         dN(1,22) = dn2xi*n2eta*n3zeta
         dN(1,23) = dn1xi*n2eta*n2zeta
         dN(1,24) = dn3xi*n2eta*n2zeta
         dN(1,25) = dn2xi*n1eta*n2zeta
         dN(1,26) = dn2xi*n3eta*n2zeta

         ! Volume
         dN(1,27) = dn2xi*n2eta*n2zeta

         ! Direction eta

         ! Bottom, corner
         dN(2, 1) = n1xi*dn1eta*n1zeta
         dN(2, 2) = n1xi*dn1eta*n3zeta
         dN(2, 3) = n3xi*dn1eta*n3zeta
         dN(2, 4) = n3xi*dn1eta*n1zeta

         ! Top, corner
         dN(2, 5) = n1xi*dn3eta*n1zeta
         dN(2, 6) = n1xi*dn3eta*n3zeta
         dN(2, 7) = n3xi*dn3eta*n3zeta
         dN(2, 8) = n3xi*dn3eta*n1zeta

         ! Bottom, edges
         dN(2, 9) = n1xi*dn1eta*n2zeta
         dN(2,10) = n2xi*dn1eta*n3zeta
         dN(2,11) = n3xi*dn1eta*n2zeta
         dN(2,12) = n2xi*dn1eta*n1zeta

         ! Top, edges
         dN(2,13) = n1xi*dn3eta*n2zeta
         dN(2,14) = n2xi*dn3eta*n3zeta
         dN(2,15) = n3xi*dn3eta*n2zeta
         dN(2,16) = n2xi*dn3eta*n1zeta

         ! Middle, edges
         dN(2,17) = n1xi*dn2eta*n1zeta
         dN(2,18) = n1xi*dn2eta*n3zeta
         dN(2,19) = n3xi*dn2eta*n3zeta
         dN(2,20) = n3xi*dn2eta*n1zeta

         ! Faces
         dN(2,21) = n2xi*dn2eta*n1zeta
         dN(2,22) = n2xi*dn2eta*n3zeta
         dN(2,23) = n1xi*dn2eta*n2zeta
         dN(2,24) = n3xi*dn2eta*n2zeta
         dN(2,25) = n2xi*dn1eta*n2zeta
         dN(2,26) = n2xi*dn3eta*n2zeta

         ! Volume
         dN(2,27) = n2xi*dn2eta*n2zeta

         ! Direction zeta

         ! Bottom, corner
         dN(3, 1) = n1xi*n1eta*dn1zeta
         dN(3, 2) = n1xi*n1eta*dn3zeta
         dN(3, 3) = n3xi*n1eta*dn3zeta
         dN(3, 4) = n3xi*n1eta*dn1zeta

         ! Top, corner
         dN(3, 5) = n1xi*n3eta*dn1zeta
         dN(3, 6) = n1xi*n3eta*dn3zeta
         dN(3, 7) = n3xi*n3eta*dn3zeta
         dN(3, 8) = n3xi*n3eta*dn1zeta

         ! Bottom, edges
         dN(3, 9) = n1xi*n1eta*dn2zeta
         dN(3,10) = n2xi*n1eta*dn3zeta
         dN(3,11) = n3xi*n1eta*dn2zeta
         dN(3,12) = n2xi*n1eta*dn1zeta

         ! Top, edges
         dN(3,13) = n1xi*n3eta*dn2zeta
         dN(3,14) = n2xi*n3eta*dn3zeta
         dN(3,15) = n3xi*n3eta*dn2zeta
         dN(3,16) = n2xi*n3eta*dn1zeta

         ! Middle, edges
         dN(3,17) = n1xi*n2eta*dn1zeta
         dN(3,18) = n1xi*n2eta*dn3zeta
         dN(3,19) = n3xi*n2eta*dn3zeta
         dN(3,20) = n3xi*n2eta*dn1zeta

         ! Faces
         dN(3,21) = n2xi*n2eta*dn1zeta
         dN(3,22) = n2xi*n2eta*dn3zeta
         dN(3,23) = n1xi*n2eta*dn2zeta
         dN(3,24) = n3xi*n2eta*dn2zeta
         dN(3,25) = n2xi*n1eta*dn2zeta
         dN(3,26) = n2xi*n3eta*dn2zeta

         ! Volume
         dN(3,27) = n2xi*n2eta*dn2zeta


      end subroutine hex27
      
      subroutine hexa_edges(ielem,nelem,nnode,npoin,ndime,connec,coord,ncorner,nedge,dist)

         implicit none
         
         integer(4), intent(in)            :: ielem, nelem, nnode, npoin, ndime
         integer(4), intent(in)            :: connec(nelem,nnode)
         real(8),    intent(in)            :: coord(npoin,ndime)
         integer(4), intent(out)           :: ncorner, nedge
         real(8),    intent(out)           :: dist(12,ndime)
         integer(4)                        :: ind(nnode)
         real(8)                           :: xp(12,ndime)
         
         ind = connec(ielem,:)
         ncorner = 8
         nedge = 12
         
         xp(1:8,1:ndime) = coord(ind(1:8),1:ndime) ! Corner coordinates
         dist(1,:) = xp(2,:)-xp(1,:)
         dist(2,:) = xp(3,:)-xp(2,:)
         dist(3,:) = xp(4,:)-xp(3,:)
         dist(4,:) = xp(1,:)-xp(4,:)
         
         dist(5,:) = xp(6,:)-xp(5,:)
         dist(6,:) = xp(7,:)-xp(6,:)
         dist(7,:) = xp(8,:)-xp(7,:)
         dist(8,:) = xp(5,:)-xp(8,:)
         
         dist(9,:) = xp(5,:)-xp(1,:)
         dist(10,:) = xp(6,:)-xp(2,:)
         dist(11,:) = xp(7,:)-xp(3,:)
         dist(12,:) = xp(8,:)-xp(4,:)
      
      end subroutine hexa_edges

end module
