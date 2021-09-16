module elem_hex

   contains  

      subroutine hex08(s,t,z,N,dN) ! HEX08 element     

         implicit none
         
         real(8), intent(in)  :: s, t, z
         real(8), intent(out) :: N(8), dN(3,8)

         N(1) = (1.0d0-s)*(1.0d0-z)*(1.0d0-t)
         N(2) = (1.0d0+s)*(1.0d0-z)*(1.0d0-t)
         N(3) = (1.0d0+s)*(1.0d0+z)*(1.0d0-t)
         N(4) = (1.0d0-s)*(1.0d0+z)*(1.0d0-t)
         N(5) = (1.0d0-s)*(1.0d0-z)*(1.0d0+t)
         N(6) = (1.0d0+s)*(1.0d0-z)*(1.0d0+t)
         N(7) = (1.0d0+s)*(1.0d0+z)*(1.0d0+t)
         N(8) = (1.0d0-s)*(1.0d0+z)*(1.0d0+t)
         N = 0.125d0*N
         
         dN(1,1) = -(1.0d0-z)*(1.0d0-t)
         dN(2,1) = -(1.0d0-s)*(1.0d0-t)
         dN(3,1) = -(1.0d0-s)*(1.0d0-z)
         dN(1,2) =  (1.0d0-z)*(1.0d0-t)
         dN(2,2) = -(1.0d0+s)*(1.0d0-t)
         dN(3,2) = -(1.0d0+s)*(1.0d0-z)
         dN(1,3) =  (1.0d0+z)*(1.0d0-t)
         dN(2,3) =  (1.0d0+s)*(1.0d0-t)
         dN(3,3) = -(1.0d0+s)*(1.0d0+z)
         dN(1,4) = -(1.0d0+z)*(1.0d0-t)
         dN(2,4) =  (1.0d0-s)*(1.0d0-t)
         dN(3,4) = -(1.0d0-s)*(1.0d0+z)
         dN(1,5) = -(1.0d0-z)*(1.0d0+t)
         dN(2,5) = -(1.0d0-s)*(1.0d0+t)
         dN(3,5) =  (1.0d0-s)*(1.0d0-z)
         dN(1,6) =  (1.0d0-z)*(1.0d0+t)
         dN(2,6) = -(1.0d0+s)*(1.0d0+t)
         dN(3,6) =  (1.0d0+s)*(1.0d0-z)
         dN(1,7) =  (1.0d0+z)*(1.0d0+t)
         dN(2,7) =  (1.0d0+s)*(1.0d0+t)
         dN(3,7) =  (1.0d0+s)*(1.0d0+z)
         dN(1,8) = -(1.0d0+z)*(1.0d0+t)
         dN(2,8) =  (1.0d0-s)*(1.0d0+t)
         dN(3,8) =  (1.0d0-s)*(1.0d0+z)
         dN = 0.125d0*dN
      
      end subroutine hex08
      
      subroutine hex27(s,t,z,N,dN) ! HEX27 element

         ! TODO: CHECK ORDERING FOR OUTPUT
         
         implicit none
         
         real(8), intent(in)  :: s, t, z
         real(8), intent(out) :: N(27), dN(3,27)
         real(8)                   :: z1,z2,z3,z4,s1,s2,s3,s4
         real(8)                   :: t1,t2,t3,t4,sq,tp,zp
         real(8)                   :: sl,tl,zl

         !
         ! Triquadratic brick
         !        
         sl=s*(s-1.0d0)
         tl=t*(t-1.0d0)
         zl=z*(z-1.0d0)
         sq=s*(s+1.0d0)
         tp=t*(t+1.0d0)
         zp=z*(z+1.0d0)
         s1= 2.0d0*s-1.0d0
         t1= 2.0d0*t-1.0d0
         z1= 2.0d0*z-1.0d0
         s2= 1.0d0-s*s
         t2= 1.0d0-t*t
         z2= 1.0d0-z*z
         s3= 1.0d0+2.0d0*s
         t3= 1.0d0+2.0d0*t
         z3= 1.0d0+2.0d0*z
         s4=-2.0d0*s
         t4=-2.0d0*t
         z4=-2.0d0*z
         N(   1) = 0.125d0*sl*tl*zl
         dN(1, 1) = 0.125d0*s1*tl*zl
         dN(2, 1) = 0.125d0*sl*t1*zl
         dN(3, 1) = 0.125d0*sl*tl*z1
         N(   2) = 0.125d0*sq*tl*zl
         dN(1, 2) = 0.125d0*s3*tl*zl
         dN(2, 2) = 0.125d0*sq*t1*zl
         dN(3, 2) = 0.125d0*sq*tl*z1
         N(   3) = 0.125d0*sq*tp*zl
         dN(1, 3) = 0.125d0*s3*tp*zl
         dN(2, 3) = 0.125d0*sq*t3*zl
         dN(3, 3) = 0.125d0*sq*tp*z1
         N(   4) = 0.125d0*sl*tp*zl
         dN(1, 4) = 0.125d0*s1*tp*zl
         dN(2, 4) = 0.125d0*sl*t3*zl
         dN(3, 4) = 0.125d0*sl*tp*z1
         N(   5) = 0.125d0*sl*tl*zp
         dN(1, 5) = 0.125d0*s1*tl*zp
         dN(2, 5) = 0.125d0*sl*t1*zp
         dN(3, 5) = 0.125d0*sl*tl*z3
         N(   6) = 0.125d0*sq*tl*zp
         dN(1, 6) = 0.125d0*s3*tl*zp
         dN(2, 6) = 0.125d0*sq*t1*zp
         dN(3, 6) = 0.125d0*sq*tl*z3
         N(   7) = 0.125d0*sq*tp*zp
         dN(1, 7) = 0.125d0*s3*tp*zp
         dN(2, 7) = 0.125d0*sq*t3*zp
         dN(3, 7) = 0.125d0*sq*tp*z3
         N(   8) = 0.125d0*sl*tp*zp
         dN(1, 8) = 0.125d0*s1*tp*zp
         dN(2, 8) = 0.125d0*sl*t3*zp
         dN(3, 8) = 0.125d0*sl*tp*z3
         N(   9) = 0.25d0*s2*tl*zl
         dN(1, 9) = 0.25d0*s4*tl*zl
         dN(2, 9) = 0.25d0*s2*t1*zl
         dN(3, 9) = 0.25d0*s2*tl*z1
         N(  10) = 0.25d0*sq*t2*zl
         dN(1,10) = 0.25d0*s3*t2*zl
         dN(2,10) = 0.25d0*sq*t4*zl
         dN(3,10) = 0.25d0*sq*t2*z1
         N(  11) = 0.25d0*s2*tp*zl
         dN(1,11) = 0.25d0*s4*tp*zl
         dN(2,11) = 0.25d0*s2*t3*zl
         dN(3,11) = 0.25d0*s2*tp*z1
         N(  12) = 0.25d0*sl*t2*zl
         dN(1,12) = 0.25d0*s1*t2*zl
         dN(2,12) = 0.25d0*sl*t4*zl
         dN(3,12) = 0.25d0*sl*t2*z1
         N(  17) = 0.25d0*sl*tl*z2
         dN(1,17) = 0.25d0*s1*tl*z2
         dN(2,17) = 0.25d0*sl*t1*z2
         dN(3,17) = 0.25d0*sl*tl*z4
         N(  18) = 0.25d0*sq*tl*z2
         dN(1,18) = 0.25d0*s3*tl*z2
         dN(2,18) = 0.25d0*sq*t1*z2
         dN(3,18) = 0.25d0*sq*tl*z4
         N(  19) = 0.25d0*sq*tp*z2
         dN(1,19) = 0.25d0*s3*tp*z2
         dN(2,19) = 0.25d0*sq*t3*z2
         dN(3,19) = 0.25d0*sq*tp*z4
         N(  20) = 0.25d0*sl*tp*z2
         dN(1,20) = 0.25d0*s1*tp*z2
         dN(2,20) = 0.25d0*sl*t3*z2
         dN(3,20) = 0.25d0*sl*tp*z4
         N(  13) = 0.25d0*s2*tl*zp
         dN(1,13) = 0.25d0*s4*tl*zp
         dN(2,13) = 0.25d0*s2*t1*zp
         dN(3,13) = 0.25d0*s2*tl*z3
         N(  14) = 0.25d0*sq*t2*zp
         dN(1,14) = 0.25d0*s3*t2*zp
         dN(2,14) = 0.25d0*sq*t4*zp
         dN(3,14) = 0.25d0*sq*t2*z3
         N(  15) = 0.25d0*s2*tp*zp
         dN(1,15) = 0.25d0*s4*tp*zp
         dN(2,15) = 0.25d0*s2*t3*zp
         dN(3,15) = 0.25d0*s2*tp*z3
         N(  16) = 0.25d0*sl*t2*zp
         dN(1,16) = 0.25d0*s1*t2*zp
         dN(2,16) = 0.25d0*sl*t4*zp
         dN(3,16) = 0.25d0*sl*t2*z3
         N(  25) = 0.5d0*s2*t2*zl
         dN(1,25) = 0.5d0*s4*t2*zl
         dN(2,25) = 0.5d0*s2*t4*zl
         dN(3,25) = 0.5d0*s2*t2*z1
         N(  23) = 0.5d0*s2*tl*z2
         dN(1,23) = 0.5d0*s4*tl*z2
         dN(2,23) = 0.5d0*s2*t1*z2
         dN(3,23) = 0.5d0*s2*tl*z4
         N(  22) = 0.5d0*sq*t2*z2
         dN(1,22) = 0.5d0*s3*t2*z2
         dN(2,22) = 0.5d0*sq*t4*z2
         dN(3,22) = 0.5d0*sq*t2*z4
         N(  24) = 0.5d0*s2*tp*z2
         dN(1,24) = 0.5d0*s4*tp*z2
         dN(2,24) = 0.5d0*s2*t3*z2
         dN(3,24) = 0.5d0*s2*tp*z4
         N(  21) = 0.5d0*sl*t2*z2
         dN(1,21) = 0.5d0*s1*t2*z2
         dN(2,21) = 0.5d0*sl*t4*z2
         dN(3,21) = 0.5d0*sl*t2*z4
         N(  26) = 0.5d0*s2*t2*zp
         dN(1,26) = 0.5d0*s4*t2*zp
         dN(2,26) = 0.5d0*s2*t4*zp
         dN(3,26) = 0.5d0*s2*t2*z3
         N(  27) = s2*t2*z2
         dN(1,27) = s4*t2*z2
         dN(2,27) = s2*t4*z2
         dN(3,27) = s2*t2*z4

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
