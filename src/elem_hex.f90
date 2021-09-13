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
                        N = 0.25d0*N

                        dN(1,1) = -1.0d0+t
                        dN(2,1) = -1.0d0+s
                        dN(1,2) =  1.0d0-t
                        dN(2,2) = -1.0d0-s
                        dN(1,3) =  1.0d0+t
                        dN(2,3) =  1.0d0+s
                        dN(1,4) = -1.0d0-t
                        dN(2,4) =  1.0d0-s
                        dN = 0.25d0*dN

                end subroutine hex08

                subroutine hex27(s,t,z,N,dN) ! QUA09 element

                        ! TODO: IMPLEMENT PROPERLY!!!!!

                        implicit none

                        real(8), intent(in)  :: s, t, z
                        real(8), intent(out) :: N(27), dN(3,27)

                        N(1) = (1.0d0-s)*(1.0d0-t)
                        N(2) = (1.0d0+s)*(1.0d0-t)
                        N(3) = (1.0d0+s)*(1.0d0+t)
                        N(4) = (1.0d0-s)*(1.0d0+t)
                        N = 0.25d0*N

                        dN(1,1) = -1.0d0+t
                        dN(2,1) = -1.0d0+s
                        dN(1,2) =  1.0d0-t
                        dN(2,2) = -1.0d0-s
                        dN(1,3) =  1.0d0+t
                        dN(2,3) =  1.0d0+s
                        dN(1,4) = -1.0d0-t
                        dN(2,4) =  1.0d0-s
                        dN = 0.25d0*dN

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
