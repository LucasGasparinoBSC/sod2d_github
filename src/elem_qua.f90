module elem_qua

        contains

                subroutine qua04(s,t,N,dN)

                        implicit none

                        real(8), intent(in)  :: s, t
                        real(8), intent(out) :: N(4), dN(2,4)

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

                end subroutine qua04

end module
