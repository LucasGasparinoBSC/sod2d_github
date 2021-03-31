module jacobian_oper

        ! TODO: make 3D compatible

        contains

                subroutine elem_jacobian(ndime,nnode,elcod,dN,Je,detJe,He)

                        implicit none

                        integer(4), intent(in)  :: ndime, nnode
                        real(8),    intent(in)  :: elcod(ndime,nnode), dN(ndime,nnode)
                        real(8),    intent(out) :: detJe, Je(ndime,ndime), He(ndime,ndime)

                        Je = matmul(elcod,transpose(dN))
                        detJe = Je(1,1)*Je(2,2)-Je(2,1)*Je(1,2)
                        He(1,1) = Je(2,2)
                        He(2,2) = Je(1,1)
                        He(1,2) = -Je(1,2)
                        He(2,1) = -Je(2,1)
                        He = (1.0d0/detJe)*He

                end subroutine

                subroutine cartesian_deriv(ndime,nnode,dN,He,dxN)

                        implicit none

                        integer(4), intent(in)  :: ndime,nnode
                        real(8)   , intent(in)  :: dN(ndime,nnode), He(ndime,ndime)
                        real(8)   , intent(out) :: dxN(ndime,nnode)

                        dxN(1,:) = He(1,1)*dN(1,:)+He(1,2)*dN(2,:)
                        dxN(2,:) = He(2,1)*dN(1,:)+He(2,2)*dN(2,:)

                end subroutine cartesian_deriv

end module jacobian_oper
