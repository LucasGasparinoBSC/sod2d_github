module jacobian_oper

        contains

                subroutine elem_jacobian(ndime,nnode,elcod,dN,detJe,He)

                        ! Computes the Jacobian transformation of an element, its determinant and
                        ! inverse. Valid for 2D and 3D elements. 3D uses cofactor method to obtain
                        ! the inverse. Dependent on element coordinates and isopar. derivatives.

                        implicit none

                        integer(4), intent(in)  :: ndime, nnode
                        real(8),    intent(in)  :: elcod(nnode,ndime), dN(ndime,nnode)
                        real(8),    intent(out) :: detJe, He(ndime,ndime)
                        real(8)                 :: Je(ndime,ndime), a(9), b(9)

                        Je = matmul(dN,elcod)
                        if (ndime == 2) then
                           detJe = Je(1,1)*Je(2,2)-Je(2,1)*Je(1,2)
                           He(1,1) = Je(2,2)
                           He(2,2) = Je(1,1)
                           He(1,2) = -Je(1,2)
                           He(2,1) = -Je(2,1)
                           He = (1.0d0/detJe)*He
                        else if (ndime == 3) then
                           detJe = Je(1,1)*Je(2,2)*Je(3,3)+Je(1,2)*Je(2,3)*Je(3,1)+Je(1,3)*Je(2,1)*Je(3,2) - &
                                   Je(3,1)*Je(2,2)*Je(1,3)-Je(3,2)*Je(2,3)*Je(1,1)-Je(3,3)*Je(2,1)*Je(1,2)
                           !
                           ! Minors for inverse
                           !
                           a(1) = Je(2,2)*Je(3,3)-Je(3,2)*Je(2,3)
                           a(2) = Je(2,1)*Je(3,3)-Je(3,1)*Je(2,3)
                           a(3) = Je(2,1)*Je(3,2)-Je(3,1)*Je(2,2)
                           a(4) = Je(1,2)*Je(3,3)-Je(3,2)*Je(1,3)
                           a(5) = Je(1,1)*Je(3,3)-Je(3,1)*Je(1,3)
                           a(6) = Je(1,1)*Je(3,2)-Je(3,1)*Je(1,2)
                           a(7) = Je(1,2)*Je(2,3)-Je(2,2)*Je(1,3)
                           a(8) = Je(1,1)*Je(2,3)-Je(2,1)*Je(1,3)
                           a(9) = Je(1,1)*Je(2,2)-Je(2,1)*Je(1,2)

                           !
                           ! Sign changes
                           !
                           a(2) = -a(2)
                           a(4) = -a(4)
                           a(6) = -a(6)
                           a(8) = -a(8)

                           !
                           ! Transpose a into b
                           !
                           b = a
                           b(2) = a(4)
                           b(3) = a(7)
                           b(4) = a(2)
                           b(6) = a(8)
                           b(7) = a(3)
                           b(8) = a(6)

                           !
                           ! Divide by detJe
                           !
                           b = (1.0d0/detJe)*b

                           !
                           ! Organize into inverse
                           !
                           He(1,1) = b(1)
                           He(1,2) = b(2)
                           He(1,3) = b(3)
                           He(2,1) = b(4)
                           He(2,2) = b(5)
                           He(2,3) = b(6)
                           He(3,1) = b(7)
                           He(3,2) = b(8)
                           He(3,3) = b(9)

                        end if

                end subroutine

                subroutine cartesian_deriv(ndime,nnode,dN,He,dxN)

                        ! Pass the isopar. derivatives to cartesian space.

                        implicit none

                        integer(4), intent(in)  :: ndime,nnode
                        real(8)   , intent(in)  :: dN(ndime,nnode), He(ndime,ndime)
                        real(8)   , intent(out) :: dxN(ndime,nnode)

                        if (ndime ==2) then
                           dxN(1,:) = He(1,1)*dN(1,:)+He(1,2)*dN(2,:)
                           dxN(2,:) = He(2,1)*dN(1,:)+He(2,2)*dN(2,:)
                        else if (ndime == 3) then
                           dxN(1,:) = He(1,1)*dN(1,:)+He(1,2)*dN(2,:)+He(1,3)*dN(3,:)
                           dxN(2,:) = He(2,1)*dN(1,:)+He(2,2)*dN(2,:)+He(2,3)*dN(3,:)
                           dxN(3,:) = He(3,1)*dN(1,:)+He(3,2)*dN(2,:)+He(3,3)*dN(3,:)
                        end if

                end subroutine cartesian_deriv

end module jacobian_oper
