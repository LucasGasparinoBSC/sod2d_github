module jacobian_oper

        contains

                subroutine elem_jacobian(ndime,nnode,nelem,npoin,ngaus,connec,coord,dNgp,wgp,gpvol,He)

                        ! Computes the Jacobian transformation of an element, its determinant and
                        ! inverse, for all Gauss points. Valid for 2D and 3D elements.
                        !3D uses cofactor method to obtain the inverse. Dependent on element coordinates and isopar. derivatives.

                        implicit none

                        integer(4), intent(in)  :: ndime, nnode, nelem, npoin, ngaus
                        integer(4), intent(in)  :: connec(nelem,nnode)
                        real(8),    intent(in)  :: coord(npoin,ndime), dNgp(ndime,nnode,ngaus), wgp(ngaus)
                        real(8),    intent(out) :: gpvol(1,ngaus,nelem), He(ndime,ndime,ngaus,nelem)
                        integer(4)              :: idime, jdime, inode, ielem, igaus
                        real(8)                 :: Je(ndime,ndime), a(9), b(9)

                        !
                        ! Initialize He and gpvol
                        !
                        !$acc kernels
                        He(:,:,:,:) = 0.0d0
                        gpvol(:,:,:) = 0.0d0
                        !$acc end kernels

                        !
                        ! Loop over elements
                        !
                        !$acc parallel loop gang private(Je,a,b) vector_length(32)
                        do ielem = 1,nelem
                           !
                           ! Loop over Gauss points
                           !
                           !$acc loop seq
                           do igaus = 1,ngaus
                              !
                              ! Compute Je at each GP for each elem
                              !
                              !$acc loop vector collapse(2)
                              do idime = 1,ndime
                                 do jdime = 1,ndime
                                    Je(idime,jdime) = 0.0d0
                                 end do
                              end do
                              !$acc loop vector collapse(2)
                              do idime = 1,ndime
                                 do jdime = 1,ndime
                                    Je(idime,jdime) = Je(idime,jdime)+dot_product(dNgp(jdime,:,igaus),coord(connec(ielem,:),idime))
                                 end do
                              end do
                              !
                              ! Compute inverse Jacobian He and gpvol = wgp*detJe
                              !
                              if (ndime == 2) then
                                 gpvol(1,igaus,ielem) = Je(1,1)*Je(2,2)-Je(2,1)*Je(1,2)
                                 He(1,1,igaus,ielem) = Je(2,2)
                                 He(2,2,igaus,ielem) = Je(1,1)
                                 He(1,2,igaus,ielem) = -Je(1,2)
                                 He(2,1,igaus,ielem) = -Je(2,1)
                                 !$acc loop vector collapse(2)
                                 do idime = 1,ndime
                                    do jdime = 1,ndime
                                       He(idime,jdime,igaus,ielem) = (1.0d0/gpvol(1,igaus,ielem))*He(idime,jdime,igaus,ielem)
                                    end do
                                 end do
                                 gpvol(1,igaus,ielem) = wgp(igaus)*gpvol(1,igaus,ielem)
                              else if (ndime == 3) then
                                 gpvol(1,igaus,ielem) = (Je(1,1)*Je(2,2)*Je(3,3)+ &
                                       Je(1,2)*Je(2,3)*Je(3,1)+Je(1,3)*Je(2,1)*Je(3,2)-Je(3,1)*Je(2,2)*Je(1,3)- &
                                       Je(3,2)*Je(2,3)*Je(1,1)-Je(3,3)*Je(2,1)*Je(1,2))
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
                                 !$acc loop seq
                                 do idime = 1,9
                                    b(idime) = a(idime)
                                 end do
                                 b(2) = a(4)
                                 b(3) = a(7)
                                 b(4) = a(2)
                                 b(6) = a(8)
                                 b(7) = a(3)
                                 b(8) = a(6)

                                 !
                                 ! Divide by detJe
                                 !
                                 !$acc loop seq
                                 do idime = 1,9
                                    b(idime) = (1.0d0/gpvol(1,igaus,ielem))*b(idime)
                                 end do

                                 !
                                 ! Organize into inverse
                                 !
                                 He(1,1,igaus,ielem) = b(1)
                                 He(1,2,igaus,ielem) = b(2)
                                 He(1,3,igaus,ielem) = b(3)
                                 He(2,1,igaus,ielem) = b(4)
                                 He(2,2,igaus,ielem) = b(5)
                                 He(2,3,igaus,ielem) = b(6)
                                 He(3,1,igaus,ielem) = b(7)
                                 He(3,2,igaus,ielem) = b(8)
                                 He(3,3,igaus,ielem) = b(9)
                                 gpvol(1,igaus,ielem) = wgp(igaus)*gpvol(1,igaus,ielem)
                              end if
                           end do
                        end do
                        !$acc end parallel loop

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
