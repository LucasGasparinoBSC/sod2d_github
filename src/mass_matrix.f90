module mass_matrix

      contains

              subroutine consistent_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Mc)

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(out) :: Mc(npoin,npoin)
                      integer(4)              :: ielem, igaus, inode, jnode, ind(nnode)
                      real(8)                 :: Me(nnode,nnode)

                      Mc = 0.0d0
                      do ielem = 1,nelem
                         ind(1:nnode) = connec(ielem,1:nnode)
                         Me = 0.0d0
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do jnode = 1,nnode
                                  Me(inode,jnode) = Me(inode,jnode) + &
                                     gpvol(1,igaus,ielem)*Ngp(igaus,inode)*Ngp(igaus,jnode)
                               end do
                            end do
                         end do
                      end do

              end subroutine consistent_mass

              subroutine lumped_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Ml)

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(out) :: Ml(npoin)
                      integer(4)              :: ielem, igaus, inode, jnode, ind(nnode)
                      real(8)                 :: Me(nnode)

                      Ml = 0.0d0
                      do ielem = 1,nelem
                         ind(1:nnode) = connec(ielem,1:nnode)
                         Me = 0.0d0
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do jnode = 1,nnode
                                  Me(inode) = Me(inode) + &
                                     (gpvol(1,igaus,ielem)/(Ngp(igaus,jnode)**2)) * &
                                     (Ngp(igaus,inode)**2)
                               end do
                            end do
                         end do
                         Ml(ind) = Ml(ind)+Me
                      end do

              end subroutine lumped_mass

end module mass_matrix
