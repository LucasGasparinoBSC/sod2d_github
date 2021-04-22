module mass_matrix

      ! TODO: create in sparse format

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
                         Mc(ind,ind) = Mc(ind,ind)+Me
                      end do

              end subroutine consistent_mass

              subroutine lumped_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Ml)

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(out) :: Ml(npoin)
                      integer(4)              :: ielem, igaus, inode, jnode, ind(nnode)
                      real(8)                 :: Me(nnode), alpha, aux1, aux2

                      Ml = 0.0d0
                      do ielem = 1,nelem
                         ind(1:nnode) = connec(ielem,1:nnode)
                         Me = 0.0d0
                         alpha = 0.0d0
                         aux1 = 0.0d0
                         aux2 = 0.0d0
                         !
                         ! tr[Mc]
                         !
                         do inode = 1,nnode
                            do igaus = 1,ngaus
                               aux1 = aux1+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*Ngp(igaus,inode)
                            end do
                         end do
                         !
                         ! elem. mass
                         !
                         do igaus = 1,ngaus
                            aux2 = aux2+gpvol(1,igaus,ielem)
                         end do
                         !
                         ! alpha
                         !
                         alpha = aux2/aux1
                         !
                         ! Me
                         !
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               Me(inode) = Me(inode)+gpvol(1,igaus,ielem)*(Ngp(igaus,inode)**2)
                            end do
                         end do
                         Me = alpha*Me
                         Ml(ind) = Ml(ind)+Me
                      end do

              end subroutine lumped_mass

end module mass_matrix
