module mass_matrix

      contains

              subroutine consistent_mass()
                      implicit none
                      Mc = 0.0d0
                      do ielem = 1,nelem
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

              subroutine lumped_mass()
              end subroutine lumped_mass

end module mass_matrix
