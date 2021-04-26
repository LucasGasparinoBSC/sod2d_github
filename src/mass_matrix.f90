module mass_matrix

      contains

              subroutine consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Forms Mc as a sparse CSR matrix, utilizing nzdom, rdom and cdom to          !
                      ! compress the elemental matrices into the full sparse assembly structure.    !
                      ! Mc is defined by A{int(Na*Nb*det(Je))}, where A{} is the assembly operator. !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus, nzdom
                      integer(4), intent(in)  :: connec(nelem,nnode), rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(out) :: Mc(nzdom)
                      integer(4)              :: ielem, igaus, inode, jnode, lnode(nnode), izdom, ipoin, jpoin
                      integer(4)              :: jzdom, rowb, rowe
                      real(8)                 :: Me(nnode,nnode)

                      !
                      ! Initialize Mc to zeros
                      !
                      Mc = 0.0d0

                      !
                      ! Loop over all elements to form Mc_e(nnode,nnode)
                      !
                      do ielem = 1,nelem
                         lnode(1:nnode) = connec(ielem,1:nnode) ! get elemental indices
                         !
                         ! Form Mc_e with Gaussian quadrature (open)
                         !
                         Me = 0.0d0
                         do igaus = 1,ngaus ! Loop over Gauss points
                            do inode = 1,nnode ! Loop over element nodes (row)
                               do jnode = 1,nnode ! Loop over element nodex (column)
                                  Me(inode,jnode) = Me(inode,jnode) + &
                                     gpvol(1,igaus,ielem)*Ngp(igaus,inode)*Ngp(igaus,jnode) ! Gaussian quad.
                               end do
                            end do
                         end do
                         !
                         ! Assemble Mc_e to CSR Mc
                         !
                         do inode = 1,nnode
                            ipoin = lnode(inode) ! Global node/Mc row index
                            rowb = rdom(ipoin)+1 ! Start of izdom for cdom
                            rowe = rdom(ipoin+1) ! end of izdom for cdom
                            do jnode = 1,nnode
                               jpoin = lnode(jnode) ! Nodes associated with ipoin on ielem
                               !
                               ! Loop over section of cdom to find out izdom
                               !
                               do jzdom = rowb,rowe
                                  if (cdom(jzdom) == jpoin) then
                                     izdom = jzdom
                                     exit
                                  end if
                               end do
                               Mc(izdom) = Mc(izdom) + Me(inode,jnode)
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
