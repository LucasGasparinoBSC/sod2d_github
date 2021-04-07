module elem_convec

      contains

              subroutine mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rmass)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), el_q(nnode,ndime)

                      Rmass = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_q(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         Rmass(ind) = Rmass(ind)+Re(1:nnode)
                      end do

              end subroutine mass_convec

              subroutine mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,q,pr,Rmom)

                      ! TODO: modify to accept 3D elements as well

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), pr(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode
                      real(8)                 :: Re(nnode,ndime), aux
                      real(8)                 :: el_q(nnode,ndime), el_u(nnode,ndime), el_pr(nnode)

                      Rmom = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_pr(1:nnode) = pr(ind)
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     !
                                     ! Compute div(qu) 
                                     !
                                     do jdime = 1,ndime
                                        aux = gpcar(jdime,jnode,igaus,ielem)*el_q(jnode,idime)*el_u(jnode,jdime)
                                     end do
                                     Re(inode,idime) = gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                             (aux+gpcar(idime,jnode,igaus,ielem)*el_pr(jnode))
                                     !Re(inode,1) = Re(inode,1) + gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                     !   (gpcar(1,jnode,igaus,ielem)*(el_q(jnode,1)*el_u(jnode,1))+ &
                                     !   gpcar(2,jnode,igaus,ielem)*(el_q(jnode,1)*el_u(jnode,2))+ &
                                     !   gpcar(1,jnode,igaus,ielem)*el_pr(jnode))
                                     !Re(inode,2) = Re(inode,2) + gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                     !   (gpcar(1,jnode,igaus,ielem)*(el_q(jnode,2)*el_u(jnode,1))+ &
                                     !   gpcar(2,jnode,igaus,ielem)*(el_q(jnode,2)*el_u(jnode,2))+ &
                                     !   gpcar(2,jnode,igaus,ielem)*el_pr(jnode))
                                  end do
                               end do
                            end do
                         end do
                         do idime = 1,ndime
                            Rmom(ind,idime) = Rmom(ind,idime)+Re(1:nnode,idime)
                         end do
                      end do

              end subroutine mom_convec

              subroutine ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,pr,E,Rener)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), pr(npoin), E(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime, ipoin
                      real(8)                 :: Re(nnode), fener(npoin,ndime)
                      real(8)                 :: el_fener(nnode,ndime)

                      do ipoin = 1,npoin
                         fener(ipoin,1:ndime) = u(ipoin,1:ndime)*(E(ipoin)+pr(ipoin))
                      end do

                      Rener = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_fener(1:nnode,1:ndime) = fener(ind,1:ndime)
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_fener(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         Rener(ind) = Rener(ind)+Re(1:nnode)
                      end do

              end subroutine ener_convec

              subroutine generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rconvec)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(out) :: Rconvec(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), el_q(nnode,ndime)

                      Rconvec = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_q(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         Rconvec(ind) = Rconvec(ind)+Re(1:nnode)
                      end do

              end subroutine generic_scalar_convec

end module elem_convec
