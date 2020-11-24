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
                      integer(4)              :: ielem, igaus, inode, jnode
                      real(8)                 :: Re(nnode), el_q(nnode,ndime)

                      Rmass = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do jnode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                              (gpcar(1,jnode,igaus,ielem)*el_q(jnode,1)+ &
                                              gpcar(2,jnode,igaus,ielem)*el_q(jnode,2))
                               end do
                            end do
                         end do
                         Rmass(ind) = Rmass(ind)+Re(1:nnode)
                      end do

              end subroutine mass_convec

              subroutine mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,q,pr,Rmom)

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
                      real(8)                 :: Re(nnode,ndime)
                      real(8)                 :: el_q(nnode,ndime), el_u(nnode,ndime), el_pr(nnode)

                      Rmom = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_pr(1:nnode) = pr(ind)
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  do inode = 1,nnode
                                     do jnode = 1,nnode
                                        Re(inode,idime) = Re(inode,idime) + &
                                           gpvol(1,igaus,ielem)*Ngp(igaus,inode) * &
                                           (gpcar(jdime,jnode,igaus,ielem) * &
                                           (el_q(jnode,idime)*el_u(jnode,jdime)) + &
                                           gpcar(idime,jnode,igaus,ielem)*el_pr(jnode))
                                     end do
                                  end do
                               end do
                            end do
                         end do
                         Rmom(ind,1:ndime) = Rmom(ind,1:ndime)+Re(1:nnode,1:ndime)
                      end do

              end subroutine mom_convec

end module elem_convec
