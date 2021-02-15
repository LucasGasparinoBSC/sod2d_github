module elem_diffu

      contains

              subroutine mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,rho,Rmass)

                      ! TODO: Add stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), el_rho(nnode)

                      Rmass = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_rho(1:nnode) = rho(ind)
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*gpcar(idime,inode,igaus,ielem)* &
                                             gpcar(idime,jnode,igaus,ielem)*el_rho(jnode)
                                  end do
                               end do
                            end do
                         end do
                         Rmass(ind) = Rmass(ind)+Re(1:nnode)
                      end do

              end subroutine mass_diffusion

              subroutine mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,Rmom)

                      ! TODO: Add. stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode
                      real(8)                 :: Re(nnode,ndime), aux(ndime,ngaus)
                      real(8)                 :: el_u(nnode,ndime), grad_u(ndime,ndime,ngaus)

                      Rmom = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         grad_u = 0.0d0
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  do jnode = 1,nnode
                                     grad_u(idime,jdime,igaus) = grad_u(idime,jdime,igaus) + &
                                             gpcar(jdime,jnode,igaus,ielem)*el_u(jnode,idime)
                                  end do
                                  do inode = 1,nnode
                                     Re(inode,idime) = Re(inode,idime) + gpvol(1,igaus,ielem) * &
                                             gpcar(jdime,inode,igaus,ielem)*grad_u(idime,jdime,igaus)
                                  end do
                               end do
                            end do
                         end do
                         Rmom(ind,1) = Rmom(ind,1)+Re(1:nnode,1)
                         Rmom(ind,2) = Rmom(ind,2)+Re(1:nnode,2)
                      end do

              end subroutine mom_diffusion

              subroutine ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,Tem,Rener)

                      ! TODO: Add stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), Tem(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode)
                      real(8)                 :: el_u(nnode,ndime), el_Tem(nnode), el_Ke(nnode)
                      real(8)                 :: grad_T(ndime,ngaus), grad_Ke(ndime,ngaus)

                      Rener = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_Tem(1:nnode) = Tem(ind)
                         !
                         ! Ke
                         !
                         do inode = 1,nnode
                            el_Ke(inode) = dot_product(el_u(inode,:),el_u(inode,:))
                         end do
                         grad_T = 0.0d0
                         grad_Ke = 0.0d0
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               do jnode = 1,nnode
                                  grad_T(idime,igaus) = grad_T(idime,igaus) + &
                                          gpcar(idime,jnode,igaus,ielem)*el_Tem(jnode)
                                  grad_Ke(idime,igaus) = grad_Ke(idime,igaus) + &
                                          gpcar(idime,jnode,igaus,ielem)*el_Ke(jnode)
                               end do
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem) * &
                                          (gpcar(idime,inode,igaus,ielem)*grad_T(idime,igaus) + &
                                          gpcar(idime,inode,igaus,ielem)*grad_Ke(idime,igaus))
                               end do
                            end do
                         end do
                         Rener(ind) = Rener(ind)+Re(1:nnode)
                      end do

              end subroutine ener_diffusion

end module elem_diffu