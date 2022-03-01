module elem_diffu

      use mod_nvtx

      ! TODO: Create unit tests for all subroutines

      contains

              subroutine mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,rho,mu_e,Rmass)

                      ! TODO: Add stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin), mu_e(nelem)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), nu_e
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Mass diffusion")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(gpcar,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         nu_e = mu_e(ielem)/maxval(abs(rho(connec(ielem,:))))
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = dot_product(gpcar(idime,:),rho(connec(ielem,:)))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                              gpcar(idime,inode)*tmp1
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+nu_e*1.0d0*Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mass_diffusion

              !
              ! Old routine
              !
              !!oldsubroutine mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,mu_e,Rmom)

              !!old        ! TODO: Add. stab. viscosity

              !!old        implicit none

              !!old        integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
              !!old        integer(4), intent(in)  :: connec(nelem,nnode)
              !!old        real(8),    intent(in)  :: Ngp(ngaus,nnode)
              !!old        real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
              !!old        real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
              !!old        real(8),    intent(in)  :: u(npoin,ndime), mu_e(nelem)
              !!old        real(8),    intent(out) :: Rmom(npoin,ndime)
              !!old        integer(4)              :: ind(nnode)
              !!old        integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
              !!old        real(8)                 :: Re(nnode,ndime)

              !!old        !$acc kernels
              !!old        Rmom(:,:) = 0.0d0
              !!old        !$acc end kernels
              !!old        call nvtxStartRange("Momentum diffusion")
              !!old        !$acc parallel loop gang private(ind,Re,gpcar) vector_length(32)
              !!old        do ielem = 1,nelem
              !!old           !$acc loop vector
              !!old           do inode = 1,nnode
              !!old              Re(inode,:) = 0.0d0
              !!old              ind(inode) = connec(ielem,inode)
              !!old           end do
              !!old           !$acc loop seq
              !!old           do igaus = 1,ngaus
              !!old              !$acc loop seq
              !!old              do idime = 1,ndime
              !!old                 !$acc loop vector
              !!old                 do inode = 1,nnode
              !!old                    gpcar(idime,inode) = dot_product(He(idime,:.igaus,ielem),dNgp(:,inode,igaus))
              !!old                 end do
              !!old              end do
              !!old           end do
              !!old           !
              !!old           ! Assembly
              !!old           !
              !!old           do idime = 1,ndime
              !!old              Rmom(ind,idime) = Rmom(ind,idime)+Re(1:nnode,idime)
              !!old           end do
              !!old        end do
              !!old        !$acc end parallel loop
              !!old        call nvtxEndRange

              !!oldend subroutine mom_diffusion

              !
              ! New routine
              !
              subroutine mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u,mu_e,Rmom)

                      ! TODO: Add. stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), mu_e(nelem)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode, kdime
                      real(8)                 :: Re(nnode,ndime), twoThirds, gpcar(ndime,nnode)
                      real(8) :: grad_1, grad_2, grad_3, grad_4, grad_5, grad_6, grad_7, grad_8, grad_9
                      real(8) :: div_1
                      real(8) :: tau_1, tau_2, tau_3, tau_4, tau_5, tau_6, tau_7, tau_8, tau_9
                      real(8) :: tmp1, tmp2, tmp3

                      twoThirds = 2.0d0/3.0d0
                      call nvtxStartRange("Momentum diffusion")
                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               Re(inode,idime) = 0.0d0
                            end do
                         end do
                         !
                         ! Gradient structure:
                         !
                         !         | u1,1 u1,2 u1,3 |
                         ! u_i,j = | u2,1 u2,2 u2,3 |
                         !         | u3,1 u3,2 u3,3 |
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            grad_1=0.0d0
                            grad_2=0.0d0
                            grad_3=0.0d0
                            grad_4=0.0d0
                            grad_5=0.0d0
                            grad_6=0.0d0
                            grad_7=0.0d0
                            grad_8=0.0d0
                            grad_9=0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !
                            ! compute dN_b,j*u_bi @ Gauss point k
                            !
                            !$acc loop vector &
                            !$acc reduction(+:grad_1,grad_2,grad_3,grad_4,grad_5,grad_6,grad_7,grad_8,grad_9)
                            do inode = 1,nnode
                               grad_1 = grad_1+gpcar(1,inode)*u(connec(ielem,inode),1)
                               grad_2 = grad_2+gpcar(2,inode)*u(connec(ielem,inode),1)
                               grad_3 = grad_3+gpcar(3,inode)*u(connec(ielem,inode),1)
                               grad_4 = grad_4+gpcar(1,inode)*u(connec(ielem,inode),2)
                               grad_5 = grad_5+gpcar(2,inode)*u(connec(ielem,inode),2)
                               grad_6 = grad_6+gpcar(3,inode)*u(connec(ielem,inode),2)
                               grad_7 = grad_7+gpcar(1,inode)*u(connec(ielem,inode),3)
                               grad_8 = grad_8+gpcar(2,inode)*u(connec(ielem,inode),3)
                               grad_9 = grad_9+gpcar(3,inode)*u(connec(ielem,inode),3)
                            end do
                            div_1 = grad_1+grad_5+grad_9 ! u_k,k = tr[grad(u)]
                            !
                            ! Compute tau_ij = mu*(u_i,j+u_j,i-(2/3)*u_k,k*d_ij)
                            !
                            tau_1 = mu_e(ielem)*(grad_1+grad_1-twoThirds*div_1)
                            tau_2 = mu_e(ielem)*(grad_2+grad_4)
                            tau_3 = mu_e(ielem)*(grad_3+grad_7)
                            tau_4 = mu_e(ielem)*(grad_4+grad_2)
                            tau_5 = mu_e(ielem)*(grad_5+grad_5-twoThirds*div_1)
                            tau_6 = mu_e(ielem)*(grad_6+grad_8)
                            tau_7 = mu_e(ielem)*(grad_7+grad_3)
                            tau_8 = mu_e(ielem)*(grad_8+grad_6)
                            tau_9 = mu_e(ielem)*(grad_9+grad_9-twoThirds*div_1)
                            !
                            ! Compute N_a,j*tau_ij @ Gauss point k
                            !
                            tmp1 = 0.0d0
                            tmp2 = 0.0d0
                            tmp3 = 0.0d0
                            !$acc loop vector reduction(+:tmp1,tmp2,tmp3)
                            do inode = 1,nnode
                               tmp1 = tmp1+(gpvol(1,igaus,ielem)*gpcar(1,inode)*tau_1)
                               tmp1 = tmp1+(gpvol(1,igaus,ielem)*gpcar(2,inode)*tau_2)
                               tmp1 = tmp1+(gpvol(1,igaus,ielem)*gpcar(3,inode)*tau_3)
                               tmp2 = tmp1+(gpvol(1,igaus,ielem)*gpcar(1,inode)*tau_4)
                               tmp2 = tmp1+(gpvol(1,igaus,ielem)*gpcar(2,inode)*tau_5)
                               tmp2 = tmp1+(gpvol(1,igaus,ielem)*gpcar(3,inode)*tau_6)
                               tmp3 = tmp1+(gpvol(1,igaus,ielem)*gpcar(1,inode)*tau_7)
                               tmp3 = tmp1+(gpvol(1,igaus,ielem)*gpcar(2,inode)*tau_8)
                               tmp3 = tmp1+(gpvol(1,igaus,ielem)*gpcar(3,inode)*tau_9)
                            end do
                            !
                            ! Gaussian quadrature
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode,1) = Re(inode,1)+tmp1
                               Re(inode,2) = Re(inode,2)+tmp2
                               Re(inode,3) = Re(inode,3)+tmp3
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               !$acc atomic update
                               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+1.0d0*Re(inode,idime)
                               !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mom_diffusion

              subroutine ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u,Tem,mu_e,Rener)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), Tem(npoin), mu_e(nelem)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), kappa_e
                      real(8)                 :: el_Ke(nnode)
                      real(8)                 :: gradT, gradKe, gpcar(ndime,nnode)

                      call nvtxStartRange("Energy diffusion")
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(el_Ke,Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         kappa_e = mu_e(ielem)*1004.0d0/0.72d0 ! Fixed Cp and Pr
                         !kappa_e = mu_e(ielem)/(1.40d0-1.0d0)
                         !$acc loop vector
                         do inode = 1,nnode
                            el_Ke(inode) = dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))/2.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               gradT = kappa_e*dot_product(gpcar(idime,:),Tem(connec(ielem,:)))
                               gradKe = mu_e(ielem)*dot_product(gpcar(idime,:),el_Ke(:))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                     gpcar(idime,inode)*(gradT+gradKe)
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+1.0d0*Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine ener_diffusion

end module elem_diffu
