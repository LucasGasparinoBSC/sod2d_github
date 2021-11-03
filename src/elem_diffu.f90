module elem_diffu

      use mod_nvtx

      ! TODO: Create unit tests for all subroutines

      contains

              subroutine mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,rho,mu_e,Rmass)

                      ! TODO: Add stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin), mu_e(nelem)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), nu_e
                      real(8)                 :: tmp1(ndime), tmp2

                      Rmass = 0.0d0
                      call nvtxStartRange("Mass diffusion")
                      !$acc parallel loop gang private(ind,tmp1,Re)
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         nu_e = mu_e(ielem)/maxval(abs(rho(ind)))
                         !!$acc loop vector collapse(2)
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop vector
                            do idime = 1,ndime
                               tmp1(idime) = dot_product(gpcar(idime,:,igaus,ielem),rho(ind))
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               tmp2 = dot_product(gpcar(:,inode,igaus,ielem),tmp1(:))
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*tmp2
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(ind(inode)) = Rmass(ind(inode))+nu_e*Re(inode)
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
              !!old        integer(4)              :: ielem, igaus, idime, jdime, inode, jnode, kdime
              !!old        real(8)                 :: Re(nnode,ndime), aux(ndime,ngaus), tau(ndime,ndime,ngaus)
              !!old        real(8)                 :: el_u(nnode,ndime), grad_u(ndime,ndime,ngaus), div_u(ndime,ndime,ngaus)

              !!old        Rmom = 0.0d0
              !!old        call nvtxStartRange("Momentum diffusion")
              !!old        do ielem = 1,nelem
              !!old           Re = 0.0d0
              !!old           ind = connec(ielem,:)
              !!old           el_u(1:nnode,1:ndime) = u(ind,1:ndime)
              !!old           grad_u = 0.0d0
              !!old           div_u = 0.0d0
              !!old           do igaus = 1,ngaus
              !!old              do idime = 1,ndime
              !!old                 !
              !!old                 ! grad_u
              !!old                 !
              !!old                 do jnode = 1,nnode
              !!old                    do jdime = 1,ndime
              !!old                       grad_u(idime,jdime,igaus) = grad_u(idime,jdime,igaus) + &
              !!old                               gpcar(jdime,jnode,igaus,ielem)*el_u(jnode,idime)
              !!old                    end do
              !!old                 end do
              !!old                 !
              !!old                 ! div_u
              !!old                 !
              !!old                 do jnode = 1,nnode
              !!old                    do kdime = 1,ndime
              !!old                       div_u(idime,idime,igaus) = div_u(idime,idime,igaus) + &
              !!old                               gpcar(kdime,jnode,igaus,ielem)*el_u(jnode,kdime)
              !!old                    end do
              !!old                 end do
              !!old                 !
              !!old                 ! tau_ij = grad_u + grad_u^T - (2/3)*div_u
              !!old                 !
              !!old                 do jdime = 1,ndime
              !!old                    tau(idime,jdime,igaus) = 1.0d0*mu_e(ielem)*(grad_u(idime,jdime,igaus) + &
              !!old                            grad_u(jdime,idime,igaus) - &
              !!old                            (2.0d0/3.0d0)*div_u(idime,jdime,igaus))
              !!old                 end do
              !!old              end do
              !!old              !
              !!old              ! div(tau)
              !!old              !
              !!old              do idime = 1,ndime
              !!old                 do inode = 1,nnode
              !!old                    do jdime = 1,ndime
              !!old                       Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem) * &
              !!old                               gpcar(jdime,inode,igaus,ielem)*tau(idime,jdime,igaus)
              !!old                    end do
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
              !!old        call nvtxEndRange

              !!oldend subroutine mom_diffusion

              !
              ! New routine
              !
              subroutine mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,mu_e,Rmom)

                      ! TODO: Add. stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), mu_e(nelem)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode, kdime
                      !real(8)                 :: tau(ndime,ndime), tmp
                      !real(8)                 :: grad_u(ndime,ndime), div_u(ndime,ndime), tmp_tau(ngaus)
                      real(8)                  :: Re(nnode,ndime), twoThirds

                      real(8) :: grad_1, grad_2, grad_3, grad_4, grad_5, grad_6, grad_7, grad_8, grad_9
                      real(8) :: div_1, div_2, div_3, div_4, div_5, div_6, div_7, div_8, div_9
                      real(8) :: tau_1, tau_2, tau_3, tau_4, tau_5, tau_6, tau_7, tau_8, tau_9
                      real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

                      twoThirds = 2.0d0/3.0d0
                      Rmom = 0.0d0
                      call nvtxStartRange("Momentum diffusion")
                      !$acc parallel loop gang private(ind,Re)
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         !$acc loop vector
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
                            div_1=0.0d0
                            div_2=0.0d0
                            div_3=0.0d0
                            div_4=0.0d0
                            div_5=0.0d0
                            div_6=0.0d0
                            div_7=0.0d0
                            div_8=0.0d0
                            div_9=0.0d0
                            tau_1=0.0d0
                            tau_2=0.0d0
                            tau_3=0.0d0
                            tau_4=0.0d0
                            tau_5=0.0d0
                            tau_6=0.0d0
                            tau_7=0.0d0
                            tau_8=0.0d0
                            tau_9=0.0d0
                            !$acc loop seq
                            do inode = 1,nnode
                               grad_1 = grad_1+gpcar(1,inode,igaus,ielem)*u(ind(inode),1)
                               grad_2 = grad_2+gpcar(2,inode,igaus,ielem)*u(ind(inode),1)
                               grad_3 = grad_3+gpcar(3,inode,igaus,ielem)*u(ind(inode),1)
                               grad_4 = grad_4+gpcar(1,inode,igaus,ielem)*u(ind(inode),2)
                               grad_5 = grad_5+gpcar(2,inode,igaus,ielem)*u(ind(inode),2)
                               grad_6 = grad_6+gpcar(3,inode,igaus,ielem)*u(ind(inode),2)
                               grad_7 = grad_7+gpcar(1,inode,igaus,ielem)*u(ind(inode),3)
                               grad_8 = grad_8+gpcar(2,inode,igaus,ielem)*u(ind(inode),3)
                               grad_9 = grad_9+gpcar(3,inode,igaus,ielem)*u(ind(inode),3)

                               div_1 = div_1+gpcar(1,inode,igaus,ielem)*u(ind(inode),1)
                               div_5 = div_5+gpcar(2,inode,igaus,ielem)*u(ind(inode),2)
                               div_9 = div_9+gpcar(3,inode,igaus,ielem)*u(ind(inode),3)
                            end do
                            tmp1 = 0.0d0
                            tmp2 = 0.0d0
                            tmp3 = 0.0d0
                            !$acc loop vector reduction(+:tmp1,tmp2,tmp3)
                            do inode = 1,nnode
                               tau_1 = mu_e(ielem)*(grad_1+grad_1-twoThirds*div_1)
                               tau_2 = mu_e(ielem)*(grad_2+grad_4-twoThirds*div_2)
                               tau_3 = mu_e(ielem)*(grad_3+grad_7-twoThirds*div_3)
                               tau_4 = mu_e(ielem)*(grad_4+grad_2-twoThirds*div_4)
                               tau_5 = mu_e(ielem)*(grad_5+grad_5-twoThirds*div_5)
                               tau_6 = mu_e(ielem)*(grad_6+grad_8-twoThirds*div_6)
                               tau_7 = mu_e(ielem)*(grad_7+grad_3-twoThirds*div_7)
                               tau_8 = mu_e(ielem)*(grad_8+grad_6-twoThirds*div_8)
                               tau_9 = mu_e(ielem)*(grad_9+grad_9-twoThirds*div_9)
                               tmp1 = tmp1+(gpvol(1,igaus,idime)*gpcar(1,inode,igaus,ielem)*tau_1)
                               tmp1 = tmp1+(gpvol(1,igaus,idime)*gpcar(2,inode,igaus,ielem)*tau_2)
                               tmp1 = tmp1+(gpvol(1,igaus,idime)*gpcar(3,inode,igaus,ielem)*tau_3)
                               tmp2 = tmp2+(gpvol(1,igaus,idime)*gpcar(1,inode,igaus,ielem)*tau_4)
                               tmp2 = tmp2+(gpvol(1,igaus,idime)*gpcar(2,inode,igaus,ielem)*tau_5)
                               tmp2 = tmp2+(gpvol(1,igaus,idime)*gpcar(3,inode,igaus,ielem)*tau_6)
                               tmp3 = tmp3+(gpvol(1,igaus,idime)*gpcar(1,inode,igaus,ielem)*tau_7)
                               tmp3 = tmp3+(gpvol(1,igaus,idime)*gpcar(2,inode,igaus,ielem)*tau_8)
                               tmp3 = tmp3+(gpvol(1,igaus,idime)*gpcar(3,inode,igaus,ielem)*tau_9)
                            end do
                            !$acc loop seq
                            do inode = 1,nnode
                               Re(inode,1) = Re(inode,1)+tmp1
                               Re(inode,2) = Re(inode,2)+tmp2
                               Re(inode,3) = Re(inode,3)+tmp3
                            end do
                         end do
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               !$acc atomic update
                               Rmom(ind(inode),idime) = Rmom(ind(inode),idime)+Re(inode,idime)
                               !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mom_diffusion

              subroutine ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,Tem,mu_e,Rener)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), Tem(npoin), mu_e(nelem)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), kappa_e
                      real(8)                 :: el_u(nnode,ndime), el_Tem(nnode), el_Ke(nnode)
                      real(8)                 :: grad_T(ndime,ngaus), grad_Ke(ndime,ngaus)
                      real(8)                 :: gradT, gradKe

                      Rener = 0.0d0
                      call nvtxStartRange("Energy diffusion")
                      !$acc parallel loop gang private(ind,grad_T,grad_Ke,el_Ke,Re)
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         kappa_e = mu_e(ielem)*1004.0d0/0.72d0 ! Fixed Cp and Pr
                         !kappa_e = mu_e(ielem)/(1.40d0-1.0d0)
                         !
                         ! Ke
                         !
                         !$acc loop vector
                         do inode = 1,nnode
                            el_Ke(inode) = dot_product(u(ind(inode),:),u(ind(inode),:))/2.0d0
                         end do
                         !grad_T = 0.0d0
                         !grad_Ke = 0.0d0
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               gradT = 0.0d0
                               gradKe = 0.0d0
                               !$acc loop vector reduction(+:gradT,gradKe)
                               do jnode = 1,nnode
                                  !grad_T(idime,igaus) = grad_T(idime,igaus) + &
                                  !        gpcar(idime,jnode,igaus,ielem)*el_Tem(jnode)
                                  !grad_Ke(idime,igaus) = grad_Ke(idime,igaus) + &
                                  !        gpcar(idime,jnode,igaus,ielem)*el_Ke(jnode)
                                  gradT = gradT + gpcar(idime,jnode,igaus,ielem)*Tem(ind(jnode))
                                  gradKe = gradKe + gpcar(idime,jnode,igaus,ielem)*el_Ke(jnode)
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem) * &
                                          (gpcar(idime,inode,igaus,ielem)*kappa_e*gradT + &
                                          gpcar(idime,inode,igaus,ielem)*mu_e(ielem)*gradKe)
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rener(ind(inode)) = Rener(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine ener_diffusion

end module elem_diffu
