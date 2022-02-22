module mod_entropy_viscosity

   use mod_nvtx

      ! TODO: Finish module and create unit tests

      contains

              subroutine residuals(nelem,ngaus,npoin,nnode,ndime,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rhok, uk, prk, qk, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                      use mass_matrix
                      use mod_solver
                      use elem_convec

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime, ppow, npoin_w
                      integer(4), intent(in)  :: connec(nelem,nnode), lpoin_w(npoin_w)
                      real(8),    intent(in)  :: Ngp(nnode,ngaus), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ml(npoin)
                      real(8),    intent(in)  :: dt, gamma_gas
                      real(8),    intent(in)  :: rhok(npoin), uk(npoin,ndime), prk(npoin), qk(npoin,ndime)     ! From substep
                      real(8),    intent(in)  :: rho(npoin,2), u(npoin,ndime,2), pr(npoin,2), q(npoin,ndime,2) ! From prediction
                      real(8),    intent(out) :: Reta(npoin), Rrho(npoin)
                      integer(4)              :: ipoin, idime
                      real(8)                 :: eta(npoin), eta_p(npoin), alpha(npoin), alpha_p(npoin)
                      real(8)                 :: f_eta(npoin,ndime), f_rho(npoin,ndime), R1(npoin), R2(npoin)
                      real(8)                 :: aux1(npoin)

                       !
                       ! Entropy function and temporal terms
                       !
                       call nvtxStartRange("Entropy transport")
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          eta(lpoin_w(ipoin)) = (rhok(lpoin_w(ipoin))/(gamma_gas-1.0d0))* &
                             log(prk(lpoin_w(ipoin))/(rhok(lpoin_w(ipoin))**gamma_gas))
                          eta_p(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)/(gamma_gas-1.0d0))* &
                             log(pr(lpoin_w(ipoin),1)/(rho(lpoin_w(ipoin),1)**gamma_gas))
                       end do
                       !$acc end parallel loop

                       !$acc parallel loop collapse(2)
                       do ipoin = 1,npoin_w
                          do idime = 1,ndime
                             f_eta(lpoin_w(ipoin),idime) = uk(lpoin_w(ipoin),idime)*eta(lpoin_w(ipoin))
                             f_rho(lpoin_w(ipoin),idime) = qk(lpoin_w(ipoin),idime)
                          end do
                       end do

                       !
                       ! Temporal eta
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          R1(lpoin_w(ipoin)) = (eta_p(lpoin_w(ipoin))-eta(lpoin_w(ipoin)))/dt  ! Temporal entropy
                          Reta(lpoin_w(ipoin)) = 0.0d0
                          alpha(lpoin_w(ipoin)) = 1.0d0
                       end do
                       !$acc end parallel loop
                       !
                       ! Entropy residual
                       !
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He, &
                                                  gpvol,f_eta,Reta,alpha)
                       !
                       ! Alter Reta with inv(Mc)
                       !
                       call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
                       !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Reta)
                       call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w,ngaus,connec,gpvol,Ngp,ppow,Ml,Reta)
                       !
                       ! Update Reta
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          Reta(lpoin_w(ipoin)) = Reta(lpoin_w(ipoin))+1.0d0*R1(lpoin_w(ipoin))
                       end do
                       !$acc end parallel loop

                       !
                       ! Temporal mass
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          R2(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)-rhok(lpoin_w(ipoin)))/dt
                          alpha(lpoin_w(ipoin)) = eta(lpoin_w(ipoin))/rhok(lpoin_w(ipoin))
                       end do
                       !$acc end parallel loop
                       !
                       ! Alter R2 with Mcw
                       !
                       call wcmass_times_vector(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,R2,aux1,alpha)
                       !
                       ! Compute weighted mass convec
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          Rrho(lpoin_w(ipoin)) = 0.0d0
                       end do
                       !$acc end parallel loop
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp, &
                                                  dNgp,He,gpvol,f_rho,Rrho,alpha)
                       !
                       ! Update Rrho with both terms
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          Rrho(lpoin_w(ipoin)) = Rrho(lpoin_w(ipoin))+1.0d0*aux1(lpoin_w(ipoin))
                       end do
                       !$acc end parallel loop
                       !
                       ! Apply solver
                       !
                       call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rrho)
                       !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rrho)
                       call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w,ngaus,connec,gpvol,Ngp,ppow,Ml,Rrho)
                       call nvtxEndRange

              end subroutine residuals

              subroutine smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho, &
                                    gamma_gas,rho,u,pr,helem,mu_e)
              
                      ! TODO: Compute element size h
              
                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, ndime, npoin, connec(nelem,nnode)
                      real(8),    intent(in)  :: Reta(npoin), Rrho(npoin), helem(nelem), gamma_gas
                      real(8),    intent(in)  :: rho(npoin), u(npoin,ndime), pr(npoin)
                      real(8),    intent(out) :: mu_e(nelem)
                      integer(4)              :: ielem, ind(nnode), inode
                      real(8)                 :: R1, R2, Ve, ue(nnode,ndime), rhoe(nnode), pre(nnode)
                      real(8)                 :: uabs, c_sound, betae
                      real(8)                 :: L3, aux

                      !$acc parallel loop gang private(ind)
                      do ielem = 1,nelem
                         !
                         ! Initialize arrays
                         !
                         ind = connec(ielem,:)             ! Element indexes
                         R1 = maxval(abs(Reta(ind)))       ! Linf norm of Reta on element
                         R2 = maxval(abs(Rrho(ind)))       ! Linf norm of Rrho on element
                         Ve = max(R1,R2)*(helem(ielem)**2) ! Normalized residual for element
                         !
                         ! Max. Wavespeed at element
                         !
                         aux = 0.0d0
                         !$acc loop vector reduction(max:aux)
                         do inode = 1,nnode
                            uabs = sqrt(dot_product(u(ind(inode),:),u(ind(inode),:))) ! Velocity mag. at element node
                            c_sound = sqrt(gamma_gas*pr(ind(inode))/rho(ind(inode)))     ! Speed of sound at node
                            L3 = abs(uabs+c_sound)                          ! L3 wavespeed
                            aux = max(aux,L3)
                         end do
                         !
                         ! Select against Upwind viscosity
                         !
                         betae = 0.5d0*helem(ielem)*aux
                         mu_e(ielem) = maxval(abs(rho(ind)))*min(Ve,betae) ! Dynamic viscosity
                         !mu_e(ielem) = betae
                      end do
                      !$acc end parallel loop

              end subroutine smart_visc

end module mod_entropy_viscosity
