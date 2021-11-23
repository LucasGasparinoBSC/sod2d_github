module mod_entropy_viscosity

   use mod_nvtx

      ! TODO: Finish module and create unit tests

      contains

              subroutine residuals(nelem,ngaus,npoin,nnode,ndime, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rhok, uk, prk, qk, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                      use mass_matrix
                      use mod_solver
                      use elem_convec

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime, ppow
                      integer(4), intent(in)  :: connec(nelem,nnode)
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
                       !$acc kernels
                       eta(:) = (rhok(:)/(gamma_gas-1.0d0))*log(prk(:)/(rhok(:)**gamma_gas))
                       eta_p(:) = (rho(:,1)/(gamma_gas-1.0d0))*log(pr(:,1)/(rho(:,1)**gamma_gas))
                       !$acc end kernels
                       do idime = 1,ndime
                          !$acc kernels
                          f_eta(:,idime) = uk(:,idime)*eta(:)
                          f_rho(:,idime) = qk(:,idime)
                          !$acc end kernels
                       end do

                       !
                       ! Temporal eta
                       !
                       !$acc kernels
                       R1(:) = (eta_p(:)-eta(:))/dt  ! Temporal entropy
                       Reta(:) = 0.0d0
                       alpha(:) = 1.0d0
                       !$acc end kernels
                       !
                       ! Entropy residual
                       !
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He, &
                                                  gpvol,f_eta,Reta,alpha)
                       !
                       ! Alter Reta with inv(Mc)
                       !
                       call lumped_solver_scal(npoin,Ml,Reta)
                       !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Reta)
                       call approx_inverse_scalar(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,ppow,Ml,Reta)
                       !
                       ! Update Reta
                       !
                       !$acc kernels
                       Reta(:) = Reta(:)+1.0d0*R1(:)
                       !$acc end kernels

                       !
                       ! Temporal mass
                       !
                       !$acc kernels
                       R2(:) = (rho(:,1)-rhok(:))/dt
                       alpha(:) = eta(:)/rhok(:)
                       !$acc end kernels
                       !
                       ! Alter R2 with Mcw
                       !
                       call wcmass_times_vector(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,R2,aux1,alpha)
                       !
                       ! Compute weighted mass convec
                       !
                       !$acc kernels
                       Rrho(:) = 0.0d0
                       !$acc end kernels
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp, &
                                                  dNgp,He,gpvol,f_rho,Rrho,alpha)
                       !
                       ! Update Rrho with both terms
                       !
                       !$acc kernels
                       Rrho(:) = Rrho(:)+1.0d0*aux1(:)
                       !$acc end kernels
                       !
                       ! Apply solver
                       !
                       call lumped_solver_scal(npoin,Ml,Rrho)
                       !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rrho)
                       call approx_inverse_scalar(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,ppow,Ml,Rrho)
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
                      real(8)                 :: L3

                      !$acc parallel loop gang private(ind)
                      do ielem = 1,nelem
                         !
                         ! Initialize arrays
                         !
                         ind = connec(ielem,:)             ! Element indexes
                         R1 = maxval(abs(Reta(ind)))       ! Linf norm of Reta on element
                         R2 = maxval(abs(Rrho(ind)))       ! Linf norm of Rrho on element
                         !Ve = max(R1,R2)*(helem(ielem)**2) ! Normalized residual for element
                         Ve = R1*(helem(ielem)**2) ! Normalized residual for element
                         !ue = u(ind,:)                     ! Element velocities
                         !rhoe = rho(ind)                   ! Element density
                         !pre = pr(ind)                     ! Element pressure
                         !
                         ! Max. Wavespeed at element
                         !
                         !$acc loop vector reduction(max:L3)
                         do inode = 1,nnode
                            uabs = sqrt(dot_product(u(ind(inode),:),u(ind(inode),:))) ! Velocity mag. at element node
                            c_sound = sqrt(gamma_gas*pr(ind(inode))/rho(ind(inode)))     ! Speed of sound at node
                            L3 = abs(uabs+c_sound)                          ! L3 wavespeed
                         end do
                         !
                         ! Select against Upwind viscosity
                         !
                         betae = 0.5d0*helem(ielem)*L3
                         !mu_e(ielem) = maxval(abs(rho(ind)))*min(Ve,betae) ! Dynamic viscosity
                         mu_e(ielem) = betae
                         !mu_e(ielem) = 0.0d0
                      end do
                      !$acc end parallel loop

              end subroutine smart_visc

end module mod_entropy_viscosity
