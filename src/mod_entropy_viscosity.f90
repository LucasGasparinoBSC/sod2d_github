module mod_entropy_viscosity

      ! TODO: Finish module and create unit tests

      contains

              subroutine residuals(nelem,ngaus,npoin,nnode,ndime, nzdom, &
                                   rdom, cdom, ppow, connec, Ngp, gpcar, gpvol, Ml, Mc, &
                                   dt, rhok, uk, prk, qk, &
                                   rho, u, pr, q, &
                                   Reta, Rrho)

                      use mod_solver
                      use elem_convec

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime, nzdom, ppow
                      integer(4), intent(in)  :: connec(nelem,nnode), rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)  :: Ngp(nnode,ngaus), gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ml(npoin), Mc(nzdom)
                      real(8),    intent(in)  :: dt
                      real(8),    intent(in)  :: rhok(npoin), uk(npoin,ndime), prk(npoin), qk(npoin,ndime)     ! From substep
                      real(8),    intent(in)  :: rho(npoin,2), u(npoin,ndime,2), pr(npoin,2), q(npoin,ndime,2) ! From prediction
                      real(8),    intent(out) :: Reta(npoin), Rrho(npoin)
                      integer(4)              :: ipoin
                      real(8)                 :: eta(npoin), eta_p(npoin), alpha(npoin), alpha_p(npoin)
                      real(8)                 :: f_eta(npoin,ndime), f_rho(npoin,ndime), R1(npoin), R2(npoin)
                      real(8)                 :: aux1(npoin)

                       !
                       ! Entropy function and temporal terms
                       !
                       do ipoin = 1,npoin
                          !
                          ! Current (substesp values)
                          !
                          eta(ipoin) = (rhok(ipoin)/0.400d0)*log(prk(ipoin)/(rhok(ipoin)**1.40d0))
                          f_eta(ipoin,1:ndime) = uk(ipoin,1:ndime)*eta(ipoin)
                          alpha(ipoin) = eta(ipoin)/rhok(ipoin)
                          f_rho(ipoin,1:ndime) = alpha(ipoin)*qk(ipoin,1:ndime)
                          !
                          ! Prediction
                          !
                          eta_p(ipoin) = (rho(ipoin,1)/0.400d0)*log(pr(ipoin,1)/(rho(ipoin,1)**1.40d0))
                          alpha_p(ipoin) = eta_p(ipoin)/rho(ipoin,1)
                          !
                          ! Temporal term
                          !
                          R1(ipoin) = (eta_p(ipoin)-eta(ipoin))/dt                              ! Temporal entropy
                          R2(ipoin) = (alpha_p(ipoin)*rho(ipoin,1)-alpha(ipoin)*rhok(ipoin))/dt ! Temporal mass
                       end do

                       !
                       ! Alter R1 and R2 with Mc
                       !
                       aux1 = R1
                       call CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,aux1,R1)
                       aux1 = R2
                       call CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,aux1,R2)

                       !
                       ! Compute both residuals
                       !
                       Reta = 0.0d0 ! Entropy residual
                       Rrho = 0.0d0 ! Mass residual
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,f_eta,Reta) ! Entropy convec
                       call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Reta)
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,f_rho,Rrho) ! Mass convec
                       call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rrho)

                       Reta = Reta+R1
                       Rrho = Rrho+R2

              end subroutine residuals

              subroutine smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho,rho,u,pr,helem,mu_e)
              
                      ! TODO: Compute element size h
              
                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, ndime, npoin, connec(nelem,nnode)
                      real(8),    intent(in)  :: Reta(npoin), Rrho(npoin), helem(nelem)
                      real(8),    intent(in)  :: rho(npoin), u(npoin,ndime), pr(npoin)
                      real(8),    intent(out) :: mu_e(nelem)
                      integer(4)              :: ielem, ind(nnode), inode
                      real(8)                 :: R1, R2, Ve, ue(nnode,ndime), rhoe(nnode), pre(nnode)
                      real(8)                 :: uabs, c_sound, L3(nnode), betae

                      do ielem = 1,nelem
                         !
                         ! Initialize arrays
                         !
                         ind = connec(ielem,:)             ! Element indexes
                         R1 = maxval(abs(Reta(ind)))       ! Linf norm of Reta on element
                         R2 = maxval(abs(Rrho(ind)))       ! Linf norm of Rrho on element
                         Ve = max(R1,R2)*(helem(ielem)**2) ! Normalized residual for element
                         ue = u(ind,:)                     ! Element velocities
                         rhoe = rho(ind)                   ! Element density
                         pre = pr(ind)                     ! Element pressure
                         !
                         ! Max. Wavespeed at element
                         !
                         do inode = 1,nnode
                            uabs = sqrt(dot_product(ue(inode,:),ue(inode,:))) ! Velocity mag. at element node
                            c_sound = sqrt(1.40d0*pre(inode)/rhoe(inode))     ! Speed of sound at node
                            L3(inode) = uabs+c_sound                          ! L3 wavespeed
                         end do
                         !
                         ! Select against Upwind viscosity
                         !
                         betae = 0.5d0*helem(ielem)*maxval(abs(L3))
                         mu_e(ielem) = maxval(abs(rhoe))*min(Ve,betae) ! Dynamic viscosity
                      end do

              end subroutine smart_visc

end module mod_entropy_viscosity
