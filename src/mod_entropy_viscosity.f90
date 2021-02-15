module mod_entropy_viscosity

      contains

              subroutine residuals(nelem,ngaus,npoin,nnode,ndime, &
                                   connec, Ngp, gpcar, gpvol, Ml, Mc, &
                                   dt, rhok, uk, prk, qk, &
                                   rho, u, pr, q, &
                                   Reta, Rrho)

                      use mod_solver
                      use elem_convec

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(nnode,ngaus), gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ml(npoin), Mc(npoin,npoin)
                      real(8),    intent(in)  :: dt
                      real(8),    intent(in)  :: rhok(npoin), uk(npoin,ndime), prk(npoin), qk(npoin,ndime)     ! From substep
                      real(8),    intent(in)  :: rho(npoin,2), u(npoin,ndime,2), pr(npoin,2), q(npoin,ndime,2) ! From prediction
                      real(8),    intent(out) :: Reta(npoin), Rrho(npoin)
                      integer(4)              :: ipoin
                      real(8)                 :: eta(npoin), eta_p(npoin), alpha(npoin), alpha_p(npoin)
                      real(8)                 :: f_eta(npoin,ndime), f_rho(npoin,ndime), R1(npoin), R2(npoin)

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
                       ! Compute both residuals
                       !
                       Reta = 0.0d0 ! Entropy residual
                       Rrho = 0.0d0 ! Mass residual
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,f_eta,Reta)    ! Entropy convec
                       call approx_inverse_scalar(npoin,Ml,Mc,Reta)
                       call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,f_rho,Rrho) ! Mass convec
                       call approx_inverse_scalar(npoin,Ml,Mc,Rrho)

                       Reta = Reta+R1
                       Rrho = Rrho+R2

              end subroutine residuals

              !!subroutine smart_visc()
              !!end subroutine smart_visc

end module mod_entropy_viscosity
