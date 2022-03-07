module time_integ

      use mod_nvtx
      use elem_convec
      use elem_diffu
      use elem_source
      use mod_solver
      use mod_entropy_viscosity

      contains

              subroutine rk_4_main(flag_predic,flag_emac,nelem,nboun,npbou,npoin,npoin_w,ndime,ngaus,nnode, &
                              ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                              rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w, &
                              ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args

                      implicit none

                      integer(4), intent(in)             :: flag_predic, flag_emac
                      integer(4), intent(in)             :: nelem, nboun, npbou, npoin, ndime, ngaus, nnode
                      integer(4), intent(in)             :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w)
                      integer(4), intent(in)             :: ppow
                      real(8),    intent(in)             :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)             :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)             :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)             :: dt, helem(nelem)
                      real(8),    intent(in)             :: Ml(npoin)
                      real(8),    intent(in)             :: Rgas, gamma_gas
                      real(8),    intent(inout)          :: rho(npoin,2)
                      real(8),    intent(inout)          :: u(npoin,ndime,2)
                      real(8),    intent(inout)          :: q(npoin,ndime,2)
                      real(8),    intent(inout)          :: pr(npoin,2)
                      real(8),    intent(inout)          :: E(npoin,2)
                      real(8),    intent(inout)          :: Tem(npoin,2)
                      real(8),    intent(inout)          :: e_int(npoin,2)
                      real(8),    intent(out)            :: mu_e(nelem)
                      integer(4), optional, intent(in)   :: ndof, nbnodes, ldof(ndof), lbnodes(nbnodes)
                      integer(4), optional, intent(in)   :: bound(nboun,npbou), bou_codes(nboun,2)
                      real(8),    optional, intent(in)   :: source_term(ndime)
                      integer(4)                         :: pos, bcode
                      integer(4)                         :: istep, ipoin, idof, idime, iboun, ipbou
                      real(8),    dimension(npoin)       :: rho_1, rho_2, rho_3, rho_4
                      real(8),    dimension(npoin,ndime) :: u_1, u_2, u_3, u_4
                      real(8),    dimension(npoin,ndime) :: q_1, q_2, q_3, q_4
                      real(8),    dimension(npoin)       :: pr_1, pr_2, pr_3, pr_4
                      real(8),    dimension(npoin)       :: E_1, E_2, E_3, E_4
                      real(8),    dimension(npoin)       :: Tem_1, Tem_2, Tem_3, Tem_4
                      real(8),    dimension(npoin)       :: e_int_1, e_int_2, e_int_3, e_int_4
                      real(8),    dimension(npoin)       :: Rmass_1, Rmass_2, Rmass_3, Rmass_4
                      real(8),    dimension(npoin)       :: Rener_1, Rener_2, Rener_3, Rener_4
                      real(8),    dimension(npoin,ndime) :: Rmom_1, Rmom_2, Rmom_3, Rmom_4
                      real(8),    dimension(npoin)       :: aux_mass, aux_ener, Reta, Rrho
                      real(8),    dimension(npoin,ndime) :: aux_mom
                      real(8)                            :: Rdiff_scal(npoin), Rdiff_vect(npoin,ndime)
                      real(8)                            :: Aemac(npoin,ndime), Femac(npoin)

                      !
                      ! Determine wheter to use prediction position or update position
                      !
                      if (flag_predic == 1) then
                         pos = 1 ! Prediction
                      else if (flag_predic == 0) then
                         pos = 2 ! Update
                      end if

                      !
                      ! Sub Step 1
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(1)'

                      !$acc kernels
                      rho_1(:) = 0.0d0
                      u_1(:,:) = 0.0d0
                      q_1(:,:) = 0.0d0
                      pr_1(:) = 0.0d0
                      E_1(:) = 0.0d0
                      Tem_1(:) = 0.0d0
                      e_int_1(:) = 0.0d0
                      !$acc end kernels

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then

                         !
                         ! Compute Reta and Rrho for selector
                         !
                         call nvtxStartRange("ENVIT")
                         call residuals(nelem,ngaus,npoin,nnode,ndime,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rho(:,2), u(:,:,2), pr(:,2), q(:,:,2), &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                         !
                         ! Compute entropy viscosity
                         !
                         call smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho, &
                                         gamma_gas,rho(:,2),u(:,:,2),pr(:,2),helem,mu_e)
                         call nvtxEndRange

                      end if

                      !
                      ! Mass
                      !
                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,q(:,:,pos),Rmass_1)
                      if (flag_predic == 0) then
                         call mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,rho(:,pos),mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rmass_1(lpoin_w(ipoin)) = Rmass_1(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass_1)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmass_1)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmass_1)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         rho_1(lpoin_w(ipoin)) = rho(lpoin_w(ipoin),pos)- &
                                                 (dt/2.0d0)*Rmass_1(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !
                      ! Momentum
                      !
                      if (flag_emac .eq. 0) then
                         !
                         ! Conservation momentum convection div(qi*uj)
                         !
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),q(:,:,pos),pr(:,pos),Rmom_1)
                      else if (flag_emac .eq. 1) then
                         !
                         ! EMAC term
                         !
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Aemac(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,pos)/sqrt(rho(lpoin_w(ipoin),pos))
                            end do
                         end do
                         !$acc end parallel loop
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Femac(lpoin_w(ipoin)) = dot_product(Aemac(lpoin_w(ipoin),:),Aemac(lpoin_w(ipoin),:))
                         end do
                         !$acc end parallel loop
                         call mom_convec_emac(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr(:,pos),Rmom_1)
                      end if
                      if(present(source_term)) then
                        call mom_source_const_vect(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rmom_1)
                      end if
                      if (flag_predic == 0) then
                         call mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),mu_e,Rdiff_vect)
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Rmom_1(lpoin_w(ipoin),idime) = Rmom_1(lpoin_w(ipoin),idime)+Rdiff_vect(lpoin_w(ipoin),idime)
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_vect(npoin,npoin_w,lpoin_w,ndime,Ml,Rmom_1)
                      !call approx_inverse_vect(ndime,npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmom_1)
                      call approx_inverse_vect(ndime,nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmom_1)
                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            q_1(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,pos)- &
                                                        (dt/2.0d0)*Rmom_1(lpoin_w(ipoin),idime)
                         end do
                      end do
                      !$acc end parallel loop

                      !
                      ! Janky boundary conditions. TODO: Fix this shite...
                      !
                      if (nboun .ne. 0) then
                         if (ndime == 2) then
                            !$acc kernels
                            q_1(lbnodes,2) = 0.0d0
                            !$acc end kernels
                         else if (ndime == 3) then
                            !
                            ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                            ! Nodes belonging to both codes will be zeroed on both directions.
                            ! Like this, there's no need to fnd intersections.
                            !
                            !$acc parallel loop gang
                            do iboun = 1,nboun
                               bcode = bou_codes(iboun,2) ! Boundary element code
                               if (bcode == 1) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),2) = 0.0d0
                                  end do
                               else if (bcode == 2) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               else if (bcode == 3) then ! non_slip wall
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),1) = 0.0d0
                                     q_1(bound(iboun,ipbou),2) = 0.0d0
                                     q_1(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               end if
                            end do
                            !$acc end parallel loop
                         end if
                      end if

                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            u_1(lpoin_w(ipoin),idime) = q_1(lpoin_w(ipoin),idime)/rho_1(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end parallel loop

                      !
                      ! Total energy
                      !
                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),pr(:,pos),E(:,pos),Rener_1)
                      if (flag_predic == 0) then
                         call ener_diffusion(nelem,ngaus,npoin,nnode,ndime, &
                                             connec,Ngp,dNgp,He,gpvol,u(:,:,pos),Tem(:,pos),mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rener_1(lpoin_w(ipoin)) = Rener_1(lpoin_w(ipoin))+Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener_1)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rener_1)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w,ngaus, &
                                                 connec,gpvol,Ngp,ppow,Ml,Rener_1)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         E_1(lpoin_w(ipoin)) = E(lpoin_w(ipoin),pos)- &
                                               (dt/2.0d0)*Rener_1(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         e_int_1(lpoin_w(ipoin)) = (E_1(lpoin_w(ipoin))/rho_1(lpoin_w(ipoin)))- &
                            0.5d0*dot_product(u_1(lpoin_w(ipoin),:),u_1(lpoin_w(ipoin),:))
                         pr_1(lpoin_w(ipoin)) = rho_1(lpoin_w(ipoin))*(gamma_gas-1.0d0)*e_int_1(lpoin_w(ipoin))
                         Tem_1(lpoin_w(ipoin)) = pr_1(lpoin_w(ipoin))/(rho_1(lpoin_w(ipoin))*Rgas)
                      end do
                      !$acc end parallel loop

                      !
                      ! Sub Step 2
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(2)'

                      !$acc kernels
                      rho_2(:) = 0.0d0
                      u_2(:,:) = 0.0d0
                      q_2(:,:) = 0.0d0
                      pr_2(:) = 0.0d0
                      E_2(:) = 0.0d0
                      Tem_2(:) = 0.0d0
                      e_int_2(:) = 0.0d0
                      !$acc end kernels

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then

                         !
                         ! Compute Reta and Rrho for selector
                         !
                         call nvtxStartRange("ENVIT")
                         call residuals(nelem,ngaus,npoin,nnode,ndime,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rho_1, u_1, pr_1, q_1, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                         !
                         ! Compute entropy viscosity
                         !
                         call smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho, &
                                         gamma_gas,rho_1,u_1,pr_1,helem,mu_e)
                         call nvtxEndRange

                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,q_1,Rmass_2)
                      if (flag_predic == 0) then
                         call mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,rho_1,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rmass_2(lpoin_w(ipoin)) = Rmass_2(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass_2)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmass_2)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmass_2)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         rho_2(lpoin_w(ipoin)) = rho(lpoin_w(ipoin),pos)- &
                                                 (dt/2.0d0)*Rmass_2(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      if (flag_emac .eq. 0) then
                         !
                         ! Conservation momentum convection div(qi*uj)
                         !
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_1,q_1,pr_1,Rmom_2)
                      else if (flag_emac .eq. 1) then
                         !
                         ! EMAC term
                         !
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Aemac(lpoin_w(ipoin),idime) = q_1(lpoin_w(ipoin),idime)/sqrt(rho_1(lpoin_w(ipoin)))
                            end do
                         end do
                         !$acc end parallel loop
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Femac(lpoin_w(ipoin)) = dot_product(Aemac(lpoin_w(ipoin),:),Aemac(lpoin_w(ipoin),:))
                         end do
                         !$acc end parallel loop
                         call mom_convec_emac(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr_1,Rmom_2)
                      end if
                      if(present(source_term)) then
                        call mom_source_const_vect(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rmom_2)
                      end if
                      if (flag_predic == 0) then
                         call mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_1,mu_e,Rdiff_vect)
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Rmom_2(lpoin_w(ipoin),idime) = Rmom_2(lpoin_w(ipoin),idime) + Rdiff_vect(lpoin_w(ipoin),idime)
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_vect(npoin,npoin_w,lpoin_w,ndime,Ml,Rmom_2)
                      !call approx_inverse_vect(ndime,npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmom_2)
                      call approx_inverse_vect(ndime,nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmom_2)
                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            q_2(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,pos)- &
                                                        (dt/2.0d0)*Rmom_2(lpoin_w(ipoin),idime)
                         end do
                      end do
                      !$acc end parallel loop

                      !
                      ! Janky boundary conditions. TODO: Fix this shite...
                      !
                      if (nboun .ne. 0) then
                         if (ndime == 2) then
                            !$acc kernels
                            q_2(lbnodes,2) = 0.0d0
                            !$acc end kernels
                         else if (ndime == 3) then
                            !
                            ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                            ! Nodes belonging to both codes will be zeroed on both directions.
                            ! Like this, there's no need to fnd intersections.
                            !
                            !$acc parallel loop gang
                            do iboun = 1,nboun
                               bcode = bou_codes(iboun,2) ! Boundary element code
                               if (bcode == 1) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_2(bound(iboun,ipbou),2) = 0.0d0
                                  end do
                               else if (bcode == 2) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_2(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               else if (bcode == 3) then ! non_slip wall
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),1) = 0.0d0
                                     q_1(bound(iboun,ipbou),2) = 0.0d0
                                     q_1(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               end if
                            end do
                            !$acc end parallel loop
                         end if
                      end if

                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            u_2(lpoin_w(ipoin),idime) = q_2(lpoin_w(ipoin),idime)/rho_2(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end parallel loop

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_1,pr_1,E_1,Rener_2)
                      if (flag_predic == 0) then
                         call ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_1,Tem_1,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rener_2(lpoin_w(ipoin)) = Rener_2(lpoin_w(ipoin))+Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener_2)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rener_2)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rener_2)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         E_2(lpoin_w(ipoin)) = E(lpoin_w(ipoin),pos)- &
                                               (dt/2.0d0)*Rener_2(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         e_int_2(lpoin_w(ipoin)) = (E_2(lpoin_w(ipoin))/rho_2(lpoin_w(ipoin)))- &
                            0.5d0*dot_product(u_2(lpoin_w(ipoin),:),u_2(lpoin_w(ipoin),:))
                         pr_2(lpoin_w(ipoin)) = rho_2(lpoin_w(ipoin))*(gamma_gas-1.0d0)*e_int_2(lpoin_w(ipoin))
                         Tem_2(lpoin_w(ipoin)) = pr_2(lpoin_w(ipoin))/(rho_2(lpoin_w(ipoin))*Rgas)
                      end do
                      !$acc end parallel loop

                      !
                      ! Sub Step 3
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(3)'

                      !$acc kernels
                      rho_3(:) = 0.0d0
                      u_3(:,:) = 0.0d0
                      q_3(:,:) = 0.0d0
                      pr_3(:) = 0.0d0
                      E_3(:) = 0.0d0
                      Tem_3(:) = 0.0d0
                      e_int_3(:) = 0.0d0
                      !$acc end kernels

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then

                         !
                         ! Compute Reta and Rrho for selector
                         !
                         call nvtxStartRange("ENVIT")
                         call residuals(nelem,ngaus,npoin,nnode,ndime,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rho_2, u_2, pr_2, q_2, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                         !
                         ! Compute entropy viscosity
                         !
                         call smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho, &
                                         gamma_gas,rho_2,u_2,pr_2,helem,mu_e)
                         call nvtxEndRange

                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,q_2,Rmass_3)
                      if (flag_predic == 0) then
                         call mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,rho_2,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rmass_3(lpoin_w(ipoin)) = Rmass_3(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass_3)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmass_3)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmass_3)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         rho_3(lpoin_w(ipoin)) = rho(lpoin_w(ipoin),pos)-(dt/1.0d0)*Rmass_3(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      if (flag_emac .eq. 0) then
                         !
                         ! Conservation momentum convection div(qi*uj)
                         !
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_2,q_2,pr_2,Rmom_3)
                      else if (flag_emac .eq. 1) then
                         !
                         ! EMAC term
                         !
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Aemac(lpoin_w(ipoin),idime) = q_2(lpoin_w(ipoin),idime)/sqrt(rho_2(lpoin_w(ipoin)))
                            end do
                         end do
                         !$acc end parallel loop
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Femac(lpoin_w(ipoin)) = dot_product(Aemac(lpoin_w(ipoin),:),Aemac(lpoin_w(ipoin),:))
                         end do
                         !$acc end parallel loop
                         call mom_convec_emac(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr_2,Rmom_3)
                      end if
                      if(present(source_term)) then
                        call mom_source_const_vect(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rmom_3)
                      end if
                      if (flag_predic == 0) then
                         call mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_2,mu_e,Rdiff_vect)
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Rmom_3(lpoin_w(ipoin),idime) = Rmom_3(lpoin_w(ipoin),idime)+Rdiff_vect(lpoin_w(ipoin),idime)
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_vect(npoin,npoin_w,lpoin_w,ndime,Ml,Rmom_3)
                      !call approx_inverse_vect(ndime,npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmom_3)
                      call approx_inverse_vect(ndime,nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmom_3)
                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            q_3(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,pos)- &
                                                        (dt/1.0d0)*Rmom_3(lpoin_w(ipoin),idime)
                         end do
                      end do
                      !$acc end parallel loop

                      !
                      ! Janky boundary conditions. TODO: Fix this shite...
                      !
                      if (nboun .ne. 0) then
                         if (ndime == 2) then
                            !$acc kernels
                            q_3(lbnodes,2) = 0.0d0
                            !$acc end kernels
                         else if (ndime == 3) then
                            !
                            ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                            ! Nodes belonging to both codes will be zeroed on both directions.
                            ! Like this, there's no need to fnd intersections.
                            !
                            !$acc parallel loop gang
                            do iboun = 1,nboun
                               bcode = bou_codes(iboun,2) ! Boundary element code
                               if (bcode == 1) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_3(bound(iboun,ipbou),2) = 0.0d0
                                  end do
                               else if (bcode == 2) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_3(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               else if (bcode == 3) then ! non_slip wall
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),1) = 0.0d0
                                     q_1(bound(iboun,ipbou),2) = 0.0d0
                                     q_1(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               end if
                            end do
                            !$acc end parallel loop
                         end if
                      end if

                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            u_3(lpoin_w(ipoin),idime) = q_3(lpoin_w(ipoin),idime)/rho_3(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end parallel loop

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_2,pr_2,E_2,Rener_3)
                      if (flag_predic == 0) then
                         call ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_2,Tem_2,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rener_3(lpoin_w(ipoin)) = Rener_3(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener_3)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rener_3)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rener_3)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         E_3(lpoin_w(ipoin)) = E(lpoin_w(ipoin),pos)- &
                                               (dt/1.0d0)*Rener_3(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         e_int_3(lpoin_w(ipoin)) = (E_3(lpoin_w(ipoin))/rho_3(lpoin_w(ipoin)))- &
                            0.5d0*dot_product(u_3(lpoin_w(ipoin),:),u_3(lpoin_w(ipoin),:))
                         pr_3(lpoin_w(ipoin)) = rho_3(lpoin_w(ipoin))*(gamma_gas-1.0d0)*e_int_3(lpoin_w(ipoin))
                         Tem_3(lpoin_w(ipoin)) = pr_3(lpoin_w(ipoin))/(rho_3(lpoin_w(ipoin))*Rgas)
                      end do
                      !$acc end parallel loop

                      !
                      ! Sub Step 4
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(4)'

                      !$acc kernels
                      rho_4(:) = 0.0d0
                      u_4(:,:) = 0.0d0
                      q_4(:,:) = 0.0d0
                      pr_4(:) = 0.0d0
                      E_4(:) = 0.0d0
                      Tem_4(:) = 0.0d0
                      e_int_4(:) = 0.0d0
                      !$acc end kernels

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then

                         !
                         ! Compute Reta and Rrho for selector
                         !
                         call nvtxStartRange("ENVIT")
                         call residuals(nelem,ngaus,npoin,nnode,ndime,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rho_3, u_3, pr_3, q_3, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                         !
                         ! Compute entropy viscosity
                         !
                         call smart_visc(nelem,nnode,ndime,npoin,connec,Reta,Rrho, &
                                         gamma_gas,rho_3,u_3,pr_3,helem,mu_e)
                         call nvtxEndRange

                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,q_3,Rmass_4)
                      if (flag_predic == 0) then
                         call mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,rho_3,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rmass_4(lpoin_w(ipoin)) = Rmass_4(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass_4)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmass_4)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmass_4)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         aux_mass(lpoin_w(ipoin)) = Rmass_1(lpoin_w(ipoin))+2.0d0*Rmass_2(lpoin_w(ipoin))+ &
                            2.0d0*Rmass_3(lpoin_w(ipoin))+Rmass_4(lpoin_w(ipoin))
                         rho_4(lpoin_w(ipoin)) = rho(lpoin_w(ipoin),pos)- &
                            (dt/6.0d0)*aux_mass(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      if (flag_emac .eq. 0) then
                         !
                         ! Conservation momentum convection div(qi*uj)
                         !
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_3,q_3,pr_3,Rmom_4)
                      else if (flag_emac .eq. 1) then
                         !
                         ! EMAC term
                         !
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Aemac(lpoin_w(ipoin),idime) = q_3(lpoin_w(ipoin),idime)/sqrt(rho_3(lpoin_w(ipoin)))
                            end do
                         end do
                         !$acc end parallel loop
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Femac(lpoin_w(ipoin)) = dot_product(Aemac(lpoin_w(ipoin),:),Aemac(lpoin_w(ipoin),:))
                         end do
                         !$acc end parallel loop
                         call mom_convec_emac(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr_3,Rmom_4)
                      end if
                      if(present(source_term)) then
                        call mom_source_const_vect(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rmom_4)
                      end if
                      if (flag_predic == 0) then
                         call mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_3,mu_e,Rdiff_vect)
                         !$acc parallel loop collapse(2)
                         do ipoin = 1,npoin_w
                            do idime = 1,ndime
                               Rmom_4(lpoin_w(ipoin),idime) = Rmom_4(lpoin_w(ipoin),idime)+Rdiff_vect(lpoin_w(ipoin),idime)
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_vect(npoin,npoin_w,lpoin_w,ndime,Ml,Rmom_4)
                      !call approx_inverse_vect(ndime,npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rmom_4)
                      call approx_inverse_vect(ndime,nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rmom_4)
                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            aux_mom(lpoin_w(ipoin),idime) = Rmom_1(lpoin_w(ipoin),idime)+2.0d0*Rmom_2(lpoin_w(ipoin),idime)+ &
                               2.0d0*Rmom_3(lpoin_w(ipoin),idime)+Rmom_4(lpoin_w(ipoin),idime)
                            q_4(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,pos)- &
                               (dt/6.0d0)*aux_mom(lpoin_w(ipoin),idime)
                         end do
                      end do
                      !$acc end parallel loop

                      !
                      ! Janky boundary conditions. TODO: Fix this shite...
                      !
                      if (nboun .ne. 0) then
                         if (ndime == 2) then
                            !$acc kernels
                            q_4(lbnodes,2) = 0.0d0
                            !$acc end kernels
                         else if (ndime == 3) then
                            !
                            ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                            ! Nodes belonging to both codes will be zeroed on both directions.
                            ! Like this, there's no need to fnd intersections.
                            !
                            !$acc parallel loop gang
                            do iboun = 1,nboun
                               bcode = bou_codes(iboun,2) ! Boundary element code
                               if (bcode == 1) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_4(bound(iboun,ipbou),2) = 0.0d0
                                  end do
                               else if (bcode == 2) then
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_4(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               else if (bcode == 3) then ! non_slip wall
                                  !$acc loop vector
                                  do ipbou = 1,npbou
                                     q_1(bound(iboun,ipbou),1) = 0.0d0
                                     q_1(bound(iboun,ipbou),2) = 0.0d0
                                     q_1(bound(iboun,ipbou),3) = 0.0d0
                                  end do
                               end if
                            end do
                            !$acc end parallel loop
                         end if
                      end if
                      
                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            u_4(lpoin_w(ipoin),idime) = q_4(lpoin_w(ipoin),idime)/rho_4(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end parallel loop

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_3,pr_3,E_3,Rener_4)
                      if (flag_predic == 0) then
                         call ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u_3,Tem_3,mu_e,Rdiff_scal)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            Rener_4(lpoin_w(ipoin)) = Rener_4(lpoin_w(ipoin)) + Rdiff_scal(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end if
                      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener_4)
                      !call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,Rener_4)
                      call approx_inverse_scalar(nelem,nnode,npoin,npoin_w,lpoin_w, &
                         ngaus,connec,gpvol,Ngp,ppow,Ml,Rener_4)
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         aux_ener(lpoin_w(ipoin)) = Rener_1(lpoin_w(ipoin))+2.0d0*Rener_2(lpoin_w(ipoin))+ &
                            2.0d0*Rener_3(lpoin_w(ipoin))+Rener_4(lpoin_w(ipoin))
                         E_4(lpoin_w(ipoin)) = E(lpoin_w(ipoin),pos)- &
                            (dt/6.0d0)*aux_ener(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         e_int_4(lpoin_w(ipoin)) = (E_4(lpoin_w(ipoin))/rho_4(lpoin_w(ipoin)))- &
                            0.5d0*dot_product(u_4(lpoin_w(ipoin),:),u_4(lpoin_w(ipoin),:))
                         pr_4(lpoin_w(ipoin)) = rho_4(lpoin_w(ipoin))*(gamma_gas-1.0d0)*e_int_4(lpoin_w(ipoin))
                         Tem_4(lpoin_w(ipoin)) = pr_4(lpoin_w(ipoin))/(rho_4(lpoin_w(ipoin))*Rgas)
                      end do
                      !$acc end parallel loop

                      !
                      ! Update
                      !

                      call nvtxStartRange("Update")
                      !$acc kernels
                      rho(:,pos) = rho_4(:)
                      u(:,:,pos) = u_4(:,:)
                      pr(:,pos) = pr_4(:)
                      E(:,pos) = E_4(:)
                      q(:,:,pos) = q_4(:,:)
                      e_int(:,pos) = e_int_4(:)
                      Tem(:,pos) = Tem_4(:)
                      !$acc end kernels
                      call nvtxEndRange

              end subroutine rk_4_main
                      
end module time_integ
