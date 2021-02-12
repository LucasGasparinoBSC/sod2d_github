module time_integ

      use elem_convec
      use elem_diffu
      use mod_solver

      contains

              subroutine rk_4_main(flag_predic,nelem,npoin,ndime,ndof,nbnodes,ngaus,nnode, &
                              ldof,lbnodes,connec,Ngp,gpcar,Ml,Mc,gpvol,dt, &
                              rho,u,q,pr,E,Tem,e_int)

                      implicit none

                      integer(4), intent(in)             :: flag_predic
                      integer(4), intent(in)             :: nelem, npoin, ndime, ngaus, nnode, ndof, nbnodes
                      integer(4), intent(in)             :: ldof(ndof), lbnodes(nbnodes), connec(nelem,nnode)
                      real(8),    intent(in)             :: Ngp(ngaus,nnode), gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)             :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)             :: dt
                      real(8),    intent(in)             :: Ml(npoin)
                      real(8),    intent(in)             :: Mc(npoin,npoin)
                      real(8),    intent(inout)          :: rho(npoin,2)
                      real(8),    intent(inout)          :: u(npoin,ndime,2)
                      real(8),    intent(inout)          :: q(npoin,ndime,2)
                      real(8),    intent(inout)          :: pr(npoin,2)
                      real(8),    intent(inout)          :: E(npoin,2)
                      real(8),    intent(inout)          :: Tem(npoin,2)
                      real(8),    intent(inout)          :: e_int(npoin,2)
                      integer(4)                         :: pos
                      integer(4)                         :: istep, ipoin, idof, idime
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
                      real(8),    dimension(npoin)       :: aux_mass, aux_ener
                      real(8),    dimension(npoin,ndime) :: aux_mom
                      real(8)                            :: Rdiff_scal(npoin), Rdiff_vect(npoin,ndime)

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

                      rho_1 = 0.0d0
                      u_1 = 0.0d0
                      q_1 = 0.0d0
                      pr_1 = 0.0d0
                      E_1 = 0.0d0
                      Tem_1 = 0.0d0
                      e_int_1 = 0.0d0

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then
                         ! call smartvisc
                      end if

                      !
                      ! Mass
                      !
                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q(:,:,pos),Rmass_1)
                      if (flag_predic == 0) then
                         call mass_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,rho(:,pos),Rdiff_scal)
                         Rmass_1 = Rmass_1 + Rdiff_scal
                      end if
                      call lumped_solver_scal(npoin,Ml,Rmass_1)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rmass_1)
                      rho_1(:) = rho(:,pos)-(dt/2.0d0)*Rmass_1(:)

                      !
                      ! Momentum
                      !
                      call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u(:,:,pos),q(:,:,pos),pr(:,pos),Rmom_1)
                      if (flag_predic == 0) then
                         call mom_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u(:,:,pos),Rdiff_vect)
                         Rmom_1 = Rmom_1 + Rdiff_vect
                      end if
                      call lumped_solver_vect(npoin,ndime,Ml,Rmom_1)
                      call approx_inverse_vect(ndime,npoin,Ml,Mc,Rmom_1)
                      q_1(:,:) = q(:,:,pos)-(dt/2.0d0)*Rmom_1(:,:)
                      q_1(lbnodes,2) = 0.0d0
                      do ipoin = 1,npoin
                         u_1(ipoin,:) = q_1(ipoin,:)/rho_1(ipoin)
                      end do

                      !
                      ! Total energy
                      !
                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u(:,:,pos),pr(:,pos),E(:,pos),Rener_1)
                      if (flag_predic == 0) then
                         call ener_diffusion(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u(:,:,pos),Tem(:,pos),Rdiff_scal)
                         Rener_1 = Rener_1 + Rdiff_scal
                      end if
                      call lumped_solver_scal(npoin,Ml,Rener_1)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rener_1)
                      E_1(:) = E(:,pos)-(dt/2.0d0)*Rener_1(:)

                      do ipoin = 1,npoin
                         e_int_1(ipoin) = (E_1(ipoin)/rho_1(ipoin))-0.5d0*dot_product(u_1(ipoin,:),u_1(ipoin,:))
                      end do
                      pr_1 = rho_1*(1.400d0-1.0d0)*e_int_1
                      Tem_1 = pr_1/(rho_1*287.0d0)

                      !
                      ! Sub Step 2
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(2)'

                      rho_2 = 0.0d0
                      u_2 = 0.0d0
                      q_2 = 0.0d0
                      pr_2 = 0.0d0
                      E_2 = 0.0d0
                      Tem_2 = 0.0d0
                      e_int_2 = 0.0d0

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then
                         ! call smartvisc
                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_1,Rmass_2)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rmass_2)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rmass_2)
                      rho_2(:) = rho(:,pos)-(dt/2.0d0)*Rmass_2(:)

                      call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_1,q_1,pr_1,Rmom_2)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_vect(npoin,ndime,Ml,Rmom_2)
                      call approx_inverse_vect(ndime,npoin,Ml,Mc,Rmom_2)
                      q_2(:,:) = q(:,:,pos)-(dt/2.0d0)*Rmom_2(:,:)
                      q_2(lbnodes,2) = 0.0d0
                      do ipoin = 1,npoin
                         u_2(ipoin,:) = q_2(ipoin,:)/rho_2(ipoin)
                      end do

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_1,pr_1,E_1,Rener_2)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rener_2)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rener_2)
                      E_2(:) = E(:,pos)-(dt/2.0d0)*Rener_2(:)

                      do ipoin = 1,npoin
                         e_int_2(ipoin) = (E_2(ipoin)/rho_2(ipoin))-0.5d0*dot_product(u_2(ipoin,:),u_2(ipoin,:))
                      end do
                      pr_2 = rho_2*(1.400d0-1.0d0)*e_int_2

                      !
                      ! Sub Step 3
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(3)'

                      rho_3 = 0.0d0
                      u_3 = 0.0d0
                      q_3 = 0.0d0
                      pr_3 = 0.0d0
                      E_3 = 0.0d0
                      Tem_3 = 0.0d0
                      e_int_3 = 0.0d0

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then
                         ! call smartvisc
                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_2,Rmass_3)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rmass_3)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rmass_3)
                      rho_3(:) = rho(:,pos)-(dt/1.0d0)*Rmass_3(:)

                      call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_2,q_2,pr_2,Rmom_3)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_vect(npoin,ndime,Ml,Rmom_3)
                      call approx_inverse_vect(ndime,npoin,Ml,Mc,Rmom_3)
                      q_3(:,:) = q(:,:,pos)-(dt/1.0d0)*Rmom_3(:,:)
                      q_3(lbnodes,2) = 0.0d0
                      do ipoin = 1,npoin
                         u_3(ipoin,:) = q_3(ipoin,:)/rho_3(ipoin)
                      end do

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_2,pr_2,E_2,Rener_3)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rener_3)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rener_3)
                      E_3(:) = E(:,pos)-(dt/1.0d0)*Rener_3(:)

                      do ipoin = 1,npoin
                         e_int_3(ipoin) = (E_3(ipoin)/rho_3(ipoin))-0.5d0*dot_product(u_3(ipoin,:),u_3(ipoin,:))
                      end do
                      pr_3 = rho_3*(1.400d0-1.0d0)*e_int_3

                      !
                      ! Sub Step 4
                      !

                      if (flag_predic == 0) write(*,*) '         SOD2D(4)'

                      rho_4 = 0.0d0
                      u_4 = 0.0d0
                      q_4 = 0.0d0
                      pr_4 = 0.0d0
                      E_4 = 0.0d0
                      Tem_4 = 0.0d0
                      e_int_4 = 0.0d0

                      !
                      ! Entropy viscosity update
                      !
                      if (flag_predic == 0) then
                         ! call smartvisc
                      end if

                      call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_3,Rmass_4)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rmass_4)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rmass_4)
                      aux_mass = Rmass_1+2.0d0*Rmass_2+2.0d0*Rmass_3+Rmass_4
                      rho_4(:) = rho(:,pos)-(dt/6.0d0)*aux_mass(:)

                      call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_3,q_3,pr_3,Rmom_4)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_vect(npoin,ndime,Ml,Rmom_4)
                      call approx_inverse_vect(ndime,npoin,Ml,Mc,Rmom_4)
                      aux_mom = Rmom_1+2.0d0*Rmom_2+2.0d0*Rmom_3+Rmom_4
                      q_4(:,:) = q(:,:,pos)-(dt/6.0d0)*aux_mom(:,:)
                      q_4(lbnodes,2) = 0.0d0
                      do ipoin = 1,npoin
                         u_4(ipoin,:) = q_4(ipoin,:)/rho_4(ipoin)
                      end do

                      call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_3,pr_3,E_3,Rener_4)
                      if (flag_predic == 0) then
                         ! call diff
                         ! Rener_1 = Rener_1 + Rdiff_sca
                      end if
                      call lumped_solver_scal(npoin,Ml,Rener_4)
                      call approx_inverse_scalar(npoin,Ml,Mc,Rmom_4)
                      aux_ener = Rener_1+2.0d0*Rener_2+2.0d0*Rener_3+Rener_4
                      E_4(:) = E(:,pos)-(dt/6.0d0)*aux_ener

                      do ipoin = 1,npoin
                         e_int_4(ipoin) = (E_4(ipoin)/rho_4(ipoin))-0.5d0*dot_product(u_4(ipoin,:),u_4(ipoin,:))
                      end do
                      pr_4 = rho_4*(1.400d0-1.0d0)*e_int_4

                      !
                      ! Update
                      !

                      rho(:,pos) = rho_4
                      u(:,:,pos) = u_4
                      pr(:,pos) = pr_4
                      E(:,pos) = E_4
                      q(:,:,pos) = q_4
                      e_int(:,pos) = e_int_4
                      Tem(:,pos) = Tem_4

              end subroutine rk_4_main
                      
end module time_integ
