module time_integ

      use elem_convec

      contains

              subroutine rk_4(nelem,npoin,ndime,ndof,nbnodes,nstep,ngaus,nnode, &
                              ldof,lbnodes,connec,Ngp,gpcar,gpvol,dt,rho,u,q,pr,E,Tem,e_int)

                      implicit none

                      integer(4), intent(in)             :: nelem, npoin, ndime, ngaus, nnode, ndof, nstep, nbnodes
                      integer(4), intent(in)             :: ldof(ndof), lbnodes(nbnodes), connec(nelem,nnode)
                      real(8),    intent(in)             :: Ngp(ngaus,nnode), gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)             :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)             :: dt
                      real(8),    intent(inout)          :: rho(npoin)
                      real(8),    intent(inout)          :: u(npoin,ndime)
                      real(8),    intent(inout)          :: q(npoin,ndime)
                      real(8),    intent(inout)          :: pr(npoin)
                      real(8),    intent(inout)          :: E(npoin)
                      real(8),    intent(inout)          :: Tem(npoin)
                      real(8),    intent(inout)          :: e_int(npoin)
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

                      write(*,*) '--| START OF TIME-INTEGRATION PROCESS...'
                      
                      do istep = 1,nstep

                         write(*,*) '   --| STEP ',istep

                         !
                         ! Sub Step 1
                         !

                         write(*,*) '         SOD2D(1)'

                         rho_1 = rho
                         u_1 = u
                         q_1 = q
                         pr_1 = pr
                         E_1 = E
                         Tem_1 = Tem
                         e_int_1 = 0.0d0

                         write(*,*) '            (MASS)'
                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rmass_1)
                         rho_1(:) = rho(:)-(dt/2.0d0)*Rmass_1(:)

                         write(*,*) '            (MOMENTUM)'
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,q,pr,Rmom_1)
                         q_1(:,:) = q(:,:)-(dt/2.0d0)*Rmom_1(:,:)
                         q_1(lbnodes,2) = 0.0d0
                         do ipoin = 1,npoin
                            u_1(ipoin,:) = q_1(ipoin,:)/rho_1(ipoin)
                         end do

                         write(*,*) '            (ENERGY)'
                         call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,pr,E,Rener_1)
                         E_1(:) = E(:)-(dt/2.0d0)*Rener_1(:)

                         do ipoin = 1,npoin
                            do idime = 1,ndime
                               e_int_1(ipoin) = e_int_1(ipoin)+u_1(ipoin,idime)**2
                            end do
                         end do
                         e_int_1 = (E_1/rho_1)-0.5d0*e_int_1
                         pr_1 = rho_1*(1.400d0-1.0d0)*e_int_1

                         !
                         ! Sub Step 2
                         !

                         write(*,*) '         SOD2D(2)'

                         rho_2 = rho
                         u_2 = u
                         q_2 = q
                         pr_2 = pr
                         E_2 = E
                         Tem_2 = Tem
                         e_int_2 = 0.0d0

                         write(*,*) '            (MASS)'
                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_1,Rmass_2)
                         rho_2(:) = rho(:)-(dt/2.0d0)*Rmass_2(:)

                         write(*,*) '            (MOMENTUM)'
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_1,q_1,pr_1,Rmom_2)
                         q_2(:,:) = q(:,:)-(dt/2.0d0)*Rmom_2(:,:)
                         q_2(lbnodes,2) = 0.0d0
                         do ipoin = 1,npoin
                            u_2(ipoin,:) = q_2(ipoin,:)/rho_2(ipoin)
                         end do

                         write(*,*) '            (ENERGY)'
                         call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_1,pr_1,E_1,Rener_2)
                         E_2(:) = E(:)-(dt/2.0d0)*Rener_2(:)

                         do ipoin = 1,npoin
                            do idime = 1,ndime
                               e_int_2(ipoin) = e_int_2(ipoin)+u_2(ipoin,idime)**2
                            end do
                         end do
                         e_int_2 = (E_2/rho_2)-0.5d0*e_int_2
                         pr_2 = rho_2*(1.400d0-1.0d0)*e_int_2

                         !
                         ! Sub Step 3
                         !

                         write(*,*) '         SOD2D(3)'

                         rho_3 = rho
                         u_3 = u
                         q_3 = q
                         pr_3 = pr
                         E_3 = E
                         Tem_3 = Tem
                         e_int_3 = e_int

                         write(*,*) '            (MASS)'
                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_2,Rmass_3)
                         rho_3(:) = rho(:)-(dt/1.0d0)*Rmass_3(:)

                         write(*,*) '            (MOMENTUM)'
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_2,q_2,pr_2,Rmom_3)
                         q_3(:,:) = q(:,:)-(dt/1.0d0)*Rmom_3(:,:)
                         q_3(lbnodes,2) = 0.0d0
                         do ipoin = 1,npoin
                            u_3(ipoin,:) = q_3(ipoin,:)/rho_3(ipoin)
                         end do

                         write(*,*) '            (ENERGY)'
                         call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_2,pr_2,E_2,Rener_3)
                         E_3(:) = E(:)-(dt/1.0d0)*Rener_3(:)

                         do ipoin = 1,npoin
                            do idime = 1,ndime
                               e_int_3(ipoin) = e_int_3(ipoin)+u_3(ipoin,idime)**2
                            end do
                         end do
                         e_int_3 = (E_3/rho_3)-0.5d0*e_int_3
                         pr_3 = rho_3*(1.400d0-1.0d0)*e_int_3

                         !
                         ! Sub Step 4
                         !

                         write(*,*) '         SOD2D(4)'

                         rho_4 = rho
                         u_4 = u
                         q_4 = q
                         pr_4 = pr
                         E_4 = E
                         Tem_4 = Tem
                         e_int_4 = 0.0d0

                         write(*,*) '            (MASS)'
                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_3,Rmass_4)
                         aux_mass = Rmass_1+2.0d0*Rmass_2+2.0d0*Rmass_3+Rmass_4
                         rho_4(:) = rho(:)-(dt/6.0d0)*aux_mass(:)

                         write(*,*) '            (MOMENTUM)'
                         call mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_3,q_3,pr_3,Rmom_4)
                         aux_mom = Rmom_1+2.0d0*Rmom_2+2.0d0*Rmom_3+Rmom_4
                         q_4(:,:) = q(:,:)-(dt/6.0d0)*aux_mom(:,:)
                         q_4(lbnodes,2) = 0.0d0
                         do ipoin = 1,npoin
                            u_4(ipoin,:) = q_4(ipoin,:)/rho_4(ipoin)
                         end do

                         write(*,*) '            (ENERGY)'
                         call ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u_3,pr_3,E_3,Rener_4)
                         aux_ener = Rener_1+2.0d0*Rener_2+2.0d0*Rener_3+Rener_4
                         E_4(:) = E(:)-(dt/6.0d0)*aux_ener

                         do ipoin = 1,npoin
                            do idime = 1,ndime
                               e_int_4(ipoin) = e_int_4(ipoin)+u_4(ipoin,idime)**2
                            end do
                         end do
                         e_int_4 = (E_4/rho_4)-0.5d0*e_int_4
                         pr_4 = rho_4*(1.400d0-1.0d0)*e_int_4

                         !
                         ! Update
                         !

                         rho = rho_4
                         u = u_4
                         pr = pr_4
                         E = E_4

                      end do

              end subroutine

end module time_integ
