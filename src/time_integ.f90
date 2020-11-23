module time_integ

      use elem_convec

      contains

              subroutine rk_4(nelem,npoin,ndime,ndof,nstep,ngaus,nnode, &
                              ldof,connec,Ngp,gpcar,gpvol,dt,rho,u,q,pr,E,Tem,e_int)

                      implicit none

                      integer(4), intent(in)             :: nelem, npoin, ndime, ngaus, nnode, ndof, nstep
                      integer(4), intent(in)             :: ldof(ndof), connec(nelem,nnode)
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
                      integer(4)                         :: istep, ipoin
                      real(8),    dimension(npoin)       :: rho_1, rho_2, rho_3, rho_4
                      real(8),    dimension(npoin,ndime) :: u_1, u_2, u_3, u_4
                      real(8),    dimension(npoin,ndime) :: q_1, q_2, q_3, q_4
                      real(8),    dimension(npoin)       :: pr_1, pr_2, pr_3, pr_4
                      real(8),    dimension(npoin)       :: E_1, E_2, E_3, E_4
                      real(8),    dimension(npoin)       :: Tem_1, Tem_2, Tem_3, Tem_4
                      real(8),    dimension(npoin)       :: e_int_1, e_int_2, e_int_3, e_int_4
                      real(8),    dimension(npoin)       :: Rmass_1, Rmass_2, Rmass_3, Rmass_4
                      real(8),    dimension(npoin)       :: aux

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
                         e_int_1 = e_int

                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rmass_1)
                         rho_1(ldof) = rho(ldof)-(dt/2.0d0)*Rmass_1(ldof)
                         do ipoin = 1,npoin
                            q_1(ipoin,:) = rho_1(ipoin)*u_1(ipoin,:)
                         end do

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
                         e_int_2 = e_int

                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_1,Rmass_2)
                         rho_2(ldof) = rho(ldof)-(dt/2.0d0)*Rmass_2(ldof)
                         do ipoin = 1,npoin
                            q_2(ipoin,:) = rho_2(ipoin)*u_2(ipoin,:)
                         end do

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

                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_2,Rmass_3)
                         rho_3(ldof) = rho(ldof)-(dt/1.0d0)*Rmass_3(ldof)
                         do ipoin = 1,npoin
                            q_3(ipoin,:) = rho_3(ipoin)*u_3(ipoin,:)
                         end do

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
                         e_int_4 = e_int

                         call mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q_3,Rmass_4)
                         aux = Rmass_1+2.0d0*Rmass_2+2.0d0*Rmass_3+Rmass_4
                         rho_4(ldof) = rho(ldof)-(dt/6.0d0)*aux(ldof)
                         do ipoin = 1,npoin
                            q_4(ipoin,:) = rho_4(ipoin)*u_4(ipoin,:)
                         end do

                         !
                         ! Update
                         !

                         rho(ldof) = rho_4(ldof)

                      end do

              end subroutine

end module time_integ
