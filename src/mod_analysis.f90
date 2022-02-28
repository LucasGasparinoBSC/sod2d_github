module mod_analysis

   contains

      subroutine volAvg_EK(nelem,npoin,ndime,nnode,ngaus,connec,gpvol,Ngp,rho0,rho,u,EK)

         implicit none

         integer(4), intent(in)  :: nelem, npoin, ndime, nnode, ngaus, connec(nelem,nnode)
         real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
         real(8),    intent(in)  :: rho0, rho(npoin), u(npoin,ndime)
         real(8),    intent(out) :: EK
         integer(4)              :: ielem, igaus, inode
         real(8)                 :: R1

         EK = 0.0d0
         do ielem = 1,nelem
            R1 = 0.0d0
            do igaus = 1,ngaus
               do inode = 1,nnode
                  R1 = R1 + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*rho(connec(ielem,inode))* &
                     dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))
               end do
            end do
            EK = EK+R1
         end do
         EK = EK/(rho0*((2.0d0*3.14159d0)**3))

      end subroutine volAvg_EK

      subroutine write_EK(time,EK)

         implicit none

         real(8), intent(in) :: time, EK

         write(666,*) time, EK

      end subroutine write_EK

end module mod_analysis
