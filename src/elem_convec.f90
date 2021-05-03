module elem_convec

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes convective term for Euler/NS equation system, as well             !
      ! as for any generic scalar transport that might occur. Based                !
      ! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
      ! This module can be passed to CUDA in order to do fine-grained parallelism. !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

              subroutine mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rmass)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(rho*u) using a standard Continuous  !
                      ! Galerkin formulation. In the above, rho is the scalar sensity and !
                      ! u is the vector of velocities in all dimensions.                  !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), el_q(nnode,ndime)

                      Rmass = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         !
                         ! Quadrature
                         !
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_q(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         Rmass(ind) = Rmass(ind)+Re(1:nnode)
                      end do

              end subroutine mass_convec

              subroutine mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,q,pr,Rmom)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(q*u) using a standard Continuous   !
                      ! Galerkin formulation. In the above, q is the momentum vector and !
                      ! u is the vector of velocities in all dimensions. The product     !
                      ! inside div() is an outer tensor product between the 2 vectors.   !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), pr(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode
                      real(8)                 :: Re(nnode,ndime), aux, divgp(ndime,ngaus), grpgp(ndime,ngaus)
                      real(8)                 :: el_q(nnode,ndime), el_u(nnode,ndime), el_pr(nnode)

                      Rmom = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_pr(1:nnode) = pr(ind)
                         do igaus = 1,ngaus
                            !
                            ! Compute divergence(qu) and grad(p) at Gauss point
                            !
                            divgp = 0.0d0
                            grpgp = 0.0d0
                            do idime = 1,ndime
                               do jnode = 1,nnode
                                  do jdime = 1,ndime
                                     aux = el_q(jnode,idime)*el_u(jnode,jdime) ! qi * uj
                                     divgp(idime,igaus) = divgp(idime,igaus) + &
                                             gpcar(jdime,jnode,igaus,ielem)*aux
                                  end do
                                  grpgp(idime,igaus) = grpgp(idime,igaus)+gpcar(idime,jnode,igaus,ielem)*el_pr(jnode)
                               end do
                               !
                               ! Quadrature
                               !
                               do inode = 1,nnode
                                  Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                          (divgp(idime,igaus)+grpgp(idime,igaus))
                               end do
                            end do
                         end do
                         !
                         ! Final assembly
                         !
                         do idime = 1,ndime
                            Rmom(ind,idime) = Rmom(ind,idime)+Re(1:nnode,idime)
                         end do
                      end do

              end subroutine mom_convec

              subroutine ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,u,pr,E,Rener)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(u*(E+p)) using a standard Continuous !
                      ! Galerkin formulation. In the above, E is the scalar total energy,  !
                      ! p is the scalar pressure and u is the vector of velocities in all  !
                      ! dimensions.                                                        !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), pr(npoin), E(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime, ipoin
                      real(8)                 :: Re(nnode), fener(npoin,ndime)
                      real(8)                 :: el_fener(nnode,ndime)

                      !
                      ! Flux f = u*(E+p)
                      !
                      do ipoin = 1,npoin
                         fener(ipoin,1:ndime) = u(ipoin,1:ndime)*(E(ipoin)+pr(ipoin))
                      end do

                      Rener = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_fener(1:nnode,1:ndime) = fener(ind,1:ndime)
                         !
                         ! Quadrature
                         !
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_fener(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         Rener(ind) = Rener(ind)+Re(1:nnode)
                      end do

              end subroutine ener_convec

              subroutine generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rconvec)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode)
                      real(8),    intent(in)  :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(out) :: Rconvec(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, jnode, idime
                      real(8)                 :: Re(nnode), el_q(nnode,ndime)

                      Rconvec = 0.0d0
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         !
                         ! Quadrature
                         !
                         do igaus = 1,ngaus
                            do inode = 1,nnode
                               do idime = 1,ndime
                                  do jnode = 1,nnode
                                     Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                 (gpcar(idime,jnode,igaus,ielem)*el_q(jnode,idime))
                                  end do
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         Rconvec(ind) = Rconvec(ind)+Re(1:nnode)
                      end do

              end subroutine generic_scalar_convec

end module elem_convec
