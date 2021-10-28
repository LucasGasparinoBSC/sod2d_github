module elem_convec

      use mod_nvtx

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
                      real(8)                 :: tmp

                      Rmass = 0.0d0
                      call nvtxStartRange("Mass Convection")
                      !$acc parallel loop gang
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         !
                         ! Quadrature
                         !
                         !$acc loop vector collapse(2)
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               tmp = dot_product(gpcar(idime,:,igaus,ielem),el_q(:,idime))
                               !$acc loop seq
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                              (tmp)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(ind(inode)) = Rmass(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      call nvtxEndRange

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
                      real(8)                 :: Re(nnode,ndime), divgp(ndime,ngaus), grpgp(ndime,ngaus)
                      real(8)                 :: el_q(nnode,ndime), el_u(nnode,ndime), el_pr(nnode)
                      real(8)                 :: tmp1,tmp2,tmp3(nnode)

                      Rmom = 0.0d0
                      call nvtxStartRange("Momentum convection")
                      !!$acc parallel loop gang
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_pr(1:nnode) = pr(ind)
                         !!$acc loop vector
                         do igaus = 1,ngaus
                            !
                            ! Compute divergence(qu) and grad(p) at Gauss point
                            !
                            divgp = 0.0d0
                            grpgp = 0.0d0
                            do idime = 1,ndime
                               tmp1 = dot_product(gpcar(idime,:,igaus,ielem),el_pr(:))
                               grpgp(idime,igaus) = grpgp(idime,igaus)+tmp1
                               !!$acc loop seq
                               do jdime = 1,ndime
                                  tmp3(1:nnode) = el_q(1:nnode,idime)*el_u(1:nnode,jdime) ! qi * uj
                                  tmp2 = dot_product(gpcar(jdime,:,igaus,ielem),tmp3(:))
                                  divgp(idime,igaus) = divgp(idime,igaus) + tmp2
                               end do
                               !
                               ! Quadrature
                               !
                               !!$acc loop seq
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
                            do inode = 1,nnode
                              !!$acc atomic update
                              Rmom(ind(inode),idime) = Rmom(ind(inode),idime)+Re(inode,idime)
                              !!$acc end atomic
                            end do
                         end do
                      end do
                      !!$acc end parallel loop
                      call nvtxEndRange

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
                      real(8)                 :: el_fener(nnode,ndime),el_u(nnode,ndime)
                      real(8)                 :: el_E(nnode),el_pr(nnode)
                      real(8)                 :: tmp

                      Rener = 0.0d0
                      call nvtxStartRange("Energy Convection")
                      !$acc parallel loop gang
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_u(1:nnode,1:ndime) = u(ind,1:ndime)
                         el_E(1:nnode) = E(ind)
                         el_pr(1:nnode) = pr(ind)
                         !
                         ! Quadrature
                         !
                         !$acc loop vector collapse(2)
                         do igaus = 1,ngaus
                            do idime = 1,ndime
                               el_fener(1:nnode,idime) = el_u(1:nnode,idime)*(el_E(1:nnode)+el_pr(1:nnode))
                               tmp = dot_product(gpcar(idime,:,igaus,ielem),el_fener(:,idime))
                               !$acc loop seq
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                              (tmp)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         do inode = 1,nnode
                            !$acc atomic update
                            Rener(ind(inode)) = Rener(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      call nvtxEndRange

              end subroutine ener_convec

              subroutine generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rconvec,alpha)

                      implicit none

                      integer(4), intent(in)           :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)           :: connec(nelem,nnode)
                      real(8),    intent(in)           :: Ngp(ngaus,nnode)
                      real(8),    intent(in)           :: gpcar(ndime,nnode,ngaus,nelem)
                      real(8),    intent(in)           :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)           :: q(npoin,ndime)
                      real(8),    intent(in), optional :: alpha(npoin)
                      real(8),    intent(out)          :: Rconvec(npoin)
                      integer(4)                       :: ind(nnode)
                      integer(4)                       :: ielem, igaus, inode, jnode, idime
                      real(8)                          :: Re(nnode), el_q(nnode,ndime), el_a(nnode)
                      real(8)                          :: tmp1, tmp2

                      Rconvec = 0.0d0
                      call nvtxStartRange("Generic Convection")
                      !$acc parallel loop gang
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         el_q(1:nnode,1:ndime) = q(ind,1:ndime)
                         if (present(alpha)) then
                            el_a(1:nnode) = alpha(ind)
                         else
                            el_a(:) = 1.0d0
                         end if
                         !
                         ! Quadrature
                         !
                         !$acc loop vector
                         do igaus = 1,ngaus
                            tmp1 = dot_product(Ngp(igaus,:),el_a(:))
                            do idime = 1,ndime
                               tmp2 = dot_product(gpcar(idime,:,igaus,ielem),el_q(:,idime))
                               !$acc loop seq
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1* &
                                              (tmp2)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         do inode = 1,nnode
                            !$acc atomic update
                            Rconvec(ind(inode)) = Rconvec(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine generic_scalar_convec

end module elem_convec
