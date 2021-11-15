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
                      real(8)                 :: Re(nnode)
                      real(8)                 :: tmp

                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Mass Convection")
                      !$acc parallel loop gang private(ind,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode)
                         end do
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp = dot_product(gpcar(idime,:,igaus,ielem),q(ind,idime))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                              (tmp)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(ind(inode)) = Rmass(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
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
                      real(8)                 :: Re(nnode,ndime), divgp, grpgp
                      real(8)                 :: tmp3(nnode), tmp2

                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Momentum convection")
                      !$acc parallel loop gang private(ind,tmp3,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            ind(inode) = connec(ielem,inode)
                            !$acc loop seq
                            do idime = 1,ndime
                               Re(inode,idime) = 0.0d0
                            end do
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !
                            ! Compute divergence(qu) and grad(p) at Gauss point
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               divgp = 0.0d0
                               grpgp = dot_product(gpcar(idime,:,igaus,ielem),pr(ind))
                               !$acc loop seq
                               do jdime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     tmp3(inode) = q(ind(inode),idime)*u(ind(inode),jdime) ! qi * uj
                                  end do
                                  tmp2 = dot_product(gpcar(jdime,:,igaus,ielem),tmp3(:))
                                  divgp = divgp+tmp2
                               end do
                               !
                               ! Quadrature
                               !
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                                    (divgp+grpgp)
                               end do
                            end do
                         end do
                         !
                         ! Final assembly
                         !
                         !$acc loop vector collapse(2)
                         do idime = 1,ndime
                            do inode = 1,nnode
                              !$acc atomic update
                              Rmom(ind(inode),idime) = Rmom(ind(inode),idime)+Re(inode,idime)
                              !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
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
                      real(8)                 :: el_fener(nnode,ndime)
                      real(8)                 :: tmp

                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Energy Convection")
                      !$acc parallel loop gang private(ind,Re,el_fener) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode  = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode)
                         end do
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               el_fener(inode,idime) = u(ind(inode),idime)*&
                                  (E(ind(inode))+pr(ind(inode)))
                            end do
                         end do
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp = dot_product(gpcar(idime,:,igaus,ielem),el_fener(:,idime))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
                                              (tmp)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rener(ind(inode)) = Rener(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
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
                      real(8)                          :: Re(nnode), el_a(nnode)
                      real(8)                          :: tmp1, tmp2

                      !$acc kernels
                      Rconvec(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Generic Convection")
                      !$acc parallel loop gang private(ind,Re,el_a) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode)
                         end do
                         if (present(alpha)) then
                            el_a(1:nnode) = alpha(ind)
                         else
                            el_a(:) = 1.0d0
                         end if
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            tmp1 = dot_product(Ngp(igaus,:),el_a(:))
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp2 = dot_product(gpcar(idime,:,igaus,ielem),q(ind,idime))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1* &
                                              (tmp2)
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         !$acc loop vector
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
