module elem_convec

      use mod_nvtx

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes convective term for Euler/NS equation system, as well             !
      ! as for any generic scalar transport that might occur. Based                !
      ! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
      ! This module can be passed to CUDA in order to do fine-grained parallelism. !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

              subroutine mass_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,q,Rmass)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(rho*u) using a standard Continuous  !
                      ! Galerkin formulation. In the above, rho is the scalar sensity and !
                      ! u is the vector of velocities in all dimensions.                  !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode)
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Mass Convection")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(ind,Re,gpcar) vector_length(32)
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
                            tmp1 = 0.0d0
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = tmp1+dot_product(gpcar(idime,:),q(ind,idime))
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                           Ngp(igaus,inode)*tmp1
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

              subroutine mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u,q,pr,Rmom)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(q*u) using a standard Continuous   !
                      ! Galerkin formulation. In the above, q is the momentum vector and !
                      ! u is the vector of velocities in all dimensions. The product     !
                      ! inside div() is an outer tensor product between the 2 vectors.   !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), pr(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
                      real(8)                 :: Re(nnode,ndime), divgp, grpgp
                      real(8)                 :: tmp1, tmp2, tmp3, gpcar(ndime,nnode)

                      call nvtxStartRange("Momentum convection")
                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(ind,Re,gpcar) vector_length(32)
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
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = dot_product(gpcar(1,:),q(ind,1)*u(ind,1)+pr(ind))
                            tmp1 = tmp1+dot_product(gpcar(2,:),q(ind,1)*u(ind,2))
                            tmp1 = tmp1+dot_product(gpcar(3,:),q(ind,1)*u(ind,3))
                            tmp2 = dot_product(gpcar(1,:),q(ind,2)*u(ind,1))
                            tmp2 = tmp2+dot_product(gpcar(2,:),q(ind,2)*u(ind,2)+pr(ind))
                            tmp2 = tmp2+dot_product(gpcar(3,:),q(ind,2)*u(ind,3))
                            tmp3 = dot_product(gpcar(1,:),q(ind,3)*u(ind,1))
                            tmp3 = tmp3+dot_product(gpcar(2,:),q(ind,3)*u(ind,2))
                            tmp3 = tmp3+dot_product(gpcar(3,:),q(ind,3)*u(ind,3)+pr(ind))
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode,1) = Re(inode,1)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1
                               Re(inode,2) = Re(inode,2)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp2
                               Re(inode,3) = Re(inode,3)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp3
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

              subroutine ener_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u,pr,E,Rener)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(u*(E+p)) using a standard Continuous !
                      ! Galerkin formulation. In the above, E is the scalar total energy,  !
                      ! p is the scalar pressure and u is the vector of velocities in all  !
                      ! dimensions.                                                        !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), pr(npoin), E(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, idime, ipoin, jdime
                      real(8)                 :: Re(nnode), fener(npoin,ndime)
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Energy Convection")
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(ind,Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode)
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = tmp1+dot_product(gpcar(idime,:),u(ind,idime)*(E(ind)+pr(ind)))
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1
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

              subroutine generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp, &
                                               dNgp,He,gpvol,q,Rconvec,alpha)

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(in)  :: alpha(npoin)
                      real(8),    intent(out) :: Rconvec(npoin)
                      integer(4)              :: ind(nnode)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: tmp1, tmp2, Re(nnode), gpcar(ndime,nnode)

                      call nvtxStartRange("Generic Convection")
                      !$acc kernels
                      Rconvec(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(ind,Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode)
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = dot_product(Ngp(igaus,:),alpha(ind))
                            tmp2 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp2 = tmp2+dot_product(gpcar(idime,:),q(ind,idime))
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                           Ngp(igaus,inode)*tmp1*tmp2
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
