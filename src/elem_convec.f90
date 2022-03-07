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
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode)
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Mass Convection")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            tmp1 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = tmp1+dot_product(gpcar(idime,:),q(connec(ielem,:),idime))
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
                            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re(inode)
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
                      integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
                      real(8)                 :: Re(nnode,ndime), divgp, grpgp
                      real(8)                 :: tmp1, tmp2, tmp3, gpcar(ndime,nnode)

                      call nvtxStartRange("Momentum convection")
                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               Re(inode,idime) = 0.0d0
                            end do
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = dot_product(gpcar(1,:),q(connec(ielem,:),1)*u(connec(ielem,:),1)+pr(connec(ielem,:)))
                            tmp1 = tmp1+dot_product(gpcar(2,:),q(connec(ielem,:),1)*u(connec(ielem,:),2))
                            tmp1 = tmp1+dot_product(gpcar(3,:),q(connec(ielem,:),1)*u(connec(ielem,:),3))
                            tmp2 = dot_product(gpcar(1,:),q(connec(ielem,:),2)*u(connec(ielem,:),1))
                            tmp2 = tmp2+dot_product(gpcar(2,:),q(connec(ielem,:),2)*u(connec(ielem,:),2)+pr(connec(ielem,:)))
                            tmp2 = tmp2+dot_product(gpcar(3,:),q(connec(ielem,:),2)*u(connec(ielem,:),3))
                            tmp3 = dot_product(gpcar(1,:),q(connec(ielem,:),3)*u(connec(ielem,:),1))
                            tmp3 = tmp3+dot_product(gpcar(2,:),q(connec(ielem,:),3)*u(connec(ielem,:),2))
                            tmp3 = tmp3+dot_product(gpcar(3,:),q(connec(ielem,:),3)*u(connec(ielem,:),3)+pr(connec(ielem,:)))
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
                              Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                              !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mom_convec

              subroutine mom_convec_emac(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr,Rmom)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = EMAC(A) using a standard Continuous    !
                      ! Galerkin model. A_i is q_i/sqrt(rho), and the EMAC term is       !
                      ! defined as:                                                      !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: Aemac(npoin,ndime), Femac(npoin), pr(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
                      real(8)                 :: Re(nnode,ndime)
                      real(8)                 :: aux(ndime), gpcar(ndime,nnode)
                      real(8)                 :: tmp1, tmp2, tmp3
                      real(8)                 :: grad_1, grad_2, grad_3
                      real(8)                 :: grad_4, grad_5, grad_6
                      real(8)                 :: grad_7, grad_8, grad_9
                      real(8)                 :: gradp_1, gradp_2, gradp_3, div_a

                      call nvtxStartRange("EMAC Momentum convection")

                      !
                      ! Start global RHS
                      !
                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      
                      !
                      ! Start elemental ops
                      !
                      !$acc parallel loop gang private(Re,gpcar,aux) vector_length(32)
                      do ielem = 1,nelem
                         !
                         ! Initialize element vector to 0
                         !
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               Re(inode,idime) = 0.0d0
                            end do
                         end do
                         !
                         ! Loop over Gauss points
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !
                            ! Create GPCAR(ndime,nnode) for each element at each Gauss point
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !
                            ! Compute grad(A) and div(A)
                            !
                            !
                            ! Gradient structure:
                            !
                            !         | u1,1 u1,2 u1,3 |
                            ! u_i,j = | u2,1 u2,2 u2,3 |
                            !         | u3,1 u3,2 u3,3 |
                            !
                            grad_1 = 0.0d0 
                            grad_2 = 0.0d0 
                            grad_3 = 0.0d0 
                            grad_4 = 0.0d0 
                            grad_5 = 0.0d0 
                            grad_6 = 0.0d0 
                            grad_7 = 0.0d0 
                            grad_8 = 0.0d0 
                            grad_9 = 0.0d0 
                            !$acc loop vector &
                            !$acc reduction(+:grad_1,grad_2,grad_3,grad_4,grad_5,grad_6,grad_7,grad_8,grad_9)
                            do inode = 1,nnode
                               grad_1 = grad_1+(gpcar(1,inode)*Aemac(connec(ielem,inode),1))
                               grad_2 = grad_2+(gpcar(2,inode)*Aemac(connec(ielem,inode),1))
                               grad_3 = grad_3+(gpcar(3,inode)*Aemac(connec(ielem,inode),1))
                               grad_4 = grad_4+(gpcar(1,inode)*Aemac(connec(ielem,inode),2))
                               grad_5 = grad_5+(gpcar(2,inode)*Aemac(connec(ielem,inode),2))
                               grad_6 = grad_6+(gpcar(3,inode)*Aemac(connec(ielem,inode),2))
                               grad_7 = grad_7+(gpcar(1,inode)*Aemac(connec(ielem,inode),3))
                               grad_8 = grad_8+(gpcar(2,inode)*Aemac(connec(ielem,inode),3))
                               grad_9 = grad_9+(gpcar(3,inode)*Aemac(connec(ielem,inode),3))
                            end do
                            !
                            ! div(A) = tr[grad(A)]
                            !
                            div_a = grad_1 + grad_5 + grad_9
                            !
                            ! Pressure gradient
                            !
                            gradp_1 = 0.0d0
                            gradp_2 = 0.0d0
                            gradp_3 = 0.0d0
                            !$acc loop vector &
                            !$acc reduction(+:gradp_1,gradp_2,gradp_3)
                            do inode = 1,nnode
                               gradp_1 = gradp_1+(gpcar(1,inode)*pr(connec(ielem,inode)))
                               gradp_2 = gradp_2+(gpcar(2,inode)*pr(connec(ielem,inode)))
                               gradp_3 = gradp_3+(gpcar(3,inode)*pr(connec(ielem,inode)))
                            end do
                            !
                            ! Interpolate A innto an auxiliary array
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               aux(idime) = dot_product(Ngp(igaus,:),Aemac(connec(ielem,:),idime))
                            end do
                            !
                            ! Compute terms grad(A)*A+grad^T(A)*A+div(A)*A+grad(P)
                            !
                            tmp1 = (grad_1+grad_1+div_a)*aux(1)+(grad_2+grad_4)*aux(2)+(grad_3+grad_7)*aux(3)+gradp_1
                            tmp2 = (grad_4+grad_2)*aux(1)+(grad_5+grad_5+div_a)*aux(2)+(grad_6+grad_8)*aux(3)+gradp_2
                            tmp3 = (grad_7+grad_3)*aux(1)+(grad_8+grad_6)*aux(2)+(grad_9+grad_9+div_a)*aux(3)+gradp_3
                            !
                            ! Subtract kinetic energy component
                            !
                            tmp1 = tmp1 - 0.5d0*dot_product(gpcar(1,:),Femac(connec(ielem,:)))
                            tmp2 = tmp2 - 0.5d0*dot_product(gpcar(2,:),Femac(connec(ielem,:)))
                            tmp3 = tmp3 - 0.5d0*dot_product(gpcar(3,:),Femac(connec(ielem,:)))
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
                              Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                              !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mom_convec_emac

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
                      integer(4)              :: ielem, igaus, inode, idime, ipoin, jdime
                      real(8)                 :: Re(nnode), fener(npoin,ndime)
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Energy Convection")
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = tmp1+dot_product(gpcar(idime,:),u(connec(ielem,:),idime)*(E(connec(ielem,:))+pr(connec(ielem,:))))
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
                            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re(inode)
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
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: tmp1, tmp2, Re(nnode), gpcar(ndime,nnode)

                      call nvtxStartRange("Generic Convection")
                      !$acc kernels
                      Rconvec(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            tmp1 = dot_product(Ngp(igaus,:),alpha(connec(ielem,:)))
                            tmp2 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp2 = tmp2+dot_product(gpcar(idime,:),q(connec(ielem,:),idime))
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
                            Rconvec(connec(ielem,inode)) = Rconvec(connec(ielem,inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine generic_scalar_convec

end module elem_convec
