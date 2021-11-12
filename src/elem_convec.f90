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
                      real(8)                 :: tmp1, tmp2

                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Mass Convection")
                      !$acc parallel loop gang private(ind,Re) vector_length(32)
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            tmp2 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tmp1 = dot_product(dNgp(jdime,:,igaus),q(ind,idime))
                                  tmp2 = tmp2+He(idime,jdime,igaus,ielem)*tmp1
                               end do
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp2
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

              !
              ! Old version
              !
              !!subroutine mom_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,dNgp,He,gpvol,u,q,pr,Rmom)

              !!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!        ! Subroutine to compute R = div(q*u) using a standard Continuous   !
              !!        ! Galerkin formulation. In the above, q is the momentum vector and !
              !!        ! u is the vector of velocities in all dimensions. The product     !
              !!        ! inside div() is an outer tensor product between the 2 vectors.   !
              !!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !!        implicit none

              !!        integer(4), intent(in)  :: nelem, ngaus, npoin, nnode, ndime
              !!        integer(4), intent(in)  :: connec(nelem,nnode)
              !!        real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
              !!        real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
              !!        real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
              !!        real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), pr(npoin)
              !!        real(8),    intent(out) :: Rmom(npoin,ndime)
              !!        integer(4)              :: ind(nnode)
              !!        integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
              !!        real(8)                 :: Re(nnode,ndime), divgp, grpgp
              !!        real(8)                 :: tmp2, tmp1

              !!        Rmom = 0.0d0
              !!        call nvtxStartRange("Momentum convection")
              !!        !$acc parallel loop gang private(ind,Re) vector_length(32)
              !!        do ielem = 1,nelem
              !!           Re = 0.0d0
              !!           ind = connec(ielem,:)
              !!           !$acc loop seq
              !!           do igaus = 1,ngaus
              !!              !
              !!              ! Compute divergence(qu) and grad(p) at Gauss point
              !!              !
              !!              !$acc loop seq
              !!              do idime = 1,ndime
              !!                 divgp = 0.0d0
              !!                 grpgp = 0.0d0
              !!                 !$acc loop seq
              !!                 do jdime = 1,ndime
              !!                    !tmp1 = dot_product(dNgp(jdime,:,igaus),pr(ind))
              !!                    !grpgp = grpgp+He(idime,jdime,igaus,ielem)*tmp1
              !!                    !$acc loop seq
              !!                    do kdime = 1,ndime
              !!                       if (jdime .eq. idime) then
              !!                          tmp2 = dot_product(dNgp(kdime,:,igaus),q(ind,idime)*u(ind,jdime)+pr(ind))
              !!                       else
              !!                          tmp2 = dot_product(dNgp(kdime,:,igaus),q(ind,idime)*u(ind,jdime))
              !!                       end if
              !!                       divgp = divgp + He(jdime,kdime,igaus,ielem)*tmp2
              !!                    end do
              !!                 end do
              !!                 !
              !!                 ! Quadrature
              !!                 !
              !!                 !$acc loop vector
              !!                 do inode = 1,nnode
              !!                    Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)* &
              !!                                      (divgp)
              !!                 end do
              !!              end do
              !!           end do
              !!           !
              !!           ! Final assembly
              !!           !
              !!           !$acc loop vector collapse(2)
              !!           do idime = 1,ndime
              !!              do inode = 1,nnode
              !!                !$acc atomic update
              !!                Rmom(ind(inode),idime) = Rmom(ind(inode),idime)+Re(inode,idime)
              !!                !$acc end atomic
              !!              end do
              !!           end do
              !!        end do
              !!        !$acc end parallel loop
              !!        call nvtxEndRange

              !!end subroutine mom_convec

              !
              ! Unrolled version
              !
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
                      integer(4)              :: ind
                      integer(4)              :: ielem, igaus, idime, jdime, inode, kdime
                      real(8)                 :: grc_1, grc_2, grc_3
                      real(8)                 :: f1,f2,f3,f4,f5,f6,f7,f8,f9
                      real(8)                 :: grad_111, grad_121,grad_131
                      real(8)                 :: grad_211, grad_221,grad_231
                      real(8)                 :: grad_311, grad_321,grad_331
                      real(8)                 :: grad_112, grad_122,grad_132
                      real(8)                 :: grad_212, grad_222,grad_232
                      real(8)                 :: grad_312, grad_322,grad_332
                      real(8)                 :: grad_113, grad_123,grad_133
                      real(8)                 :: grad_213, grad_223,grad_233
                      real(8)                 :: grad_313, grad_323,grad_333

                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Momentum convection")
                      !$acc parallel loop gang vector_length(32)
                      do ielem = 1,nelem
                         !
                         ! Compute Fij = qi*uj+p*dij
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            grad_111 = 0.0d0
                            grad_121 = 0.0d0
                            grad_131 = 0.0d0
                            grad_211 = 0.0d0
                            grad_221 = 0.0d0
                            grad_231 = 0.0d0
                            grad_311 = 0.0d0
                            grad_321 = 0.0d0
                            grad_331 = 0.0d0
                            grad_112 = 0.0d0
                            grad_122 = 0.0d0
                            grad_132 = 0.0d0
                            grad_212 = 0.0d0
                            grad_222 = 0.0d0
                            grad_232 = 0.0d0
                            grad_312 = 0.0d0
                            grad_322 = 0.0d0
                            grad_332 = 0.0d0
                            grad_113 = 0.0d0
                            grad_123 = 0.0d0
                            grad_133 = 0.0d0
                            grad_213 = 0.0d0
                            grad_223 = 0.0d0
                            grad_233 = 0.0d0
                            grad_313 = 0.0d0
                            grad_323 = 0.0d0
                            grad_333 = 0.0d0
                            !$acc loop vector &
                            !$acc reduction(+:grad_111,grad_121,grad_131) &
                            !$acc reduction(+:grad_211,grad_221,grad_231) &
                            !$acc reduction(+:grad_311,grad_321,grad_331) &
                            !$acc reduction(+:grad_112,grad_122,grad_132) &
                            !$acc reduction(+:grad_212,grad_222,grad_232) &
                            !$acc reduction(+:grad_312,grad_322,grad_332) &
                            !$acc reduction(+:grad_113,grad_123,grad_133) &
                            !$acc reduction(+:grad_213,grad_223,grad_233) &
                            !$acc reduction(+:grad_313,grad_323,grad_333)
                            do inode = 1,nnode
                               ind = connec(ielem,inode)
                               f1 = q(ind,1)*u(ind,1)+pr(ind)
                               f2 = q(ind,1)*u(ind,2)
                               f3 = q(ind,1)*u(ind,3)
                               f4 = q(ind,2)*u(ind,1)
                               f5 = q(ind,2)*u(ind,2)+pr(ind)
                               f6 = q(ind,2)*u(ind,3)
                               f7 = q(ind,3)*u(ind,1)
                               f8 = q(ind,3)*u(ind,2)
                               f9 = q(ind,3)*u(ind,3)+pr(ind)
                               grad_111 = grad_111+dNgp(1,inode,igaus)*f1
                               grad_121 = grad_121+dNgp(1,inode,igaus)*f2
                               grad_131 = grad_131+dNgp(1,inode,igaus)*f3
                               grad_211 = grad_211+dNgp(1,inode,igaus)*f4
                               grad_221 = grad_221+dNgp(1,inode,igaus)*f5
                               grad_231 = grad_231+dNgp(1,inode,igaus)*f6
                               grad_311 = grad_311+dNgp(1,inode,igaus)*f7
                               grad_321 = grad_321+dNgp(1,inode,igaus)*f8
                               grad_331 = grad_331+dNgp(1,inode,igaus)*f9
                               grad_112 = grad_112+dNgp(2,inode,igaus)*f1
                               grad_122 = grad_122+dNgp(2,inode,igaus)*f2
                               grad_132 = grad_132+dNgp(2,inode,igaus)*f3
                               grad_212 = grad_212+dNgp(2,inode,igaus)*f4
                               grad_222 = grad_222+dNgp(2,inode,igaus)*f5
                               grad_232 = grad_232+dNgp(2,inode,igaus)*f6
                               grad_312 = grad_312+dNgp(2,inode,igaus)*f7
                               grad_322 = grad_322+dNgp(2,inode,igaus)*f8
                               grad_332 = grad_332+dNgp(2,inode,igaus)*f9
                               grad_113 = grad_113+dNgp(3,inode,igaus)*f1
                               grad_123 = grad_123+dNgp(3,inode,igaus)*f2
                               grad_133 = grad_133+dNgp(3,inode,igaus)*f3
                               grad_213 = grad_213+dNgp(3,inode,igaus)*f4
                               grad_223 = grad_223+dNgp(3,inode,igaus)*f5
                               grad_233 = grad_233+dNgp(3,inode,igaus)*f6
                               grad_313 = grad_313+dNgp(3,inode,igaus)*f7
                               grad_323 = grad_323+dNgp(3,inode,igaus)*f8
                               grad_333 = grad_333+dNgp(3,inode,igaus)*f9
                            end do
                            grc_1 = He(1,1,igaus,ielem)*grad_111+He(1,2,igaus,ielem)*grad_112+He(1,3,igaus,ielem)*grad_113
                            grc_1 = grc_1+He(2,1,igaus,ielem)*grad_121+He(2,2,igaus,ielem)*grad_122+He(2,3,igaus,ielem)*grad_123
                            grc_1 = grc_1+He(3,1,igaus,ielem)*grad_131+He(3,2,igaus,ielem)*grad_132+He(3,3,igaus,ielem)*grad_133
                            grc_2 = He(1,1,igaus,ielem)*grad_211+He(1,2,igaus,ielem)*grad_212+He(1,3,igaus,ielem)*grad_213
                            grc_2 = grc_2+He(2,1,igaus,ielem)*grad_221+He(2,2,igaus,ielem)*grad_222+He(2,3,igaus,ielem)*grad_223
                            grc_2 = grc_2+He(3,1,igaus,ielem)*grad_231+He(3,2,igaus,ielem)*grad_232+He(3,3,igaus,ielem)*grad_233
                            grc_3 = He(1,1,igaus,ielem)*grad_311+He(1,2,igaus,ielem)*grad_312+He(1,3,igaus,ielem)*grad_313
                            grc_3 = grc_3+He(2,1,igaus,ielem)*grad_321+He(2,2,igaus,ielem)*grad_322+He(2,3,igaus,ielem)*grad_323
                            grc_3 = grc_3+He(3,1,igaus,ielem)*grad_331+He(3,2,igaus,ielem)*grad_332+He(3,3,igaus,ielem)*grad_333
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            ind = connec(ielem,inode)
                            !$acc atomic update
                            Rmom(ind,1) = Rmom(ind,1)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*grc_1
                            !$acc end atomic
                            !$acc atomic update
                            Rmom(ind,2) = Rmom(ind,2)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*grc_2
                            !$acc end atomic
                            !$acc atomic update
                            Rmom(ind,3) = Rmom(ind,3)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*grc_3
                            !$acc end atomic
                         end do
                         !!!
                         !!! Final assembly
                         !!!
                         !!!$acc loop vector collapse(2)
                         !!do idime = 1,ndime
                         !!   do inode = 1,nnode
                         !!     ind = connec(ielem,inode)
                         !!     !$acc atomic update
                         !!     Rmom(ind,idime) = Rmom(ind,idime)+Re(inode,idime)
                         !!     !$acc end atomic
                         !!   end do
                         !!end do
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
                      real(8)                 :: el_fener(nnode,ndime)
                      real(8)                 :: tmp1, tmp2

                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Energy Convection")
                      !$acc parallel loop gang private(ind,Re,el_fener) vector_length(32)
                      do ielem = 1,nelem
                         Re = 0.0d0
                         ind = connec(ielem,:)
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
                            tmp2 = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tmp1 = dot_product(dNgp(jdime,:,igaus),el_fener(:,idime))
                                  tmp2 = tmp2+He(idime,jdime,igaus,ielem)*tmp1
                               end do
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp2
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
                      integer(4)              :: ind
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: tmp1, tmp2, tmp3

                      !$acc kernels
                      Rconvec(:) = 0.0d0
                      !$acc end kernels
                      call nvtxStartRange("Generic Convection")
                      !$acc parallel loop gang vector_length(32)
                      do ielem = 1,nelem
                         !
                         ! Quadrature
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            tmp3 = 0.0d0
                            tmp1 = dot_product(Ngp(igaus,:),alpha(connec(ielem,:)))
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tmp2 = dot_product(dNgp(jdime,:,igaus),q(connec(ielem,:),idime))
                                  tmp3 = tmp3+He(idime,jdime,igaus,ielem)*tmp2
                               end do
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               ind = connec(ielem,inode)
                               !$acc atomic update
                               Rconvec(ind) = Rconvec(ind)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1*tmp2
                               !$acc end atomic
                            end do
                         end do
                         !!
                         !! Assembly
                         !!
                         !!$acc loop vector
                         !do inode = 1,nnode
                         !   !$acc atomic update
                         !   Rconvec(ind(inode)) = Rconvec(ind(inode))+Re(inode)
                         !   !$acc end atomic
                         !end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine generic_scalar_convec

end module elem_convec
