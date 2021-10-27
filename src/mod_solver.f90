module mod_solver

      use mod_nvtx

      contains

              subroutine lumped_solver_scal(npoin,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin
                      real(8),    intent(in)    :: Ml(npoin)
                      real(8),    intent(inout) :: R(npoin)

                      R = R/Ml

              end subroutine lumped_solver_scal

              subroutine lumped_solver_vect(npoin,ndime,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin, ndime
                      real(8),    intent(in)    :: Ml(npoin)
                      real(8),    intent(inout) :: R(npoin,ndime)
                      integer(4)                :: idime

                      do idime = 1,ndime
                         R(:,idime) = R(:,idime)/Ml
                      end do

              end subroutine lumped_solver_vect

              subroutine approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,R)

                      implicit none

                      integer(4), intent(in)       :: npoin, nzdom, ppow
                      integer(4), intent(in)       :: rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)       :: Ml(npoin), Mc(nzdom)
                      real(8),    intent(inout)    :: R(npoin)
                      integer(4)                   :: ipoin, jpoin, ipow, izdom, rowb, rowe
                      real(8),    dimension(npoin) :: b, v, x
                      real(8)                      :: Ar(nzdom)

                      !
                      ! Compute Ar
                      !
                      !
                      call nvtxStartRange("Scalar APINV")

                      Ar = Mc

                      do ipoin = 1,npoin
                         rowb = rdom(ipoin)+1
                         rowe = rdom(ipoin+1)
                         do izdom = rowb,rowe
                            if(cdom(izdom) == ipoin) then
                               Ar(izdom) = Ml(ipoin)-Ar(izdom)
                            else
                               Ar(izdom) = -Ar(izdom)
                            end if
                            Ar(izdom) = Ar(izdom)/Ml(ipoin)
                         end do
                      end do

                      !
                      ! Initialize series at k=0
                      !
                      b = R
                      v = b
                      x = b

                      !
                      ! Step over sucessive powers
                      !
                      do ipow = 1,ppow
                         call CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Ar,v,b)
                         v = b
                         x = x+v
                      end do
                      R = x
                      call nvtxEndRange

              end subroutine approx_inverse_scalar

              subroutine approx_inverse_vect(ndime,npoin,nzdom,rdom,cdom,ppow,Ml,Mc,R)

                      implicit none

                      integer(4), intent(in)       :: ndime, npoin, nzdom, ppow
                      integer(4), intent(in)       :: rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)       :: Ml(npoin), Mc(nzdom)
                      real(8),    intent(inout)    :: R(npoin,ndime)
                      integer(4)                   :: idime, ipoin, jpoin, ipow, izdom, rowb, rowe
                      real(8),    dimension(npoin) :: b, v, x
                      real(8)                      :: Ar(nzdom)

                      !
                      ! Compute Ar
                      !
                      call nvtxStartRange("Vector APINV")
                      Ar = Mc
                      do ipoin = 1,npoin
                         rowb = rdom(ipoin)+1
                         rowe = rdom(ipoin+1)
                         do izdom = rowb,rowe
                            if(cdom(izdom) == ipoin) then
                               Ar(izdom) = Ml(ipoin)-Ar(izdom)
                            else
                               Ar(izdom) = -Ar(izdom)
                            end if
                            Ar(izdom) = Ar(izdom)/Ml(ipoin)
                         end do
                      end do

                      do idime = 1,ndime
                         !
                         ! Initialize series at k=0
                         !
                         b = R(:,idime)
                         v = b
                         x = b

                         !
                         ! Step over sucessive powers
                         !
                         do ipow = 1,ppow
                            call CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Ar,v,b)
                            v = b
                            x = x+v
                         end do
                         R(:,idime) = x
                      end do
                      call nvtxEndRange

              end subroutine approx_inverse_vect

              subroutine CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,v,u)

                      implicit none

                      integer(4), intent(in)  :: npoin, nzdom
                      integer(4), intent(in)  :: rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)  :: Mc(nzdom), v(npoin)
                      real(8),    intent(out) :: u(npoin)
                      integer(4)              :: ipoin, izdom, jpoin, rowb, rowe

                      call nvtxStartRange("SPMV")
                      u = 0.0d0
                      do ipoin = 1,npoin
                         !
                         ! Get CSR section for row ipoin
                         !
                         rowb = rdom(ipoin)+1
                         rowe = rdom(ipoin+1)
                         !
                         ! Loop inside CSR section
                         !
                         do izdom = rowb,rowe
                            jpoin = cdom(izdom) ! Col. index
                            u(ipoin) = u(ipoin)+Mc(izdom)*v(jpoin) ! Dot product
                         end do
                      end do
                      call nvtxEndRange

              end subroutine CSR_SpMV_scal

end module mod_solver
