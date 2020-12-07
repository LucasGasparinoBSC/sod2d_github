module mod_solver

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

              subroutine approx_inverse_scalar(npoin,Ml,Mc,R)

                      implicit none

                      integer(4), intent(in)       :: npoin
                      real(8),    intent(in)       :: Ml(npoin), Mc(npoin,npoin)
                      real(8),    intent(inout)    :: R(npoin)
                      integer(4)                   :: ipoin, jpoin, ipow, ppow
                      real(8),    dimension(npoin) :: b, v, x
                      real(8)                      :: Ar(npoin,npoin)

                      !
                      ! Compute Ar
                      !
                      Ar = Mc
                      do ipoin = 1,npoin
                         do jpoin = 1,npoin
                            if (jpoin == ipoin) then
                               Ar(ipoin,jpoin) = 1.0d0-(Ar(ipoin,jpoin)/Ml(ipoin))
                            else
                               Ar(ipoin,jpoin) = -(Ar(ipoin,jpoin)/Ml(ipoin))
                            end if
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
                      ppow = 10
                      do ipow = 1,ppow
                         v = matmul(Ar,v)
                         x = x+v
                      end do
                      R = x

              end subroutine approx_inverse_scalar

              subroutine approx_inverse_vect(ndime,npoin,Ml,Mc,R)

                      implicit none

                      integer(4), intent(in)       :: ndime, npoin
                      real(8),    intent(in)       :: Ml(npoin), Mc(npoin,npoin)
                      real(8),    intent(inout)    :: R(npoin,ndime)
                      integer(4)                   :: idime, ipoin, jpoin, ipow, ppow
                      real(8),    dimension(npoin) :: b, v, x
                      real(8)                      :: Ar(npoin,npoin)

                      !
                      ! Compute Ar
                      !
                      Ar = Mc
                      do ipoin = 1,npoin
                         do jpoin = 1,npoin
                            if (jpoin == ipoin) then
                               Ar(ipoin,jpoin) = 1.0d0-(Ar(ipoin,jpoin)/Ml(ipoin))
                            else
                               Ar(ipoin,jpoin) = -(Ar(ipoin,jpoin)/Ml(ipoin))
                            end if
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
                         ppow = 10
                         do ipow = 1,ppow
                            v = matmul(Ar,v)
                            x = x+v
                         end do
                         R(:,idime) = x
                      end do

              end subroutine approx_inverse_vect

end module mod_solver
