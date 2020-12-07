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

              subroutine approx_inverse_scalar()
                      implicit none       
              end subroutine approx_inverse_scalar

end module mod_solver
