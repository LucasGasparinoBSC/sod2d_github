module quadrature_rules

        contains

                subroutine gll_qua(ndime,ngaus,xgp,wgp)

                        implicit none

                        integer(4), intent(in)  :: ndime, ngaus
                        real(8),    intent(out) :: xgp(ngaus,ndime), wgp(ngaus)
                        real(8)                 :: q

                        if (ngaus == 1) then
                           xgp(1,1) = 0.0d0
                           xgp(1,2) = 0.0d0
                           wgp(1) = 4.0d0
                        else if (ngaus == 4) then
                           q = 1.0d0/sqrt(3.0d0)
                           xgp(1,1) = -q
                           xgp(1,2) = -q
                           xgp(2,1) =  q
                           xgp(2,2) = -q
                           xgp(3,1) =  q
                           xgp(3,2) =  q
                           xgp(4,1) = -q
                           xgp(4,2) =  q
                           wgp(1) = 1.0d0
                           wgp(2) = 1.0d0
                           wgp(3) = 1.0d0
                           wgp(4) = 1.0d0
                        else if (ngaus == 9) then
                           write(*,*) 'NOT CODED YET!'
                           xgp = 0.0d0
                           wgp = 0.0d0
                        end if

                end subroutine gll_qua

                !subroutine gll_hex()
                !        implicit none
                !        if (ngaus == 1) then
                !        else if (ngaus == 8) then
                !        else if (ngaus == 27) then
                !        end if
                !end subroutine gll_hex

end module quadrature_rules
