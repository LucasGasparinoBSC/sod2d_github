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

                subroutine gll_hex(ndime,ngaus,xgp,wgp)

                        implicit none

                        integer(4), intent(in)  :: ndime, ngaus
                        real(8),    intent(out) :: xgp(ngaus,ndime), wgp(ngaus)
                        real(8)                 :: q, sl, s0, sr, ws, w0

                        if (ngaus == 1) then
                           xgp(1,1) = 0.0d0
                           xgp(1,2) = 0.0d0
                           xgp(1,3) = 0.0d0
                           wgp(1) = 4.0d0
                        else if (ngaus == 8) then
                           q = 1.0d0/sqrt(3.0d0)
                           xgp(1,1) = -q ! xi
                           xgp(1,2) = -q ! eta
                           xgp(1,3) = -q ! zeta
                           xgp(2,1) = -q
                           xgp(2,2) = -q
                           xgp(2,3) =  q
                           xgp(3,1) =  q
                           xgp(3,2) = -q
                           xgp(3,3) =  q
                           xgp(4,1) =  q
                           xgp(4,2) = -q
                           xgp(4,3) = -q
                           xgp(5,1) = -q ! xi
                           xgp(5,2) =  q ! eta
                           xgp(5,3) = -q ! zeta
                           xgp(6,1) = -q
                           xgp(6,2) =  q
                           xgp(6,3) =  q
                           xgp(7,1) =  q
                           xgp(7,2) =  q
                           xgp(7,3) =  q
                           xgp(8,1) =  q
                           xgp(8,2) =  q
                           xgp(8,3) = -q
                           wgp(1) = 1.0d0
                           wgp(2) = 1.0d0
                           wgp(3) = 1.0d0
                           wgp(4) = 1.0d0
                           wgp(5) = 1.0d0
                           wgp(6) = 1.0d0
                           wgp(7) = 1.0d0
                           wgp(8) = 1.0d0
                        else if (ngaus == 27) then
                           sl = -sqrt(0.6d0)
                           s0 = 0.0d0
                           sr = sqrt(0.6d0)
                           xgp(1,1:3) = [sl,sl,sl]
                           xgp(2,1:3) = [sl,sl,s0]
                           xgp(3,1:3) = [sl,sl,sr]
                           xgp(4,1:3) = [s0,sl,sl]
                           xgp(5,1:3) = [s0,sl,s0]
                           xgp(6,1:3) = [s0,sl,sr]
                           xgp(7,1:3) = [sr,sl,sl]
                           xgp(8,1:3) = [sr,sl,s0]
                           xgp(9,1:3) = [sr,sl,sr]
                           xgp(10,1:3) = [sl,s0,sl]
                           xgp(11,1:3) = [sl,s0,s0]
                           xgp(12,1:3) = [sl,s0,sr]
                           xgp(13,1:3) = [s0,s0,sl]
                           xgp(14,1:3) = [s0,s0,s0]
                           xgp(15,1:3) = [s0,s0,sr]
                           xgp(16,1:3) = [sr,s0,sl]
                           xgp(17,1:3) = [sr,s0,s0]
                           xgp(18,1:3) = [sr,s0,sr]
                           xgp(19,1:3) = [sl,sr,sl]
                           xgp(20,1:3) = [sl,sr,s0]
                           xgp(21,1:3) = [sl,sr,sr]
                           xgp(22,1:3) = [s0,sr,sl]
                           xgp(23,1:3) = [s0,sr,s0]
                           xgp(24,1:3) = [s0,sr,sr]
                           xgp(25,1:3) = [sr,sr,sl]
                           xgp(26,1:3) = [sr,sr,s0]
                           xgp(27,1:3) = [sr,sr,sr]

                           ws = 5.0d0/9.0d0
                           w0 = 8.0d0/9.0d0
                           wgp(1) = ws*ws*ws
                           wgp(2) = ws*ws*w0
                           wgp(3) = ws*ws*ws
                           wgp(4) = w0*ws*ws
                           wgp(5) = w0*ws*w0
                           wgp(6) = w0*ws*ws
                           wgp(7) = ws*ws*ws
                           wgp(8) = ws*ws*w0
                           wgp(9) = ws*ws*ws
                           wgp(10) = ws*w0*ws
                           wgp(11) = ws*w0*w0
                           wgp(12) = ws*w0*ws
                           wgp(13) = w0*w0*ws
                           wgp(14) = w0*w0*w0
                           wgp(15) = w0*w0*ws
                           wgp(16) = ws*w0*ws
                           wgp(17) = ws*w0*w0
                           wgp(18) = ws*w0*ws
                           wgp(19) = ws*ws*ws
                           wgp(20) = ws*ws*w0
                           wgp(21) = ws*ws*ws
                           wgp(22) = w0*ws*ws
                           wgp(23) = w0*ws*w0
                           wgp(24) = w0*ws*ws
                           wgp(25) = ws*ws*ws
                           wgp(26) = ws*ws*w0
                           wgp(27) = ws*ws*ws
                        else if (ngaus == 64) then
                           write(*,*) 'NOT CODED YET!'
                           xgp = 0.0d0
                           wgp = 0.0d0
                        end if

                end subroutine gll_hex


end module quadrature_rules
