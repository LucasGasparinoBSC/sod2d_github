program unitt_CSR_apinv

        use mod_solver

        implicit none

        integer(4)              :: npoin, nzdom, ipoin, ppow
        integer(4), allocatable :: rdom(:), cdom(:)
        real(8), allocatable    :: Ml(:), Mc(:), R(:)

        !
        ! Basic data
        !
        npoin = 6   ! Nodes on mesh
        nzdom = 28  ! Non-zero entries
        ppow = 1000 ! APINV order (iterations)

        !
        ! Create CSR matrix
        !
        allocate(Mc(nzdom))
        allocate(Ml(npoin))
        allocate(cdom(nzdom))
        allocate(rdom(npoin+1))

        !
        ! Lumped mass
        !
        Ml(1:6) = (/9.0d0, 18.0d0, 9.0d0, 9.0d0, 18.0d0, 9.0d0/)

        !
        ! Consistent mass
        !
        Mc(1:10)  = (/4.0d0, 2.0d0, 1.0d0, 2.0d0, 2.0d0, 8.0d0, 4.0d0, 1.0d0, 2.0d0, 1.0d0/)
        Mc(11:20) = (/2.0d0, 4.0d0, 2.0d0, 1.0d0, 1.0d0, 2.0d0, 4.0d0, 2.0d0, 1.0d0, 4.0d0/)
        Mc(21:28) = (/8.0d0, 2.0d0, 1.0d0, 2.0d0, 2.0d0, 1.0d0, 2.0d0, 4.0d0/)

        cdom(1:10)  = (/1, 2, 5, 6, 1, 2, 5, 6, 3, 4/)
        cdom(11:20) = (/2, 3, 4, 5, 2, 3, 4, 5, 1, 2/)
        cdom(21:28) = (/5, 6, 3, 4, 1, 2, 5, 6/)

        rdom(1:7) = (/0, 4, 10, 14, 18, 24, 28/)

        !
        ! Create vector
        !
        allocate(R(npoin))
        R(1:6) = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0/)

        !
        ! Pre-multiply R with inv(Ml)
        !
        do ipoin = 1,npoin
           R(ipoin) = R(ipoin)/Ml(ipoin)
        end do

        !
        ! Call CSR apinv solver
        !
        call approx_inverse_scalar(npoin,nzdom,rdom,cdom,ppow,Ml,Mc,R)
        write(*,*) '--| sum(R) = ',sum(R)

        if (sum(R) .gt. 2.333333333d0+0.001d0) then
                write(*,*) '--| sum(R) = ',sum(R)
                write(*,*) '--| APINV SOLVER FAILED!'
                STOP 1
        else if (sum(R) .lt. 2.333333333d0-0.001d0) then
                write(*,*) '--| sum(R) = ',sum(R)
                write(*,*) '--| APINV SOLVER FAILED!'
                STOP 1
        end if

end program unitt_CSR_apinv
