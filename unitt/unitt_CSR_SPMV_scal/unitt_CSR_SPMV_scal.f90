program unitt_CSR_SPMV_scal

        use mod_solver

        implicit none

        integer(4)              :: npoin, nzdom
        integer(4), allocatable :: rdom(:), cdom(:)
        real(8), allocatable    :: Mc(:), u(:), v(:)

        !
        ! Basic data
        !
        npoin = 6
        nzdom = 28

        !
        ! Create CSR matrix
        !
        allocate(Mc(nzdom))
        allocate(cdom(nzdom))
        allocate(rdom(npoin+1))

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
        allocate(v(npoin))
        v(1:6) = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0/)

        !
        ! Call SpMV subroutine
        !
        allocate(u(npoin))
        call CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,v,u)

        if (sum(u) .gt. 252.0d0+0.0000000001d0) then
                write(*,*) '--| sum(u) = ',sum(u)
                write(*,*) '--| SPMV FAILED!'
        else if (sum(u) .lt. 252.0d0-0.0000000001d0) then
                write(*,*) '--| sum(u) = ',sum(u)
                write(*,*) '--| SPMV FAILED!'
        end if

end program unitt_CSR_SPMV_scal
