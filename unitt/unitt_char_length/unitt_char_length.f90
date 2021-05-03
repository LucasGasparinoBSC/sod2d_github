program unitt_char_length

        use mod_geom

        implicit none

        integer(4)              :: ndime, nelem, npoin, nnode
        integer(4)              :: ielem
        integer(4), allocatable :: connec(:,:)
        real(8)                 :: he
        real(8),    allocatable :: coord(:,:), helem(:)

        !
        ! Basic data
        !
        ndime = 2
        nelem = 2
        nnode = 4
        npoin = 6

        !
        ! Mesh table
        !
        allocate(connec(nelem,nnode))
        connec(1,:) = (/1, 2, 5, 6/)
        connec(2,:) = (/2, 3, 4, 5/)

        !
        ! Coordinates
        !
        allocate(coord(npoin,ndime))
        coord(1,:) = (/1.0d0, 0.0d0/)
        coord(2,:) = (/2.0d0, 0.0d0/)
        coord(3,:) = (/4.0d0, 0.0d0/)
        coord(4,:) = (/4.0d0, 4.0d0/)
        coord(5,:) = (/2.0d0, 4.0d0/)
        coord(6,:) = (/0.0d0, 4.0d0/)

        !
        ! Compute characteristic lengths for each element
        !
        allocate(helem(nelem))
        do ielem = 1,nelem
           call char_length(ielem,nelem,nnode,npoin,ndime,connec,coord,he)
           print*, 'ielem, he || ', ielem, he
           helem(ielem) = he
        end do

        if (helem(1) .ne. 1.0d0 .or. helem(2) .ne. 2.0d0) then
           write(*,*) '--| TEST FAILED!'
           STOP 1
        end if

end program unitt_char_length
