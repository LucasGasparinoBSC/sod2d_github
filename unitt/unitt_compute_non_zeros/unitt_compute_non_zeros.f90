program unitt_compute_non_zeros

      use mod_graph

      implicit none

      integer(4) :: nelem, npoin, nnode, nzdom
      integer(4), allocatable :: connec(:,:)

      !
      ! Basic data
      !
      nelem = 9  ! Number of elements
      npoin = 16 ! Number of nodes
      nnode = 4  ! Nodes per element

      !
      ! Connectivity table
      !
      allocate(connec(nelem,nnode))
      connec(1,:) = (/1, 2, 13, 12/)
      connec(2,:) = (/2, 3, 14, 13/)
      connec(3,:) = (/3, 4, 5, 14/)
      connec(4,:) = (/12, 13, 16, 11/)
      connec(5,:) = (/13, 14, 15, 16/)
      connec(6,:) = (/14, 5, 6, 15/)
      connec(7,:) = (/11, 16, 9, 10/)
      connec(8,:) = (/16, 15, 8, 9/)
      connec(9,:) = (/15, 6, 7, 8/)

      !
      ! Call subroutine to compute nzdom
      !
      call compute_nzdom(npoin,nnode,nelem,connec,nzdom)
      write(*,*) 'Number of non-zero entries is := ',nzdom
      if (nzdom .ne. 100) then
         write(*,*) '--| SUBROUTINE RETURNED WRONG NUMBER OF NONZERO ENTRIES!'
         stop 1
      end if

end program unitt_compute_non_zeros
