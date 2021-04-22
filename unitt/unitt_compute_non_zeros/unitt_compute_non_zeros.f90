program unitt_compute_non_zeros

      use mod_graph

      implicit none

      integer(4) :: nelem, npoin, nnode, nzdom, izdom
      integer(4), allocatable :: connec(:,:)
      integer(4), allocatable :: rdom(:), cdom(:), aux_cdom(:)

      !
      ! Basic data
      !
      nelem = 3  ! Number of elements
      npoin = 8 ! Number of nodes
      nnode = 4  ! Nodes per element

      !
      ! Connectivity table
      !
      allocate(connec(nelem,nnode))
      connec(1,:) = (/1, 2, 5, 8/)
      connec(2,:) = (/2, 3, 4, 5/)
      connec(3,:) = (/8, 5, 6, 7/)

      !!connec(1,:) = (/1, 2, 13, 12/)
      !!connec(2,:) = (/2, 3, 14, 13/)
      !!connec(3,:) = (/3, 4, 5, 14/)
      !!connec(4,:) = (/12, 13, 16, 11/)
      !!connec(5,:) = (/13, 14, 15, 16/)
      !!connec(6,:) = (/14, 5, 6, 15/)
      !!connec(7,:) = (/11, 16, 9, 10/)
      !!connec(8,:) = (/16, 15, 8, 9/)
      !!connec(9,:) = (/15, 6, 7, 8/)

      !
      ! Allocate graph data
      !
      allocate(rdom(npoin+1))
      allocate(aux_cdom(nelem*nnode*nnode))

      !
      ! Call subroutine to compute nzdom
      !
      call compute_nzdom(npoin,nnode,nelem,connec,nzdom,rdom,aux_cdom)

      !
      ! Transfer aux_cdom to cdom
      !
      allocate(cdom(nzdom))
      do izdom = 1,nzdom
         cdom(izdom) = aux_cdom(izdom)
      end do
      deallocate(aux_cdom)

      !
      ! Check
      !
      write(*,*) 'Number of non-zero entries is := ',nzdom
      if (nzdom .ne. 40) then
         write(*,*) '--| SUBROUTINE RETURNED WRONG NUMBER OF NONZERO ENTRIES!'
         stop 1
      else if (rdom(npoin+1) .ne. nzdom) then
         write(*,*) '--| BAD RDOM ARRAY!'
         STOP 1
      else if (cdom(1) .ne. 1 .or. cdom(10) .ne. 4 .or. cdom(nzdom) .ne. 7) then
         write(*,*) '--| BAD CDOM ARRAY!'
         STOP 1
      end if

end program unitt_compute_non_zeros
