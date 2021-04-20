module mod_graph

      contains

              subroutine compute_nzdom(npoin,nnode,nelem,connec,nzdom)

                      implicit none

                      integer(4), intent(in)  :: npoin, nnode, nelem
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      integer(4), intent(out) :: nzdom
                      integer(4)              :: ielem, inode, jnode, ipoin, jpoin, count_repeated
                      integer(4)              :: aux_inc, lnode(nnode), aux1(nnode)
                      integer(4)              :: list_used(npoin)

                      nzdom = 0 ! non-zero entries initialized to zero

                      !
                      ! Loop over line entries
                      !
                      do ipoin = 1,npoin
                         print*, 'ipoin = ',ipoin
                         list_used = 0 ! Initialize used nodes for this entry to 0
                         !
                         ! Loop over elements
                         !
                         do ielem = 1,nelem
                            print*, '   ielem = ',ielem
                            !
                            ! Verify if node is on connectivity of ielem; if no entries match, procveed to next element
                            !
                            aux1(1:nnode) = connec(ielem,1:nnode)
                            do inode = 1,nnode
                               if (aux1(inode) == ipoin) then ! ipoin matches entry on connectivity
                                  print*, '      inode = ',inode
                                  lnode = aux1
                                  !
                                  ! loop over element nodes
                                  !
                                  count_repeated = 0
                                  do jnode = 1,nnode
                                     jpoin = lnode(jnode)
                                     !
                                     ! Check if node has already been used on previous elements (connections)
                                     !
                                     if (list_used(jpoin) .ne. 0) then
                                        print*, 'Found repeated node = ',jpoin
                                        count_repeated = count_repeated+1
                                        print*, count_repeated
                                     else if (list_used(jpoin) == 0) then
                                        count_repeated = count_repeated
                                     end if
                                  end do
                                  aux_inc = count_repeated
                                  nzdom = nzdom+(nnode-aux_inc)
                                  print*, '   nzdom = ',nzdom
                                  list_used(lnode) = 1
                                  exit ! Move to next element
                               end if
                            end do
                         end do
                      end do

              end subroutine compute_nzdom

end module mod_graph
