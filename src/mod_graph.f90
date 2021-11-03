module mod_graph

      use mod_nvtx

      contains

              subroutine compute_nzdom(npoin,nnode,nelem,connec,nzdom,rdom,aux_cdom)

                      implicit none

                      integer(4), intent(in)  :: npoin, nnode, nelem
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      integer(4), intent(out) :: nzdom, rdom(npoin+1), aux_cdom(nelem*nnode*nnode)
                      integer(4)              :: ielem, inode, jnode, ipoin, jpoin, count_repeated
                      integer(4)              :: aux_inc, lnode(nnode), aux1(nnode)
                      integer(4)              :: list_used(npoin)
                      integer(4)              :: nnz_loc

                      call nvtxStartRange("Compute graph")
                      nzdom = 0    ! non-zero entries initialized to zero
                      rdom = 0     ! Row idexes for CSR
                      aux_cdom = 0 ! Auxiliary col. indexes when nzdom is still unknown

                      !
                      ! Loop over line entries
                      !
                      do ipoin = 1,npoin ! Row index
                         list_used = 0 ! Initialize used nodes for this entry to 0
                         !
                         ! Loop over elements
                         !
                         do ielem = 1,nelem
                            !
                            ! Verify if node is on connectivity of ielem; if no entries match, procveed to next element
                            !
                            aux1(1:nnode) = connec(ielem,1:nnode)
                            do inode = 1,nnode
                               if (aux1(inode) == ipoin) then ! ipoin matches entry on connectivity
                                  lnode = aux1
                                  !
                                  ! loop over element nodes
                                  !
                                  count_repeated = 0
                                  nnz_loc = 0 ! Counter for non-zeros on line, used to advance in aux_cdom
                                  do jnode = 1,nnode
                                     jpoin = lnode(jnode) ! Column index
                                     !
                                     ! Check if node has already been used on previous elements (connections)
                                     !
                                     if (list_used(jpoin) .ne. 0) then
                                        count_repeated = count_repeated+1
                                     else if (list_used(jpoin) == 0) then
                                        count_repeated = count_repeated
                                        nnz_loc = nnz_loc+1
                                        aux_cdom(nzdom+nnz_loc) = jpoin
                                     end if
                                  end do
                                  aux_inc = count_repeated
                                  nzdom = nzdom+(nnode-aux_inc)
                                  rdom(ipoin+1) = nzdom
                                  list_used(lnode) = 1
                                  exit ! Move to next element
                               end if
                            end do
                         end do
                      end do
                      call nvtxEndRange

              end subroutine compute_nzdom

end module mod_graph
