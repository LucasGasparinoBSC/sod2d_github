module mod_period

   contains

      subroutine periodic_ops(nelem,npoin,nboun,npbou,npoin_w,nnode,nper, &
                              lpoin_w,connec,connec_orig,masSla,bound,bound_orig)

         implicit none

         integer(4), intent(in)               :: nelem, npoin, nnode, nboun, npbou, nper, masSla(nper,2)
         integer(4), intent(out)              :: npoin_w, connec_orig(nelem,nnode)
         integer(4), intent(inout)            :: connec(nelem,nnode)
         integer(4), optional, intent(out)    :: bound_orig(nboun,npbou)
         integer(4), optional, intent(inout)  :: bound(nboun,npbou)
         integer(4), allocatable, intent(out) :: lpoin_w(:)
         integer(4)                           :: ielem, inode, iboun, ipbou, ipoin, iper, counter
         integer(4), allocatable              :: aux1(:)

         !
         ! Copy connec and bound to dummy
         !
         !$acc kernels
         connec_orig(:,:) = connec(:,:)
         !$acc end kernels
         
         if(present(bound) .and. present(bound_orig)) then
            !$acc kernels
            bound_orig(:,:) = bound(:,:)
            !$acc end kernels
         end if

         !
         ! Modify connec with Master/Slave relations
         !
         !$acc parallel loop gang vector_length(32)
         do ielem = 1,nelem
            !$acc loop vector
            do inode = 1,nnode
               !$acc loop seq
               do iper = 1,nper
                  if (connec(ielem,inode) .eq. masSla(iper,2)) then
                     connec(ielem,inode) = masSla(iper,1)
                  end if
               end do
            end do
         end do
         !$acc end parallel loop

         !
         ! Modify bound with Master/Slave relations
         !
         if(present(bound)) then
            !$acc parallel loop gang
            do iboun = 1,nboun
               !$acc loop vector
               do ipbou = 1,npbou
                  !$acc loop seq
                  do iper = 1,nper
                     if (bound(iboun,ipbou) .eq. masSla(iper,2)) then
                        bound(iboun,ipbou) = masSla(iper,1)
                     end if
                  end do
               end do
            end do
            !$acc end parallel loop
         end if

         !
         ! Compute npoin_w (npoin without periodic nodes)
         !
         npoin_w = npoin-nper

         !
         ! Create list of working nodes (exclude slaves from list)
         !
         allocate(aux1(npoin))
         allocate(lpoin_w(npoin_w))

         !$acc parallel loop
         do ipoin = 1,npoin
            aux1(ipoin) = ipoin
         end do
         !$acc end parallel loop
         
         do iper = 1,nper
            do ipoin = 1,npoin
               if (masSla(iper,2) .eq. ipoin) then
                  aux1(ipoin) = 0
                  exit
               end if
            end do
         end do

         counter = 0
         do ipoin = 1,npoin
            if (aux1(ipoin) .ne. 0) then
               counter = counter+1
               lpoin_w(counter) = aux1(ipoin)
            end if
         end do
         deallocate(aux1)

      end subroutine periodic_ops

end module mod_period
