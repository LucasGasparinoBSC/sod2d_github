module mod_period

   contains

      subroutine periodic_ops(nelem,npoin,npoin_orig,nnode,nper,lpoin,connec,connec_orig,masSla)

         implicit none

         integer(4), intent(in)               :: nelem, nnode, nper, masSla(nper,2)
         integer(4), intent(out)              :: npoin_orig, connec_orig(nelem,nnode)
         integer(4), intent(inout)            :: npoin, connec(nelem,nnode)
         integer(4), allocatable, intent(out) :: lpoin(:)
         integer(4)                           :: ielem, inode, ipoin, iper, counter
         integer(4), allocatable              :: aux1(:)

         !
         ! Copy connec to dummy
         !
         !$acc kernels
         connec_orig(:,:) = connec(:,:)
         !$acc end kernels

         !
         ! Modify connec with Master/Slave relations
         !
         !$acc parallel loop gang
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
         ! Recompute npoin without periodic nodes
         !
         npoin_orig = npoin
         npoin = npoin-nper

         !
         ! Create list of working nodes (exclude slaves from list)
         !
         allocate(aux1(npoin_orig))
         allocate(lpoin(npoin))

         !$acc parallel loop
         do ipoin = 1,npoin_orig
            aux1(ipoin) = ipoin
         end do
         !$acc end parallel loop
         
         do iper = 1,nper
            do ipoin = 1,npoin_orig
               if (masSla(iper,2) .eq. ipoin) then
                  aux1(ipoin) = 0
                  exit
               end if
            end do
         end do

         counter = 0
         do ipoin = 1,npoin_orig
            if (aux1(ipoin) .ne. 0) then
               counter = counter+1
               lpoin(counter) = aux1(ipoin)
            end if
         end do
         deallocate(aux1)

      end subroutine periodic_ops

end module mod_period
