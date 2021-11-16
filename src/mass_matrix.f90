module mass_matrix

   use mod_nvtx

      contains

              subroutine consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc,weight)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Forms Mc as a sparse CSR matrix, utilizing nzdom, rdom and cdom to          !
                      ! compress the elemental matrices into the full sparse assembly structure.    !
                      ! Mc is defined by A{int(Na*Nb*det(Je))}, where A{} is the assembly operator. !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)           :: nelem, nnode, npoin, ngaus, nzdom
                      integer(4), intent(in)           :: connec(nelem,nnode), rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)           :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(in), optional :: weight(npoin)
                      real(8),    intent(out)          :: Mc(nzdom)
                      integer(4)                       :: ielem, igaus, inode, jnode, lnode(nnode), izdom, ipoin, jpoin
                      integer(4)                       :: jzdom, rowb, rowe
                      real(8)                          :: Me(nnode,nnode), el_w(nnode), tmp

                      !
                      ! Initialize Mc to zeros
                      !
                      call nvtxStartRange("Mass Matrix")
                      Mc = 0.0d0

                      !
                      ! Loop over all elements to form Mc_e(nnode,nnode)
                      !
                      do ielem = 1,nelem
                         lnode(1:nnode) = connec(ielem,1:nnode) ! get elemental indices
                         if(present(weight))then
                            el_w(1:nnode) = weight(lnode)
                         else
                            el_w(:) = 1.0d0
                         end if

                         !
                         ! Form Mc_e with Gaussian quadrature (open)
                         !
                         call nvtxStartRange("Elemental Matrix")
                         Me = 0.0d0
                         do igaus = 1,ngaus ! Loop over Gauss points
                            tmp = dot_product(Ngp(igaus,:),el_w(:))
                            do inode = 1,nnode ! Loop over element nodes (row)
                               do jnode = 1,nnode ! Loop over element nodex (column)
                                  Me(inode,jnode) = Me(inode,jnode) + &
                                     gpvol(1,igaus,ielem)*(tmp)* &
                                     Ngp(igaus,inode)*Ngp(igaus,jnode) ! Gaussian quad.
                               end do
                            end do
                         end do
                         call nvtxEndRange
                         !
                         ! Assemble Mc_e to CSR Mc
                         !
                         call nvtxStartRange("Assembly")
                         do inode = 1,nnode
                            ipoin = lnode(inode) ! Global node/Mc row index
                            rowb = rdom(ipoin)+1 ! Start of izdom for cdom
                            rowe = rdom(ipoin+1) ! end of izdom for cdom
                            do jnode = 1,nnode
                               jpoin = lnode(jnode) ! Nodes associated with ipoin on ielem
                               !
                               ! Loop over section of cdom to find out izdom
                               jzdom = rowb
                               do while ((cdom(jzdom) .ne. jpoin) .and. (jzdom .le. rowe))
                                  jzdom = jzdom+1
                               end do
                               Mc(jzdom) = Mc(jzdom) + Me(inode,jnode)
                            end do
                         end do
                         call nvtxEndRange
                      end do
                      call nvtxEndRange

              end subroutine consistent_mass

              subroutine lumped_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Ml,weight)

                      implicit none

                      integer(4), intent(in)           :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)           :: connec(nelem,nnode)
                      real(8),    intent(in)           :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                      real(8),    intent(in), optional :: weight(npoin)
                      real(8),    intent(out)          :: Ml(npoin)
                      integer(4)                       :: ielem, igaus, inode, jnode, ind(nnode)
                      real(8)                          :: Me(nnode), alpha, aux1, aux2, el_w(nnode)

                      !$acc kernels
                      Ml(:) = 0.0d0
                      !$acc end kernels

                      !$acc parallel loop gang private(ind,Me,el_w)
                      do ielem = 1,nelem
                         ind(1:nnode) = connec(ielem,1:nnode)
                         if (present(weight)) then
                            el_w(1:nnode) = weight(ind)
                         else
                            el_w(:) = 1.0d0
                         end if
                         Me = 0.0d0
                         aux1 = 0.0d0
                         !
                         ! tr[Mc]
                         !
                         !$acc loop vector reduction(+:aux1)
                         do inode = 1,nnode
                            !$acc loop seq
                            do igaus = 1,ngaus
                            aux1 = aux1+gpvol(1,igaus,ielem)*(dot_product(Ngp(igaus,:),el_w(:)))* &
                                   Ngp(igaus,inode)*Ngp(igaus,inode)
                            end do
                         end do
                         !
                         ! elem. mass
                         !
                         !!!$acc loop seq
                         !!do igaus = 1,ngaus
                         !!   aux2 = aux2+gpvol(1,igaus,ielem)
                         !!end do
                         !
                         ! Me
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop vector
                            do inode = 1,nnode
                               Me(inode) = Me(inode)+gpvol(1,igaus,ielem)*(Ngp(igaus,inode)**2)
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            aux2 = 0.0d0
                            !$acc loop seq
                            do igaus = 1,ngaus
                               aux2 = aux2+gpvol(1,igaus,ielem)
                            end do
                            !$acc atomic update
                            Ml(ind(inode)) = Ml(ind(inode))+((aux2/aux1)*Me(inode))
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop

              end subroutine lumped_mass

              subroutine cmass_times_vector(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,v,Rmc)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Forms Mc as a sparse CSR matrix, utilizing nzdom, rdom and cdom to          !
                      ! compress the elemental matrices into the full sparse assembly structure.    !
                      ! Mc is defined by A{int(Na*Nb*det(Je))}, where A{} is the assembly operator. !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode), v(npoin)
                      real(8),    intent(out) :: Rmc(npoin)
                      integer(4)              :: ielem, igaus, inode, ind(nnode)
                      real(8)                 :: Re(nnode), tmp2

                      !
                      ! Initialize Mc to zeros
                      !
                      call nvtxStartRange("Cmass times vector")
                      !$acc kernels
                      Rmc(:) = 0.0d0
                      !$acc end kernels

                      !$acc parallel loop gang private(ind,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode) ! get elemental indices
                         end do
                         !
                         ! Form Re with Gaussian quadrature (open)
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus ! Loop over Gauss points
                            tmp2 = dot_product(Ngp(igaus,:),v(ind))
                            !$acc loop vector
                            do inode = 1,nnode ! Loop over element nodes (row)
                               Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                           Ngp(igaus,inode)*tmp2
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmc(ind(inode)) = Rmc(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine cmass_times_vector

              subroutine wcmass_times_vector(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,v,Rmc,weight)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Forms Mc as a sparse CSR matrix, utilizing nzdom, rdom and cdom to          !
                      ! compress the elemental matrices into the full sparse assembly structure.    !
                      ! Mc is defined by A{int(Na*Nb*det(Je))}, where A{} is the assembly operator. !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, nnode, npoin, ngaus
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode), v(npoin)
                      real(8),    intent(in)  :: weight(npoin)
                      real(8),    intent(out) :: Rmc(npoin)
                      integer(4)              :: ielem, igaus, inode, ind(nnode)
                      real(8)                 :: Re(nnode), tmp1, tmp2

                      !
                      ! Initialize Mc to zeros
                      !
                      call nvtxStartRange("Cmass times vector")
                      !$acc kernels
                      Rmc(:) = 0.0d0
                      !$acc end kernels

                      !
                      ! Loop over all elements to form Mc_e(nnode,nnode)
                      !
                      !$acc parallel loop gang private(ind,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                            ind(inode) = connec(ielem,inode) ! get elemental indices
                         end do
                         !
                         ! Form Re with Gaussian quadrature (open)
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus ! Loop over Gauss points
                            tmp1 = dot_product(Ngp(igaus,:),weight(ind))
                            tmp2 = dot_product(Ngp(igaus,:),v(ind))
                            !$acc loop vector
                            do inode = 1,nnode ! Loop over element nodes (row)
                               Re(inode) = Re(inode) + &
                                  gpvol(1,igaus,ielem)*(tmp1)* &
                                  Ngp(igaus,inode)*tmp2
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmc(ind(inode)) = Rmc(ind(inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine wcmass_times_vector

end module mass_matrix
