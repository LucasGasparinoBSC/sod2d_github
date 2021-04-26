program unitt_CSR_mc

        use elem_qua
        use mod_graph
        use quadrature_rules
        use mass_matrix

        implicit none

        integer(4)              :: ndime, nelem, npoin, nnode, ngaus, nzdom
        integer(4)              :: idime, ielem, ipoin, inode, igaus, izdom
        integer(4), allocatable :: connec(:,:), rdom(:), cdom(:), aux_cdom(:)
        real(8)                 :: s, t
        real(8), allocatable    :: wgp(:), xgp(:,:)
        real(8), allocatable    :: N(:), dN(:,:)
        real(8), allocatable    :: Ngp(:,:), dNgp(:,:,:)
        real(8), allocatable    :: gpvol(:,:,:)
        real(8), allocatable    :: Mc(:)

        !
        ! Mesh data
        !
        ndime = 2
        nelem = 2
        npoin = 6
        nnode = 4
        ngaus = 4

        !
        ! Element connectivity
        !
        allocate(connec(nelem,npoin))
        connec(1,1:4) = (/1, 2, 5, 6/)
        connec(2,1:4) = (/2, 3, 4, 5/)

        !
        ! Compute mesh graph (nzdom, rdom, cdom)
        !
        allocate(rdom(npoin+1))
        allocate(aux_cdom(nelem*nnode*nnode))
        call compute_nzdom(npoin,nnode,nelem,connec,nzdom,rdom,aux_cdom)
        allocate(cdom(nzdom))
        do izdom = 1,nzdom
           cdom(izdom) = aux_cdom(izdom)
        end do
        deallocate(aux_cdom)

        !
        ! Gaussian table
        !
        allocate(wgp(ngaus))
        allocate(xgp(ngaus,ndime))
        call gll_qua(ndime,ngaus,xgp,wgp) ! QUA_XX quadraturea

        !
        ! Isopar. data
        !
        allocate(N(nnode),dN(ndime,nnode))                ! dN is dummy
        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus)) ! dNgp is dummy
        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           call qua04(s,t,N,dN)
           Ngp(igaus,:) = N
           dNgp(:,:,igaus) = dN
        end do

        !
        ! Compute dummy gpvol
        !
        allocate(gpvol(1,ngaus,nelem))
        do ielem = 1,nelem
           do igaus = 1,ngaus
              gpvol(1,igaus,ielem) = wgp(igaus) ! Assumes det(Je) = 1.0d0
           end do
        end do

        !
        ! Compute Mc CSR
        !
        allocate(Mc(nzdom))
        call consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc)

        if (sum(Mc) .gt. 8.0d0+0.00000000001d0) then
                write(*,*) '--| ERROR IN FORMING Mc CSR!'
                stop 1
        else if (sum(Mc) .lt. 8.0d0-0.00000000001d0) then
                write(*,*) '--| ERROR IN FORMING Mc CSR!'
                stop 1
        end if

end program unitt_CSR_mc
