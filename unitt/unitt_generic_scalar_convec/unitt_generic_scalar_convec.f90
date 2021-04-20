program unitt_generic_scalar_convec

        !*********************************************************************!
        ! Unit test for generic scalar convection subroutine. Uses hard-coded !
        ! information to be independent from other subroutines.               !
        !*********************************************************************!

        use elem_convec

        implicit none

        integer(4)                 :: ndime, nnode, ngaus
        integer(4)                 :: idime, inode, igaus
        integer(4)                 :: nelem, npoin
        integer(4)                 :: ielem, ipoin
        integer(4), allocatable    :: connec(:,:)
        real(8),    allocatable    :: coord(:,:), elcod(:,:)
        real(8),    allocatable    :: xgp(:,:), wgp(:)
        real(8),    allocatable    :: N(:), dN(:,:), Ngp(:,:), dNgp(:,:,:)
        real(8),    allocatable    :: Je(:,:), He(:,:)
        real(8),    allocatable    :: struc_J(:,:,:,:), struc_H(:,:,:,:), struc_detJ(:,:,:)
        real(8),    allocatable    :: dxN(:,:), gpcar(:,:,:,:), gpvol(:,:,:)
        real(8),    allocatable    :: q(:,:), Rconvec(:)
        real(8)                    :: s, t, detJe

        !*********************************************************************!
        ! Basic data                                                          !
        !*********************************************************************!

        ndime = 2
        nnode = 4
        nelem = 2
        npoin = 6
        ngaus = 4

        allocate(connec(nelem,nnode))
        allocate(coord(npoin,ndime))

        connec(1,:) = (/1,2,5,6/)
        connec(2,:) = (/2,3,4,5/)

        coord(1,:) = (/0.0d0, 0.0d0/)
        coord(2,:) = (/2.0d0, 0.0d0/)
        coord(3,:) = (/3.0d0, 0.0d0/)
        coord(4,:) = (/3.5d0, 3.0d0/)
        coord(5,:) = (/2.0d0, 2.0d0/)
        coord(6,:) = (/0.0d0, 2.0d0/)

        !*********************************************************************!
        ! Allocate variables                                                  !
        !*********************************************************************!

        !
        ! Last rank is for prediction-advance related to entropy viscosity,
        ! where 1 is prediction, 2 is final value
        !
        allocate(q(npoin,ndime))  ! momentum
        q(1,:) = (/1.0d0, 1.0d0/)
        q(2,:) = (/1.5d0, 1.5d0/)
        q(3,:) = (/1.0d0, 1.0d0/)
        q(4,:) = (/1.0d0, 1.0d0/)
        q(5,:) = (/1.5d0, 1.5d0/)
        q(6,:) = (/1.0d0, 1.0d0/)

        !*********************************************************************!
        ! Generate GLL table                                                  !
        !*********************************************************************!

        allocate(xgp(ngaus,ndime))
        allocate(wgp(ngaus))

        xgp(1,:) = (/-1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0)/)
        xgp(2,:) = (/1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0)/)
        xgp(3,:) = (/1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0)/)
        xgp(4,:) = (/-1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0)/)
        wgp(1) = 1.0d0
        wgp(2) = 1.0d0
        wgp(3) = 1.0d0
        wgp(4) = 1.0d0

        !*********************************************************************!
        ! Generate N and dN for all GP                                        !
        !*********************************************************************!

        allocate(N(nnode),dN(ndime,nnode))
        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))

        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           N(1) = (1.0d0-s)*(1.0d0-t)
           N(2) = (1.0d0+s)*(1.0d0-t)
           N(3) = (1.0d0+s)*(1.0d0+t)
           N(4) = (1.0d0-s)*(1.0d0+t)
           N = 0.25d0*N
           dN(1,1) = (-1.0d0+t)
           dN(1,2) = ( 1.0d0-t)
           dN(1,3) = ( 1.0d0+t)
           dN(1,4) = (-1.0d0-t)
           dN(2,1) = (-1.0d0+s)
           dN(2,2) = (-1.0d0-s)
           dN(2,3) = ( 1.0d0+s)
           dN(2,4) = ( 1.0d0-s)
           dN = 0.25d0*dN
           Ngp(igaus,:) = N
           dNgp(:,:,igaus) = dN
        end do

        !*********************************************************************!
        ! Generate Jacobian related information                               !
        !*********************************************************************!

        allocate(elcod(ndime,nnode))
        allocate(Je(ndime,ndime))
        allocate(He(ndime,ndime))
        allocate(struc_J(ndime,ndime,ngaus,nelem))
        allocate(struc_H(ndime,ndime,ngaus,nelem))
        allocate(struc_detJ(1,ngaus,nelem))
        allocate(gpvol(1,ngaus,nelem))
        allocate(dxN(ndime,nnode))
        allocate(gpcar(ndime,nnode,ngaus,nelem))

        do ielem = 1,nelem
           do idime = 1,ndime
              elcod(idime,1:nnode) = coord(connec(ielem,1:nnode),idime)
           end do
           do igaus = 1,ngaus
              dN = dNgp(:,:,igaus)
              Je = matmul(elcod,transpose(dN))
              detJe = Je(1,1)*Je(2,2)-Je(2,1)*Je(1,2)
              He(1,1) = Je(2,2)
              He(2,2) = Je(1,1)
              He(1,2) = -Je(1,2)
              He(2,1) = -Je(2,1)
              He = (1.0d0/detJe)*He
              struc_J(:,:,igaus,ielem) = Je
              struc_detJ(1,igaus,ielem) = detJe
              struc_H(:,:,igaus,ielem) = He
              dxN(1,:) = He(1,1)*dN(1,:)+He(1,2)*dN(2,:)
              dxN(2,:) = He(2,1)*dN(1,:)+He(2,2)*dN(2,:)
              gpcar(:,:,igaus,ielem) = dxN(:,:)
              gpvol(1,igaus,ielem) = wgp(igaus)*detJe
           end do
        end do

        !*********************************************************************!
        ! Call convection subroutine                                          !
        !*********************************************************************!

        allocate(Rconvec(npoin))
        call generic_scalar_convec(nelem,ngaus,npoin,nnode,ndime,connec,Ngp,gpcar,gpvol,q,Rconvec)

        do ipoin = 1,npoin
           write(*,*) ipoin, Rconvec(ipoin)
        end do

        if (abs(sum(Rconvec)) .gt. 0.00000000000001d0) then
                write(*,*) "|sum(Rconvec)| = ", abs(sum(Rconvec))
                stop 1
        else
                write(*,*) "|sum(Rconvec)| = ", abs(sum(Rconvec))
                write(*,*) "TEST PASSED!"
        end if

end program unitt_generic_scalar_convec
