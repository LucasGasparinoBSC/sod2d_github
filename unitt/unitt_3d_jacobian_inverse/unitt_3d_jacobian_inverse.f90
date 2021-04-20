program unitt_3d_jacobian_inverse

        use jacobian_oper

        implicit none

        integer(4)           :: ndime, nnode
        real(8)              :: detJe, detHe
        real(8), allocatable :: elcod(:,:), dN(:,:), Je(:,:), He(:,:)

        ndime = 3
        nnode = 4

        allocate(elcod(ndime,nnode))
        allocate(dN(ndime,nnode))
        allocate(Je(ndime,ndime))
        allocate(He(ndime,ndime))

        elcod(1,1:4) = (/0.0d0, 2.0d0, 4.0d0, 2.0d0/)
        elcod(2,1:4) = (/0.0d0, 0.0d0, 0.0d0, 2.0d0/)
        elcod(3,1:4) = (/0.0d0, 2.0d0, 0.0d0, 0.0d0/)

        dN(1,1:4) = (/-1.0d0, 1.0d0, 0.0d0, 0.0d0/)
        dN(2,1:4) = (/-1.0d0, 0.0d0, 1.0d0, 0.0d0/)
        dN(3,1:4) = (/-1.0d0, 0.0d0, 0.0d0, 1.0d0/)
        call elem_jacobian(ndime,nnode,elcod,dN,Je,detJe,He)

        detHe = He(1,1)*He(2,2)*He(3,3)+He(1,2)*He(2,3)*He(3,1)+He(1,3)*He(2,1)*He(3,2) - &
                He(3,1)*He(2,2)*He(1,3)-He(3,2)*He(2,3)*He(1,1)-He(3,3)*He(2,1)*He(1,2)

        if (detHe .gt. 0.0)

end program unitt_3d_jacobian_inverse
