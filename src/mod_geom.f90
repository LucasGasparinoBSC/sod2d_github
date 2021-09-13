module mod_geom

      use elem_qua
      use elem_hex

      contains

         subroutine char_length(ielem,nelem,nnode,npoin,ndime,connec,coord,he)

                 implicit none

                 integer(4), intent(in)  :: ielem, nelem, nnode, npoin, ndime, connec(nelem,nnode)
                 real(8),    intent(in)  :: coord(npoin,ndime)
                 real(8),    intent(out) :: he
                 integer(4)              :: iedge, ncorner, nedge
                 real(8)                 :: dist(12,ndime), dist2, aux

                 !
                 ! Compute r = x2-x1 for all element edges
                 !
                 if (ndime == 2) then
                    if (nnode == 3 .or. nnode == 6) then ! TRI_XX
                       write(*,*) "ELEMENT TRI_XX NOT CODED!"
                    else if (nnode == 4 .or. nnode == 9) then ! QUA_XX
                       call quad_edges(ielem,nelem,nnode,npoin,ndime,connec,coord,ncorner,nedge,dist(1:4,1:ndime))
                    else
                       write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                    end if
                 else if (ndime == 3) then
                    if (nnode == 4 .or. nnode == 10) then ! TET_XX
                       write(*,*) "ELEMENT TET_XX NOT CODED!"
                    else if (nnode == 8 .or. nnode == 27 .or. nnode == 64) then ! HEX_XX
                       call hexa_edges(ielem,nelem,nnode,npoin,ndime,connec,coord,ncorner,nedge,dist(1:12,1:ndime))
                    else
                       write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                    end if
                 else
                    write(*,*) "BY SIGMAR NO!"
                 end if

                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 dist2 = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
                 do iedge = 1,nedge
                    aux = sqrt(dot_product(dist(iedge,:),dist(iedge,:)))
                    dist2 = min(dist2,aux)
                 end do
                 he = dist2

         end subroutine char_length

end module mod_geom
