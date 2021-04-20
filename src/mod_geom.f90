module mod_geom

      ! TODO: Finish it ffs...

      use elem_qua

      contains

         subroutine char_length(ielem,nelem,nnode,npoin,ndime,connec,coord,he)

                 implicit none

                 integer(4), intent(in)  :: ielem, nelem, nnode, npoin, ndime, connec(nelem,nnode)
                 real(8),    intent(in)  :: coord(npoin,ndime)
                 real(8),    intent(out) :: he
                 integer(4)              :: iedge, ncorner, nedge
                 real(8)                 :: dist(12,ndime), dist2(12)

                 dist = 0.0d0
                 if (ndime == 2) then
                         if (nnode == 3 .or. nnode == 6) then ! TRI_XX
                                 write(*,*) "ELEMENT TRI_XX NOT CODED!"
                         else if (nnode == 4 .or. nnode == 9) then ! QUA_XX
                                 call quad_edges(ielem,nelem,nnode,npoin,ndime,connec,coord,ncorner,nedge,dist)
                         else
                                 write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                         end if
                 else if (ndime == 3) then
                         write(*,*) "NOT CODED YET!"
                 else
                         write(*,*) "BY SIGMAR NO!"
                 end if

                 do iedge = 1,nedge
                    dist2(iedge) = sqrt(dot_product(dist(iedge,:),dist(iedge,:)))
                 end do
                 he = minval(dist2)

         end subroutine char_length

end module mod_geom
