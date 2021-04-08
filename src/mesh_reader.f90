module mesh_reader

      contains

              subroutine read_dims(file_path,file_name,npoin,nelem,nboun)

                      implicit none

                      character(500), intent(in)  :: file_path, file_name
                      integer(4)    , intent(out) :: npoin, nelem, nboun
                      character(500)              :: file_type, line

                      write(file_type,*) ".dims.dat"

                      write(*,*) "--| READING DIMS FILE..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      read(99,*) line, npoin
                      write(*,*) "--| NODES ON MESH : ",npoin
                      read(99,*) line, nelem
                      write(*,*) "--| ELEMENTS ON MESH : ",nelem
                      read(99,*) line, nboun
                      write(*,*) "--| BOUNDARY ELEMENTS ON MESH : ",nboun
                      write(*,*) "--| END OF DIMS FILE!"
                      close(99)

              end subroutine read_dims

              subroutine read_geo_dat(file_path,file_name,npoin,nelem,nboun,nnode,ndime,npbou,connec,bound,coord)

                      implicit none

                      character(500), intent(in)  :: file_path, file_name
                      integer(4)    , intent(in)  :: npoin, nelem, nboun, nnode, ndime, npbou
                      integer(4)    , intent(out) :: connec(nelem,nnode), bound(nboun,npbou)
                      real(8)       , intent(out) :: coord(npoin,ndime)
                      integer(4)                  :: iline, int1, inode, idime, aux(nnode+1), bou_aux(npbou+1)
                      character(500)              :: file_type, line

                      write(file_type,*) ".geo.dat"

                      write(*,*) "--| READING GEO.DAT FILE..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      !
                      ! Nodes/element section, blank for now
                      !
                      read(99,*) ! Section header
                      do iline = 1,nelem
                         read(99,*)
                      end do
                      read(99,*) ! Section ender
                      !
                      ! Connectivity table section
                      !
                      write(*,*) "--| READING ELEMENT TABLE..."
                      read(99,*) ! Section header
                      do iline = 1,nelem
                         read(99,'(a)') line
                         read(line,*) (aux(inode), inode=1,nnode+1)
                         connec(iline,1:nnode) = aux(2:nnode+1)
                      end do
                      read(99,*) ! Section ender
                      !
                      ! Nodal coordinates section
                      !
                      write(*,*) "--| READING COORDINATES..."
                      read(99,*) ! Section header
                      do iline = 1,npoin
                         if (ndime == 2) then
                            read(99,*) int1, coord(iline,1), coord(iline,2)
                         else if (ndime == 3) then
                            read(99,*) int1, coord(iline,1), coord(iline,2), coord(iline,3)
                         end if
                      end do
                      read(99,*) ! Section ender
                      !
                      ! Boundary nodes section
                      !
                      write(*,*) "--| READING BOUNDARIES..."
                      read(99,*) line! Section header
                      do iline = 1,nboun
                         read(99,'(a)') line
                         read(line,*) (bou_aux(inode), inode=1,npbou+1)
                         bound(iline,1:npbou) = bou_aux(2:npbou+1)
                      end do
                      close(99)
                      write(*,*) "--| END OF GEO.DAT FILE!"

              end subroutine read_geo_dat

end module mesh_reader
