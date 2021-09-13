module inicond_reader

      contains

              subroutine read_veloc(ndime,npoin,file_path,u)

                      implicit none

                      character(500), intent(in)  :: file_path
                      integer(4)    , intent(in)  :: ndime, npoin
                      real(8)       , intent(out) :: u(npoin,ndime)
                      character(500)              :: file_type, file_name
                      integer(4)                  :: ipoin, ind

                      write(file_type,*) ".alya"
                      write(file_name,*) "VELOC"
                      write(*,*) "--| READING FILE VELOC.alya..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      do ipoin = 1,npoin
                        if (ndime == 2) then
                           read(99,*) ind, u(ipoin,1), u(ipoin,2)
                        else if (ndime == 3) then
                           read(99,*) ind, u(ipoin,1), u(ipoin,2), u(ipoin,3)
                        end if
                      end do
                      close(99)

              end subroutine read_veloc

              subroutine read_densi(npoin,file_path,rho)

                      implicit none

                      character(500), intent(in)  :: file_path
                      integer(4)    , intent(in)  :: npoin
                      real(8)       , intent(out) :: rho(npoin)
                      character(500)              :: file_type, file_name
                      integer(4)                  :: ipoin, ind

                      write(file_type,*) ".alya"
                      write(file_name,*) "DENSI"
                      write(*,*) "--| READING FILE DENSI.alya..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      do ipoin = 1,npoin
                         read(99,*) ind, rho(ipoin)
                      end do
                      close(99)

              end subroutine read_densi

              subroutine read_press(npoin,file_path,pr)

                      implicit none

                      character(500), intent(in)  :: file_path
                      integer(4)    , intent(in)  :: npoin
                      real(8)       , intent(out) :: pr(npoin)
                      character(500)              :: file_type, file_name
                      integer(4)                  :: ipoin, ind

                      write(file_type,*) ".alya"
                      write(file_name,*) "PRESS"
                      write(*,*) "--| READING FILE PRESS.alya..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      do ipoin = 1,npoin
                         read(99,*) ind, pr(ipoin)
                      end do
                      close(99)

              end subroutine read_press

              subroutine read_temper(npoin,file_path,Tem)

                      implicit none

                      character(500), intent(in)  :: file_path
                      integer(4)    , intent(in)  :: npoin
                      real(8)       , intent(out) :: Tem(npoin)
                      character(500)              :: file_type, file_name
                      integer(4)                  :: ipoin, ind

                      write(file_type,*) ".alya"
                      write(file_name,*) "TEMPE"
                      write(*,*) "--| READING FILE TEMPE.alya..."
                      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
                      do ipoin = 1,npoin
                         read(99,*) ind, Tem(ipoin)
                      end do
                      close(99)

              end subroutine read_temper

end module inicond_reader
