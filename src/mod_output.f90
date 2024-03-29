module mod_output
        contains
                subroutine write_vtk_ascii(istep,ndime,npoin,nelem,nnode,coord,connec, &
                                           rho,u,pr,E,mu_e)
                   implicit none

                   integer(4), intent(in)                           :: istep, ndime, npoin, nelem, nnode
                   integer(4), intent(in)                           :: connec(nelem,nnode)
                   real(8)   , intent(in)                           :: coord(npoin,ndime)
                   real(8)   , intent(in), dimension(npoin)         :: rho, pr, E
                   real(8)   , intent(in), dimension(npoin,ndime)   :: u
                   real(8)   , intent(in), dimension(nelem)         :: mu_e
                   integer(4)                                       :: i, ivtk=9
                   integer(4)            , dimension(nelem,nnode+1) :: cells
                   integer(4)            , dimension(nelem)         :: cellTypes
                   real(8)               , dimension(npoin,3)       :: points, u3d
                   character(500)                                   :: filename
                   character(8)                                     :: str1, str2

                   !
                   ! Pass coordinates to a suitable 3D generic format
                   !
                   points = 0.0d0
                   points(:,1:ndime) = coord(:,1:ndime)

                   !
                   ! Pass vector data to a suitable 3D generic format
                   !
                   u3d = 0.0d0
                   u3d(:,1:ndime) = u(:,1:ndime)

                   !
                   ! Pass cell list to VTK format
                   !
                   cells(:,1) = nnode
                   cells(:,2:nnode+1) = connec(:,1:nnode)
                   cells(:,2:nnode+1) = cells(:,2:nnode+1)-1

                   !
                   ! Define cell types
                   !
                   if (ndime .eq. 2) then
                      if (nnode .eq. 4) then ! QUA04
                         cellTypes = 9
                      end if
                   else if (ndime .eq. 3) then
                      if (nnode .eq. 8) then ! HEX08
                         cellTypes = 12
                      else if (nnode .eq. 27) then ! HEX27
                         cellTypes = 29
                      end if
                   end if

                   !
                   ! Open file with ascii input
                   !
                   write(filename,'("vtkTstep_",i0,".vtk")') istep
                   open(unit=ivtk,file=filename,status='replace') ! Binary file access with stream
                   
                   !
                   ! Write header in ascii format
                   !
                   write(ivtk,'(a)') '# vtk DataFile Version 3.0'
                   write(ivtk,'(a)') 'unstr_grid'
                   write(ivtk,'(a)') 'ASCII' 
                   write(ivtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
                   
                   !
                   ! Write points
                   !
                   write(str1(1:8),'(i8)') npoin
                   write(ivtk,'(a)') 'POINTS '//trim(str1)//'  double'
                   do i = 1,npoin
                      write(ivtk,*) points(i,:)
                   end do
                   
                   !
                   ! Write cells
                   !
                   write(str1(1:8),'(i8)') nelem
                   write(str2(1:8),'(i8)') nelem*(nnode+1)
                   write(ivtk,'(a)') 'CELLS '//trim(str1)//trim(str2)
                   do i = 1,nelem
                      write(ivtk,*) cells(i,:)
                   end do
                   
                   
                   ! Write cell types
                   !
                   write(str1(1:8),'(i8)') nelem
                   write(ivtk,'(a)') 'CELL_TYPES '//trim(str1)
                   do i = 1,nelem
                      write(ivtk,*) cellTypes(i)
                   end do
                   
                   !
                   ! Write point scalar data
                   !
                   write(str1(1:8),'(i8)') npoin
                   write(ivtk,'(a)') 'POINT_DATA '//trim(str1)
                   write(ivtk,'(a)') 'SCALARS '//' DENSI '//' double '//' 1'
                   write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
                   do i = 1,npoin
                      write(ivtk,*) rho(i)
                   end do
                   write(ivtk,'(a)') 'SCALARS '//' PRESS '//' double '//' 1'
                   write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
                   do i = 1,npoin
                      write(ivtk,*) pr(i)
                   end do
                   write(ivtk,'(a)') 'SCALARS '//' TENER '//' double '//' 1'
                   write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
                   do i = 1,npoin
                      write(ivtk,*) E(i)
                   end do
                   
                   !
                   ! Write point vector data
                   !
                   write(str1(1:8),'(i8)') npoin
                   write(ivtk,'(a)') 'VECTORS '//' VELOC '//' double'
                   do i = 1,npoin
                      write(ivtk,*) u3d(i,:)
                   end do
                   
                   !!
                   !! Write cell scalar data
                   !!
                   write(str1(1:8),'(i8)') nelem
                   write(ivtk,'(a)') 'CELL_DATA '//trim(str1)
                   write(ivtk,'(a)') 'SCALARS '//' ENVIT '//' double'
                   write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
                   do i = 1,nelem
                      write(ivtk,*) mu_e(i)
                   end do
                   
                   close(ivtk)

                end subroutine

end module mod_output
