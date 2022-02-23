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
            write(ivtk,*) u3d(i,1:3)
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

      subroutine write_vtk_binary(isPeriodic,istep,ndime,npoin,nelem,nnode,coord,connec, &
                                 rho,u,pr,E,mu_e,nper,masSla)
         implicit none
      
         integer(4), intent(in)                            :: isPeriodic, nper
         integer(4), intent(in)                            :: istep, ndime, npoin, nelem, nnode
         integer(4), intent(in)                            :: connec(nelem,nnode)
         integer(4), intent(in), optional                  :: masSla(nper,2)
         real(8)   , intent(in)                            :: coord(npoin,ndime)
         real(8)   , intent(inout), dimension(npoin)       :: rho, pr, E
         real(8)   , intent(inout), dimension(npoin,ndime) :: u
         real(8)   , intent(inout), dimension(nelem)       :: mu_e
         integer(4)                                        :: i, iper, ivtk=9
         integer(4)            , dimension(nelem,nnode+1)  :: cells
         integer(4)            , dimension(nelem)          :: cellTypes
         real(8)               , dimension(npoin,3)        :: points, u3d
         character(500)                                    :: filename
         character(80)                                     :: buffer
         character(8)                                      :: str1, str2
         character(1)                                      :: lf

         lf = achar(10)
      
         !
         ! Pass coordinates to a suitable 3D generic format
         !
         !$acc kernels
         points(:,:) = 0.0d0
         points(:,1:ndime) = coord(:,1:ndime)
         !$acc end kernels
         print*, 'kernel1: ok!'
      
         !
         ! If case is periodic, adjust slave nodes
         !
         if (isPeriodic .eq.1 .and. present(masSla)) then
            !!$acc parallel loop
            do iper = 1,nper
               u(masSla(iper,2),1:ndime) = u(masSla(iper,1),1:ndime)
               rho(masSla(iper,2)) = rho(masSla(iper,1))
               pr(masSla(iper,2)) = pr(masSla(iper,1))
               E(masSla(iper,2)) = E(masSla(iper,1))
            end do
            !!$acc end parallel loop
         end if
         print*, 'kernel2: ok!'

         !
         ! Pass vector data to a suitable 3D generic format
         !
         !$acc kernels
         u3d(:,:) = 0.0d0
         u3d(:,1:ndime) = u(:,1:ndime)
         !$acc end kernels
      
         !
         ! Pass cell list to VTK format
         !
         !$acc kernels
         cells(:,1) = nnode
         cells(:,2:nnode+1) = connec(:,1:nnode)
         cells(:,2:nnode+1) = cells(:,2:nnode+1)-1 ! maybe can be removed?
         !$acc end kernels
      
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
         open(unit=ivtk,file=filename,status='replace',access='stream',convert='BIG_ENDIAN') ! Binary file access with stream
         
         !
         ! Write header in ascii format
         !
         write(ivtk) '# vtk DataFile Version 3.0'//lf
         write(ivtk) 'unstr_grid'//lf
         write(ivtk) 'BINARY'//lf
         write(ivtk) 'DATASET UNSTRUCTURED_GRID'//lf//lf
         
         !
         ! Write points
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) 'POINTS '//str1//'  double'//lf
         do i = 1,npoin
            write(ivtk) points(i,:)
         end do
         
         !
         ! Write cells
         !
         write(str1(1:8),'(i8)') nelem
         write(str2(1:8),'(i8)') nelem*(nnode+1)
         write(ivtk) lf//lf//'CELLS '//str1//' '//str2//lf
         do i = 1,nelem
            write(ivtk) cells(i,:)
         end do
         
         !
         ! Write cell types
         !
         write(str1(1:8),'(i8)') nelem
         write(ivtk) lf//lf//'CELL_TYPES '//str1//lf
         do i = 1,nelem
            write(ivtk) cellTypes(i)
         end do
         
         !
         ! Write point scalar data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'POINT_DATA '//str1//lf
         write(ivtk) 'SCALARS DENSI double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) rho(i)
         end do
         write(ivtk) lf//lf//'SCALARS PRESS double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) pr(i)
         end do
         write(ivtk) lf//lf//'SCALARS TENER double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) E(i)
         end do
         
         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS VELOC double'//lf
         do i = 1,npoin
            write(ivtk) u3d(i,:)
         end do
         
         !
         ! Write cell scalar data
         !
         write(str1(1:8),'(i8)') nelem
         write(ivtk) lf//lf//'CELL_DATA '//str1//lf
         write(ivtk) 'SCALARS ENVIT double'//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,nelem
            write(ivtk) mu_e(i)
         end do
         
         close(ivtk)
      
      end subroutine write_vtk_binary

end module mod_output
