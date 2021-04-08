program sod2d

        !*********************************************************************!
        ! Computes the sod shock tube problem over a 2D domain using FEM.     !
        ! Stabilized through Entropy viscosity method.                        !
        !*********************************************************************!

        use elem_qua
        use jacobian_oper
        use quadrature_rules
        use mesh_reader
        use inicond_reader
        use mass_matrix

        use time_integ

        implicit none

        integer(4)                 :: ndime, nnode, ngaus, nstep
        integer(4)                 :: idime, inode, igaus, istep
        integer(4)                 :: nelem, npoin, nboun, npbou
        integer(4)                 :: ielem, ipoin, iboun, ipbou
        integer(4)                 :: idof, ndof, nbnodes, ibnodes
        integer(4)                 :: flag_predic
        integer(4), allocatable    :: connec(:,:), bound(:,:), ldof(:), lbnodes(:)
        integer(4), allocatable    :: aux1(:)
        real(8),    allocatable    :: coord(:,:), elcod(:,:)
        real(8),    allocatable    :: xgp(:,:), wgp(:)
        real(8),    allocatable    :: N(:), dN(:,:), Ngp(:,:), dNgp(:,:,:)
        real(8),    allocatable    :: Je(:,:), He(:,:)
        real(8),    allocatable    :: struc_J(:,:,:,:), struc_H(:,:,:,:), struc_detJ(:,:,:)
        real(8),    allocatable    :: dxN(:,:), gpcar(:,:,:,:), gpvol(:,:,:)
        real(8),    allocatable    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:)
        real(8),    allocatable    :: Mc(:,:), Ml(:)
        real(8)                    :: s, t, detJe
        real(8)                    :: Rgas, gamma_gas, Cp, Cv
        real(8)                    :: dt, cfl
        character(500)             :: file_path
        character(500)             :: file_name, dumpfile
        character(5)               :: matrix_type, solver_type

        !*********************************************************************!
        ! Basic data, hardcoded for now                                       !
        !*********************************************************************!

        write(*,*) "--| ENTER PROBLEM DIMENSION (2 OR 3) :"
        read(*,*) ndime
        nnode = 4 ! TODO: need to allow for mixed elements...
        npbou = 2 ! TODO: Need to get his from somewhere...
        nstep = 100 ! TODO: Needs to be input...
        Rgas = 287.00d0
        Cp = 1004.00d0
        gamma_gas = 1.40d0
        Cv = Cp/gamma_gas
        dt = 0.0001d0 ! TODO: make it adaptive...

        !*********************************************************************!
        ! Read mesh in Alya format                                            !
        !*********************************************************************!

        write(file_path,*) "./mesh/"
        write(*,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(file_path))," !"
        write(*,*) "--| ENTER NAME OF MESH RELATED FILES :"
        read(*,*) file_name
        call read_dims(file_path,file_name,npoin,nelem,nboun)
        allocate(connec(nelem,nnode))
        allocate(bound(nboun,npbou))
        allocate(coord(npoin,ndime))
        call read_geo_dat(file_path,file_name,npoin,nelem,nboun,nnode,ndime,npbou,connec,bound,coord)

        !*********************************************************************!
        ! Generate list of "free" nodes                                       !
        !*********************************************************************!

        allocate(aux1(npoin))

        do ipoin = 1,npoin
           aux1(ipoin) = ipoin
        end do

        ndof = 0
        do ipoin = 1,npoin
           do iboun = 1,nboun
              do ipbou = 1,npbou
                 if (bound(iboun,ipbou) == ipoin) then
                    aux1(ipoin) = 0
                    ndof = ndof+1
                    exit
                 end if
              end do
           end do
        end do
        nbnodes = ndof
        ndof = npoin-ndof
        write(*,*) '--| TOTAL FREE NODES := ',ndof

        allocate(ldof(ndof))
        allocate(lbnodes(nbnodes))

        idof = 0
        ibnodes = 0
        do ipoin = 1,npoin
           if (aux1(ipoin) .ne. 0) then
              idof = idof+1
              ldof(idof) = aux1(ipoin)
           else
              ibnodes = ibnodes+1
              lbnodes(ibnodes) = ipoin
           end if
        end do

        !*********************************************************************!
        ! Allocate variables                                                  !
        !*********************************************************************!

        !
        ! Last rank is for prediction-advance related to entropy viscosity,
        ! where 1 is prediction, 2 is final value
        !
        allocate(u(npoin,ndime,2))  ! Velocity
        allocate(q(npoin,ndime,2))  ! momentum
        allocate(rho(npoin,2))      ! Density
        allocate(pr(npoin,2))       ! Pressure
        allocate(E(npoin,2))        ! Total Energy
        allocate(Tem(npoin,2))      ! Temperature
        allocate(e_int(npoin,2))    ! Internal Energy

        !*********************************************************************!
        ! Read initial conditions                                             !
        !*********************************************************************!

        call read_veloc(ndime,npoin,file_path,u(:,:,2))
        call read_densi(npoin,file_path,rho(:,2))
        call read_press(npoin,file_path,pr(:,2)) ! Can be switched for TEMPE

        !
        ! File dump
        !
        write(dumpfile,'("tstep_",i0,".dat")') 0
        open(unit = 99+1,file = dumpfile,form="formatted",status="replace",action="write")
        do ipoin = 1,npoin
           if (coord(ipoin,2) > -0.045d0 .and. coord(ipoin,2) < 0.045d0) then
              write(99+1,"(f8.4, f16.8, f16.8, f16.8)") coord(ipoin,1), rho(ipoin,2), u(ipoin,1,2), pr(ipoin,2)
           end if
        end do
        close(unit=99+1)

        !*********************************************************************!
        ! Generate complementary info                                         !
        !*********************************************************************!

        ! Assuming u, rho, p as IC:

        do ipoin = 1,npoin
           e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(gamma_gas-1.0d0))
           Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*Rgas)
           E(ipoin,2) = rho(ipoin,2)*(0.5d0*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
           q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
        end do

        !*********************************************************************!
        ! Generate GLL table                                                  !
        !*********************************************************************!

        write(*,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

        if (nnode == 4) then
                ngaus = 4
        else if (nnode == 9) then
                ngaus = 9
        end if
        !ngaus = 1 ! Test value

        allocate(xgp(ngaus,ndime))
        allocate(wgp(ngaus))

        call gll_qua(ndime,ngaus,xgp,wgp)

        !*********************************************************************!
        ! Generate N and dN for all GP                                        !
        !*********************************************************************!

        write(*,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."

        allocate(N(nnode),dN(ndime,nnode))
        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))

        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           call qua04(s,t,N,dN)
           Ngp(igaus,:) = N
           dNgp(:,:,igaus) = dN
        end do

        !*********************************************************************!
        ! Generate Jacobian related information                               !
        !*********************************************************************!

        write(*,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

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
           if (ndime == 2) then
              elcod(1,1:nnode) = coord(connec(ielem,1:nnode),1)
              elcod(2,1:nnode) = coord(connec(ielem,1:nnode),2)
           else if (ndime == 3) then
              elcod(1,1:nnode) = coord(connec(ielem,1:nnode),1)
              elcod(2,1:nnode) = coord(connec(ielem,1:nnode),2)
              elcod(3,1:nnode) = coord(connec(ielem,1:nnode),3)
           end if
           do igaus = 1,ngaus
              dN = dNgp(:,:,igaus)
              call elem_jacobian(ndime,nnode,elcod,dN,Je,detJe,He)
              struc_J(:,:,igaus,ielem) = Je
              struc_detJ(1,igaus,ielem) = detJe
              struc_H(:,:,igaus,ielem) = He
              call cartesian_deriv(ndime,nnode,dN,He,dxN)
              gpcar(:,:,igaus,ielem) = dxN(:,:)
              gpvol(1,igaus,ielem) = wgp(igaus)*detJe
           end do
        end do

        !*********************************************************************!
        ! Compute mass matrix (Lumped and Consistent) and set solver type      !
        !*********************************************************************!

        write(*,*) '--| COMPUTING LUMPED MASS MATRIX...'
        allocate(Ml(npoin))
        call lumped_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Ml)
        solver_type = 'LUMSO'

        write(*,*) '--| COMPUTING CONSISTENT MASS MATRIX...'
        allocate(Mc(npoin,npoin))
        call consistent_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Mc)
        write(*,*) '--| ENTER REQUIRED SOLVER FOR MASS MATRIX:'
        write(*,*) '--| AVAILABLE SOLVERS ARE: LUMSO, APINV, CONGR:'
        read(*,*) solver_type
        write(*,*) '--| USING SOLVER ',solver_type,' FOR MASS MATRIX'

        !*********************************************************************!
        ! Start of time stepping                                              !
        !*********************************************************************!

        do istep = 1,nstep

           write(*,*) '   --| STEP: ', istep

           !
           ! Prediction
           !
           flag_predic = 1
           rho(:,1) = rho(:,2)
           u(:,:,1) = u(:,:,2)
           q(:,:,1) = q(:,:,2)
           pr(:,1) = pr(:,2)
           E(:,1) = E(:,2)
           Tem(:,1) = Tem(:,2)
           e_int(:,1) = e_int(:,2)

           call rk_4_main(flag_predic,nelem,npoin,ndime,ndof,nbnodes,ngaus,nnode, &
                     ldof,lbnodes,connec,Ngp,gpcar,Ml,Mc,gpvol,dt, &
                     rho,u,q,pr,E,Tem,e_int)

           !
           ! Advance with entropy viscosity
           !
           flag_predic = 0
           call rk_4_main(flag_predic,nelem,npoin,ndime,ndof,nbnodes,ngaus,nnode, &
                     ldof,lbnodes,connec,Ngp,gpcar,Ml,Mc,gpvol,dt, &
                     rho,u,q,pr,E,Tem,e_int)

        end do

        !
        ! File dump
        !
        write(dumpfile,'("tstep_",i0,".dat")') 20
        open(unit = 99+1,file = dumpfile,form="formatted",status="replace",action="write")
        do ipoin = 1,npoin
           if (coord(ipoin,2) > -0.045d0 .and. coord(ipoin,2) < 0.045d0) then
              write(99+1,"(f8.4, f16.8, f16.8, f16.8)") coord(ipoin,1), rho(ipoin,2), u(ipoin,1,2), pr(ipoin,2)
           end if
        end do
        close(unit=99+1)

end program sod2d
