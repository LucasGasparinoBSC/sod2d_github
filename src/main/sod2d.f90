program sod2d

        !*********************************************************************!
        ! Computes the sod shock tube problem over a 2D domain using FEM.     !
        ! Stabilized through Entropy viscosity method.                        !
        !*********************************************************************!

        use mod_nvtx
#ifdef GPU
        use cudafor
        use mod_gpu_vars
#endif

        use omp_lib
        use elem_qua
        use elem_hex
        use jacobian_oper
        use quadrature_rules
        use mesh_reader
        use inicond_reader
        use mass_matrix
        use mod_graph
        use mod_geom
        use mod_output

        use time_integ

        implicit none

        integer(4)                 :: ndime, nnode, ngaus, nstep, nzdom
        integer(4)                 :: idime, inode, igaus, istep, izdom
        integer(4)                 :: nelem, npoin, nboun, npbou, nbcodes
        integer(4)                 :: ielem, ipoin, iboun, ipbou
        integer(4)                 :: idof, ndof, nbnodes, ibnodes
        integer(4)                 :: ppow, porder
        integer(4)                 :: flag_predic
        integer(4), allocatable    :: rdom(:), cdom(:), aux_cdom(:)
        integer(4), allocatable    :: connec(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
        integer(4), allocatable    :: aux1(:)
        real(8),    allocatable    :: coord(:,:), elcod(:,:), helem(:)
        real(8),    allocatable    :: xgp(:,:), wgp(:)
        real(8),    allocatable    :: N(:), dN(:,:), Ngp(:,:), dNgp(:,:,:)
        real(8),    allocatable    :: Je(:,:), He(:,:)
        real(8),    allocatable    :: struc_J(:,:,:,:), struc_H(:,:,:,:), struc_detJ(:,:,:)
        real(8),    allocatable    :: dxN(:,:), gpcar(:,:,:,:), gpvol(:,:,:)
        real(8),    allocatable    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:)
        real(8),    allocatable    :: Mc(:), Ml(:)
        real(8),    allocatable    :: mu_e(:)
        real(8)                    :: s, t, z, detJe
        real(8)                    :: Rgas, gamma_gas, Cp, Cv
        real(8)                    :: dt, cfl, he_aux
        character(500)             :: file_path
        character(500)             :: file_name, dumpfile
        character(5)               :: matrix_type, solver_type
        character(4)               :: timeStep

        integer(4) :: counter

        !*********************************************************************!
        ! Basic data, hardcoded for now                                       !
        !*********************************************************************!

        write(*,*) "--| ENTER PROBLEM DIMENSION (2 OR 3) :"
        !read(*,*) ndime
        ndime = 3 ! NVVP
        nnode = 27 ! TODO: need to allow for mixed elements...
        porder = 2 ! Element order
        npbou = 9 ! TODO: Need to get his from somewhere...
        nstep = 100 ! TODO: Needs to be input...
        Rgas = 287.00d0
        !Rgas = 1.00d0
        Cp = 1004.00d0
        gamma_gas = 1.40d0
        Cv = Cp/gamma_gas
        dt = 0.0025d0/2.0d0 ! TODO: make it adaptive...

        !*********************************************************************!
        ! Read mesh in Alya format                                            !
        !*********************************************************************!

        write(file_path,*) "./mesh/"
        write(*,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(file_path))," !"
        write(*,*) "--| ENTER NAME OF MESH RELATED FILES :"
        !read(*,*) file_name
        write(file_name,*) "shock_tube" ! NVVP
        call read_dims(file_path,file_name,npoin,nelem,nboun)
        allocate(connec(nelem,nnode))
        allocate(bound(nboun,npbou))
        allocate(bou_codes(nboun,2))
        allocate(coord(npoin,ndime))
        call read_geo_dat(file_path,file_name,npoin,nelem,nboun,nnode,ndime,npbou,connec,bound,coord)
        call read_fixbou(file_path,file_name,nboun,nbcodes,bou_codes)

        !*********************************************************************!
        ! Compute characteristic size of elements                             !
        !*********************************************************************!

        allocate(helem(nelem))
        do ielem = 1,nelem
           call char_length(ielem,nelem,nnode,npoin,ndime,connec,coord,he_aux)
           helem(ielem) = he_aux
        end do

        !*********************************************************************!
        ! Create mesh graph for CSR matrices                                  !
        !*********************************************************************!

        write(*,*) "--| PERFORMING GRAPH OPERATIONS..."
        allocate(rdom(npoin+1))                                          ! Implicit row indexing
        allocate(aux_cdom(nelem*nnode*nnode))                            ! Preliminary cdom for subroutine
        call compute_nzdom(npoin,nnode,nelem,connec,nzdom,rdom,aux_cdom) ! Computes nzdom, rdom and aux_cdom
        allocate(cdom(nzdom))                                            ! Row indexes with proper array size
        do izdom = 1,nzdom
           cdom(izdom) = aux_cdom(izdom)
        end do
        deallocate(aux_cdom)
        write(*,*) "--| END OF GRAPH OPERATIONS!"

        !*********************************************************************!
        ! Generate list of "free" nodes                                       !
        !*********************************************************************!

        write(*,*) "--| SPLITTING BOUNDARY NODES FROM DOFs..."
        allocate(aux1(npoin))

        !
        ! Fill aux1 with all nodes in order
        !
        do ipoin = 1,npoin
           aux1(ipoin) = ipoin
        end do

        !
        ! Zero aux1 entries that belong to a boundary
        !
        do ipoin = 1,npoin       ! Loop over all nodes to fill aux()
           do iboun = 1,nboun    ! Loop over element edges belonging to a boundary
              do ipbou = 1,npbou ! Loop over nodes on an element face belonging to a boundary
                 if (bound(iboun,ipbou) == ipoin) then
                    aux1(ipoin) = 0
                    exit
                 end if
              end do
           end do
        end do

        !
        ! Determine how many nodes are boundary nodes
        !
        ndof = 0
        do ipoin = 1,npoin
           if (aux1(ipoin) == 0) then
              ndof = ndof+1
           end if
        end do

        nbnodes = ndof
        ndof = npoin-ndof
        write(*,*) '--| TOTAL FREE NODES := ',ndof
        write(*,*) '--| TOTAL BOUNDARY NODES := ',nbnodes

        allocate(ldof(ndof))
        allocate(lbnodes(nbnodes))

        !
        ! Split aux1 into the 2 lists
        !
        idof = 0    ! Counter for free nodes
        ibnodes = 0 ! Counter for boundary nodes
        do ipoin = 1,npoin
           if (aux1(ipoin) == 0) then
              ibnodes = ibnodes+1
              lbnodes(ibnodes) = ipoin
           else
              idof = idof+1
              ldof(idof) = aux1(ipoin)
           end if
        end do

        !*********************************************************************!
        ! Allocate variables                                                  !
        !*********************************************************************!

        WRITE(*,*) "--| ALLOCATING MAIN VARIABLES"
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
        allocate(mu_e(nelem))       ! Elemental viscosity

        !*********************************************************************!
        ! Read initial conditions                                             !
        !*********************************************************************!

        call read_veloc(ndime,npoin,file_path,u(:,:,2))
        call read_densi(npoin,file_path,rho(:,2))
        call read_press(npoin,file_path,pr(:,2)) ! Can be switched for TEMPE

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
        mu_e = 0.0d0

        !
        ! Call VTK output
        !
        call write_vtk_ascii(0,ndime,npoin,nelem,nnode,coord,connec, &
                             rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e)

        !*********************************************************************!
        ! Generate GLL table                                                  !
        !*********************************************************************!

        write(*,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

        ! TODO: allow for more element types...

        if (ndime == 2) then ! 2D elements
           if (nnode == 3) then ! TRI03
               ngaus = 3
               write(*,*) 'NOT CODED YET!'
               STOP 1
           else if (nnode == 6) then ! TRI06
               ngaus = 7
               write(*,*) 'NOT CODED YET!'
               STOP 1
           else if (nnode == 4) then ! QUA04
               ngaus = 4
           else if (nnode == 9) then ! QUA09
               ngaus = 9
           end if
        else if (ndime == 3) then ! 3D elements
           if (nnode == 4) then ! TET04
              ngaus = 4
              write(*,*) 'NOT CODED YET!'
              stop 1
           else if (nnode == 8) then ! HEX08
              ngaus = 8
           else if (nnode == 27) then ! HEX27
              ngaus = 27
              write(*,*) '--| USING 27 GAUSS NODES PER ELEMENT!'
           else if (nnode == 64) then ! HEX64
              ngaus = 64
           else
              write(*,*) 'ELEMENT DOES NOT EXIST, OR IS NOT CODED YET!!'
              stop 1
           end if
        else
           write(*,*) 'ONLY 2D AND 3D ELEMMENTS SUPPORTED!'
           stop 1
        end if

        ! Option to run with 1 Gauss point per element, for testing and debugging. Will override the previous section.

#ifdef TGAUS
        ngaus = 1 ! Test value
#endif

        allocate(xgp(ngaus,ndime))
        allocate(wgp(ngaus))

        if (ndime == 2) then
           if (nnode == (porder+1)**2) then ! QUA_XX of order porder
              call gll_qua(ndime,ngaus,xgp,wgp)
           else if (nnode == 3 .or. nnode == 6 .or. nnode == 10) then ! TRI_XX
              write(*,*) '--| NOT CODED YET!'
              STOP 1
           end if
        else if (ndime == 3) then
           if (nnode == (porder+1)**3) then ! HEX_XX
              call gll_hex(ndime,ngaus,xgp,wgp)
           else if (nnode == 4 .or. nnode == 10 .or. nnode == 20) then ! TET_XX
              write(*,*) '--| NOT CODED YET!'
              STOP 1
           end if
        end if

        !*********************************************************************!
        ! Generate N and dN for all GP                                        !
        !*********************************************************************!

        ! TODO: Allow for more element types

        write(*,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."

        allocate(N(nnode),dN(ndime,nnode))
        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))

        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           if (ndime == 2) then
              if (nnode == 4) then
                 call qua04(s,t,N,dN)
              else if (nnode == 9) then
                 !call qua09(s,t,N,dN)
                 write(*,*) '--| NOT CODED YET!'
                 STOP 1
              end if
           else if (ndime == 3) then
              z = xgp(igaus,3)
              if (nnode == 8) then
                 call hex08(s,t,z,N,dN)
              else if (nnode == 27) then
                 call hex27(s,t,z,N,dN)
              else if (nnode == 64) then
                 !call hex64(s,t,z,N,dN)
              else
              end if
           end if
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
        ! Compute mass matrix (Lumped and Consistent) and set solver type     !
        !*********************************************************************!

        write(*,*) '--| COMPUTING LUMPED MASS MATRIX...'
        allocate(Ml(npoin))
        call lumped_mass(nelem,nnode,npoin,ngaus,connec,gpvol,Ngp,Ml)
        solver_type = 'LUMSO'

        write(*,*) '--| COMPUTING CONSISTENT MASS MATRIX...'
        allocate(Mc(nzdom))
        call consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc)
        write(*,*) '--| ENTER REQUIRED SOLVER FOR MASS MATRIX:'
        write(*,*) '--| AVAILABLE SOLVERS ARE: LUMSO, APINV:'
        !read(*,*) solver_type
        write(solver_type,'(a)') "APINV" ! NVVP
        if (solver_type == 'APINV') then
                write(*,*) '--| ENTER NUMBER OF ITERATIONS FOR APINV SOLVER:'
                !read(*,*) ppow
                ppow = 4 ! NVVP
        end if
        write(*,*) '--| USING SOLVER ',solver_type,' FOR MASS MATRIX'

        !*********************************************************************!
        ! Create variables on GPU                                             !
        !*********************************************************************!

#ifdef GPU

        ! Range with standard color
        call nvtxStartRange("Memory Management")

        ! Mesh info

        allocate(connec_d(nelem,nnode))
        allocate(lbnodes_d(nbnodes))

        connec_d = connec
        lbnodes_d = lbnodes

        ! Primary vars.

        allocate(rho_d(npoin,2))      ! Density
        allocate(u_d(npoin,ndime,2))  ! Velocity
        allocate(q_d(npoin,ndime,2))  ! momentum
        allocate(pr_d(npoin,2))       ! Pressure
        allocate(E_d(npoin,2))        ! Total Energy
        allocate(Tem_d(npoin,2))      ! Temperature
        allocate(e_int_d(npoin,2))    ! Internal Energy
        allocate(mu_e_d(nelem))       ! Elemental viscosity

        rho_d = rho
        u_d = u
        q_d = q
        pr_d = pr
        E_d = E
        Tem_d = Tem
        e_int_d = e_int

        ! Mass matrices

        allocate(Ml_d(npoin))
        allocate(Mc_d(nzdom))

        Ml_d = Ml
        Mc_d = Mc

        ! Elemental info

        allocate(Ngp_d(ngaus,npoin))
        allocate(gpvol_d(1,ngaus,nelem))
        allocate(gpcar_d(ndime,nnode,ngaus,nelem))

        Ngp_d = Ngp
        gpvol_d = gpvol
        gpcar_d = gpcar_d

        ! End nvtx range
        call nvtxEndRange

#endif

        !*********************************************************************!
        ! Start of time stepping                                              !
        !*********************************************************************!

        counter = 1

        call nvtxStartRange("Start RK4")
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

           ! nvtx range for full RK
           write(timeStep,'(i4)') istep
           call nvtxStartRange("RK4 step "//timeStep,istep)

           call rk_4_main(flag_predic,nelem,nboun,npbou,npoin,ndime,ndof,nbnodes,ngaus,nnode, &
                     ppow, nzdom,rdom,cdom,ldof,lbnodes,connec,bound,bou_codes, &
                     Ngp,gpcar,Ml,Mc,gpvol,dt,helem,Rgas,gamma_gas, &
                     rho,u,q,pr,E,Tem,e_int,mu_e)

           !
           ! Advance with entropy viscosity
           !
           flag_predic = 0
           call rk_4_main(flag_predic,nelem,nboun,npbou,npoin,ndime,ndof,nbnodes,ngaus,nnode, &
                     ppow, nzdom,rdom,cdom,ldof,lbnodes,connec,bound,bou_codes, &
                     Ngp,gpcar,Ml,Mc,gpvol,dt,helem,Rgas,gamma_gas, &
                     rho,u,q,pr,E,Tem,e_int,mu_e)
           call nvtxEndRange

           !
           ! Call VTK output
           !
           call write_vtk_ascii(counter,ndime,npoin,nelem,nnode,coord,connec, &
                                rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e)

           counter = counter+1

        end do
        call nvtxEndRange

end program sod2d
