module mod_gpu_vars

   implicit none

#ifdef GPU

   ! Mesh information

   integer(4), device, allocatable :: connec_d(:,:)
   integer(4), device, allocatable :: lbnodes_d(:)

   ! Primary variables

   real(8), device, allocatable    :: mu_e_d(:)
   real(8), device, allocatable    :: rho_d(:,:)
   real(8), device, allocatable    :: pr_d(:,:)
   real(8), device, allocatable    :: E_d(:,:)
   real(8), device, allocatable    :: Tem_d(:,:)
   real(8), device, allocatable    :: e_int_d(:,:)
   real(8), device, allocatable    :: u_d(:,:,:)
   real(8), device, allocatable    :: q_d(:,:,:)

   ! Mass matrices

   real(8), device, allocatable    :: Ml_d(:)
   real(8), device, allocatable    :: Mc_d(:)

   ! Elemental info

   real(8), device, allocatable    :: gpvol_d(:,:,:)
   real(8), device, allocatable    :: gpcar_d(:,:,:,:)
   real(8), device, allocatable    :: Ngp_d(:,:)

#endif

end module mod_gpu_vars
