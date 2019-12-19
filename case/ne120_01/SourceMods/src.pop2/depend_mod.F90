module depend_mod
   use kinds_mod, only: r4, r8, int_kind, char_len, log_kind, rtavg
   use blocks, only: nx_block, ny_block, block, get_block
   use time_management, only: max_blocks_clinic, km, nt, mix_pass, c2dtt
public
#ifndef NOATH
   real (r8), dimension (nx_block,ny_block,2,max_blocks_clinic),public :: kxyu 
   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) :: & 
      KXU,KYU
#else
   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) :: & 
      KXU,KYU
#endif
   integer (int_kind),public :: &
      tavg_WVEL,         &! Vertical Velocity
      tavg_WVEL_2,       &! Vertical Velocity		
      tavg_WVEL2,        &! Vertical Velocity Squared
      tavg_UEU,          &! flux of zonal momentum across east  face
      tavg_VNU,          &! flux of zonal momentum across north face
      tavg_WTU,          &! flux of zonal momentum across top   face
      tavg_UEV,          &! flux of merid momentum across east  face
      tavg_VNV,          &! flux of merid momentum across north face
      tavg_WTV,          &! flux of merid momentum across top   face
      tavg_PV,           &! potential vorticity
      tavg_Q,            &! z-derivative of pot density
      tavg_PD,           &! potential density 
      tavg_RHOU,         &! pot density times U velocity
      tavg_RHOV,         &! pot density times V velocity
      tavg_PVWM,         &! pot vorticity flux through bottom
      tavg_PVWP,         &! pot vorticity flux through top
      tavg_UPV,          &! pot vorticity flux through east  face
      tavg_VPV,          &! pot vorticity flux through north face
      tavg_URHO,         &! pot density   flux through east  face
      tavg_VRHO,         &! pot density   flux through north face
      tavg_WRHO,         &! pot density   flux through top   face
      tavg_UQ,           &! advection of Q across east  face
      tavg_VQ             ! advection of Q across north face
end module
