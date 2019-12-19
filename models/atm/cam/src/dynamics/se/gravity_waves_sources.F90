module gravity_waves_sources
  use derivative_mod, only : derivative_t
  use dimensions_mod, only : np,nlev
  use edgetype_mod, only       : EdgeBuffer_t
  use element_mod, only    : element_t
  use hybrid_mod, only     : hybrid_t
  use kinds, only          : real_kind
  use shr_kind_mod, only   : r8 => shr_kind_r8
  use thread_mod, only     : horz_num_threads
  use dyn_comp, only       : dyn_import_t

  implicit none
  private
  save
  
  !! gravity_waves_sources created by S Santos, 10 Aug 2011
  !! 
  !! gws_src_fnct starts parallel environment and computes frontogenesis
  !!   for use by WACCM (via dp_coupling)
 
  public  :: gws_src_fnct
  public  :: gws_init
  private :: compute_frontogenesis
  
  type (EdgeBuffer_t) :: edge3
  type (derivative_t), allocatable   :: deriv(:)
  real(r8) :: psurf_ref

!----------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------

  subroutine gws_init(dyn_in)
    use parallel_mod, only   : par
    use edge_mod, only       : initEdgeBuffer
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    implicit none

    type (dyn_import_t), intent(inout) :: dyn_in
    
    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    call initEdgeBuffer(par, edge3,dyn_in%elem, 3*nlev)
!JMD Not quite sure this is needed
    allocate(deriv(0:horz_num_threads-1))
    
    psurf_ref = hypi(plev+1)

  end subroutine gws_init

  subroutine gws_src_fnct(elem, tl, frontgf, frontga)
    use derivative_mod, only  : derivinit
    use dimensions_mod, only  : npsq, nelemd
    use dof_mod, only         : UniquePoints
    use hybrid_mod, only      : config_thread_region, init_loop_ranges, get_loop_ranges
    use parallel_mod, only    : par
    use ppgrid, only          : pver
    use thread_mod, only      : omp_get_thread_num
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    integer, intent(in)          :: tl
    real (kind=real_kind), intent(out) :: frontgf(npsq,pver,nelemd)
    real (kind=real_kind), intent(out) :: frontga(npsq,pver,nelemd)
    
    ! Local variables
    type (hybrid_t) :: hybrid
    integer :: nets, nete, ithr, ncols, ie
    real(kind=real_kind), allocatable  ::  frontgf_thr(:,:,:,:)
    real(kind=real_kind), allocatable  ::  frontga_thr(:,:,:,:)
    
    call init_loop_ranges(nelemd)


    !$OMP PARALLEL NUM_THREADS(horz_num_threads),  DEFAULT(SHARED), PRIVATE(nets,nete,hybrid,ie,ncols,frontgf_thr,frontga_thr) 
    if (horz_num_threads == 1) then
      hybrid = config_thread_region(par,'serial')
    else
      hybrid = config_thread_region(par,'horizontal')
    endif
    call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
    write(*,200) nets, nete
    200 format(2x, 2i4, ' nets, nete - gravity')

!JMD call to derivinit... don't know what it needs to be a thread private call
    call derivinit(deriv(hybrid%ithr))
    allocate(frontgf_thr(np,np,nlev,nets:nete))
    allocate(frontga_thr(np,np,nlev,nets:nete))
    call compute_frontogenesis(frontgf_thr,frontga_thr,tl,elem,deriv(hybrid%ithr),hybrid,nets,nete)
    do ie=nets,nete
       ncols = elem(ie)%idxP%NumUniquePts
       call UniquePoints(elem(ie)%idxP, nlev, frontgf_thr(:,:,:,ie), frontgf(1:ncols,:,ie))
       call UniquePoints(elem(ie)%idxP, nlev, frontga_thr(:,:,:,ie), frontga(1:ncols,:,ie))
    end do
    deallocate(frontga_thr)
    deallocate(frontgf_thr)
    !$OMP END PARALLEL

  end subroutine gws_src_fnct
  
  subroutine compute_frontogenesis(frontgf,frontga,tl,elem,ederiv,hybrid,nets,nete)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenesis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  ! 
  ! Original by Mark Taylor, July 2011
  ! Change by Santos, 10 Aug 2011:
  ! Integrated into gravity_waves_sources module, several arguments made global
  !  to prevent repeated allocation/initialization
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physical_constants, only : kappa
    use derivative_mod, only : gradient_sphere, ugradv_sphere
    use edge_mod, only : edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use dyn_comp, only : hvcoord
    implicit none
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in) :: ederiv
    integer, intent(in) :: nets,nete
    integer, intent(in) :: tl ! timelevel to use
    real(kind=real_kind), intent(out) ::  frontgf(np,np,nlev,nets:nete)
    real(kind=real_kind), intent(out) ::  frontga(np,np,nlev,nets:nete)
  
    ! local
    integer :: k,kptr,i,j,ie,component
    real(kind=real_kind)  ::  gradth(np,np,2,nlev,nets:nete)  ! grad(theta)
    real(kind=real_kind)  :: p(np,np)        ! pressure at mid points
    real(kind=real_kind)  :: theta(np,np)    ! potential temperature at mid points
    real(kind=real_kind)  ::  C(np,np,2)     
    
    do ie=nets,nete
       do k=1,nlev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! potential temperature: theta = T (p/p0)^kappa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          
          ! pressure at mid points
          p(:,:)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
          theta(:,:) = elem(ie)%state%T(:,:,k,tl)*(psurf_ref / p(:,:))**kappa
          ! gradth(:,:,:,k,ie) = gradient_sphere(theta,ederiv,elem(ie)%Dinv)
          call gradient_sphere(theta,ederiv,elem(ie)%Dinv,gradth(:,:,:,k,ie))
          
          ! compute C = (grad(theta) dot grad ) u
          ! C(:,:,:) = ugradv_sphere(gradth(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),ederiv,elem(ie))
          call ugradv_sphere(gradth(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),ederiv,elem(ie),C(:,:,:))
          
          ! gradth dot C
          frontgf(:,:,k,ie) = -( C(:,:,1)*gradth(:,:,1,k,ie) +  C(:,:,2)*gradth(:,:,2,k,ie)  )
          
          ! apply mass matrix
          gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
          gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
          frontgf(:,:,k,ie)=frontgf(:,:,k,ie)*elem(ie)%spheremp(:,:)
          
       enddo
       ! pack
       call edgeVpack(edge3, frontgf(:,:,:,ie),nlev,0,ie)
       call edgeVpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
    enddo
    call bndry_exchangeV(hybrid,edge3)
    do ie=nets,nete
       call edgeVunpack(edge3, frontgf(:,:,:,ie),nlev,0,ie)
       call edgeVunpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
       ! apply inverse mass matrix,
       do k=1,nlev
          gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
          gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
          frontgf(:,:,k,ie)=frontgf(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          
          ! Frontogenesis angle
          frontga(:,:,k,ie) = atan2 ( gradth(:,:,2,k,ie) , gradth(:,:,1,k,ie) + 1.e-10_real_kind )
       enddo
    enddo
  end subroutine compute_frontogenesis


end module gravity_waves_sources
