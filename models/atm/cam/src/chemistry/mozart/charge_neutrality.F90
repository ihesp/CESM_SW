
      module charge_neutrality

      use shr_kind_mod,      only : r8 => shr_kind_r8
      use cam_logfile,       only : iulog

      implicit none

      private
      public :: charge_balance
      public :: charge_fix     ! temporary, for fixing charge balance after vertical diffusion
                               ! without converting mass mixing ratios to volume
                               ! mean mass assumed to be mwdry

      contains

      subroutine charge_balance( ncol, conc )
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

        use ppgrid,       only : pver
        use mo_chem_utls, only : get_spc_ndx
        use chem_mods,    only : gas_pcnst

        implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      integer,  intent(in)          :: ncol
      real(r8), intent(inout)       :: conc(ncol,pver,gas_pcnst)         ! concentration

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n
      integer  :: elec_ndx
      real(r8) :: wrk(ncol,pver)

      elec_ndx = get_spc_ndx('e')
#ifdef CB_DEBUG
      write(iulog,*) ' '
      write(iulog,*) '------------------------------------------------------------------'
      write(iulog,*) 'charge_balance: e ndx,offset = ',elec_ndx,offset
      write(iulog,*) 'charge_balance: size of conc = ',size(conc,dim=1),' x ',size(conc,dim=2),' x ',size(conc,dim=3)
#endif
      if( elec_ndx > 0 ) then
	 wrk(:,:) = 0._r8
         n = get_spc_ndx('Np')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('N2p')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('Op')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('O2p')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('NOp')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
#ifdef CB_DEBUG
         write(iulog,*) 'charge_balance: electron concentration before balance'
         write(iulog,'(1p,5g15.7)') conc(1,:,elec_ndx)
         write(iulog,*) 'charge_balance: electron concentration after  balance'
         write(iulog,'(1p,5g15.7)') wrk(1,:)
         write(iulog,*) '------------------------------------------------------------------'
         write(iulog,*) ' '
#endif
         conc(:ncol,:,elec_ndx) = wrk(:ncol,:)
      end if

      end subroutine charge_balance


      subroutine charge_fix(state, pbuf)
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

      use ppgrid,              only : pcols, pver
      use constituents,        only : cnst_get_ind, cnst_mw
      use physconst,           only : mwdry                   ! molecular weight of dry air
      use physconst,           only : mbarv                       ! Constituent dependent mbar
      use phys_control,        only : waccmx_is
      use short_lived_species, only : slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species in pbuf
      use mo_chem_utls,        only : get_spc_ndx
      use chem_mods,           only : adv_mass
      use physics_buffer,      only : pbuf_get_field,physics_buffer_desc ! Needed to get variables from physics buffer
      use physics_types,       only : physics_state
      
      implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      type(physics_state), intent(inout), target :: state
      type(physics_buffer_desc), pointer :: pbuf(:)    ! physics buffer

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n, ns, nc
      integer  :: elec_ndx, elec_sndx
      integer  :: lchnk                 !Chunk number from state structure
      integer  :: ncol                  !Number of columns in this chunk from state structure

      real(r8) :: wrk(pcols,pver)
      real(r8) :: mbar(pcols,pver)  ! mean mass (=mwdry) used to fake out optimizer to get
                                    ! identical answers to old code
				   
      real(r8), dimension(:,:,:), pointer   :: q         ! model mass mixing ratios
      real(r8), dimension(:,:),   pointer   :: qs        ! Pointer to access fields in pbuf
      real(r8), dimension(:,:),   pointer   :: qse       ! Pointer to access electrons in pbuf

!-----------------------------------------------------------------------      
      lchnk = state%lchnk
      ncol  = state%ncol
      q => state%q

     !------------------------------------------------
     ! assume that mbar = mwdry except for WACCM-X
     !------------------------------------------------
      mbar(:ncol,:) = mwdry
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) mbar(:ncol,:) = mbarv(:ncol,:,lchnk)

     !-----------------------------------------------------
     ! Get index to access electron mass mixing ratio
     !-----------------------------------------------------     
      call cnst_get_ind( 'e', elec_ndx, abort=.false. )  
      if (elec_ndx < 0) elec_sndx = slvd_index( 'e' )
    
      !--------------------------------------------------------------------
      ! If electrons are in state%q or pbuf, add up ions to get electrons
      !--------------------------------------------------------------------  
      if( elec_ndx > 0 .or. elec_sndx > 0) then
	 wrk(:,:) = 0._r8
         call cnst_get_ind( 'Np', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'Np' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('Np')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'N2p', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'N2p' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('N2p')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'Op', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'Op' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('Op')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'O2p', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'O2p' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('O2p')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'NOp', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'NOp' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('NOp')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
	 
         !--------------------------------------------------------------------------------------
         !  Total ions now in wrk array so determine electrons.  qse is a pointer to pbuf and
	 !  q is a pointer to state%q 
         !--------------------------------------------------------------------------------------
         if (elec_ndx < 0) then 
            ns = slvd_index( 'e' )
            call pbuf_get_field(pbuf, slvd_pbf_ndx, qse, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('e')
      	    do k = 1,pver
	      qse(:ncol,k) = adv_mass(nc) * wrk(:ncol,k) / mbar(:ncol,k)
	    end do	   
         else
           do k = 1,pver
	     q(:ncol,k,elec_ndx) = cnst_mw(elec_ndx) * wrk(:ncol,k) / mbar(:ncol,k)
	   end do
	 end if  
      end if

      end subroutine charge_fix

      end module charge_neutrality
