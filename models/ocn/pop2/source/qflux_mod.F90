!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module qflux_mod

!BOP
! !MODULE: qflux_mod
!
! !DESCRIPTION:
!  This module supports time-averaged qflux computations
!
! !REVISION HISTORY:
! SVN:$Id: qflux_mod.F90 87926 2017-12-12 00:31:20Z nanr $
!

! !USES:

   use kinds_mod
   use constants
   use exit_mod
   use time_management
   use tavg

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_qflux
 
!EOP
!BOC

   integer (int_kind) :: &
      tavg_QFLUX,        &! tavg id for QFLUX 
      tavg_QFLUX_2        ! tavg id for QFLUX  daily

!EOC
!***********************************************************************

      contains

!***********************************************************************

   subroutine init_qflux

   character (char_len) :: string
 
   
   string = 'Internal Ocean Heat Flux Due to Ice Formation; ' /&
             &/ 'heat of fusion > 0 or ice-melting potential < 0 '
 
   call define_tavg_field(tavg_QFLUX,'QFLUX',2,               &
                          long_name=trim(string),             &
                          tavg_method=tavg_method_qflux,      &
                          units='Watts/meter^2',              &
                          coordinates  ='TLONG TLAT time',    &
                          grid_loc ='2111')
   call define_tavg_field(tavg_QFLUX_2,'QFLUX_2',2,               &
                          long_name=trim(string),             &
                          tavg_method=tavg_method_qflux,      &
                          units='Watts/meter^2',              &
                          coordinates  ='TLONG TLAT time',    &
                          grid_loc ='2111')

   end subroutine init_qflux
 
 
!***********************************************************************

 end module qflux_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
