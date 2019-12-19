
 module POP_SolversMod

!BOP
! !MODULE: POP_SolversMod
!
! !DESCRIPTION:
!  This module contains routines and operators for solving the elliptic
!  system for surface pressure in the barotropic mode.
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use blocks, only: block, get_block
   use POP_ErrorMod
   use POP_CommMod
   use POP_ConfigMod
   use POP_IOUnitsMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_GridHorzMod
   use POP_ReductionsMod
   use POP_BroadcastMod
   use POP_RedistributeMod
   use POP_HaloMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_DomainSizeMod
   use time_management,only:matsuno_ts
   use domain
   use grid
   use time_management

#ifndef  NOATH
#define ZYU_ATH

!!#define ZYU_SAMEDISTRIBUTE
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_SolversInit,     &
             POP_SolversPrep,     &
             POP_SolversRun,      &
             POP_SolversDiagonal, &
             POP_SolversGetDiagnostics

! !PUBLIC DATA MEMBERS:

!-----------------------------------------------------------------------
!
!  other operator and preconditioning weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (POP_r8), dimension (:,:,:), allocatable, public :: & 
      mMaskTropic,     &! land mask in barotropic distribution 
      btropWgtCenter,  &! barotropic operater center coefficient
      btropWgtNorth,   &! barotropic operater north  coefficient
      btropWgtEast,    &! barotropic operater east   coefficient
      btropWgtNE        ! barotropic operater northeast coefficient

   real (POP_r8), dimension (:,:,:), allocatable, public :: & 
      centerWgtClinicIndep, &! time indep  center wgt on clinic distrb
      centerWgtClinic        ! time depend center wgt on clinic distrb  

!-----------------------------------------------------------------------
!
!  EVP preconditioner influence matrix and block coefficient matrix,
!  Saved in EVP sub-block format. Dimension switched:
!  (X,Y,block_ID) ==> (subblock_X, subblock_Y, subblock_ID, block_ID) 
!-----------------------------------------------------------------------

   real (POP_r8), dimension (:,:,:,:), allocatable, public :: & 
      EvpRinv               ! EVP influence matrix

   real (POP_r8), dimension (:,:,:,:), allocatable, public :: & 
      EvpCenterWgt,        &! reshape coefficients into blocks 
      InvEvpCenterWgt,     &! save inverse for computational efficiency
      EvpNeWgt,            &! 
      InvEvpNeWgt          

   integer (POP_i4),dimension(:,:,:),allocatable :: &
      landIndx              ! index of land EVP blocks 

!EOP
!BOC
   !*** preconditioning operator coefficients (ninept operator)

   real (POP_r8), dimension (:,:,:), allocatable :: & 
      precondCenter,              &
      precondNorth, precondSouth, &
      precondEast,  precondWest,  &
      precondNE,    precondSE,    &
      precondNW,    precondSW

!-----------------------------------------------------------------------
!
!  supported solvers and preconditioners
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      solverChoice

   character (3), parameter :: &
      solverChoicePCG = 'pcg'

   character (4), parameter :: & 
      solverChoicePCSI = 'PCSI'
   character (9), parameter :: &
      solverChoiceChronGear = 'ChronGear'

   character (POP_charLength) :: &
      preconditionerChoice

   character (8), parameter :: &
      precondChoiceDiag = 'diagonal'
   character (4), parameter :: &
      precondChoiceFile = 'file'
   character (3), parameter :: &
      precondChoiceEvp = 'evp' 

   logical (POP_logical) :: &
      usePreconditioner

!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::      &
      maxIterations,        &! max number of solver iterations
      convergenceCheckStart,&! start checking from given steps
      convergenceCheckFreq   ! check convergence every freq steps

   real (POP_r8) ::         &
      convergenceCriterion, &! convergence error criterion
      LanczosconvergenceCriterion, &! Lanczos convergence criterion
      residualNorm           ! residual normalization

   !*** convergence diagnostics

   integer (POP_i4), public ::      &
      numIterations          ! accumulated no of iterations (diagnostic)

   real (POP_r8) ::         &
      rmsResidual            ! residual (also a diagnostic)

   real (POP_r8) ::        & 
      PcsiMaxEigs,         &! Largest eigenvalues for PCSI method
      PcsiMinEigs           ! smallest eigenvalues for PCSI method

   integer(POP_i4), parameter  :: &
      EvpXbs  = 8,        &! EVP block size along longitude
      EvpYbs  = 8          ! EVP block size along latitude

   integer (POP_i4) :: &
      MaxLanczosStep        ! Max lanczos steps to get eigenvalues 

   integer (POP_i4) :: &
      EvpXnb,          &    ! sub block number along longitude
      EvpYnb                ! sub block number along latitude

   integer (POP_i4),dimension(:),allocatable :: &
      EvpXbidx,          &  ! index of sub blocks start point 
      EvpYbidx              ! along longtitude and latitude

   integer,save :: SAMEDISTRIBUTE =0

!EOC
!***********************************************************************

 contains
!***********************************************************************
!BOP
! !SUBROUTINE: POP_SolversPrep
! !INTERFACE:

 subroutine POP_SolversPrep(errorCode)

! !DESCRIPTION:
! solver preprocessing : 
! 1. compute eigenvalues if PCSI is used 
! 2. Preprocessing EVP influence matrix if EVP preconditioning is used 
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu et al., Tsinghua University
!BOP
! !INPUT PARAMETERS:


! !INPUT/OUTPUT PARAMETERS:


! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code
   integer (POP_i4) :: &
      numProcs           ! number of processors in barotropic distrib

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      numBlocks,          &! number of blocks in barotropic distrib
      bid,                &! local counters
      nx1,ny1             ! POP_nxBlock -1


   errorCode = POP_Success
   if (POP_myTask == POP_masterTask) then
      write(POP_stdout,POP_blankFormat)
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,'(a42)') ' Solver Preprocessing (barotropic solver)'
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,POP_blankFormat)
   endif

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numProcs = numProcs)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error getting num procs')
      return
   endif

   call POP_RedistributeBlocks(btropWgtCenter,  POP_distrbTropic, &
                               centerWgtClinic, POP_distrbClinic, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing center operator weight')
      return
   endif

   if (POP_myTask < numProcs) then
      call POP_DistributionGet(POP_distrbTropic, errorCode, &
                               numLocalBlocks = numBlocks)

!-----------------------------------------------------------------------
!
!   EVP preprocessing
!
!-----------------------------------------------------------------------
      if (trim(preconditionerChoice) == precondChoiceEvp) then

          if (POP_myTask == POP_masterTask) then
            write(POP_stdout,'(a50)')  & 
              'Solver EVP preconditioning START '
          endif

          EvpCenterWgt(:,:,:,:) = 0.0
          InvEvpCenterWgt(:,:,:,:) = 0.0
          EvpNeWgt(:,:,:,:) = 0.0
          InvEvpNeWgt(:,:,:,:) = 0.0

 
  !call getlooptime(start_POP_SolversMod_loop1)
          do bid=1,numBlocks
            nx1 = POP_nxBlock-1
            ny1 = POP_nyBlock-1
            call EvpPre(btropWgtCenter(2:nx1,2:ny1,bid),btropWgtNE(2:nx1,2:ny1,bid),&
                        EvpCenterWgt(:,:,:,bid),EvpNeWgt(:,:,:,bid),EvpRinv(:,:,:,bid),&
                        landIndx(:,:,bid), nx1-1,ny1-1,EvpXbs,EvpYbs,&
                        EvpXbidx,EvpYbidx,EvpXnb,EvpYnb,bid)
            

            ! Save inverse of coefficents for efficiency
            where(EvpCenterWgt(:,:,:,bid) .ne. 0.0) 
              InvEvpCenterWgt(:,:,:,bid) = &
                          1.0_POP_r8/EvpCenterWgt(:,:,:,bid)
            end where
            where(EvpNeWgt(:,:,:,bid) .ne. 0.0)
              InvEvpNeWgt(:,:,:,bid) = 1.0_POP_r8/EvpNeWgt(:,:,:,bid)
            end where

          end do 
  !call getlooptime(stop_POP_SolversMod_loop1)
  !if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop1 costs w zhuangy: ',stop_POP_SolversMod_loop1-start_POP_SolversMod_loop1,' (us);'
  !endif
 

          if (POP_myTask == POP_masterTask) then
            write(POP_stdout,'(a35)') ' Solver EVP preconditioning END'
          endif

      endif

   
!-----------------------------------------------------------------------
!
! call lanczos function to compute eigenvalues for PCSI
!
!-----------------------------------------------------------------------

      if(trim(solverChoice) .eq. solverChoicePCSI) then
        call PcsiLanczos(POP_nxBlock, POP_nyBlock, numBlocks, &
                    MaxLanczosStep, PcsiMaxEigs,PcsiMinEigs, errorCode)

        if (POP_myTask == POP_masterTask) then 
            write(POP_stdout,*) "LANCZOS EIGS: ",PcsiMinEigs, PcsiMaxEigs
        endif 
      endif 

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversRun: error computing eigenvalues in Lanczos')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversPrep

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversRun
! !INTERFACE:

 subroutine POP_SolversRun(sfcPressure, rhsClinic, errorCode)

! !DESCRIPTION:
!  Solves the elliptic equation for surface pressure by calling
!  the requested solver routine.  Also redistributes necessary
!  array to the barotropic distribution of blocks for better performance
!  of the solver.
!  The elliptic equation is
!  \begin{equation}
!     AF = B
!  \end{equation}
!  where $F$ is a field (eg surface pressure), $B$ is the right hand side
!  and $A$ is the operator defined as
!  \begin{equation}
!     AF = a \nabla\cdot(H \nabla F)
!  \end{equation}
!  where $a$ is the cell area.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      rhsClinic         ! right-hand-side of linear system
                        !  for blocks in baroclinic distribution

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
     sfcPressure              ! on input,  initial guess in baroclinic distrb
                        ! on output, final solution for sfc pressure

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      numProcs           ! number of processors in barotropic distrib

   real (POP_r8), dimension(size(sfcPressure,dim=1), &
                            size(sfcPressure,dim=2), &
                            size(sfcPressure,dim=3)) :: &
      pressTropic,     &! surface pressure in barotropic distribution
      rhsTropic         ! right hand side  in barotropic distribution
  !!  integer :: st=0,ed=0
!-----------------------------------------------------------------------
!
!  switch to the barotropic distribution for iterative solvers
!
!-----------------------------------------------------------------------



   errorCode = POP_Success

   if(SAMEDISTRIBUTE==1) then

     btropWgtCenter  =  centerWgtClinic
     pressTropic     =  sfcPressure
     rhsTropic       =  rhsClinic

   else

   call POP_RedistributeBlocks(btropWgtCenter,  POP_distrbTropic, &
                               centerWgtClinic, POP_distrbClinic, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing center operator weight')
      return
   endif

   call POP_RedistributeBlocks(pressTropic, POP_distrbTropic, &
                               sfcPressure, POP_distrbClinic, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing pressure')
      return
   endif

   call POP_RedistributeBlocks(rhsTropic, POP_distrbTropic, &
                               rhsClinic, POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing right hand side')
      return
   endif


   endif
!-----------------------------------------------------------------------
!
!  call proper routine based on user choice of solver
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numProcs = numProcs)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error getting num procs')
      return
   endif

   if (POP_myTask < numProcs) then

!-----------------------------------------------------------------------
!
!  solver preprocessing if PCSI or EVP preconditioner is used 
!
!-----------------------------------------------------------------------

      select case(trim(solverChoice))
!!280000-150000us
      case (solverChoicePCSI)  
     !!label
     !! call   getlooptime(st)
         call PCSI(pressTropic, rhsTropic, errorCode)
     !! call   getlooptime(ed)
     !!  if (POP_myTask == POP_masterTask) then
      !!           WRITE (*,*) 'PCSI cost=',ed-st,'us'
      !! end if

      
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolverRun: error in PCSI')
            return
         endif

      case (solverChoicePCG)

         call pcg(pressTropic, rhsTropic, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolverRun: error in PCG')
            return
         endif

      case (solverChoiceChronGear)

         call ChronGear(pressTropic, rhsTropic, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolverRun: error in ChronGear')
            return
         endif

      end select
   endif

!-----------------------------------------------------------------------
!
!  switch solution back to the baroclinic distribution
!
!-----------------------------------------------------------------------



   if(SAMEDISTRIBUTE==1) then
      sfcPressure = pressTropic
   

   else


   call POP_RedistributeBlocks(sfcPressure, POP_distrbClinic, &
                               pressTropic, POP_distrbTropic, errorCode)



   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing pressure back')
      return
   endif

   endif
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversRun

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversInit
! !INTERFACE:

 subroutine POP_SolversInit(errorCode)

! !DESCRIPTION:
!  This routine initializes choice of solver, calculates the 
!  coefficients of the 9-point stencils for the barotropic operator and
!  reads in a preconditioning matrix if requested.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
!         from {x,y} components of divergence
!       HU = depth at U points
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,lEvp,nbEvp,  &! dummy loop counters
      configUnit,      &! unit for configuration file
      numBlocksClinic, &! num local blocks in baroclinic distribution
      numBlocksTropic, &! num local blocks in barotropic distribution
      iblock,          &! block counter
      istat             ! status flag for allocates

   character (POP_charLength) :: &
      preconditionerFile  ! file containing preconditioner

   real (POP_r8) ::    &
      xne,xse,xnw,xsw, &! contribution to coefficients from x,y
      yne,yse,ynw,ysw, &!   components of divergence
      ase,anw,asw

   real (POP_r8), dimension(:,:,:), allocatable :: &
      work0,           &! temp space for computing barotropic
      workNorth,       &!
      workEast,        &!
      workNE,          &!
      mMaskTmp          ! land mask in barotropic distribution

!    real(POP_r8) :: const_para
!    integer,save :: init_flag=0 
!-----------------------------------------------------------------------
!
!  read solver choice and solver constants from configuration file
!
!-----------------------------------------------------------------------
!   if(init_flag.eq.0) then
!        call  athread_init()
!          init_flag=1
!   end if

   errorCode = POP_Success

   if (POP_myTask == POP_masterTask) then
      write(POP_stdout,POP_blankFormat)
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,'(a35)') ' Solver options (barotropic solver)'
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,POP_blankFormat)
   endif

   call POP_ConfigOpen(configUnit, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error opening config file')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'solverChoice',    &
   call POP_ConfigRead(configUnit, 'solvers', 'solverchoice',    &
                       solverChoice, solverChoicePCG, errorCode, &
                       outStringBefore = 'Solver choice: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver choice')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'convergenceCriterion', &
   call POP_ConfigRead(configUnit, 'solvers', 'convergencecriterion', &
                   convergenceCriterion, 1.e-12_POP_r8, errorCode,    &
                   outStringBefore = 'Solver converged for err < ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver convergence criterion')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'maxIterations',    &
   call POP_ConfigRead(configUnit, 'solvers', 'maxiterations',    &
                       maxIterations, 1000, errorCode, &
                       outStringBefore = 'Solver maximum iterations: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver max iterations')
      return
   endif
!-----------------------------------------------------------------------
!
! Read in Lanczos parameters for PCSI
!
!-----------------------------------------------------------------------

!  call POP_ConfigRead(configUnit, 'solvers', 'LanczosconvergenceCriterion', &
   call POP_ConfigRead(configUnit, 'solvers', 'lanczosconvergencecriterion', &
                   LanczosconvergenceCriterion, 0.1_POP_r8, errorCode,    &
                   outStringBefore = 'Lanczos converged for err < ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading Lanczos convergence criterion')
      return
   endif

   call POP_ConfigRead(configUnit, 'solvers', 'maxlanczosstep', &
                       MaxLanczosStep, 20, errorCode,           &
                       outStringBefore = 'Maximum Lanczos step ',  &
                       outStringAfter  = ' iterations')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading lanczos step')
      return
   endif
!  call POP_ConfigRead(configUnit, 'solvers', 'convergenceCheckstart', &
   call POP_ConfigRead(configUnit, 'solvers', 'convergencecheckstart', &
                       convergenceCheckStart, 60, errorCode,           &
                       outStringBefore = 'Check convergence start from ',  &
                       outStringAfter  = ' iterations')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading convergence check start steps')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'convergenceCheckFreq', &
   call POP_ConfigRead(configUnit, 'solvers', 'convergencecheckfreq', &
                       convergenceCheckFreq, 10, errorCode,           &
                       outStringBefore = 'Check convergence every ',  &
                       outStringAfter  = ' iterations')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading convergence check frequency')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'preconditionerChoice', &
   call POP_ConfigRead(configUnit, 'solvers', 'preconditionerchoice', &
                       preconditionerChoice, precondChoiceDiag,       &
                       errorCode,                                     &
                       outStringBefore = 'Preconditioner choice: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading preconditioner choice')
      return
   endif

   if (trim(preconditionerChoice) == precondChoiceFile) then
!     call POP_ConfigRead(configUnit, 'solvers', 'preconditionerFile', &
      call POP_ConfigRead(configUnit, 'solvers', 'preconditionerfile', &
                 preconditionerFile, 'UnknownPrecondFile', errorCode,  &
                 outStringBefore = 'Reading preconditioner from file: ')

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversInit: error reading preconditioner file name')
         return
      endif
   endif

   call POP_ConfigClose(configUnit, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error closing config file')
      return
   endif

!-----------------------------------------------------------------------
!
!  check inputs
!
!-----------------------------------------------------------------------

   select case(trim(solverChoice))
   case(solverChoicePCG)
   case(solverChoicePCSI) 
   case(solverChoiceChronGear)
   case default
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: unknown solver - must be pcg, ChronGear, PCSI')
      return
   end select


   select case (trim(preconditionerChoice))
   case(precondChoiceDiag)
      usePreconditioner = .false.   ! default is diagonal
   case(precondChoiceEvp)
      usePreconditioner = .true.
   case(precondChoiceFile)
      usePreconditioner = .true.
   case default
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: unknown preconditioner choice')
      return
   end select

!-----------------------------------------------------------------------
!
!   Matsnuo timestep is not supported in PCSI or EVP. Because it  uses 
!   a half timestep, which introduces a different coefficient matrix A. 
!
!-----------------------------------------------------------------------
   if (matsuno_ts .and. (trim(preconditionerChoice) == precondChoiceEvp &
                   .or.  trim(solverChoice) == solverChoicePCSI)) then
      call POP_ErrorSet(errorCode, &
      'POP_SolversInit: Matsnuo timestep is not supported in &
       PCSI solver or EVP preconditioning')
      return
   endif

!-----------------------------------------------------------------------
!
!  compute nine point operator coefficients: compute on baroclinic
!  decomposition first where grid info defined and redistribute
!  to barotropic distribution
!  leave center coefficients in baroclinic distribution to facilitate 
!  easy time-dependent changes in barotropic routine
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbClinic, errorCode, &
                            numLocalBlocks = numBlocksClinic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error retrieving clinic local block count')
      return
   endif

   allocate(work0     (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workNorth (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workEast  (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workNE    (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            mMaskTmp  (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
       centerWgtClinicIndep(POP_nxBlock,POP_nyBlock,numBlocksClinic), &
       centerWgtClinic     (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error allocating temporary arrays')
      return
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,i,j,xne,xse,xnw,xsw,yne,yse,ynw,ysw, &
   !$OMP                     ase,anw,asw)

!!attention
    
   !! const_para = 1.0_POP_r8
   !! call loop_1(const_para,POP_nxblock,POP_nyblock,numBlocksClinic) 
      
      do iblock = 1,numBlocksClinic
      work0                (:,:,iblock) = 0.0_POP_r8
      workNorth            (:,:,iblock) = 0.0_POP_r8
      workEast             (:,:,iblock) = 0.0_POP_r8
      workNE               (:,:,iblock) = 0.0_POP_r8
      mMaskTmp             (:,:,iblock) = 0.0_POP_r8
      centerWgtClinicIndep (:,:,iblock) = 0.0_POP_r8
      do j=2,POP_nyBlock
      do i=2,POP_nxBlock
         xne = 0.25_POP_r8*HU(i  ,j  ,iblock)*DXUR(i  ,j  ,iblock)* &
                                              DYU (i  ,j  ,iblock)
         xse = 0.25_POP_r8*HU(i  ,j-1,iblock)*DXUR(i  ,j-1,iblock)* &
                                              DYU (i  ,j-1,iblock)
         xnw = 0.25_POP_r8*HU(i-1,j  ,iblock)*DXUR(i-1,j  ,iblock)* &
                                              DYU (i-1,j  ,iblock)
         xsw = 0.25_POP_r8*HU(i-1,j-1,iblock)*DXUR(i-1,j-1,iblock)* &
                                              DYU (i-1,j-1,iblock)

         yne = 0.25_POP_r8*HU(i  ,j  ,iblock)*DYUR(i  ,j  ,iblock)* &
                                              DXU (i  ,j  ,iblock)
         yse = 0.25_POP_r8*HU(i  ,j-1,iblock)*DYUR(i  ,j-1,iblock)* &
                                              DXU (i  ,j-1,iblock)
         ynw = 0.25_POP_r8*HU(i-1,j  ,iblock)*DYUR(i-1,j  ,iblock)* &
                                              DXU (i-1,j  ,iblock)
         ysw = 0.25_POP_r8*HU(i-1,j-1,iblock)*DYUR(i-1,j-1,iblock)* &
                                              DXU (i-1,j-1,iblock)
         workNE(i,j,iblock) = xne + yne
         ase                = xse + yse
         anw                = xnw + ynw
         asw                = xsw + ysw
         workEast (i,j,iblock)  = xne + xse - yne - yse
         workNorth(i,j,iblock)  = yne + ynw - xne - xnw

         centerWgtClinicIndep(i,j,iblock) = &
                        -(workNE(i,j,iblock) + ase + anw + asw)

         work0   (i,j,iblock) = TAREA (i,j,iblock)**2
         mMaskTmp(i,j,iblock) = RCALCT(i,j,iblock)

      end do
      end do
  end do
 

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  redistribute operator weights and mask to barotropic distribution
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocksTropic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error retrieving tropic local block count')
      return
   endif

   allocate(btropWgtCenter (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtNorth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtEast   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtNE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            mMaskTropic    (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error allocating operator weights')
      return
   endif

   call POP_RedistributeBlocks(btropWgtNorth, POP_distrbTropic, &
                               workNorth,     POP_distrbClinic, errorCode)

   if(errorCode == 0) then
        SAMEDISTRIBUTE = 1
   endif     

  ! if (errorCode /= POP_Success) then
  !    call POP_ErrorSet(errorCode, &
  !       'POP_SolversInit: error redistributing north operator weight')
  !    return
  ! endif

   call POP_RedistributeBlocks(btropWgtEast, POP_distrbTropic, &
                               workEast,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing east operator weight')
      return
   endif

   call POP_RedistributeBlocks(btropWgtNE, POP_distrbTropic, &
                               workNE,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing NE operator weight')
      return
   endif

   call POP_RedistributeBlocks(mMaskTropic, POP_distrbTropic, &
                               mMaskTmp,    POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing mask')
      return
   endif

!-----------------------------------------------------------------------
!
!  calculate normalization constant (darea,darea) for rmsResidual
!  in cgr routine.
!
!-----------------------------------------------------------------------

   residualNorm = 1.0_POP_r8/POP_GlobalSum(work0, POP_distrbClinic, &
                                           POP_gridHorzLocCenter,   &
                                           errorCode,               &
                                           mMask = mMaskTmp)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error computing normalization factor')
      return
   endif

   convergenceCriterion = convergenceCriterion**2/residualNorm

   deallocate(mMaskTmp, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error deallocating temp mask')
      return
   endif

!-----------------------------------------------------------------------
!
!  setup preconditioner if required
!
!-----------------------------------------------------------------------

   if (usePreconditioner) then

     if (trim(preconditionerChoice) == precondChoiceFile) then
       call POP_ErrorSet(errorCode, &
          'POP_SolversInit: preconditioner not supported')
       return

       allocate(                                                   &
          precondCenter (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondNorth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondSouth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondEast   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondWest   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondNE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondSE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondNW     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          precondSW     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
          stat = istat)

       if (istat > 0) then
          call POP_ErrorSet(errorCode, &
             'POP_SolversInit: error allocating preconditioner')
          return
       endif

     else if (trim(preconditionerChoice) == precondChoiceEvp) then


         call EvpBlockPartition(POP_nxBlock-2,EvpXbs,EvpXnb,EvpXbidx,errorCode)
         call EvpBlockPartition(POP_nyBlock-2,EvpYbs,EvpYnb,EvpYbidx,errorCode)
         
         if (POP_myTask == POP_masterTask) then
           write(POP_stdout,'(a35)') 'EVP block distribution'
           write(POP_stdout,*) EvpXbidx(:)
           write(POP_stdout,*) EvpYbidx(:)
         endif

         lEvp  = EvpXbs + EvpYbs -1
         nbEvp = EvpXnb*EvpYnb
         allocate(landIndx   (EvpXnb,EvpYnb,numBlocksTropic)          ,&         
                  EvpRinv    (levp,levp,EvpXnb*EvpYnb,numBlocksTropic),&
                  EvpNeWgt   (EvpXbs+2,EvpYbs+2,nbEvp,numBlocksTropic),&
                  InvEvpNeWgt(EvpXbs+2,EvpYbs+2,nbEvp,numBlocksTropic),&
                 EvpCenterWgt(EvpXbs+2,EvpYbs+2,nbEvp,numBlocksTropic),&
              InvEvpCenterWgt(EvpXbs+2,EvpYbs+2,nbEvp,numBlocksTropic),&
           stat = istat )

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversInit: error allocating EVP arrays')
            return
         endif
     endif 
   endif

!-----------------------------------------------------------------------
!
!  clean up temporary arrays
!
!-----------------------------------------------------------------------

   deallocate(work0, workNorth, workEast, workNE, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error deallocating temp mask')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversInit

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversDiagonal
! !INTERFACE:

 subroutine POP_SolversDiagonal(diagonalCorrection, blockIndx, errorCode)

! !DESCRIPTION:
!  This routine corrects the center barotropic operator by 
!  subtracting off the time dependent diagnonal term computed in
!  the barotropic driver.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:), intent(in) :: &
      diagonalCorrection   ! time dependent diagonal term

   integer (POP_i4), intent(in) :: &
      blockIndx            ! local block index

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  subtract the time dependent diagonal term from the center operator
!  coefficient
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
  !!call myrpcc(start_POP_SolversMod_loop17)
!  call loop_diag(blockIndx,centerWgtClinic,centerWgtClinicIndep,diagonalCorrection)//slower ?
   
  
      centerWgtClinic(:,:,blockIndx) =                               &
                            centerWgtClinicIndep(:,:,blockIndx) - &
                             diagonalCorrection(:,:)

 !   if (POP_myTask == POP_masterTask) then
 !   WRITE (*,*) '#####START######',POP_nxblock,POP_nyblock,loc(diagonalCorrection(1,2)),loc(diagonalCorrection(1,1)),'##### END ########'
!!  WRITE (*,*) 'POP_SolversMod_sigal-org costszhuangy:',stop_POP_SolversMod_loop1-start_POP_SolversMod_loop1,' (us);'
 !   endif
                  
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversDiagonal

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversGetDiagnostics
! !INTERFACE:

 subroutine POP_SolversGetDiagnostics(iterationCount, residual, &
                                      errorCode)

! !DESCRIPTION:
!  This routine returns the latest iteration count and residual
!  for diagnostic purposes.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
      residual          ! latest residual from solver

   integer (POP_i4), intent(out) :: &
      iterationCount,  &! latest iteration count from solver
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  return the diagnostic quantities
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   residual = rmsResidual
   iterationCount = numIterations

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversGetDiagnostics

!***********************************************************************
!BOP
! !IROUTINE: pcg
! !INTERFACE:

 subroutine pcg(X,B,errorCode)

! !DESCRIPTION:
!  This routine uses a preconditioned conjugate-gradient solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,m,           &! local iteration counter
      nx, ny,          &! horizontal extents
      numBlocks,       &! number of local blocks
      iblock            ! local block     counter

   real (POP_r8) ::          &
      eta0,eta1,rr        ! scalar inner product results

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2), &
                                          size(X,dim=3)) :: &
      R,                 &! residual (b-Ax)
      S,                 &! conjugate direction vector
      Q,work0,work1       ! various cg intermediate results

!-----------------------------------------------------------------------
!
!  compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocks)

   nx = size(X,dim=1)
   ny = size(X,dim=2)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversPCG: error retrieving local block count')
      return
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

  ! call getlooptime(start_POP_SolversMod_loop3)
   do iblock=1,numBlocks
      call btropOperator(S,X,iblock)
          
      do j=1,ny
      do i=1,nx
         R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
         S(i,j,iblock) = 0.0_POP_r8
      end do
      end do
   end do ! block loop
  !call getlooptime(stop_POP_SolversMod_loop3)
  !if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop3 costs w zhuangy: ',stop_POP_SolversMod_loop3-start_POP_SolversMod_loop3,' (us);'
  !endif
 

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  initialize fields and scalars
!
!-----------------------------------------------------------------------

   call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                       POP_fieldKindScalar, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversPCG: error updating initial residual halo')
      return
   endif

   eta0 =1.0_POP_r8 
   numIterations = maxIterations
 
!-----------------------------------------------------------------------
!
!  iterate
!
!-----------------------------------------------------------------------

!!attention  not run in 5 hours 
  !call getlooptime(start_POP_SolversMod_loop4)
   iterationLoop: do m = 1, maxIterations

!-----------------------------------------------------------------------
!
!     calculate (PC)r 
!     diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         if (usePreconditioner) then
            call preconditioner(work1,R,iblock)
         else
            do j=1,ny
            do i=1,nx

               if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
                  work1(i,j,iblock) = R(i,j,iblock)/ &
                                      btropWgtCenter(i,j,iblock)
               else
                  work1(i,j,iblock) = 0.0_POP_r8
               endif
            end do
            end do
         endif

         do j=1,ny
         do i=1,nx
            work0(i,j,iblock) = R(i,j,iblock)*work1(i,j,iblock)
         end do
         end do
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     update conjugate direction vector s
!
!-----------------------------------------------------------------------

      if (usePreconditioner) then
         call POP_HaloUpdate(work1, POP_haloTropic, &
                             POP_gridHorzLocCenter, &
                             POP_fieldKindScalar, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error updating halo for preconditioner')
            return
         endif
      endif

      !*** (r,(PC)r)
      eta1 = POP_GlobalSum(work0, POP_distrbTropic, &
                           POP_gridHorzLocCenter,   &
                           errorCode, mMask = mMaskTropic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error in initial dot product')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            S(i,j,iblock) = work1(i,j,iblock) + &
                                S(i,j,iblock)*(eta1/eta0) 
         end do
         end do

!-----------------------------------------------------------------------
!
!        compute As
!
!-----------------------------------------------------------------------

         call btropOperator(Q,S,iblock)
         do j=1,ny
         do i=1,nx
            work0(i,j,iblock) = Q(i,j,iblock)*S(i,j,iblock)
         end do
         end do

      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     compute next solution and residual
!
!-----------------------------------------------------------------------

      call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error updating Q halo')
         return
      endif

      eta0 = eta1
      eta1 = eta0/POP_GlobalSum(work0, POP_distrbTropic,          &
                                POP_gridHorzLocCenter, errorCode, &
                                mMask = mMaskTropic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error in dot product')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            X(i,j,iblock) = X(i,j,iblock) + eta1*S(i,j,iblock)
            R(i,j,iblock) = R(i,j,iblock) - eta1*Q(i,j,iblock)
         end do
         end do

         if (mod(m,convergenceCheckFreq) == 0) then

            call btropOperator(R,X,iblock)
            do j=1,ny
            do i=1,nx
               R(i,j,iblock) = B(i,j,iblock) - R(i,j,iblock)
               work0(i,j,iblock) = R(i,j,iblock)*R(i,j,iblock)
            end do
            end do
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,convergenceCheckFreq) == 0) then

         call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                             POP_fieldKindScalar, errorCode)


         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error updating residual halo in convrg')
            return
         endif

         rr = POP_GlobalSum(work0, POP_distrbTropic, &
                            POP_gridHorzLocCenter,   &
                            errorCode, mMask = mMaskTropic)   ! (r,r)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error computing convergence dot prod')
            return
         endif

         if (rr < convergenceCriterion) then
            numIterations = m
            exit iterationLoop
         endif

      endif

   enddo iterationLoop
  !call getlooptime(stop_POP_SolversMod_loop4)
  !if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop4 costs w zhuangy: ',stop_POP_SolversMod_loop4-start_POP_SolversMod_loop4,' (us);'
  !endif
 

   rmsResidual = sqrt(rr*residualNorm)

   if (numIterations == maxIterations) then
      if (convergenceCriterion /= 0.0_POP_r8) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: solver not converged')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg

!***********************************************************************
!BOP
! !IROUTINE: PCSI
! !INTERFACE:

 subroutine PCSI(X,B,errorCode)  

! !DESCRIPTION:
!  This routine implements the Preconditioned Classical Stiefel Iteration 
!  (PCSI)  solver for solving the linear system $Ax=b$.
!  It uses the two extreme eigenvalues of $A$ instead of the norm of
!  residual. Thus, PCSI eliminates global reductions in each iteration.
!  The eigenvalues are estimated by routine PcsiLanczos.
!  PCSI supports all kinds of preconditioners supported by PCG/ChronGear, 
!  including diagonal and EVP preconditioning.
!
!  References:
!     Stiefel, E. L. (1958). Kernel polynomial in linear algebra and their
!        numerical applications, in: Further contributions to the 
!        determination of eigenvalues. NBS Applied Math. Ser., 49, 1-22.
!     Hu, Y., Huang, X., Wang, X., Fu, H., Xu, S., Ruan, H., Xue, W. and Yang, G. (2013). 
!        A scalable barotropic mode solver for the parallel ocean program. 
!        In Euro-Par 2013  Parallel Processing (pp. 739-750) Springer Berlin Heidelberg.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu et al., Tsinghua University

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!!!!   real(8)  :: st,st1,st2,st3,ed,ed1,ed2,ed3,st01,st02,st03
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,m,           &! local iteration counter
      nx, ny,          &! horizontal dimensions
      numBlocks,       &! number of local blocks
      iblock            ! local block counter


   real (POP_r8) ::          &
      csalpha,csbeta,csy,csomga,rr ! scalar inner product results
   
   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2), &
                                          size(X,dim=3)) :: &
      R,                 &! residual (b-Ax)
      S,                 &! conjugate direction vector
      Q,work0,work1,     &! various cg intermediate results
      A0R

    real (POP_r8), dimension(EvpXbs+2,EvpYbs+2,EvpXnb*EvpYnb):: &
      f               ! reshaped X  
 
    integer (POP_i4) :: &
      js,je,is,ie,lm,ln,l,ib,ierr               ! dummy counters

 
   errorCode = POP_Success
    
   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocks)
   nx = size(X,dim=1)
   ny = size(X,dim=2)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversPCSI: error retrieving local block count')
      return
   endif

   if (.not. usePreconditioner) then
   
      !--- diagonal preconditioner if preconditioner not specified
      !$OMP PARALLEL DO PRIVATE(iblock)
  !call getlooptime(start_POP_SolversMod_loop5)
      do iblock=1,numBlocks
         do j=1,ny
         do i=1,nx
            if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
               A0R(i,j,iblock) = 1.0_POP_r8 /btropWgtCenter(i,j,iblock)
            else
               A0R(i,j,iblock) = 0.0_POP_r8
            endif
         end do
         end do
      end do 
  !call getlooptime(stop_POP_SolversMod_loop5)
  !if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop5 costs w zhuangy: ',stop_POP_SolversMod_loop5-start_POP_SolversMod_loop5,' (us);'
!  endif
 
      !$OMP END PARALLEL DO
   endif


!-----------------------------------------------------------------------
!
!  step 1 : compute iteration parameters by eigenvalues
!  $\alpha =\frac{2}{\mu -\nu}$, $ \beta = \frac{\mu +\nu}{\mu -\nu}$,
!  $\gamma = \frac{\beta}{\alpha}$, $\omega_0 =\frac{ 2}{\gamma}$
!
!-----------------------------------------------------------------------
    
    csalpha = 2.0_POP_r8/(PcsiMaxEigs-PcsiMinEigs)
    csbeta = (PcsiMaxEigs+PcsiMinEigs)/(PcsiMaxEigs-PcsiMinEigs)
    csy = csbeta/csalpha
    csomga = 2.0_POP_r8/csy
    
!-----------------------------------------------------------------------
!
! step 2 : compute initial residual and initialize X
! $\textbf{r}_0 = \textbf{b}-\textbf{B}\textbf{x}_0$; 
! $\textbf{x}_1 =\textbf{x}_0 -\gamma^{-1}\textbf{M}^{-1}\textbf{r}_0$; 
! $\textbf{r}_1 =\textbf{b} -\textbf{B}\textbf{x}_1$; 
!
!-----------------------------------------------------------------------
!!   call getlooptime(st1)
      !$OMP PARALLEL DO PRIVATE(iblock,i,j)
 
    do iblock=1,numBlocks

        call btropOperator(S,X,iblock)
  !     call loop_6(iblock,B,S,R)
        do j=1,ny
        do i=1,nx  
          R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
        end do 
        end do 
    end do ! block loop
   

   !$OMP END PARALLEL DO
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)
 

        if (usePreconditioner) then
        do iblock=1,numBlocks
#if defined(ZYU_ATH_YY)

            work1(:,:,iblock) = 0.0_POP_r8  

            call loop_preconditioner(nx_block,ny_block,EvpYnb,EvpXnb,EvpYbidx(1),EvpXbidx(1), &
                 landIndx(1,1,iblock),R(1,1,iblock),f(1,1,1),work1(1,1,iblock),InvEvpCenterWgt(1,1,1,iblock), &
                 EvpCenterWgt(1,1,1,iblock),EvpNeWgt(1,1,1,iblock),EvpRinv(1,1,1,iblock),InvEvpNeWgt(1,1,1,iblock))


            do j=1,ny
            do i=1,nx
                Q(i,j,iblock) = (1.0_POP_r8/csy)*R(i,j,iblock)
           end do
           end do

#else

      
      call preconditioner(work1,R,iblock)

            do j=1,ny
            do i=1,nx
                R(i,j,iblock) = work1(i,j,iblock)
                Q(i,j,iblock) = (1.0_POP_r8/csy)*R(i,j,iblock)
           end do
           end do

#endif
         end do

        else
        do iblock=1,numBlocks
            do j=1,ny
            do i=1,nx
                R(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
                Q(i,j,iblock) = (1.0_POP_r8/csy)*R(i,j,iblock)
           end do
           end do


        end do
        endif
      !  call loop_7(iblock,csy,R,Q)
  

 
   !$OMP END PARALLEL DO    
    call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
      POP_fieldKindScalar, errorCode)

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, &
          'POP_SolversPCSI: error updating Q halo')
       return
    endif

    !$OMP PARALLEL DO PRIVATE(iblock,i,j)
 

    do iblock=1,numBlocks
     ! call loop_8_1(iblock)
      do j=1,ny
      do i=1,nx
       X(i,j,iblock) = X(i,j,iblock) + Q(i,j,iblock)
      end do
      end do
    
      call btropOperator(S,X,iblock)      
!      call loop_8_2(iblock,R,B,S)
      do j=1,ny
      do i=1,nx
          R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
      end do
      end do    
    end do ! block loop
 
   !$OMP END PARALLEL DO
  !! call getlooptime(st1)     !! 500 us
    
    call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
         POP_fieldKindScalar, errorCode)
  !!  call getlooptime(ed1)
    if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
        'POP_SolversPCSI: error updating initial residual halo')
        return
    endif
    
    numIterations = maxIterations
    

  !!   call getlooptime(ed1)
    !-----------------------------------------------------------------------
    !
    !  iterate
    !
    !-----------------------------------------------------------------------
    
!!!attention    cost lots of time  210000-140000 us



    iterationLoop: do m = 1, maxIterations


!-----------------------------------------------------------------------
!
! step 3 : update iteration parameter 
! $\omega_k = 1/(\gamma - \frac{1}{4\alpha^2}\omega_{k-1})$
!
!-----------------------------------------------------------------------
    csomga = 1.0_POP_r8/(csy-csomga/(4.0_POP_r8*csalpha*csalpha))

!-----------------------------------------------------------------------
!
! step 4 : preconditioning 
! $\textbf{r}'_{k-1} =\textbf{M}^{-1}\textbf{r}_{k-1}$
!
!-----------------------------------------------------------------------
!----PCSI-ATH---!



      if (usePreconditioner) then
      do iblock=1,numBlocks
#if defined(ZYU_ATH)


#if defined(ZYU_TIME1)   
         call getlooptime(st1) 
#endif
if (iblock.eq.1)  then
         call loop_preconditioner(nx_block,ny_block,EvpYnb,EvpXnb,EvpYbidx(1),EvpXbidx(1),                      &
              landIndx(1,1,iblock),R(1,1,iblock),InvEvpCenterWgt(1,1,1,iblock),      &
              EvpCenterWgt(1,1,1,iblock),EvpNeWgt(1,1,1,iblock),EvpRinv(1,1,1,iblock),InvEvpNeWgt(1,1,1,iblock),m)

#if defined(ZYU_TIME1)
        call getlooptime(ed1)
         if(POP_myTask ==0)then
              write(*,*)"preconditioner-ath cost=",ed1-st1,"us"
         endif
#endif

else
        call preconditioner(work1,R,iblock)


          do j=1,ny
          do i=1,nx
              R(i,j,iblock) = work1(i,j,iblock)
          end do
          end do



endif

         
#else   


          call preconditioner(work1,R,iblock) 


          do j=1,ny
          do i=1,nx
              R(i,j,iblock) = work1(i,j,iblock)
          end do 
          end do
#endif
      end do
      else
      do iblock=1,numBlocks
          do j=1,ny
          do i=1,nx
                R(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
          end do 
          end do
      end do
      endif 
 


    !$OMP END PARALLEL DO

    !!! assume that halosize = 2
    call POP_HaloUpdate_zy(R, POP_haloTropic, POP_gridHorzLocCenter, &
                           POP_fieldKindScalar, errorCode)



    if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
        'POP_SolversPCSI: error updating R halo')
        return
    endif


!-----------------------------------------------------------------------
!
! step 5 : compute X increment and update X and R
! $\Delta \textbf{x}_{k} =\omega_k\textbf{r}'_{k-1}+(\gamma \omega_k-1)\Delta \textbf{x}_{k-1}$; 
! $\textbf{x}_{k} =\textbf{x}_{k-1}+\Delta \textbf{x}_{k-1}$;  
! $\textbf{r}_{k} =b- \textbf{B}\textbf{x}_{k}$;
!
!-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(iblock)

  if ((mod(m,convergenceCheckFreq) == 0) .and.(m.ge.convergenceCheckStart)) then
    do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
        Q(i,j,iblock) = csomga*R(i,j,iblock)+ (csy*csomga-1.0_POP_r8)*Q(i,j,iblock)
        X(i,j,iblock) = X(i,j,iblock) + Q(i,j,iblock)
      end do
      end do

     call btropOperator(S,X,iblock)

      do j=1,ny
      do i=1,nx
          R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
      end do 
      end do

          do j=1,ny
          do i=1,nx
          WORK0(i,j,iblock) = R(i,j,iblock)*R(i,j,iblock)
          end do
          end do
    
   end do ! block loop
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !
    !     test for convergence
    !
    !-----------------------------------------------------------------------
    
    
        rr = POP_GlobalSum(work0, POP_distrbTropic, &
                           POP_gridHorzLocCenter,   &
                           errorCode, mMask = mMaskTropic)   ! (r,r)
    
        if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_SolversPCSI: error computing convergence dot prod')
            return
        endif
    
        if (rr < convergenceCriterion) then
            numIterations = m
            exit iterationLoop
        endif
    

   else

    do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
        Q(i,j,iblock) = csomga*R(i,j,iblock)+(csy*csomga-1.0_POP_r8)*Q(i,j,iblock)
        X(i,j,iblock) = X(i,j,iblock) + Q(i,j,iblock)
      end do
      end do

     call btropOperator(S,X,iblock)

      do j=1,ny
      do i=1,nx
          R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
      end do
      end do
  
    enddo

   endif
    
    enddo iterationLoop

#if defined(ZYU_ATH_YY)
    
         if(POP_myTask ==0)then
              write(*,*)"PCSI-m=",m
         endif
#endif

 
    rmsResidual = sqrt(rr*residualNorm)


!-----------------------------------------------------------------------
!EOC

 end subroutine PCSI
!***********************************************************************
!BOP
! !IROUTINE: ChronGear
! !INTERFACE:

 subroutine ChronGear(X,B,errorCode)

! !DESCRIPTION:
!  This routine implements the Chronopoulos-Gear conjugate-gradient 
!  solver with preconditioner for solving the linear system $Ax=b$.
!  It is a rearranged conjugate gradient solver that reduces the 
!  number of inner products per iteration from two to one (not 
!  counting convergence check). Both the operator $A$ and 
!  preconditioner are nine-point stencils. If no preconditioner has 
!  been supplied, a diagonal preconditioner is applied.  Convergence 
!  is checked every {\em ncheck} steps.
!
!
!  References:
!     Dongarra, J. and V. Eijkhout. LAPACK Working Note 159.
!        Finite-choice algorithm optimization in conjugate gradients.
!        Tech. Rep. ut-cs-03-502. Computer Science Department.
!        University of Tennessee, Knoxville. 2003.
!
!     D Azevedo, E.F., V.L. Eijkhout, and C.H. Romine. LAPACK Working
!        Note 56. Conjugate gradient algorithms with reduced
!        synchronization overhead on distributed memory multiprocessors.
!        Tech. Rep. CS-93-185.  Computer Science Department.
!        University of Tennessee, Knoxville. 1993.
!
! !REVISION HISTORY:
!  this routine implemented by Frank Bryan et al., NCAR

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,m,           &! local iteration counter
      nx, ny,          &! horizontal extents
      numBlocks,       &! number of local blocks
      iblock            ! local block counter

   real (POP_r8) :: & ! scalar results
      cgAlpha, cgBeta, cgSigma, cgDelta, cgRhoOld, cgRho, rr

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2),    &
                                          size(X,dim=3)) :: &
      R,                  &! residual (b-Ax)
      S,                  &! conjugate direction vector
      Q,Z,AZ,WORK0,       &! various cg intermediate results
      A0R                  ! diagonal preconditioner

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2),    &
                                        2,size(X,dim=3)) :: & 
      WORKN              ! WORK array 

   real (POP_r8), dimension(2) :: &
      sumN               ! global sum results for multiple arrays

!-----------------------------------------------------------------------
!
!  initialize some scalars
!
!-----------------------------------------------------------------------

   errorCode = POP_Success


   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocks)


   nx = size(X,dim=1)
   ny = size(X,dim=2)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error retrieving local block count')
      return
   endif

   cgRho = 1.0_POP_r8
   numIterations = maxIterations


!-----------------------------------------------------------------------
!
!  compute initial residual and initialize other arrays
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

 
 ! call getlooptime(start_POP_SolversMod_loop10)
   do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
         R    (i,j,iblock)   = 0.0_POP_r8
         S    (i,j,iblock)   = 0.0_POP_r8
         Z    (i,j,iblock)   = 0.0_POP_r8
         Q    (i,j,iblock)   = 0.0_POP_r8
         AZ   (i,j,iblock)   = 0.0_POP_r8
         WORK0(i,j,iblock)   = 0.0_POP_r8
         WORKN(i,j,1,iblock) = 0.0_POP_r8
         WORKN(i,j,2,iblock) = 0.0_POP_r8
      end do
      end do

      if (.not. usePreconditioner) then
     
         !--- diagonal preconditioner if preconditioner not specified
         do j=1,ny
         do i=1,nx
            if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
               A0R(i,j,iblock) = 1.0_POP_r8/btropWgtCenter(i,j,iblock)
            else
               A0R(i,j,iblock) = 0.0_POP_r8
            endif
         end do
         end do
      endif

     

      ! use S as a temp here for Ax

      call btropOperator(S,X,iblock)
      do j=1,ny
      do i=1,nx
         R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock) ! b-Ax
      end do
      end do
   end do ! block loop
 ! call getlooptime(stop_POP_SolversMod_loop10)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop10 costs w zhuangy: ',stop_POP_SolversMod_loop10-start_POP_SolversMod_loop10,' (us);'
 ! endif
 

   !$OMP END PARALLEL DO

   call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error updating initial residual halo')
      return
   endif

!-----------------------------------------------------------------------
!
!    take one pass of standard algorithm
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

 
 ! call getlooptime(start_POP_SolversMod_loop11)
   do iblock=1,numBlocks

      !---- calculate (PC)r store in Z
      if (usePreconditioner) then
         call preconditioner(Z,R,iblock)
      else   ! use diagonal preconditioner
         do j=1,ny
         do i=1,nx
            Z(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
         end do
         end do
      endif

   end do
 ! call getlooptime(stop_POP_SolversMod_loop11)
 ! if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop11 costs w zhuangy: ',stop_POP_SolversMod_loop11-start_POP_SolversMod_loop11,' (us);'
  !endif
 
   !$OMP END PARALLEL DO


   if (usePreconditioner) then
     call POP_HaloUpdate(Z, POP_haloTropic, POP_gridHorzLocCenter, &
                            POP_fieldKindScalar, errorCode)

     if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
           'POP_SolversChronGear: error updating Z halo')
        return
     endif
   endif


   !$OMP PARALLEL DO PRIVATE(iblock,i,j)
 
 ! call getlooptime(start_POP_SolversMod_loop12)
   do iblock=1,numBlocks
      !---- Compute intermediate result for dot product
      !---- update conjugate direction vector S
      do j=1,ny
      do i=1,nx
         WORKN(i,j,1,iblock) = R(i,j,iblock)*Z(i,j,iblock)
         S(i,j,iblock) =  Z(i,j,iblock)
      end do
      end do

      !---- compute Q = A * S
      call btropOperator(Q,S,iblock)

      !---- compute intermediate result for dot product
      do j=1,ny
      do i=1,nx
         WORKN(i,j,2,iblock) = S(i,j,iblock)*Q(i,j,iblock)
      end do
      end do

   end do
 ! call getlooptime(stop_POP_SolversMod_loop12)
  !if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop12 costs w zhuangy: ',stop_POP_SolversMod_loop12-start_POP_SolversMod_loop12,' (us);'
  !endif
 
   !$OMP END PARALLEL DO

   call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error updating Q halo')
      return
   endif


   !---- Form dot products
   sumN = POP_GlobalSum(WORKN, POP_distrbTropic,        &
                               POP_gridHorzLocCenter,   &
                               errorCode, mMask = mMaskTropic)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error in initial dot products')
      return
   endif

   cgRhoOld = sumN(1) !(r,PCr)
   cgSigma  = sumN(2) !(s,As)
   cgAlpha  = cgRhoOld/cgSigma

   !---- compute first solution and residual
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)
 
 ! call getlooptime(start_POP_SolversMod_loop13)
   do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
         X(i,j,iblock) = X(i,j,iblock) + cgAlpha*S(i,j,iblock)
         R(i,j,iblock) = R(i,j,iblock) - cgAlpha*Q(i,j,iblock)
      end do
      end do

   end do
 ! call getlooptime(stop_POP_SolversMod_loop13)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop13 costs w zhuangy: ',stop_POP_SolversMod_loop13-start_POP_SolversMod_loop13,' (us);'
 ! endif
 
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     iterate
!
!-----------------------------------------------------------------------

 
 ! call getlooptime(start_POP_SolversMod_loop14)
   iterationLoop: do m = 1, maxIterations

!-----------------------------------------------------------------------
!
!     calculate (PC)r and A*(Pc)r
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)
      do iblock=1,numBlocks

         if (usePreconditioner) then
            call preconditioner(Z,R,iblock)
         else
            do j=1,ny
            do i=1,nx
               Z(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
            end do
            end do
         endif
      end do
      !$OMP END PARALLEL DO


      ! Here we move the haloupdating forward in order to use EVP preconditioning
      ! Based on the assumption that halo_size = 2. Otherwise, its not correct
      call POP_HaloUpdate(Z, POP_haloTropic, POP_gridHorzLocCenter, &
                             POP_fieldKindScalar, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: error updating Z halo')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)
      do iblock=1,numBlocks

        call btropOperator(AZ,Z,iblock)

         !--- intermediate results for inner products

         do j=1,ny
         do i=1,nx
            WORKN(i,j,1,iblock) =  R(i,j,iblock)*Z(i,j,iblock)
            WORKN(i,j,2,iblock) = AZ(i,j,iblock)*Z(i,j,iblock)
         end do
         end do

      end do
      !$OMP END PARALLEL DO


      sumN = POP_GlobalSum(WORKN, POP_distrbTropic,        &
                                  POP_gridHorzLocCenter,   &
                                  errorCode, mMask = mMaskTropic)

;
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: error in dot products')
         return
      endif

      cgRho    = sumN(1)     ! (r,(PC)r)
      cgDelta  = sumN(2)     ! (A (PC)r,(PC)r)
      cgBeta   = cgRho/cgRhoOld
      cgSigma  = cgDelta - (cgBeta**2)*cgSigma
      cgAlpha  = cgRho/cgSigma
      cgRhoOld = cgRho

!-----------------------------------------------------------------------
!
!     compute S and Q
!     compute next solution and residual
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock, i,j)
      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            S(i,j,iblock) =  Z(i,j,iblock) + cgBeta *S(i,j,iblock)
            Q(i,j,iblock) = AZ(i,j,iblock) + cgBeta *Q(i,j,iblock)
            X(i,j,iblock) =  X(i,j,iblock) + cgAlpha*S(i,j,iblock)
            R(i,j,iblock) =  R(i,j,iblock) - cgAlpha*Q(i,j,iblock)
         end do
         end do

         !--- recompute residual as b-Ax for convergence check
         if (mod(m,convergenceCheckFreq) == 0) then

            !--- Reset residual using r = b - Ax
            !--- (r,r) for norm of residual
            call btropOperator(Z,X,iblock)
            do j=1,ny
            do i=1,nx
               R(i,j,iblock) = B(i,j,iblock) - Z(i,j,iblock)
               WORK0(i,j,iblock) = R(i,j,iblock)*R(i,j,iblock)
            end do
            end do

         endif
      end do
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     test for convergence if it is time
!
!-----------------------------------------------------------------------

      if (mod(m,convergenceCheckFreq) == 0) then

         !--- update ghost cells for next iteration
         call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                                POP_fieldKindScalar, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversChronGear: error updating residual halo in convrg')
            return
         endif

         !--- residual norm for convergence 
         rr = POP_GlobalSum(work0, POP_distrbTropic,        &! (r,r)
                                   POP_gridHorzLocCenter,   &
                                   errorCode, mMask = mMaskTropic)


         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversChronGear: error computing convergence dot prod')
            return
         endif

         if (rr < convergenceCriterion) then
            numIterations = m
            exit iterationLoop
         endif

      endif

   end do iterationLoop
 ! call getlooptime(stop_POP_SolversMod_loop14)
 ! if (POP_myTask == POP_masterTask) then
  !   WRITE (*,*) 'POP_SolversMod_loop14 costs w zhuangy: ',stop_POP_SolversMod_loop14-start_POP_SolversMod_loop14,' (us);'
 ! endif
 

   rmsResidual = sqrt(rr*residualNorm)

   if (numIterations == maxIterations) then
      if (convergenceCriterion /= 0.0_POP_r8) then

         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: solver not converged')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ChronGear

!***********************************************************************
!BOP
! !IROUTINE: preconditioner
! !INTERFACE:

 subroutine preconditioner(PX,X,bid)

! !DESCRIPTION:
!  This function applies a precomputed preconditioner as a nine-point
!  stencil operator.
!
! !REVISION HISTORY:
!  add EVP preconditioning interface by Yong Hu et al., Tsinghua University

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: & 
      X                     ! array to be operated on 

   integer (POP_i4), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      PX                  ! nine point operator result

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (POP_r8), dimension(EvpXbs+2,EvpYbs+2,EvpXnb*EvpYnb):: &
      f               ! reshaped X

   integer (POP_i4) :: &
      i,j,js,je,is,ie,lm,ln,l,ib               ! dummy counters

!-----------------------------------------------------------------------

   PX(:,:,bid) = 0.0_POP_r8
 
!-----------------------------------------------------------------------
!
!  EVP preconditioning 
!
!-----------------------------------------------------------------------

  if (trim(preconditionerChoice) == precondChoiceEvp) then
       f = 0.0_POP_r8 
 
!! call getlooptime(start_POP_SolversMod_loop16)

!!  call loop_pop_solvers16(EvpYnb,EvpXnb,EvpCenterWgt,EvpNeWgt,InvEvpNeWgt,EvpRinv,PX,f)
!    call Evp_loop(EvpYnb,EvpXnb,bid,    EvpYbidx,EvpXbidx,landIndx,X,f,
!                  PX,InvEvpCenterWgt,EvpCenterWgt,EvpNeWgt,EvpRinv)
       do j = 1, EvpYnb
         js = EvpYbidx(j)
         je = EvpYbidx(j+1) +1
         lm = (je-js) +1
         do i = 1, EvpXnb
           is = EvpXbidx(i) 
           ie = EvpXbidx(i+1) +1
           ln = (ie-is) +1
           l  = ln + lm -5
           ib = (j-1)*EvpXnb+i
           f(2:ln-1,2:lm-1,ib) = X(is+1:ie-1,js+1:je-1,bid)
           if (landIndx(i,j,bid) == 1 ) then 

             ! diagonal preconditioning for blocks containing land potins
             PX(is+1:ie-1,js+1:je-1,bid) = & 
             f(2:ln-1,2:lm-1,ib)*InvEvpCenterWgt(2:ln-1,2:lm-1,ib,bid)

           else if (landIndx(i,j,bid) == 0 ) then 
             !!call evp()  
             ! EVP solver on sub-blocks
             call ExplicitEvp(EvpCenterWgt(1:ln,1:lm,ib,bid),&
                  EvpNeWgt(1:ln,1:lm,ib,bid),InvEvpNeWgt(1:ln,1:lm,ib,bid),&
                  EvpRinv(1:l,1:l,ib,bid),PX(is+1:ie-1,js+1:je-1,bid),f(1:ln,1:lm,ib),ln,lm)

         !  else 
         !    write(POP_stdout,'(a35,3I5.3)') 'EVP Error: unpreconditioned block ',&
         !                                    i,j,landIndx(i,j,bid)
           endif 

         end do 
       end do 
!!  call getlooptime(stop_POP_SolversMod_loop16)
!  if (POP_myTask == POP_masterTask) then
   !  WRITE (*,*) 'POP_SolversMod_loop16:',EvpYnb,EvpXnb
   !  WRITE (*,*) 'POP_SolversMod_loop16 costs w zhuangy: ',stop_POP_SolversMod_loop16-start_POP_SolversMod_loop16,' (us);',EvpYnb,EvpXnb
!  endif
 
   else  if (trim(preconditionerChoice) == precondChoiceFile) then
!! not run now! 
 ! call getlooptime(start_POP_SolversMod_loop15)
       do j=2,size(X,dim=2)-1
       do i=2,size(X,dim=1)-1
          PX(i,j,bid) = precondNE    (i,j,bid)*X(i+1,j+1,bid) + &
                        precondNW    (i,j,bid)*X(i-1,j+1,bid) + &
                        precondSE    (i,j,bid)*X(i+1,j-1,bid) + &
                        precondSW    (i,j,bid)*X(i-1,j-1,bid) + &
                        precondNorth (i,j,bid)*X(i  ,j+1,bid) + &
                        precondSouth (i,j,bid)*X(i  ,j-1,bid) + &
                        precondEast  (i,j,bid)*X(i+1,j  ,bid) + &
                        precondWest  (i,j,bid)*X(i-1,j  ,bid) + &
                        precondCenter(i,j,bid)*X(i  ,j  ,bid)
       end do
       end do
 ! call getlooptime(stop_POP_SolversMod_loop15)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop15 costs w zhuangy: ',stop_POP_SolversMod_loop15-start_POP_SolversMod_loop15,' (us);'
 ! endif
   end if
!-----------------------------------------------------------------------
!EOC

 end subroutine preconditioner

!***********************************************************************
!BOP
! !IROUTINE: btropOperator
! !INTERFACE:

 subroutine btropOperator(AX,X,bid)

! !DESCRIPTION:
!  This routine applies the nine-point stencil operator for the
!  barotropic solver.  It takes advantage of some 9pt weights being 
!  shifted versions of others.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: & 
      X                  ! array to be operated on 

   integer (POP_i4), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      AX            ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j                ! dummy counters
#if defined(ZYU_TIME) 
    real(8) :: st,ed,st1,ed1;
#endif
#if defined(ZYU_CHECK)
    real(8) :: maxdiff=0.0_POP_r8 ,diff=0.0_POP_r8
    real(8),dimension(128,128,4) :: check_AX
#endif
!-----------------------------------------------------------------------

#if defined(ZYU_ATH_SMALL)

#if defined(ZYU_TIME1)
   call getlooptime(st)
#endif
  AX(:,1,bid) = 0.0_POP_r8
  AX(:,POP_nyblock,bid) = 0.0_POP_r8
  call loop_pop_solvers_operator(btropWgtCenter(1,1,bid),btropWgtNorth(1,1,bid),btropWgtEast(1,1,bid) &
                      ,btropWgtNE(1,1,bid),AX(1,1,bid),X(1,1,bid),POP_nxblock,POP_nyblock)

#else

#if defined(ZYU_TIME1)
   call getlooptime(ed)
#endif
#if defined(ZYU_CHECK)
!!  call  myrpcc(ed)
   do j=1,size(X,dim=1)-1
   do i=1,size(X,dim=1)-1
         check_AX(i,j,bid)  = AX(i,j,bid)
   end do
   end do
#endif
#if defined(ZYU_TIME1)
  call getlooptime(st1)
#endif
    AX(:,:,bid) = 0.0_POP_r8
    do j=2,size(X,dim=2)-1
    do i=2,size(X,dim=1)-1
      AX(i,j,bid) = btropWgtCenter(i  ,j  ,bid)*X(i  ,j  ,bid) + &
                    btropWgtNorth (i  ,j  ,bid)*X(i  ,j+1,bid) + &
                    btropWgtNorth (i  ,j-1,bid)*X(i  ,j-1,bid) + &
                    btropWgtEast  (i  ,j  ,bid)*X(i+1,j  ,bid) + &
                    btropWgtEast  (i-1,j  ,bid)*X(i-1,j  ,bid) + &
                    btropWgtNE    (i  ,j  ,bid)*X(i+1,j+1,bid) + &
                    btropWgtNE    (i  ,j-1,bid)*X(i+1,j-1,bid) + &
                    btropWgtNE    (i-1,j  ,bid)*X(i-1,j+1,bid) + &
                    btropWgtNE    (i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do
#if defined(ZYU_TIME1)
  call getlooptime(ed1)
  if (POP_myTask == POP_masterTask) then   
         WRITE (*,*) 'ath-cost= ',ed-st,' (us) org-cost ',ed1-st1,' (us)'
  endif
#endif 

#if defined(ZYU_CHECK)
   do j=1,size(X,dim=1)-1
   do i=1,size(X,dim=1)-1
           diff   = AX(i,j,bid) -  check_AX(i,j,bid)
           maxdiff = MAX(maxdiff , ABS(diff))       
   end do
   end do 
  if (POP_myTask ==  POP_masterTask) then
       if (POP_myTask.lt.3) then
         WRITE (*,*) '####diff:',maxdiff,'(2,2)=',AX(2,2,bid)
       endif
  endif
#endif

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine btropOperator


 subroutine EvpPre(cc,ne,evpcc,evpne,rinv,landIndx,n,m,nn,mm,ndi,mdi,nb,mb,bid)

! !DESCRIPTION:
! !prepare EVP preconditioning and save reshaped array for efficiency 
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University
! !INPUT/OUTPUT PARAMETERS:
   integer(POP_i4),intent(in) :: n,m  ! total block size
   integer(POP_i4),intent(in) :: nn,mm  ! small block ideal size
   real(POP_r8),dimension(n,m),intent(in) :: cc,ne
   integer(POP_i4),intent(in) :: nb, mb  ! blocks on x and y direction 
   integer(POP_i4),intent(in) :: ndi(nb+1),mdi(mb+1) ! bound index 

! !INPUT/OUTPUT PARAMETERS:
   real(POP_r8),dimension(nn+2,mm+2,nb*mb),intent(inout) :: evpcc,evpne
   real(POP_r8),dimension(nn+mm-1,nn+mm-1,nb*mb),intent(inout):: rinv
   integer(POP_i4),dimension(nb,mb),intent(inout):: landIndx
   integer(POP_i4),intent(in):: bid  ! bartropic blocks id

   ! LOCAL  VARIABLES
   integer (POP_i4) :: &
      is,js,ie,je,          &! start and end index for EVP sub block
      lm,ln,l,              &! EVP sub block size
      i,j,ib,               &! index of sub block
      iland                 ! index of sub block

   type (block) :: &
      this_block   ! block information for this block

   this_block = get_block(blocks_tropic(bid),bid)
 
!  call getlooptime(start_POP_SolversMod_loop18)
   do j = 1, mb
     js = mdi(j)-1
     je = mdi(j+1) 
     lm = (je-js) +1
     do i = 1, nb
       is = ndi(i)-1
       ie = ndi(i+1)
       ln = (ie-is) +1
       ib = (j-1)*nb+i
       ! Reshape coefficents into EVP sub-blocks  for efficiency
       evpcc(1:ln,1:lm,ib) =  cc(is:ie,js:je)
       evpne(1:ln,1:lm,ib) =  ne(is:ie,js:je)

       iland = 0 
       if (minval(abs(ne(is+1:ie-1,js+1:je-1))) == 0.0_POP_r8) then 
         iland = 1
       endif 
       if (is+2 <this_block%ib .or. ie-1 > this_block%ie) then
         iland = 2
       endif 
       if (js+2 <this_block%jb .or. je-1 > this_block%je) then
         iland = 3
       endif 

       if( iland > 0 ) then 
         ! mark land sub-blocks
         landIndx(i,j) = 1
         rinv(:,:,ib) = 0.0_POP_r8
       else 
         ! EVP preprocessing ocean sub-blocks
         landIndx(i,j) = 0
         l  = ln + lm -5
         call ExplicitBlockEvpPre(cc(is:ie,js:je),ne(is:ie,js:je), &
                                  rinv(1:l,1:l,ib),ln,lm)
       endif 
     end do 
   end do 
!  call getlooptime(stop_POP_SolversMod_loop18)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop18 costs w zhuangy: ',stop_POP_SolversMod_loop18-start_POP_SolversMod_loop18,' (us);'
!  endif
 
   write(POP_stdout,'(a8,I5,a15,I5,a15,I5)') 'PROC',POP_myTask, & 
                    ' EVP blocks :',nb*mb, ' land blocks :', sum(landIndx)

  end subroutine EvpPre
  
  subroutine ExplicitBlockEvpPre(cc,ne,rinv,n,m)
! !DESCRIPTION:
!  This routine implements the preprocessing of explicit EVP method for
!  solve a nine point elliptic equation on sub blocks.
!  Four coefficients related to north, south, east and west are small compared 
!  with the other five coefficients, ignored here to reduce computation.
!
!  References:
!     Roache, P. J. (1995). Elliptic marching methods and domain
!        decomposition (Vol.5). CRC press.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University
 
 ! !INPUT PARAMETERS:
   integer,intent(in) :: n,m                      !sub-block size [n,m]
   real(POP_r8),dimension(n,m),intent(in) :: cc   !center coefficients
   real(POP_r8),dimension(n,m),intent(in) :: ne   !four corners coefficients

 ! !OUTPUT PARAMETERS:
   real(POP_r8),dimension(n+m-5,n+m-5),intent(inout) :: rinv !preconditioning matrix

 ! LOCAL VARIABLES
   integer :: i,j,k,ii
   integer :: nm                          ! length of error and final vectors
   integer (POP_i4)::  errorCode          ! returned error code
   real(POP_r8),dimension(n,m) :: y       ! temporary array

   real(POP_r8),dimension(n +m -5,n +m -5) :: rin,WORK
   real(POP_r8):: maxvalr
   real(POP_r8),dimension(n +m -5) :: r

   nm = n +m -5
   y(:,:) = 0.0_POP_r8

   ! five points marching from initial error vectors (west)
 
!  call getlooptime(start_POP_SolversMod_oop21)
   do ii = 1,m-2
     y(2,m-ii) = 1.0_POP_r8
     do j = 2, m-1
       do i = 2, n-1
         y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
                -  ne(i,j-1)   * y(i+1,j-1)  &
                -  ne(i-1,j)   * y(i-1,j+1)  & 
                -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
       end do
     end do

     ! get the final error vector f
     do i = 1,n-2
       rin(ii,i) = -y(i+2,m)
     end do 

     do j = 1,m-3
       rin(ii,n-2+j) = -y(n,m-j)
     end do 

     y(2,m-ii) = 0.0_POP_r8
   end do 
!  call getlooptime(stop_POP_SolversMod_loop19)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop19 costs w zhuangy: ',stop_POP_SolversMod_loop19-start_POP_SolversMod_loop19,' (us);'
!  endif
 
   
   ! five points marching from inital error vectors (south)
 
!  call getlooptime(start_POP_SolversMod_loop20)
   do ii = 1,n-3
     y(ii+2,2) = 1.0_POP_r8
     do j = 2, m-1
       do i = 2, n-1
         y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
                -  ne(i,j-1)   * y(i+1,j-1)  &
                -  ne(i-1,j)   * y(i-1,j+1)  & 
                -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
       end do
     end do
     
     ! get the final error vector f
     do i = 1,n-2
       rin(m-2+ii,i) = -y(i+2,m)
     end do 

     do j = 1,m-3
       rin(m-2+ii,n-2+j) = -y(n,m-j)
     end do 

     y(ii+2,2) = 0.0_POP_r8
   end do 
!  call getlooptime(stop_POP_SolversMod_loop20)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop20 costs w zhuangy: ',stop_POP_SolversMod_loop20-start_POP_SolversMod_loop20,' (us);'
!  endif
 

   WORK(:,:) = rin(:,:)
   call inverse(WORK,rinv,nm)

   !! check pre rinv 
   WORK(:,:) = 0.0_POP_r8
 
!  call getlooptime(start_POP_SolversMod_loop21)
   do j = 1,nm
     do i = 1,nm
       do k = 1,nm
         WORK(i,j) = WORK(i,j) + rinv(i,k)*rin(k,j)
       end do 
       if (i == j ) then 
         WORK(i,j) = WORK(i,j) -1.0
       endif 
     end do 
   end do 
!  call getlooptime(stop_POP_SolversMod_loop21)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop21 costs w zhuangy: ',stop_POP_SolversMod_loop21-start_POP_SolversMod_loop21,' (us);'
!  endif
 
   maxvalr = maxval(abs(WORK))
   !write(POP_stdout,*) "MAXVAL(RINV*RIN) = ", maxvalr
   if (maxvalr > 1.0e-8 ) then 
     write(POP_stdout,*) 'maxvalr ', maxvalr
     call POP_ErrorSet(errorCode, &
          'POP_EXPLICITPRE: error in computing the inverse, error > 1.0e-8 ;&
           Check EVP sub-block size!')
     return
   endif 

  end subroutine 

  subroutine ExplicitEvp(cc,ne,ine,rinv,tu,f,n,m)
! !DESCRIPTION:
!  This routine implements the EVP method to explicitly solve 
!  a nine point elliptic equation on sub blocks.
!  Four coefficients related to north, south, east and west are 
!  small compared with the rest five coefficients, thus ignored 
!  here to reduce computation
!
!  References:
!     Roache, P. J. (1995). Elliptic marching methods and domain 
!        decomposition (Vol.5). CRC press.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University
 
 ! !INPUT PARAMETERS:
   integer(POP_i4),intent(in) :: n,m                    !sub-block size [n,m]
   real(POP_r8),dimension(n,m),intent(in) :: cc,&! center weight
                                             ne,&! northeast weight
                                             ine ! inverse of ne
   real(POP_r8),dimension(n,m),intent(in) :: f   ! residual 
   real(POP_r8),dimension(n+m-5,n+m-5),intent(in) :: rinv
                                                 ! preconditioning matrx

 ! !OUTPUT PARAMETERS:
   real(POP_r8),dimension(n-2,m-2),intent(inout) :: tu ! solution
 
 
   !LOCAL VARIABLES
   integer(POP_i4) :: i,j,k       ! local counters
   integer(POP_i4) :: nm          ! length of initial and final error
   real(POP_r8),dimension(n,m) :: y ! temporary array
   real(POP_r8),dimension(n+m-5) :: r !final error vector


 !  call getlooptime(start_POP_SolversMod_loop25)
 
   nm = n+m-5
   y(:,:) = 0.0_POP_r8
   y(2:n-1,2:m-1) = tu(:,:) 
   
   ! marching from intial error vectors (west and south)
   do j = 2, m-1
     do i = 2, n-1
         y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
                -  ne(i,j-1)   * y(i+1,j-1)  &
                -  ne(i-1,j)   * y(i-1,j+1)  & 
                -  ne(i-1,j-1) * y(i-1,j-1) ) *ine(i,j) 
     end do
   end do

   ! get final error vectors (east and north)
   r(1:n-2) = y(3:n,m)
   r(n-1:n+m-5) = y(n,m-1:3:-1) 
 
   !compute intial error 
   do j = 1,m-2
       do k = 1,nm
         y(2,m-j)  = y(2,m-j) + rinv(k,j)*r(k)
       end do 
   end do 
   do i = 1,n-3
       do k = 1,nm
         y(i+2,2)  = y(i+2,2) + rinv(k,m-2+i)*r(k)
       end do 
   end do 
   !marching again based on the adjusted intial vector
 
!  call getlooptime(start_POP_SolversMod_loop25)
   do j = 2, m-2
     do i = 2, n-2
         y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
                -  ne(i,j-1)   * y(i+1,j-1)  &
                -  ne(i-1,j)   * y(i-1,j+1)  & 
                -  ne(i-1,j-1) * y(i-1,j-1) ) *ine(i,j) 
       end do
   end do
   !copy results 
   tu(1:n-2,1:m-2) = y(2:n-1,2:m-1) 

!   call getlooptime(stop_POP_SolversMod_loop25)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop25 costs w zhuangy:',stop_POP_SolversMod_loop25-start_POP_SolversMod_loop25,' (us);'
!  endif

   
  end subroutine 
 
 !***********************************************************************
   subroutine PcsiLanczos(nx, ny, nb, maxlanczosstep, u, v, errorCode)
 
! !DESCRIPTION:
!  This routine implements the Lanczos based method to solve  
!  the maximum and minimum eigenvalues of the coefficient matrix $A$.
!  Or $AM$ if a preconditioner $M^{-1}$ is provided.
!
!  References:
!     Hu, Y., Huang, X., Wang, X., Fu, H., Xu, S., Ruan, H., ... & Yang, G. (2013). A
!        scalable barotropic mode solver for the parallel ocean program. In Euro-Par 2013
!        Parallel Processing (pp. 739-750). Springer Berlin Heidelberg.
!     Paige, C. C. (1980). Accuracy and effectiveness of the Lanczos algorithm for the
!        symmetric eigenproblem. Linear algebra and its applications, 34, 235-258.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University
 
 ! !INPUT PARAMETERS:
 
   integer (POP_i4),intent(in) :: &
      nx, ny, nb         ! block size and numblocktropic
   integer (POP_i4),intent(in) :: &
      maxlanczosstep     ! max lanczos steps
 
   real (POP_r8), intent(inout) :: &
      u, v             !  max and min eigenvalues
 
   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code
 
 ! LOCAL VARIABLES
   real (POP_r8), dimension(nx, ny,nb) :: &
      R,                  &! r 
      S,                  &! preconditioned r : Mr
      Q,                  &! normalized r : r/||r||_A
      Q1,                 &! previous step q
      P,                  &! preconditioned q : Mq
      A0R,                &! inverse of diagonal of A
      WORK,WORK1           ! temporary vectors
 
   real (POP_r8), dimension(maxlanczosstep) :: &
      vcsa, vcsb                ! diagonal and off-diagonal 
                                ! elements of tridiagonal matrix $T$
 
   real (POP_r8), dimension(:), allocatable :: &
      mcsa, mcsb                ! temporary vectors saving 

 
   integer (POP_i4) :: &
      i, j, m,                 &! local iteration counter
      info,istat,              &! local flags 
      iblock                    ! local block counter
 
   real (POP_r8) :: &
        csa,                   &! inner product alpha = (p,r)
        csb,                   &! inner product beta = (r,s)
        csc                     ! vector magnitude  ||r0||_A

   real (POP_r8) :: & 
        mineig   

   real(POP_r8) :: const1,const2

   integer :: jend,iend
     
   errorCode = POP_Success
   
   ! compute diagonal inverse for diagonal preconditioning
!!attention
 !  const1=0.0_POP_r8
 !  const2=1.0_POP_r8
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) '#### test ####' !! nb=2
 ! endif 
 
 !  call getlooptime(start26)
!!   call get_timer(timer26,'LOOP26',nblocks_clinic, distrb_clinic%nprocs) 
!!   call timer_start(timer26)
 !  call loop_pop_solvers26(btropWgtCenter,A0R,const1,const2,nb,R,Q,Q1) 
!!   call timer_stop(timer26)
 !  call getlooptime(end26) 
 !  call getlooptime(start_POP_SolversMod_loop26)

   do iblock=1,nb
      do j=1,ny
      do i=1,nx
        if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
           A0R(i,j,iblock) = 1.0_POP_r8 /btropWgtCenter(i,j,iblock)
        else
           A0R(i,j,iblock) = 0.0_POP_r8
        endif
      end do
      end do
   end do ! block loop
 ! call getlooptime(stop_POP_SolversMod_loop26)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop26 costs w zhuangy:',stop_POP_SolversMod_loop26-start_POP_SolversMod_loop26,' (us);','athread:',end26-start26
 ! endif
! 
 
   ! set random vector R, intialize temporary vectors P, Q
!!! loop27+loop26 
!  call getlooptime(start_POP_SolversMod_loop27)
   do iblock=1,nb
      do j=1,ny
      do i=1,nx
        R(i,j,iblock) = 1.0_POP_r8
        Q(i,j,iblock) = 0.0_POP_r8
        Q1(i,j,iblock) = 0.0_POP_r8
      end do
      end do
   end do ! block loop
 ! call getlooptime(stop_POP_SolversMod_loop27)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop27 costs w zhuangy: ',stop_POP_SolversMod_loop27-start_POP_SolversMod_loop27,' (us);'
!  endif
 
   
   ! preconditioning s = M^{-1}r 
 
!  call getlooptime(start_POP_SolversMod_loop28)
   if (usePreconditioner) then
   do iblock=1,nb
    
        call preconditioner(S,R,iblock) 
        do j=1,ny
        do i=1,nx
            WORK(i,j,iblock) = S(i,j,iblock)*R(i,j,iblock)
        end do
        end do
   end do ! block loop
    else
    do iblock=1,nb   
            do j=1,ny
            do i=1,nx
                 S(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
                 WORK(i,j,iblock) = S(i,j,iblock)*R(i,j,iblock)
            end do
            end do
   end do
   endif




!  call getlooptime(stop_POP_SolversMod_loop28)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop28 costs w zhuangy: ',stop_POP_SolversMod_loop28-start_POP_SolversMod_loop28,' (us);'
!  endif
 
 
   ! (R,S)  $A$ and $M$  is negative defined 
   csc =- POP_GlobalSum(WORK, POP_distrbTropic, &
              POP_gridHorzLocCenter,   &
              errorCode, mMask = mMaskTropic)

   ! q = R/||R||_{AM}  normalize R based on measurement of $AM^{-1}$
   if (csc .gt. 0.0_POP_r8) then 
 
 ! call getlooptime(start_POP_SolversMod_loop29)
      do iblock=1,nb
      do j=1,ny
      do i=1,nx
            Q(i,j,iblock) = (1/sqrt(csc))*R(i,j,iblock)
      end do 
      end do 
      end do ! block loop
!  call getlooptime(stop_POP_SolversMod_loop29)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop29 costs w zhuangy: ',stop_POP_SolversMod_loop29-start_POP_SolversMod_loop29,' (us);'
!  endif
 
   else 
      ! either random vector R == 0 or A is sigular 
      call POP_ErrorSet(errorCode, &
         'POP_PcsiLanczos: error estimating in lanczos: b == 0 !!!')
      return
   endif 
 
   call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)
 
   csb= 0.0_POP_r8
   u= 0.0_POP_r8
   v= 0.0_POP_r8
   mineig= 1.0_POP_r8
 
!-----------------------------------------------------------------------
!
!  lanczos iterations
!
!-----------------------------------------------------------------------
 
!  call getlooptime(start_POP_SolversMod_loop30)
   lanczos_iter: do m = 1, maxlanczosstep

   !  preconditioning q
!!  call getlooptime(start_POP_SolversMod17)
     if (usePreconditioner) then
     do iblock=1,nb
            call preconditioner(P,Q,iblock)
     end do    
     else
     do iblock=1,nb
            do j=1,ny
            do i=1,nx
                P(i,j,iblock) = Q(i,j,iblock)*A0R(i,j,iblock)
            end do 
            end do 
     end do ! block loop
     end if
!! call getlooptime(stop_POP_SolversMod17)

    call POP_HaloUpdate(P, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)

   ! r = Mp - beta_{j-1} q_{j-1}

!  call getlooptime(start_POP_SolversMod_loop30_2)
   do iblock=1,nb
      call btropOperator(WORK1,P,iblock) 
        
      !  jend=size(P,dim=2)-1
      !  iend=size(P,dim=1)-1  
      !  WORK1(:,1,iblock) = 0.0_POP_r8
      !  WORK1(:,jend+1,iblock) = 0.0_POP_r8
     
       !! call getlooptime(start_POP_SolversMod_loop30_2)
      !  call loop_pop_solvers30_2(jend,iend,iblock,btropWgtCenter,btropWgtNorth,btropWgtEast &
      !                ,btropWgtNE,WORK1,P,R,csb,Q1,WORK)
      !!  call getlooptime(stop_POP_SolversMod_loop30_2)
  
    do j=1,ny
      do i=1,nx
         R(i,j,iblock) =  WORK1(i,j,iblock)-csb*Q1(i,j,iblock)
          WORK(i,j,iblock) = P(i,j,iblock)*R(i,j,iblock)
      end do
      end do
   end do ! block loop
!  call getlooptime(stop_POP_SolversMod_loop30_2)
    csa = -POP_GlobalSum(WORK, POP_distrbTropic, &
           POP_gridHorzLocCenter,   &
           errorCode, mMask = mMaskTropic)   ! (r,r)
  
!!attention

  
 ! call getlooptime(start_POP_SolversMod_loop30_3)
!   call loop_pop_solvers30_3(R,Q,csa) 
 !  call getlooptime(stop_POP_SolversMod_loop29)
   ! r = r - alpha*q 
 !  call getlooptime(start_POP_SolversMod_loop30_3)
   do iblock=1,nb
      do j=1,ny
     do i=1,nx
        R(i,j,iblock) = R(i,j,iblock)-csa*Q(i,j,iblock)
      end do
      end do 
   end do 
 ! call getlooptime(stop_POP_SolversMod_loop30_3)
   call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)
   ! preconditioning R 

 ! call getlooptime(start_POP_SolversMod_loop30_4)
   do iblock=1,nb
        if (usePreconditioner) then
            call preconditioner(S,R,iblock)
        else
            do j=1,ny
            do i=1,nx
                S(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
            end do 
            end do 
        endif 
   !     call loop_pop_solvers30_4(iblock,S,R,WORK)
        do j=1,ny
        do i=1,nx
          WORK(i,j,iblock) = S(i,j,iblock)*R(i,j,iblock)
        end do
        end do 
   end do ! block loop
!  call getlooptime(stop_POP_SolversMod_loop30_4)
   csc = -POP_GlobalSum(WORK, POP_distrbTropic, &
                POP_gridHorzLocCenter,   &
                errorCode, mMask = mMaskTropic)   ! (r,s)
   csb = sqrt(csc)
 
   vcsa(m) = csa
   vcsb(m) = csb
   
   ! Gershgorin circle theorem to estimate the largest eigenvaule 
   if (m ==1 ) then
    u = vcsa(1) + vcsb(1)
   else 
    u = max(u, vcsa(m) + vcsb(m) + vcsb(m-1))
   endif

   ! normalize r  : q = r/||r||_A

    
   if (csb /= 0.0_POP_r8) then 

 !    call getlooptime(start_POP_SolversMod_loop30_5)
 !     call loop_pop_solvers30_5(Q1,Q,R,csb)
       
   !   call getlooptime(start_POP_SolversMod_loop30_5)
      do iblock=1,nb
         do j=1,ny
         do i=1,nx
          Q1(i,j,iblock) = Q(i,j,iblock)
          Q(i,j,iblock) = (1/csb)*R(i,j,iblock)
          end do 
          end do 
       end do ! block loop
   !
 !      call getlooptime(stop_POP_SolversMod_loop30_5)
   else 
      call POP_ErrorSet(errorCode, &
         'POP_PcsiLanczos: error estimating in lanczos: b == 0 !!!')
      return
   endif 
       
   
!-----------------------------------------------------------------------
!
!  compute eigenvalues of the tridiagonal matrix T every 10 steps
!
!-----------------------------------------------------------------------
    
   if ((mod(m,10) == 0) .or. (m ==maxlanczosstep)) then
      if (POP_myTask == POP_masterTask ) then

        allocate(mcsa(m), mcsb(m), stat = istat)
        if (istat > 0) then
           call POP_ErrorSet(errorCode, &
              'POP_PcsiLanczos: error allocating Lanczos Tridiagonal')
           return
        endif

        do i = 1, m-1
            mcsa(i) = vcsa(i) 
            mcsb(i+1) = vcsb(i)
        end do
        mcsa(m) = vcsa(m)
        mcsb(1) = 0.0
            
        info = 0
        ! compute smallest eigenvalue
        call ratqr(m, 1.0e-8_POP_r8, mcsa,mcsb,v,info)
        if (info /= 0) then 
           call POP_ErrorSet(errorCode, &
              'POP_PcsiLanczos: error estimating smallest eigenvalue !!!')
           return
        endif

        write(POP_stdout,'(a18,i3,a8,2e15.7)') "Lanczos steps ", m, " eigs:", v, u

        deallocate(mcsa, mcsb,stat = istat)
        if (istat > 0) then
           call POP_ErrorSet(errorCode, &
              'POP_PcsiLanczos: error deallocating Lanczos Tridiagonal')
           return
        endif
      endif 

      call POP_Broadcast(v,POP_masterTask,errorCode)

      ! check convergence of the smallest eigenvalue
      if (abs(1-v/mineig ) <  LanczosconvergenceCriterion) then 
          exit lanczos_iter
      endif 

      mineig = v  ! save previous value 
   endif 
    
   enddo lanczos_iter
!  call getlooptime(stop_POP_SolversMod_loop30)
!  if (POP_myTask ==  POP_masterTask) then
! WRITE (*,*) 'POP_SolversMod_loop30 costs w zhuangy: ',stop_POP_SolversMod_loop30-start_POP_SolversMod_loop30,' (us);itor-index:',maxlanczosstep
! WRITE (*,*) 'POP_SolversMod_loop30_1 costs w zhuangy:',stop_POP_SolversMod_loop30_1-start_POP_SolversMod_loop30_1,' (us);'
! WRITE (*,*) 'POP_SolversMod_loop30_2 costs w zhuangy:',stop_POP_SolversMod_loop30_2-start_POP_SolversMod_loop30_2,' (us);'!athread_30_2:',stop_POP_SolversMod_loop10-start_POP_SolversMod_loop10,'us'
! WRITE (*,*) 'POP_SolversMod_loop30_3 costs w zhuangy:',stop_POP_SolversMod_loop30_3-start_POP_SolversMod_loop30_3,' (us);'!athread_30_3:',stop_POP_SolversMod_loop29-start_POP_SolversMod_loop29,'us'
! WRITE (*,*) 'POP_SolversMod_loop30_4 costs w zhuangy:',stop_POP_SolversMod_loop30_4-start_POP_SolversMod_loop30_4,' (us);'
! WRITE (*,*) 'POP_SolversMod_loop30_5 costs w zhuangy:',stop_POP_SolversMod_loop30_5-start_POP_SolversMod_loop30_5,' (us);'!athread_30_5:',stop_POP_SolversMod_loop15-start_POP_SolversMod_loop15,'us'
!endif
 
 

!-----------------------------------------------------------------------
!EOC

 end subroutine PcsiLanczos

 subroutine EvpBlockPartition(m,mm,mb,mdi,errorCode)
! !DESCRIPTION:
!  This routine implements the strategy to divide the local barotropic block 
!  into smaller blocks in one direction.
!  EVP can not handle block larger than [12x12] while a small block size has
!  an adverse effect on the convergence rate.
!
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University

 ! !INPUT PARAMETERS:
   integer(POP_i4), intent(in) :: m              ! size of original block
   integer(POP_i4), intent(in) :: mm             ! user-defined size

 ! !OUTPUT PARAMETERS:
   integer(POP_i4), dimension(:),allocatable, intent(inout) :: &
                                  mdi  !start index of each small block
   integer(POP_i4), intent(inout) :: mb          ! number of blocks
   integer(POP_i4), intent(out) :: &
      errorCode          ! returned error code
  
 ! !LOCAL  VARIABLE
   integer :: i, istat    

   mb = (m-3)/mm +1 
   allocate(mdi(mb+1), stat = istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_EvpBlockPartition: error allocating EVP subBlock Index')
      return
   endif

   mdi(1) = 2 

   if ( mb == 1 ) then 
     ! block size is not larger than sub-block size 
     mdi(mb+1) = m  
   else 
     ! more than one sub block
 
 ! call getlooptime(start_POP_SolversMod_loop31)
     do i = 1, mb-2
       mdi(i+1) = 2+i*mm
     end do 
 ! call getlooptime(stop_POP_SolversMod_loop31)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop31 costs w zhuangy: ',stop_POP_SolversMod_loop31-start_POP_SolversMod_loop31,' (us);'
 ! endif
 

     ! in case mod(m-2, mm) /= 0, make sure the last block not too small
     mdi(mb) = (mdi(mb-1) +m )/2 
     mdi(mb+1) = m
   endif 
 end subroutine EvpBlockPartition

 subroutine inverse(a,c,n)
! !DESCRIPTION:
! This routine implemetns the Inverse matrix method 
! based on Doolittle LU factorization for Ax=b.
! Follow : Alex G. December 2009
!-----------------------------------------------------------

! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University

 ! !INPUT PARAMETERS:
  integer(POP_i4),intent(in) ::  n    !dimension
  real(POP_r8),intent(inout) :: a(n,n)!array of coefficients for matrix A

 ! !OUTPUT PARAMETERS:
  real(POP_r8),intent(inout) :: c(n,n) !inverse matrix of A

  
 ! !LOCAL VARIABLES
  real(POP_r8) :: L(n,n), U(n,n), b(n), d(n), x(n)
  real(POP_r8) :: coeff
  integer::  i, j, k

  ! step 0: initialization for matrices L and U and b
  L(:,:)=0.0
  U(:,:)=0.0
  b(:)=0.0

  ! step 1: forward elimination
 
!!  call getlooptime(start_POP_SolversMod_loop35)
  do k=1, n-1
    do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
        a(i,j) = a(i,j)-coeff*a(k,j)
      end do
    end do
  end do
!  call getlooptime(stop_POP_SolversMod_loop32)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop32 costs w zhuangy: ',stop_POP_SolversMod_loop32-start_POP_SolversMod_loop32,' (us);'
!  endif
 

  ! Step 2: prepare L and U matrices.
 
!  call getlooptime(start_POP_SolversMod_loop33)
  do i=1,n
    L(i,i) = 1.0
  end do
!  call getlooptime(stop_POP_SolversMod_loop33)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop33 costs w zhuangy: ',stop_POP_SolversMod_loop33-start_POP_SolversMod_loop33,' (us);'
!  endif
 

 
!  call getlooptime(start_POP_SolversMod_loop34)
  do j=1,n
    do i=1,j
      U(i,j) = a(i,j)
    end do
  end do
!  call getlooptime(stop_POP_SolversMod_loop34)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop34 costs w zhuangy: ',stop_POP_SolversMod_loop34-start_POP_SolversMod_loop34,' (us);'
!  endif
 

  ! Step 3: compute columns of the inverse matrix C
 
!  call getlooptime(start_POP_SolversMod_loop35)
  do k=1,n
    b(k)=1.0
    d(1) = b(1)

    do i=2,n
      d(i)=b(i)
      do j=1,i-1
        d(i)=d(i)-L(i,j)*d(j)
      end do
    end do

    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
    end do

    ! fill the solutions  x(n) into column  k  of  C
    do  i=1,n
      c(i,k) = x(i)
    end do
    b(k)=0.0
  end do 
!  call getlooptime(stop_POP_SolversMod_loop35)
!  if (POP_myTask == POP_masterTask) then
!     WRITE (*,*) 'POP_SolversMod_loop35 costs w zhuangy: ',stop_POP_SolversMod_loop35-start_POP_SolversMod_loop35,' (us);'
!  endif
 

 end subroutine inverse

 subroutine  ratqr(n, eps1, d, e, mineig,ierr )
! !DESCRIPTION:
!  This subroutine finds the algebraically smallest 
!  eigenvalue of a positive definite symmetric tridiagonal matrix by the
!  rational QR method with Newton corrections.
!  Follow EISPACK lib subroutine RATQR 
!
! ! INPUT VARIABLES
   integer(POP_i4), intent(in) ::  n    ! matrix order
   real(POP_r8),intent(in) :: d(n)      ! diagonal elements
   real(POP_r8),intent(in) :: e(n)      ! off-diagonal elements, e(1) is arbitrary
   real(POP_r8),intent(in) :: eps1      ! convergence tolerance

!  ! OUTPUT VARIABLES
   integer(POP_i4), intent(inout) ::  ierr  ! status flag 
   real(POP_r8), intent(inout) ::  mineig   !smallest eigevalue

!  ! LOCAL  VARIABLES
   real(POP_r8), dimension(n) ::  e2, bd, w
   real(POP_r8) :: f, ep, delta, err
   integer(POP_i4) i, ii ! local counters
   real(POP_r8) p, q, qp,r,s,tot

   ierr = 0
   w(1:n) = d(1:n)
   err = 0.0_POP_r8 
   s = 0.0_POP_r8 
   tot = w(1)
   q = 0.0_POP_r8 

 
 ! call getlooptime(start_POP_SolversMod_loop36)
   do i = 1, n
      p = q
      bd(i) = e(i)*e(i)
      q = 0.0_POP_r8 

      if ( i /= n ) then
        q = abs ( e(i+1) )
      endif

      tot = min ( w(i) - p - q, tot )
   end do
 ! call getlooptime(stop_POP_SolversMod_loop36)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop36 costs w zhuangy: ',stop_POP_SolversMod_loop36-start_POP_SolversMod_loop36,' (us);'
 ! endif
 
   bd(1) = 0.0_POP_r8 

   if ( tot < 0.0_POP_r8  ) then
     tot = 0.0_POP_r8 
   else 
     w(1:n) = w(1:n) - tot
   endif

!  QR transformation.

 
 ! call getlooptime(start_POP_SolversMod_loop37)
   do while (.true.)
     tot = tot + s
     delta = w(n) - s
     i = n
     if ( delta <= eps1 ) then
       exit
     endif 

     f = bd(n) / delta
     qp = delta + f
     p = 1.0

     do ii = 1, n - 1
       i = n - ii
       q = w(i) - s - f
       r = q / qp
       p = p * r + 1.0
       ep = f * r
       w(i+1) = qp + ep
       delta = q - ep

       ! check convergence 
       if ( delta <= eps1 ) exit

       f = bd(i) / q
       qp = delta + f
       bd(i+1) = qp * ep
     end do
     
     ! check convergence 
     if ( delta <= eps1 ) exit

     w(1) = qp
     s = qp / p

     !  Set error: irregular end of iteration.
     if ( tot + s <= tot ) then
       ierr = 1
       return
     endif

   end do 
 ! call getlooptime(stop_POP_SolversMod_loop37)
 ! if (POP_myTask == POP_masterTask) then
 !    WRITE (*,*) 'POP_SolversMod_loop37 costs w zhuangy: ',stop_POP_SolversMod_loop37-start_POP_SolversMod_loop37,' (us);'
 ! endif
 

   w(1) = tot
   err = err + abs ( delta)
   bd(1) = err
   mineig = w(1)

   return
 end subroutine ratqr

 end module POP_SolversMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!----total 37 loops found!----
