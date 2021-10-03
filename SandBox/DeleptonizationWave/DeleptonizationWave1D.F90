PROGRAM DeleptonizationWave1D

  USE KindModule, ONLY: &
    DP, SqrtTiny, Half, Zero, &
    Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_Initialize, &
    Timer_InputOutput, &
    Timer_Evolve, &
    Timer_PositivityLimiter, &
    Timer_PL_In, &
    Timer_PL_Points, &
    Timer_PL_CellAverage, &
    Timer_PL_Theta_1, &
    Timer_PL_Theta_2, &
    Timer_PL_Out
  USE ReferenceElementModuleZ, ONLY: &
    InitializeReferenceElementZ
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T,  iAF_Ye, iAF_S,  iAF_E, iAF_Me, &
    iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE, &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ApplyEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE InitializationModule, ONLY: &
    InitializeFields_DeleptonizationWave
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
#ifdef TWOMOMENT_ORDER_1
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    InitializeNonlinearSolverTally, &
    FinalizeNonlinearSolverTally, &
    WriteNonlinearSolverTally
  USE TimeSteppingModule_Flash, ONLY: &
    Update_IMEX_PDARS
#elif TWOMOMENT_ORDER_V
  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
    InitializeTroubledCellIndicator_TwoMoment, &
    FinalizeTroubledCellIndicator_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Initialize_IMEX_RK, &
    Update_IMEX_RK
#endif

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(64) :: ProfileName
  CHARACTER(32) :: TimeSteppingScheme
  LOGICAL       :: wrt
  LOGICAL       :: UseSlopeLimiter, UsePositivityLimiter , UseEnergyLimiter
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: EvolveEuler, EvolveTwoMoment
  INTEGER       :: iCycle, iCycleD, iCycleT
  INTEGER       :: nE, bcE, nX(3), nNodes, nSpecies
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt
  REAL(DP)      :: eL, eR, ZoomE
  REAL(DP)      :: xL(3), xR(3), ZoomX(3)

  nNodes   = 2
  nSpecies = 2

  TimeSteppingScheme = 'IMEX_PDARS'

  nX = [ 128, 1, 1 ]
  xL = [ 0.0_DP           , 0.0_DP, 0.0_DP ]
  xR = [ 3.0d2 * Kilometer, Pi    , TwoPi  ]
  ZoomX = [ 1.011986923647337_DP, 1.0_DP, 1.0_DP ]

  nE = 16
  eL = 0.0d0 * MeV
  eR = 3.0d2 * MeV
  ZoomE = 1.158291374972257_DP
  bcE   = 10

  UseSlopeLimiter          = .FALSE.
  UsePositivityLimiter     = .TRUE.
  UseEnergyLimiter         = .FALSE.
  UseTroubledCellIndicator = .FALSE.
  EvolveEuler              = .FALSE.
  EvolveTwoMoment          = .TRUE.

  ProfileName = 'input_thornado_VX_100ms.dat'

  t       = 0.0_DP
  t_end   = 1.0d0 * Millisecond
  dt_wrt  = 1.0d-1 * Millisecond
  iCycleD = 1
  iCycleT = 10

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'DeleptonizationWave1D', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 01, 0, 0 ], &
           bcX_Option &
             = [ 30, 0, 0 ], &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nE_Option &
             = nE, &
           swE_Option &
             = 1, &
           bcE_Option &
             = bcE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           ZoomE_Option &
             = ZoomE, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Initialize Timers ---

  CALL InitializeTimers

  CALL TimersStart( Timer_Initialize )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElementZ

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option = 'EquationOfStateTable.h5' )

  ! --- Initialize Opacities ---

  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_EmAb_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &
           OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &
           OpacityTableName_Pair_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
           Verbose_Option = .TRUE. )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

#ifdef TWOMOMENT_ORDER_V
  ! --- Initialize Troubled Cell Indicator ---

  CALL InitializeTroubledCellIndicator_TwoMoment &
         ( UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           C_TCI_Option &
             = Zero, &
           Verbose_Option &
             = .TRUE. )

  ! --- Initialize Limiter ---

    call InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter, &
             Verbose_Option &
               = .TRUE. )

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = 0.0d0 + SqrtTiny, &
           Max_1_Option = 1.0d0 - EPSILON(1.0d0), &
           Min_2_Option = 0.0d0 + SqrtTiny, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           UseEnergyLimiter_Option &
             = UseEnergyLimiter, &
           Verbose_Option &
             = .TRUE. )

#elif TWOMOMENT_ORDER_1
  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = 0.0d0 + SqrtTiny, &
           Max_1_Option = 1.0d0 - EPSILON(1.0d0), &
           Min_2_Option = 0.0d0 + SqrtTiny, &
           UsePositivityLimiter_Option &
             = UsePositivityLimiter, &
           Verbose_Option &
             = .TRUE. )

#endif

  ! --- Set Initial Condition ---

  CALL InitializeFields_DeleptonizationWave( ProfileName )

  ! --- Initialize Time Stepper ---

#ifdef TWOMOMENT_ORDER_V
   CALL Initialize_IMEX_RK &
          ( TRIM( TimeSteppingScheme ), &
            EvolveEuler_Option = EvolveEuler, &
            EvolveTwoMoment_Option = EvolveTwoMoment )
#endif

  ! --- Write Initial Condition Before Limiter ---

  CALL TimersStart( Timer_InputOutput )

  CALL ComputeFromConserved_Fluid

  CALL ComputeFromConserved_Radiation

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

#ifdef TWOMOMENT_ORDER_1
  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

#elif TWOMOMENT_ORDER_V
  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )
#endif

  ! Reset these timers
  Timer_PositivityLimiter = Zero
  Timer_PL_In             = Zero
  Timer_PL_Points         = Zero
  Timer_PL_CellAverage    = Zero
  Timer_PL_Theta_1        = Zero
  Timer_PL_Theta_2        = Zero
  Timer_PL_Out            = Zero

  ! --- Write Initial Condition After Limiter ---

  CALL TimersStart( Timer_InputOutput )

  CALL ComputeFromConserved_Fluid

  CALL ComputeFromConserved_Radiation

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  CALL TimersStop( Timer_Initialize )

#ifdef TWOMOMENT_ORDER_1
  CALL InitializeNonlinearSolverTally
#endif

  ! --- Evolve ---

  t_wrt = dt_wrt
  wrt   = .FALSE.

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, &
    ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    dt = Half * MINVAL( MeshX(1) % Width(iX_B0(1):iX_E0(1)) ) &
           / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t / Millisecond, &
          '', 'dt = ', dt / Millisecond

    END IF

    CALL TimersStart( Timer_Evolve )

#ifdef TWOMOMENT_ORDER_V 
    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit )

#elif TWOMOMENT_ORDER_1
    CALL Update_IMEX_PDARS &
           ( dt, uCF, uCR, &
             Explicit_Option = .TRUE., &
             Implicit_Option = .TRUE., &
             SingleStage_Option = .FALSE., &
             CallFromThornado_Option = .TRUE. )
#endif

    t = t + dt

    CALL TimersStop( Timer_Evolve )

    IF( wrt )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( uCF, uCR )
#endif

      CALL TimersStart( Timer_InputOutput )

      CALL ComputeFromConserved_Fluid

      CALL ComputeFromConserved_Radiation

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL TimersStop( Timer_InputOutput )

      wrt = .FALSE.

    END IF

#ifdef TWOMOMENT_ORDER_1
    IF( MOD( iCycle, iCycleT ) == 0 )THEN

      CALL WriteNonlinearSolverTally( t )

    END IF
#endif

  END DO

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( uCF, uCR )
#endif

  ! --- Write Final Solution ---

  CALL TimersStart( Timer_InputOutput )

  CALL ComputeFromConserved_Fluid

  CALL ComputeFromConserved_Radiation

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  ! --- Finalize ---

  CALL FinalizeTimers

#ifdef TWOMOMENT_ORDER_1
  CALL FinalizeNonlinearSolverTally
#endif

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

#ifdef TWOMOMENT_ORDER_V
  CALL FinalizeTroubledCellIndicator_TwoMoment

  CALL FinalizeSlopeLimiter_TwoMoment
#endif

  CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

CONTAINS


  SUBROUTINE ComputeFromConserved_Fluid

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputePrimitive_Euler &
             ( uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeThermodynamicStates_Auxiliary_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye) )

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), &
               uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), &
               uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), &
               uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), &
               uAF(:,iX1,iX2,iX3,iAF_Gm) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_Fluid


  SUBROUTINE ComputeFromConserved_Radiation

  END SUBROUTINE ComputeFromConserved_Radiation


END PROGRAM DeleptonizationWave1D
