PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
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
  USE TimeSteppingModule_IMEX_RK, ONLY: &
    Initialize_IMEX_RK, &
    Finalize_IMEX_RK, &
    Update_IMEX_RK
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE InitializationModule, ONLY: &
    InitializeFields, &
    ComputeError
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE dgDiscretizationModule_Collisions, ONLY: &
    InitializeCollisions, &
    FinalizeCollisions, &
    ComputeIncrement_M1_DG_Implicit, &
    ComputeCorrection_M1_DG_Implicit

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(8)  :: Direction
  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER       :: nE, nX(3), bcX(3), nNodes
  REAL(DP)      :: t, dt, t_end, wTime
  REAL(DP)      :: xL(3), xR(3)
  REAL(DP)      :: eL,    eR
  REAL(DP)      :: N0, SigmaA, SigmaS
  REAL(DP)      :: Radius = 1.0d16
  REAL(DP)      :: Min_1, Max_1, Min_2

  ProgramName = 'SineWaveDiffusion'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'SineWaveStreaming' )

      ! --- Minerbo Closure Only ---

      Direction = 'X'

      nX = [ 16, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'SSPRK3'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+1
      iCycleD   = 10
      iCycleW   = 10
      maxCycles = 10000

    CASE( 'SineWaveDamping' )

      ! --- Minerbo Closure Only ---

      nX = [ 32, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_P_A2'

      N0     = 0.0_DP
      SigmaA = 1.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+1
      iCycleD   = 10
      iCycleW   = 100
      maxCycles = 100000

    CASE( 'SineWaveDiffusion' )

      nX = [ 16, 1, 1 ]
      xL = [ - 3.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ + 3.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 1, 1, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PC2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 1.0d+2

      UsePositivityLimiter = .FALSE.

      Min_1 = - HUGE( One ) ! --- Min Density
      Max_1 = + HUGE( One ) ! --- Max Density
      Min_2 = - HUGE( One ) ! --- Min "Gamma"

      t_end     = 1.0d+2
      iCycleD   = 10
      iCycleW   = 2000
      maxCycles = 100000

    CASE( 'PackedBeam' )

      nX = [ 400, 1, 1 ]
      xL = [ - 1.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ + 1.0_DP, 1.0_DP, 1.0_DP ]

      bcX = [ 2, 0, 0 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 3

      TimeSteppingScheme = 'IMEX_PC2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero ! --- Min Density
      Max_1 = One  ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 8.5d-1
      iCycleD   = 10
      iCycleW   = 10
      maxCycles = 10000

    CASE( 'LineSource' )

      nX = [ 128, 128, 1 ]
      xL = [   0.00_DP,   0.00_DP, 0.0_DP ]
      xR = [ + 1.25_DP, + 1.25_DP, 1.0_DP ]

      bcX = [ 32, 32, 1 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'SSPRK2'

      N0     = 0.0_DP
      SigmaA = 0.0_DP
      SigmaS = 0.0_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero ! --- Min Density
      Max_1 = One  ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 1.0d+0
      iCycleD   = 10
      iCycleW   = 10
      maxCycles = 10000

    CASE( 'HomogeneousSphere' )

      nX = [ 64, 64, 64 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, 2.0_DP, 2.0_DP ]

      bcX = [ 32, 32, 32 ]

      nE = 1
      eL = 0.0_DP
      eR = 1.0_DP

      nNodes = 2

      TimeSteppingScheme = 'IMEX_P_A2'

      N0     = 1.00_DP
      SigmaA = 10.0_DP
      SigmaS = 0.00_DP
      Radius = 1.00_DP

      UsePositivityLimiter = .TRUE.

      Min_1 = Zero ! --- Min Density
      Max_1 = One  ! --- Max Density
      Min_2 = Zero ! --- Min "Gamma"

      t_end     = 5.0d-0
      iCycleD   = 10
      iCycleW   = 50
      maxCycles = 100000

  END SELECT

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           BasicInitialization_Option &
             = .TRUE. )

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

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Time Stepper ---

  CALL Initialize_IMEX_RK( TRIM( TimeSteppingScheme ) )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Initialize Implicit Solver ---

  CALL InitializeCollisions &
         ( N0_Option = N0, SigmaA0_Option = SigmaA, &
           SigmaS0_Option = SigmaS, Radius_Option = Radius )

  ! --- Set Initial Condition ---

  CALL InitializeFields &
         ( Direction_Option = TRIM( Direction ) )

  ! --- Initialize Positivity Limiter ---

  CALL InitializePositivityLimiter_TwoMoment &
         ( Min_1_Option = Min_1, Max_1_Option = Max_1, Min_2_Option = Min_2, &
           UsePositivityLimiter_Option = UsePositivityLimiter )

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! --- Write Initial Condition ---

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  wTime = MPI_WTIME( )

  t  = 0.0d-0
  dt = 0.5_DP * MINVAL( (xR-xL) / DBLE( nX ) ) &
       / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t, ' to t = ', t_end
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCR, &
             ComputeIncrement_TwoMoment_Explicit, &
             ComputeIncrement_M1_DG_Implicit, &
             ComputeCorrection_M1_DG_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

    END IF

  END DO

  CALL WriteFieldsHDF( Time = t, WriteRF_Option = .TRUE. )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', wTime, ' s'
  WRITE(*,*)

  CALL ComputeError( Time = t, SigmaA = SigmaA, SigmaS = SigmaS )

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL Finalize_IMEX_RK

  CALL FinalizeCollisions

  CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
