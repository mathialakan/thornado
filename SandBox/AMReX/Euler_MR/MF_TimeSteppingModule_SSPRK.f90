MODULE MF_TimeSteppingModule_SSPRK

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_regrid, &
    amrex_get_numlevels

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
  USE MF_Euler_TallyModule, ONLY: &
    IncrementOffGridTally_Euler_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_FieldsModule, ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uDF, &
    MF_OffGridFlux_Euler
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    nMaxLevels, &
    swX, &
    CFL, &
    nNodes, &
    nStages, &
    do_reflux, &
    t_new, &
    dt, &
    UseAMR
  USE RefluxModule_Euler, ONLY: &
    Reflux_Euler_MF
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_UpdateFluid, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_CopyMultiFab

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFluid_SSPRK_MF
  PUBLIC :: UpdateFluid_SSPRK_MF
  PUBLIC :: FinalizeFluid_SSPRK_MF

  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  LOGICAL :: Verbose

CONTAINS


  SUBROUTINE InitializeFluid_SSPRK_MF( Verbose_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL InitializeSSPRK

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Timestepper'
      WRITE(*,'(A)') &
        '    -----------------'
      WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages
      WRITE(*,'(A5,A,ES10.3E3)') '', 'CFL:           ', &
        CFL * ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Butcher Table:'
      WRITE(*,'(A5,A)') '', '--------------'
      DO iS = 1, nStages
        WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(iS), a_SSPRK(iS,1:nStages)
      END DO
      WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)
      WRITE(*,*)

    END IF

  END SUBROUTINE InitializeFluid_SSPRK_MF


  SUBROUTINE FinalizeFluid_SSPRK_MF

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

  END SUBROUTINE FinalizeFluid_SSPRK_MF


  SUBROUTINE UpdateFluid_SSPRK_MF

    TYPE(amrex_multifab) :: MF_U(1:nStages,0:nMaxLevels-1)
    TYPE(amrex_multifab) :: MF_D(1:nStages,0:nMaxLevels-1)

    INTEGER :: iS, jS, nCompCF
    INTEGER :: iLevel

    REAL(DP) :: dM_OffGrid_Euler(1:nCF,0:nMaxLevels-1)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

    dM_OffGrid_Euler = Zero

    nCompCF = nDOFX * nCF

    DO iS = 1, nStages

      IF( iS .EQ. 1 )THEN

        IF( UseAMR )THEN

          DO iLevel = 0, nLevels-1

            IF( iLevel .LT. nLevels-1 ) &
              CALL amrex_regrid( iLevel, t_new(iLevel) )

          END DO

          nLevels = amrex_get_numlevels()

          CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

        END IF

      END IF

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_U(iS,iLevel), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL amrex_multifab_build &
               ( MF_D(iS,iLevel), MF_uCF(iLevel) % BA, &
                 MF_uCF(iLevel) % DM, nCompCF, swX )

        CALL MF_U(iS,iLevel) % COPY( MF_uCF(iLevel), 1, 1, nCompCF, swX )

      END DO ! iLevel

      DO jS = 1, iS-1

        DO iLevel = 0, nLevels-1

          IF( a_SSPRK(iS,jS) .NE. Zero ) &
            CALL MF_U(iS,iLevel) &
                   % LinComb( One, MF_U(iS,iLevel), 1, &
                              dt(iLevel) * a_SSPRK(iS,jS), MF_D(jS,iLevel), 1, &
                              1, nCompCF, 0 )

        END DO ! iLevel

      END DO ! jS

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        CALL ApplySlopeLimiter_Euler_MF &
               ( t_new, MF_uGF, MF_U(iS,:), MF_uDF )

        CALL ApplyPositivityLimiter_Euler_MF &
               ( MF_uGF, MF_U(iS,:), MF_uDF )

        CALL ComputeIncrement_Euler_MF &
               ( t_new, MF_uGF, MF_U(iS,:), MF_uDF, MF_D(iS,:) )

        DO iLevel = 0, nLevels-1

          dM_OffGrid_Euler(:,iLevel) &
            = dM_OffGrid_Euler(:,iLevel) &
                + dt(iLevel) * w_SSPRK(iS) * MF_OffGridFlux_Euler(:,iLevel)

        END DO

        IF( do_reflux ) CALL Reflux_Euler_MF( MF_uGF, MF_D(iS,:) )

      END IF ! a(:,iS) .NE. Zero OR w(iS) .NE. Zero

    END DO ! iS

    DO iS = 1, nStages

      DO iLevel = 0, nLevels-1

        IF( w_SSPRK(iS) .NE. Zero ) &
          CALL MF_uCF(iLevel) &
                 % LinComb( One, MF_uCF(iLevel), 1, &
                            dt(iLevel) * w_SSPRK(iS), MF_D(iS,iLevel), 1, 1, &
                            nCompCF, 0 )

      END DO ! iLevel

    END DO ! iS

    DO iLevel = 0, nLevels-1

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_U(iS,iLevel) )
        CALL amrex_multifab_destroy( MF_D(iS,iLevel) )

      END DO ! iS

    END DO ! iLevel

    CALL ApplySlopeLimiter_Euler_MF( t_new, MF_uGF, MF_uCF, MF_uDF )
    CALL ApplyPositivityLimiter_Euler_MF( MF_uGF, MF_uCF, MF_uDF )

    CALL IncrementOffGridTally_Euler_MF( dM_OffGrid_Euler )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_UpdateFluid )

  END SUBROUTINE UpdateFluid_SSPRK_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE InitializeSSPRK

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK

    SELECT CASE ( nStages )

      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_DP
        w_SSPRK(1)   = 1.0_DP

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSPRK(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSPRK(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSPRK(1:3)   = [ 1.0_DP / 6.0_DP, &
                           1.0_DP / 6.0_DP, &
                           2.0_DP / 3.0_DP ]

    END SELECT

    DO iS = 1, nStages

      c_SSPRK(iS) = SUM( a_SSPRK(iS,1:iS-1) )

    END DO

  END SUBROUTINE InitializeSSPRK


  SUBROUTINE AllocateButcherTables_SSPRK

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


END MODULE MF_TimeSteppingModule_SSPRK
