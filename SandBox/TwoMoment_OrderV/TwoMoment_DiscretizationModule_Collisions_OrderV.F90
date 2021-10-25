MODULE TwoMoment_DiscretizationModule_Collisions_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One, Half, Three, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFE, &
    nDOFZ
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions, &
    Timer_Collisions_Zero, &
    Timer_Collisions_Permute, &
    Timer_Collisions_PrimitiveFluid, &
    Timer_Collisions_Solve
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeEddingtonTensorComponents_dd
  USE TwoMoment_OpacityModule_OrderV, ONLY: &
    uOP, iOP_D0, iOP_Chi, iOP_Sigma, nOP
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit

  INTEGER :: iE_B0,    iE_E0
  INTEGER :: iX_B0(3), iX_E0(3)
  INTEGER :: nZ(4), nE, nX(3), nE_G, nX_G

  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    N_PTR, G1_PTR, G2_PTR, G3_PTR
  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    D0_PTR, Chi_PTR, Sigma_PTR
  REAL(DP), POINTER, CONTIGUOUS, DIMENSION(:) :: &
    dN_PTR, dG1_PTR, dG2_PTR, dG3_PTR

  REAL(DP), ALLOCATABLE, TARGET :: GX_N(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: PF_N(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: CF_N(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: CR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: dCR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: OP_N(:,:,:,:)

  INTEGER,  ALLOCATABLE :: PositionIndexZ(:)
  INTEGER,  ALLOCATABLE :: nIterations_Simple(:)

CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U_F, dU_F, U_R, dU_R )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      dt
    REAL(DP), INTENT(in)    :: &
      GE  (1:nDOFE, &
           iZ_B1(1):iZ_E1(1), &
           1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX  (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F (1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R (1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR,1:nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS, iGF, iCF, iOP
    INTEGER :: iX1, iX2, iX3, iE
    INTEGER :: iNodeZ, iNodeX, iNodeE, iN_X, iN_E

    CALL TimersStart( Timer_Collisions )

    CALL InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0 ) &
    !$OMP MAP( alloc: dU_F, dU_R )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( GX, U_F, U_R, uOP, iZ_B1, iZ_E1, iX_B0 ) &
    !$ACC CREATE( dU_F, dU_R )
#endif

    CALL TimersStart( Timer_Collisions_Zero )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( dU_F, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iCF = 1, nCF
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)

      DO iNodeX = 1, nDOFX

        dU_F(iNodeX,iZ2,iZ3,iZ4,iCF) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO



#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRESENT( dU_R, iZ_B1, iZ_E1 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(7)
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        dU_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO


    CALL TimersStop( Timer_Collisions_Zero )





    CALL TimersStart( Timer_Collisions_Permute )

    ! --- Arrange Geometry Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, GX_N, GX )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iX1, iX2, iX3, iNodeX )
#endif
    DO iGF  = 1, nGF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      GX_N(iN_X,iGF) = GX(iNodeX,iX1,iX2,iX3,iGF)

    END DO
    END DO


    ! --- Arrange Fluid Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iNodeX, iX1, iX2, iX3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iNodeX, iX1, iX2, iX3 ) &
    !$ACC PRESENT( nX, iX_B0, CF_N, U_F )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iX1, iX2, iX3, iNodeX )
#endif
    DO iCF  = 1, nCF
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      CF_N(iN_X,iCF) = U_F(iNodeX,iX1,iX2,iX3,iCF)

    END DO
    END DO


    ! --- Arrange Radiation Fields ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeE, iNodeZ, iNodeX, iX1, iX2, iX3, iE )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iX1, iX2, iX3, iE ) &
    !$ACC PRESENT( nZ, nX, iX_B0, CR_N, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iE, iX1, iX2, iX3, iNodeE, iNodeX, iNodeZ )
#endif
    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      CR_N(iN_E,iN_X,iS,iCR) = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO


    ! --- Arrange Opacities ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iNodeE, iNodeZ, iNodeX, iX1, iX2, iX3, iE )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRIVATE( iNodeZ, iNodeX, iNodeE, iX1, iX2, iX3, iE ) &
    !$ACC PRESENT( nZ, nX, iX_B0, OP_N, uOP )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( iE, iX1, iX2, iX3, iNodeE, iNodeX, iNodeZ )
#endif
    DO iS   = 1, nSpecies
    DO iOP  = 1, nOP
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

      OP_N(iN_E,iN_X,iS,iOP) = uOP(iNodeZ,iE,iX1,iX2,iX3,iOP,iS)

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Collisions_Permute )


    CALL TimersStart( Timer_Collisions_PrimitiveFluid )



#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( CF_N, PF_N, GX_N )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( CF_N(iN_X,iCF_D ), &
               CF_N(iN_X,iCF_S1), &
               CF_N(iN_X,iCF_S2), &
               CF_N(iN_X,iCF_S3), &
               CF_N(iN_X,iCF_E ), &
               CF_N(iN_X,iCF_Ne), &
               PF_N(iN_X,iPF_D ), &
               PF_N(iN_X,iPF_V1), &
               PF_N(iN_X,iPF_V2), &
               PF_N(iN_X,iPF_V3), &
               PF_N(iN_X,iPF_E ), &
               PF_N(iN_X,iPF_Ne), &
               GX_N(iN_X,iGF_Gm_dd_11), &
               GX_N(iN_X,iGF_Gm_dd_22), &
               GX_N(iN_X,iGF_Gm_dd_33) )

    END DO

    CALL TimersStop( Timer_Collisions_PrimitiveFluid )


    CALL TimersStart( Timer_Collisions_Solve )


! #if   defined( THORNADO_OMP_OL )
!     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
! #elif defined( THORNADO_OACC   )
!     !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
!     !$ACC PRESENT( CR_N, PF_N, GX_N, OP_N, dCR_N )
! #elif defined( THORNADO_OMP    )
!     !$OMP PARALLEL DO SIMD COLLAPSE(3)
! #endif
    ! DO iN_X = 1, nX_G
    ! DO iN_E = 1, nE_G
    ! DO iS   = 1, nSpecies
    !
    !   CALL ComputeIncrement_FixedPoint_Richardson &
    !          ( dt, &
    !            CR_N (iN_E,iN_X,iS,iCR_N  ), &
    !            CR_N (iN_E,iN_X,iS,iCR_G1 ), &
    !            CR_N (iN_E,iN_X,iS,iCR_G2 ), &
    !            CR_N (iN_E,iN_X,iS,iCR_G3 ), &
    !            PF_N(iN_X,iPF_V1), &
    !            PF_N(iN_X,iPF_V2), &
    !            PF_N(iN_X,iPF_V3), &
    !            GX_N(iN_X,iGF_Gm_dd_11), &
    !            GX_N(iN_X,iGF_Gm_dd_22), &
    !            GX_N(iN_X,iGF_Gm_dd_33), &
    !            OP_N (iN_E,iN_X,iS,iOP_D0    ), &
    !            OP_N (iN_E,iN_X,iS,iOP_Chi   ), &
    !            OP_N (iN_E,iN_X,iS,iOP_Sigma ), &
    !            dCR_N(iN_E,iN_X,iS,iCR_N     ), &
    !            dCR_N(iN_E,iN_X,iS,iCR_G1    ), &
    !            dCR_N(iN_E,iN_X,iS,iCR_G2    ), &
    !            dCR_N(iN_E,iN_X,iS,iCR_G3    ) )
    !
    ! END DO
    ! END DO
    ! END DO


    CALL ComputeIncrement_FixedPoint_Vector_Richardson &
           ( dt, N_PTR, G1_PTR, G2_PTR, G3_PTR, &
             PF_N(:,iPF_V1), PF_N(:,iPF_V2), PF_N(:,iPF_V3), &
             GX_N(:,iGF_Gm_dd_11), GX_N(:,iGF_Gm_dd_22), GX_N(:,iGF_Gm_dd_33), &
             D0_PTR, Chi_PTR, Sigma_PTR, &
             dN_PTR, dG1_PTR, dG2_PTR, dG3_PTR, &
             PositionIndexZ, nIterations_Simple)

    CALL TimersStop( Timer_Collisions_Solve )

    CALL TimersStart( Timer_Collisions_Permute )

    ! --- Revert Radiation Increment ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeE, iNodeX, iN_X, iN_E )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) &
    !$ACC PRIVATE( iNodeX, iNodeE, iN_X, iN_E ) &
    !$ACC PRESENT( nX, iX_B0, iX_E0, dU_R, dCR_N )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP PRIVATE( iNodeE, iNodeX, iN_E, iN_X )
#endif
    DO iS  = 1, nSpecies
    DO iCR = 1, nCR
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        iN_X = iNodeX &
                 + (iX1-iX_B0(1)) * nDOFX &
                 + (iX2-iX_B0(2)) * nDOFX * nX(1) &
                 + (iX3-iX_B0(3)) * nDOFX * nX(1) * nX(2)

        iN_E = iNodeE &
                 + (iE-iE_B0) * nDOFE

        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) = dCR_N(iN_E,iN_X,iS,iCR)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Collisions_Permute )


#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: dU_F, dU_R ) &
    !$OMP MAP( release: GX, U_F, U_R, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( dU_F, dU_R ) &
    !$ACC DELETE( GX, U_F, U_R, iZ_B1, iZ_E1 )
#endif

    CALL FinalizeCollisions

    CALL TimersStop( Timer_Collisions )

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit


  SUBROUTINE ComputeIncrement_FixedPoint &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      dN, dG_d_1, dG_d_2, dG_d_3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, j, k, mk, INFO
    REAL(DP) :: D, I_d_1, I_d_2, I_d_3, Kappa
    REAL(DP) ::    I_u_1, I_u_2, I_u_3
    REAL(DP) :: A_d_1, A_d_2, A_d_3, k_dd(3,3)
    ! REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: D_00, D_ii
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: LMAT(4,4), DET, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)

    Kappa = Chi + Sigma

    D_00 = One + dt * Chi
    D_ii = One + dt * Kappa

    ! --- Constant Vector ---

    CVEC = [ N + dt * Chi * D_0, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      ! CALL ComputeEddingtonTensorComponents_dd &
      !        ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      A_d_1 = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
      A_d_2 = V_u_1 * k_dd(1,2) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
      A_d_3 = V_u_1 * k_dd(1,3) + V_u_2 * k_dd(2,3) + V_u_3 * k_dd(3,3)

      DET = ( D_00 * D_ii &
              - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 + V_u_3 * A_d_3 ) ) * D_ii**2

      LMAT(1,1) = D_ii**3
      LMAT(2,1) = - A_d_1 * D_ii**2
      LMAT(3,1) = - A_d_2 * D_ii**2
      LMAT(4,1) = - A_d_3 * D_ii**2

      LMAT(1,2) = - V_u_1 * D_ii**2
      LMAT(2,2) = D_00 * D_ii**2 &
                    - ( V_u_2 * A_d_2 + V_u_3 * A_d_3 ) * D_ii
      LMAT(3,2) = V_u_1 * A_d_2 * D_ii
      LMAT(4,2) = V_u_1 * A_d_3 * D_ii

      LMAT(1,3) = - V_u_2 * D_ii**2
      LMAT(2,3) = V_u_2 * A_d_1 * D_ii
      LMAT(3,3) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_3 * A_d_3 ) * D_ii
      LMAT(4,3) = V_u_2 * A_d_3 * D_ii

      LMAT(1,4) = - V_u_3 * D_ii**2
      LMAT(2,4) = V_u_3 * A_d_1 * D_ii
      LMAT(3,4) = V_u_3 * A_d_2 * D_ii
      LMAT(4,4) = D_00 * D_ii**2 &
                    - ( V_u_1 * A_d_1 + V_u_2 * A_d_2 ) * D_ii

      LMAT = LMAT / DET

      ! PRINT*, "Doing DGEMV"
      ! CALL DGEMV( 'N', 4, 4, One, LMAT, 4, CVEC, 1, Zero, GVEC(:,mk), 1 )

      PRINT*, "Doing GPU multiplication"
      GVEC(:,mk) = Zero

      DO j = 1, 4
      DO i = 1, 4

        GVEC(i,mk) = GVEC(i,mk) + LMAT(i,j) * CVEC(j)

      END DO
      END DO

      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        ! --- Anderson Accelerated Fixed-Point ---

        ! BVEC = - FVEC(:,mk)
        !
        ! AMAT(:,1:mk-1) &
        !  = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )
        !
        ! CALL DGELS( 'N', 4, mk-1, 1, AMAT(:,1:mk-1), 4, BVEC, 4, &
        !            WORK, LWORK, INFO )
        !
        ! Alpha(1:mk-1) = BVEC(1:mk-1)
        ! Alpha(mk)     = One - SUM( Alpha(1:mk-1) )

        !CALL SolveAlpha_LS( M, mk, FVEC, Alpha )

        Alpha = Alpha_LS( M, mk, FVEC )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        !GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        !FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

        !CALL ShiftVectors( M, mk, FVEC, GVEC )

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    dN     = Chi * ( D_0 - D )
    dG_d_1 = - Kappa * I_d_1
    dG_d_2 = - Kappa * I_d_2
    dG_d_3 = - Kappa * I_d_3

    ! IF( k == MaxIterations )THEN

    !   PRINT*
    !   PRINT*, "ComputeIncrement_FixedPoint"
    !   PRINT*
    !   PRINT*, "  N     = ", N
    !   PRINT*, "  G_d_1 = ", G_d_1
    !   PRINT*, "  G_d_2 = ", G_d_2
    !   PRINT*, "  G_d_3 = ", G_d_3
    !   PRINT*
    !   PRINT*, "  V_u_1 = ", V_u_1
    !   PRINT*, "  V_u_2 = ", V_u_2
    !   PRINT*, "  V_u_3 = ", V_u_3
    !   PRINT*

    !   PRINT*, "  Converged with k = ", k

    !   PRINT*
    !   PRINT*, "  FVECm = ", FVECm
    !   PRINT*

    !   PRINT*
    !   PRINT*, "  D     = ", D
    !   PRINT*, "  I_u_1 = ", I_u_1
    !   PRINT*, "  I_u_2 = ", I_u_2
    !   PRINT*, "  I_u_3 = ", I_u_3
    !   PRINT*

    ! END IF

  END SUBROUTINE ComputeIncrement_FixedPoint

  SUBROUTINE ComputeIncrement_FixedPoint_Richardson &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      dN, dG_d_1, dG_d_2, dG_d_3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: dt
    REAL(DP), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    LOGICAL  :: CONVERGED
    INTEGER  :: i, j, k, mk, INFO
    REAL(DP) :: D, I_d_1, I_d_2, I_d_3, Kappa
    REAL(DP) ::    I_u_1, I_u_2, I_u_3
    REAL(DP) :: k_dd(3,3)
    ! REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: D_00, D_ii
    REAL(DP) :: UVEC(4), CVEC(4)
    REAL(DP) :: GVEC(4,M), GVECm(4)
    REAL(DP) :: FVEC(4,M), FVECm(4)
    REAL(DP) :: vMag, Omega, vI, vK, Alpha(M)
    REAL(DP) :: BVEC(4), AMAT(4,M), WORK(LWORK)

    Kappa = Chi + Sigma

    D_00 = One + dt * Chi
    D_ii = One + dt * Kappa

    ! --- Constant Vector ---

    CVEC = [ N, G_d_1, G_d_2, G_d_3 ]

    ! --- Initial Guess ---

    D     = N
    I_u_1 = Zero
    I_u_2 = Zero
    I_u_3 = Zero

    I_d_1 = Gm_dd_11 * I_u_1
    I_d_2 = Gm_dd_22 * I_u_2
    I_d_3 = Gm_dd_33 * I_u_3

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIterations )

      k  = k + 1
      mk = MIN( M, k )

      UVEC = [ D, I_d_1, I_d_2, I_d_3 ]

      ! CALL ComputeEddingtonTensorComponents_dd &
      !        ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

      k_dd = EddingtonTensorComponents_dd &
               ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )


      vMag = SQRT(V_u_1 * Gm_dd_11 * V_u_1 &
                + V_u_2 * Gm_dd_22 * V_u_2 &
                + V_u_3 * Gm_dd_33 * V_u_3)

      Omega = One / (One + vMag)

      vI = V_u_1 * UVEC(2) + V_u_2 *  UVEC(3) + V_u_3 *  UVEC(4)
      GVEC(1,mk) = (One - Omega) * UVEC(1) + &
                   Omega / D_00 * (CVEC(1) + dt * Chi * D_0 - vI)

      DO j = 1, 3
        vK = V_u_1 * k_dd(j,1) + V_u_2 * k_dd(j,2) + V_u_3 * k_dd(j,3)
        GVEC(j+1,mk) = (One - Omega) * UVEC(j+1) + &
                     Omega / D_ii * (CVEC(j+1) - vK * UVEC(1))
      END DO


      FVEC(:,mk) = GVEC(:,mk) - UVEC

      IF( mk == 1 )THEN

        ! --- Picard Iteration ---

        GVECm = GVEC(:,mk)

      ELSE

        Alpha = Alpha_LS( M, mk, FVEC )

        GVECm = Zero
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm = GVECm - UVEC

      IF( ALL( ABS( FVECm ) <= Rtol * ABS( CVEC ) ) )THEN

        CONVERGED = .TRUE.

      END IF

      UVEC = GVECm

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        FVEC = ShiftVec( M, mk, FVEC )
        GVEC = ShiftVec( M, mk, GVEC )

      END IF

      D     = UVEC(1)
      I_d_1 = UVEC(2); I_u_1 = I_d_1 / Gm_dd_11
      I_d_2 = UVEC(3); I_u_2 = I_d_2 / Gm_dd_22
      I_d_3 = UVEC(4); I_u_3 = I_d_3 / Gm_dd_33

    END DO

    dN     = Chi * ( D_0 - D )
    dG_d_1 = - Kappa * I_d_1
    dG_d_2 = - Kappa * I_d_2
    dG_d_3 = - Kappa * I_d_3

    ! IF( k == MaxIterations )THEN

    !   PRINT*
    !   PRINT*, "ComputeIncrement_FixedPoint"
    !   PRINT*
    !   PRINT*, "  N     = ", N
    !   PRINT*, "  G_d_1 = ", G_d_1
    !   PRINT*, "  G_d_2 = ", G_d_2
    !   PRINT*, "  G_d_3 = ", G_d_3
    !   PRINT*
    !   PRINT*, "  V_u_1 = ", V_u_1
    !   PRINT*, "  V_u_2 = ", V_u_2
    !   PRINT*, "  V_u_3 = ", V_u_3
    !   PRINT*

    !   PRINT*, "  Converged with k = ", k

    !   PRINT*
    !   PRINT*, "  FVECm = ", FVECm
    !   PRINT*

    !   PRINT*
    !   PRINT*, "  D     = ", D
    !   PRINT*, "  I_u_1 = ", I_u_1
    !   PRINT*, "  I_u_2 = ", I_u_2
    !   PRINT*, "  I_u_3 = ", I_u_3
    !   PRINT*

    ! END IF

  END SUBROUTINE ComputeIncrement_FixedPoint_Richardson

  SUBROUTINE ComputeIncrement_FixedPoint_Vector_Richardson &
    ( dt, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, D_0, Chi, Sigma, &
      dN, dG_d_1, dG_d_2, dG_d_3, PositionIndexZ, nIterations_Option )

    REAL(DP),               INTENT(in)  :: dt
    REAL(DP), DIMENSION(:), INTENT(in)  :: N, G_d_1, G_d_2, G_d_3
    REAL(DP), DIMENSION(:), INTENT(in)  ::    V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:), INTENT(in)  :: D_0, Chi, Sigma
    REAL(DP), DIMENSION(:), INTENT(out) :: dN, dG_d_1, dG_d_2, dG_d_3
    INTEGER,  DIMENSION(:), INTENT(in)  :: PositionIndexZ
    INTEGER,  DIMENSION(:), INTENT(out), OPTIONAL :: nIterations_Option

    ! --- Parameters ---

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: MaxIterations = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08

    ! --- Local Variables ---

    INTEGER  :: nZ
    INTEGER  :: iX, iZ
    INTEGER  :: k, Mk, iM, i, j

    REAL(DP) :: FTMP(4,M), GTMP(4,M)
    REAL(DP) :: SUM1, k_dd(3,3)
    REAL(DP) :: vMag, Omega, vI, vK
    LOGICAL  :: CONVERGED

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FVEC, GVEC
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: CVEC, UVEC, FVECm, GVECm, Alpha
    REAL(DP), DIMENSION(:),     ALLOCATABLE :: D, I_u_1, I_u_2, I_u_3
    REAL(DP), DIMENSION(:),     ALLOCATABLE :: D_00, D_ii, Kappa
    LOGICAL,  DIMENSION(:),     ALLOCATABLE :: ITERATE
    INTEGER,  DIMENSION(:),     ALLOCATABLE :: nIterations

    nZ = SIZE( N, 1 )


    ALLOCATE( FVEC(4,M,nZ) )
    ALLOCATE( GVEC(4,M,nZ) )

    ALLOCATE( CVEC (4,nZ) )
    ALLOCATE( UVEC (4,nZ) )
    ALLOCATE( FVECm(4,nZ) )
    ALLOCATE( GVECm(4,nZ) )
    ALLOCATE( Alpha(M,nZ) )

    ALLOCATE( D    (nZ) )
    ALLOCATE( I_u_1(nZ) )
    ALLOCATE( I_u_2(nZ) )
    ALLOCATE( I_u_3(nZ) )

    ALLOCATE( D_00 (nZ) )
    ALLOCATE( D_ii (nZ) )
    ALLOCATE( Kappa(nZ) )

    ALLOCATE( ITERATE(nZ) )
    ALLOCATE( nIterations(nZ) )

    ITERATE = .TRUE.

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: ITERATE ) &
  !$OMP MAP( alloc: FVEC, GVEC, CVEC, UVEC, &
  !$OMP             FVECm, GVECm, Alpha, nIterations, &
  !$OMP             D, I_u_1, I_u_2, I_u_3, &
  !$OMP             Kappa, D_00, D_ii )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA ASYNC &
  !$ACC COPYIN( ITERATE ) &
  !$ACC CREATE( FVEC, GVEC, CVEC, UVEC, &
  !$ACC         FVECm, GVECm, Alpha, nIterations, &
  !$ACC         D, I_u_1, I_u_2, I_u_3, &
  !$ACC         Kappa, D_00, D_ii )
#endif

    ! --- Initial Guess ---

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
  !$ACC PRESENT( CVEC, N, G_d_1, G_d_2, G_d_3, &
  !$ACC          D, I_u_1, I_u_2, I_u_3, &
  !$ACC          Chi, Sigma, Kappa, D_00, D_ii )
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO
#endif
    DO iZ = 1, nZ
      CVEC(iCR_N ,iZ) = N    (iZ)
      CVEC(iCR_G1,iZ) = G_d_1(iZ)
      CVEC(iCR_G2,iZ) = G_d_2(iZ)
      CVEC(iCR_G3,iZ) = G_d_3(iZ)

      D    (iZ) = N(iZ)
      I_u_1(iZ) = Zero
      I_u_2(iZ) = Zero
      I_u_3(iZ) = Zero

      Kappa(iZ) = Chi(iZ) + Sigma(iZ)

      D_00(iZ) = One + dt * Chi(iZ)
      D_ii(iZ) = One + dt * Kappa(iZ)
    END DO


    k = 0
    DO WHILE( ANY( ITERATE ) .AND. k < MaxIterations )

      k = k + 1
      Mk = MIN( M, k )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX, k_dd, vMag, Omega, vI, vK )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
    !$ACC PRIVATE( iX, k_dd, vMag, Omega, vI, vK ) &
    !$ACC PRESENT( ITERATE, UVEC, CVEC, GVEC, FVEC, GVECm, FVECm, &
    !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, V_u_1, V_u_2, V_u_3, &
    !$ACC          Chi, D_0, D_00, D_ii )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX, k_dd, vMag, Omega, vI, vK )
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          UVEC(iPR_D ,iZ) = D    (iZ)
          UVEC(iPR_I1,iZ) = I_u_1(iZ) * Gm_dd_11(iX)
          UVEC(iPR_I2,iZ) = I_u_2(iZ) * Gm_dd_22(iX)
          UVEC(iPR_I3,iZ) = I_u_3(iZ) * Gm_dd_33(iX)

          k_dd = EddingtonTensorComponents_dd &
                   ( D(iZ), I_u_1(iZ), I_u_2(iZ), I_u_3(iZ), &
                     Gm_dd_11(iX), Gm_dd_22(iX), Gm_dd_33(iX) )

          vMag = SQRT(   V_u_1(iX) * Gm_dd_11(iX) * V_u_1(iX) &
                       + V_u_2(iX) * Gm_dd_22(iX) * V_u_2(iX) &
                       + V_u_3(iX) * Gm_dd_33(iX) * V_u_3(iX) )

          Omega = One / ( One + vMag )

          vI =   V_u_1(iX) * UVEC(iPR_I1,iZ) &
               + V_u_2(iX) * UVEC(iPR_I2,iZ) &
               + V_u_3(iX) * UVEC(iPR_I3,iZ)

          GVECm(1,iZ) = (One - Omega) * UVEC(iPR_D,iZ) &
                        + Omega / D_00(iZ) * &
                        ( CVEC(iCR_N,iZ) + dt * Chi(iZ) * D_0(iZ) - vI )

          DO j = 1, 3

            vK =   V_u_1(iX) * k_dd(j,1) &
                 + V_u_2(iX) * k_dd(j,2) &
                 + V_u_3(iX) * k_dd(j,3)

            GVECm(j+1,iZ) = (One - Omega) * UVEC(j+1,iZ) &
                            + Omega / D_ii(iZ) * &
                            ( CVEC(j+1,iZ) - vK * UVEC(iPR_D,iZ) )

          END DO

          DO i = 1, 4

            FVECm(i,iZ) = GVECm(i,iZ) - UVEC(i,iZ)

            GVEC(i,Mk,iZ) = GVECm(i,iZ)
            FVEC(i,Mk,iZ) = FVECm(i,iZ)

          END DO

        END IF
      END DO



      IF ( Mk > 1 ) THEN

        CALL Alpha_LS_Vector &
               ( ITERATE, nZ, M, Mk, FVECm, FVEC, Alpha )

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC &
      !$ACC PRIVATE( SUM1 ) &
      !$ACC PRESENT( ITERATE, GVECm, FVECm, GVEC, UVEC, Alpha )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#endif
        DO iZ = 1, nZ
          DO i = 1, 4
            IF ( ITERATE(iZ) ) THEN
              SUM1 = Zero
              DO iM = 1, Mk
                SUM1 = SUM1 + GVEC(i,iM,iZ) * Alpha(iM,iZ)
              END DO
              GVECm(i,iZ) = SUM1
              FVECm(i,iZ) = GVECm(i,iZ) - UVEC(i,iZ)
            END IF
          END DO
        END DO
      END IF



#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
    !$ACC PRIVATE( iX, CONVERGED, FTMP, GTMP ) &
    !$ACC PRESENT( ITERATE, UVEC, CVEC, GVECm, FVECm, GVEC, FVEC, &
    !$ACC          PositionIndexZ, D, I_u_1, I_u_2, I_u_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX, CONVERGED, FTMP, GTMP )
#endif
      DO iZ = 1, nZ
        IF ( ITERATE(iZ) ) THEN

          iX = PositionIndexZ(iZ)

          D    (iZ) = GVECm(iPR_D ,iZ)
          I_u_1(iZ) = GVECm(iPR_I1,iZ) / Gm_dd_11(iX)
          I_u_2(iZ) = GVECm(iPR_I2,iZ) / Gm_dd_22(iX)
          I_u_3(iZ) = GVECm(iPR_I3,iZ) / Gm_dd_33(iX)

          ! CONVERGED = SQRT( SUM( FVECm(:,iZ)**2 ) ) <= &
          !                        Rtol * SQRT( SUM( CVEC(:,iZ)**2 ) )
          CONVERGED = ALL( ABS( FVECm(:,iZ) ) <= Rtol * ABS( CVEC(:,iZ )) )

          IF ( CONVERGED ) THEN
            ITERATE(iZ) = .FALSE.
            nIterations(iZ) = k
          ELSE IF ( Mk == M ) THEN
            DO j = 1, Mk - 1
              DO i = 1, 4
                FTMP(i,j) = FVEC(i,j+1,iZ)
                GTMP(i,j) = GVEC(i,j+1,iZ)
              END DO
            END DO
            DO j = 1, Mk - 1
              DO i = 1, 4
                FVEC(i,j,iZ) = FTMP(i,j)
                GVEC(i,j,iZ) = GTMP(i,j)
              END DO
            END DO
          END IF
        END IF
      END DO
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( ITERATE )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( ITERATE ) ASYNC
    !$ACC WAIT
#endif


    END DO


    IF( PRESENT( nIterations_Option ) ) THEN
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
    !$ACC PRESENT( nIterations, nIterations_Option )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iZ = 1, nZ
        nIterations_Option(iZ) = nIterations(iZ)
      END DO
    END IF


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
    !$ACC PRESENT( dN, dG_d_1, dG_d_2, dG_d_3, &
    !$ACC          GVECm, D_0, Chi, Kappa )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
        DO iZ = 1, nZ
          dN(iZ)     = Chi(iZ) * ( D_0(iZ) - GVECm(iPR_D ,iZ) )
          dG_d_1(iZ) = - Kappa(iZ) * GVECm(iPR_I1 ,iZ)
          dG_d_2(iZ) = - Kappa(iZ) * GVECm(iPR_I2 ,iZ)
          dG_d_3(iZ) = - Kappa(iZ) * GVECm(iPR_I3 ,iZ)
        END DO

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: FVEC, GVEC, CVEC, UVEC, &
  !$OMP               FVECm, GVECm, Alpha, ITERATE, nIterations, &
  !$OMP               D, I_u_1, I_u_2, I_u_3, &
  !$OMP               Kappa, D_00, D_ii )
#elif defined( THORNADO_OACC   )
  !$ACC EXIT DATA WAIT &
  !$ACC DELETE( FVEC, GVEC, CVEC, UVEC, &
  !$ACC         FVECm, GVECm, Alpha, ITERATE, nIterations, &
  !$ACC         D, I_u_1, I_u_2, I_u_3, &
  !$ACC         Kappa, D_00, D_ii )
#endif

    DEALLOCATE( FVEC, GVEC )
    DEALLOCATE( CVEC, UVEC, FVECm, GVECm, Alpha )
    DEALLOCATE( D, I_u_1, I_u_2, I_u_3 )
    DEALLOCATE( D_00, D_ii, Kappa )
    DEALLOCATE( ITERATE, nIterations )


  END SUBROUTINE ComputeIncrement_FixedPoint_Vector_Richardson


  SUBROUTINE Alpha_LS_Vector &
    ( MASK, nZ, M, Mk, Fm, F, Alpha )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: nZ, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA
    REAL(DP) :: A1, A2, B
    INTEGER  :: iZ, iP

    IF ( Mk > 1 ) THEN

      IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AB1, B )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
        !$ACC PRIVATE( AA11, AB1, B ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AB1, B )
#endif
        DO iZ = 1, nZ
          IF ( MASK(iZ) ) THEN

            AA11 = Zero
            AB1 = Zero

            DO iP = 1, 4

              A1 = F(iP,1,iZ) - Fm(iP,iZ)
              B  = - Fm(iP,iZ)

              AA11 = AA11 + A1 * A1
              AB1  = AB1  + A1 * B

            END DO

            Alpha(1,iZ) = AB1 / ( AA11 + SqrtTiny )
            Alpha(2,iZ) = One - Alpha(1,iZ)

          END IF
        END DO

      ELSE IF ( Mk == 3 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR ASYNC &
        !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B ) &
        !$ACC PRESENT( MASK, Alpha, F, Fm )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA, A1, A2, B )
#endif
        DO iZ = 1, nZ
          IF ( MASK(iZ) ) THEN

            AA11 = Zero
            AA12 = Zero
            AA22 = Zero
            AB1  = Zero
            AB2  = Zero

            DO iP = 1, 4

              A1 = F(iP,1,iZ) - Fm(iP,iZ)
              A2 = F(iP,2,iZ) - Fm(iP,iZ)
              B  = - Fm(iP,iZ)

              AA11 = AA11 + A1 * A1
              AA12 = AA12 + A1 * A2
              AA22 = AA22 + A2 * A2

              AB1  = AB1  + A1 * B
              AB2  = AB2  + A2 * B

            END DO

            DET_AA = AA11*AA22 - AA12*AA12

            Alpha(1,iZ) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
            Alpha(2,iZ) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA
            Alpha(3,iZ) = One - Alpha(1,iZ) - Alpha(2,iZ)

          END IF
        END DO

      ELSE IF ( Mk > 3 ) THEN

        ! --- Not Implemented ---

      END IF

    END IF

  END SUBROUTINE Alpha_LS_Vector


  FUNCTION Alpha_LS( M, mk, FVEC )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: FVEC(4,M)
    REAL(DP)             :: Alpha_LS(M)

    INTEGER  :: i
    REAL(DP) :: BVEC(4), AMAT(4,M)
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

    BVEC = - FVEC(:,mk)

    DO i = 1, mk - 1

      AMAT(:,i) = FVEC(:,i) - FVEC(:,mk)

    END DO

    IF( mk == 2 )THEN

      AA11 = Zero
      AB1  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)

      END DO

      BVEC(1) = AB1 / AA11

    ELSEIF( mk == 3 )THEN

      AA11 = Zero
      AA12 = Zero
      AA22 = Zero
      AB1  = Zero
      AB2  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AA12 = AA12 + AMAT(i,1) * AMAT(i,2)
        AA22 = AA22 + AMAT(i,2) * AMAT(i,2)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)
        AB2  = AB2  + AMAT(i,2) * BVEC(i)

      END DO

      DET_AA = AA11 * AA22 - AA12 * AA12

      BVEC(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
      BVEC(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

    ELSEIF( mk > 3 )THEN

      ! --- Not Implemented ---

    END IF

    SUM1 = Zero
    DO i = 1, mk - 1

      Alpha_LS(i) = BVEC(i)

      SUM1 = SUM1 + BVEC(i)

    END DO

    Alpha_LS(mk) = One - SUM1

    RETURN
  END FUNCTION Alpha_LS


  FUNCTION ShiftVec( M, mk, Vec )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: M, mk
    REAL(DP), INTENT(in) :: Vec(4,M)
    REAL(DP)             :: ShiftVec(4,M)

    INTEGER  :: i, j
    REAL(DP) :: VecTMP(4,M)

    DO j = 1, mk - 1
    DO i = 1, 4

      VecTMP(i,j) = Vec(i,j+1)

    END DO
    END DO

    DO j = 1, mk - 1
    DO i = 1, 4

      ShiftVec(i,j) = VecTMP(i,j)

    END DO
    END DO

    RETURN
  END FUNCTION ShiftVec




  SUBROUTINE SolveAlpha_LS( M, mk, FVEC, Alpha )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: M, mk
    REAL(DP), INTENT(inout) :: FVEC(4,M), Alpha(M)

    INTEGER  :: i
    REAL(DP) :: BVEC(4), AMAT(4,M)
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

    BVEC = - FVEC(:,mk)

    DO i = 1, mk - 1

      AMAT(:,i) = FVEC(:,i) - FVEC(:,mk)

    END DO

    IF( mk == 2 )THEN

      AA11 = Zero
      AB1  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)

      END DO

      BVEC(1) = AB1 / AA11

    ELSEIF( mk == 3 )THEN

      AA11 = Zero
      AA12 = Zero
      AA22 = Zero
      AB1  = Zero
      AB2  = Zero

      DO i = 1, 4

        AA11 = AA11 + AMAT(i,1) * AMAT(i,1)
        AA12 = AA12 + AMAT(i,1) * AMAT(i,2)
        AA22 = AA22 + AMAT(i,2) * AMAT(i,2)
        AB1  = AB1  + AMAT(i,1) * BVEC(i)
        AB2  = AB2  + AMAT(i,2) * BVEC(i)

      END DO

      DET_AA = AA11 * AA22 - AA12 * AA12

      BVEC(1) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
      BVEC(2) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

    ELSEIF( mk > 3 )THEN

      PRINT*, "mk > 3"

      STOP

    END IF

    SUM1 = Zero
    DO i = 1, mk - 1

      Alpha(i) = BVEC(i)

      SUM1 = SUM1 + BVEC(i)

    END DO

    Alpha(mk) = One - SUM1

  END SUBROUTINE SolveAlpha_LS


  SUBROUTINE ShiftVectors( M, mk, FVEC, GVEC )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: M, mk
    REAL(DP), INTENT(inout) :: FVEC(4,M), GVEC(4,M)

    INTEGER  :: i, j
    REAL(DP) :: FTMP(4,M), GTMP(4,M)

    DO j = 1, mk - 1
    DO i = 1, 4

      FTMP(i,j) = FVEC(i,j+1)
      GTMP(i,j) = GVEC(i,j+1)

    END DO
    END DO

    DO j = 1, mk - 1
    DO i = 1, 4

      FVEC(i,j) = FTMP(i,j)
      GVEC(i,j) = GTMP(i,j)

    END DO
    END DO

  END SUBROUTINE ShiftVectors

  FUNCTION EddingtonTensorComponents_dd &
    ( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: &
      D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP)              :: &
      EddingtonTensorComponents_dd(3,3)

    INTEGER  :: i, j
    REAL(DP) :: FF, EF, a, b
    REAL(DP) :: h_d(3), Gm_dd(3,3)

    FF = FluxFactor( D, I_u_1, I_u_2, I_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    EF = EddingtonFactor( D, FF )

    a = Half * ( One - EF )
    b = Half * ( Three * EF - One )

    h_d(1) = Gm_dd_11 * I_u_1 / ( FF * D )
    h_d(2) = Gm_dd_22 * I_u_2 / ( FF * D )
    h_d(3) = Gm_dd_33 * I_u_3 / ( FF * D )

    Gm_dd = Zero
    Gm_dd(1,1) = Gm_dd_11
    Gm_dd(2,2) = Gm_dd_22
    Gm_dd(3,3) = Gm_dd_33

    DO j = 1, 3
    DO i = 1, 3

      EddingtonTensorComponents_dd(i,j) &
        = a * Gm_dd(i,j) + b * h_d(i) * h_d(j)

    END DO
    END DO

    RETURN
  END FUNCTION EddingtonTensorComponents_dd

  SUBROUTINE InitializeCollisions( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iS, iN_E, iN_X, iZ, nZ_G

    iE_B0 = iZ_B0(1);   iE_E0 = iZ_E0(1)
    iX_B0 = iZ_B0(2:4); iX_E0 = iZ_E0(2:4)

    nZ = iZ_E0 - iZ_B0 + 1
    nE = iE_E0 - iE_B0 + 1
    nX = iX_E0 - iX_B0 + 1

    nE_G = nDOFE * nE
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( GX_N(nX_G,nGF) )
    ALLOCATE( PF_N(nX_G,nPF) )
    ALLOCATE( CF_N(nX_G,nCF) )
    ALLOCATE( CR_N (nE_G,nX_G,nSpecies,nCR) )
    ALLOCATE( dCR_N(nE_G,nX_G,nSpecies,nCR) )
    ALLOCATE( OP_N (nE_G,nX_G,nSpecies,nOP) )

    nZ_G = nE_G * nX_G * nSpecies

    ALLOCATE( PositionIndexZ     (nZ_G) )
    ALLOCATE( nIterations_Simple (nZ_G) )

    iZ = 0
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      iZ = iZ + 1

      PositionIndexZ(iZ) = iN_X

    END DO
    END DO
    END DO


    N_PTR    (1:nZ_G) => CR_N  (:,:,:,iCR_N     )
    G1_PTR   (1:nZ_G) => CR_N  (:,:,:,iCR_G1    )
    G2_PTR   (1:nZ_G) => CR_N  (:,:,:,iCR_G2    )
    G3_PTR   (1:nZ_G) => CR_N  (:,:,:,iCR_G3    )
    D0_PTR   (1:nZ_G) => OP_N  (:,:,:,iOP_D0    )
    Chi_PTR  (1:nZ_G) => OP_N  (:,:,:,iOP_Chi   )
    Sigma_PTR(1:nZ_G) => OP_N  (:,:,:,iOP_Sigma )
    dN_PTR   (1:nZ_G) => dCR_N (:,:,:,iCR_N     )
    dG1_PTR  (1:nZ_G) => dCR_N (:,:,:,iCR_G1    )
    dG2_PTR  (1:nZ_G) => dCR_N (:,:,:,iCR_G2    )
    dG3_PTR  (1:nZ_G) => dCR_N (:,:,:,iCR_G3    )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, nZ, nX, PositionIndexZ ) &
    !$OMP MAP( alloc: GX_N, CF_N, PF_N, CR_N, dCR_N, OP_N, nIterations_Simple )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, nZ, nX, PositionIndexZ ) &
    !$ACC CREATE( GX_N, CF_N, PF_N, CR_N, dCR_N, OP_N, nIterations_Simple )
#endif

  END SUBROUTINE InitializeCollisions


  SUBROUTINE FinalizeCollisions

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iX_B0, iX_E0, nZ, nX, PositionIndexZ, &
    !$OMP               GX_N, CF_N, PF_N, CR_N, dCR_N, OP_N, nIterations_Simple )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( iX_B0, iX_E0, nZ, nX, PositionIndexZ, &
    !$ACC         GX_N, CF_N, PF_N, CR_N, dCR_N, OP_N, nIterations_Simple )
#endif

  DEALLOCATE( GX_N, PF_N, CF_N, CR_N, dCR_N, OP_N )
  DEALLOCATE( PositionIndexZ, nIterations_Simple )


  NULLIFY( N_PTR, G1_PTR, G2_PTR, G3_PTR )
  NULLIFY( D0_PTR, Chi_PTR, Sigma_PTR )
  NULLIFY( dN_PTR, dG1_PTR, dG2_PTR, dG3_PTR )




  END SUBROUTINE FinalizeCollisions


END MODULE TwoMoment_DiscretizationModule_Collisions_OrderV
