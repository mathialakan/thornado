MODULE MF_TwoMoment_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1
  USE RadiationFieldsModule,            ONLY: &
    nCR, nPR, iCR_G1, iCR_G2, iCR_G3, iCR_N, &
    iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE FluidFieldsModule,            ONLY: &
    nCF, nPF, iCF_D, iCF_E, iCF_Ne, iCF_S1, iCF_S2, iCF_S3, &
    iPF_D, iPF_E, iPF_Ne, iPF_V1, iPF_V2, iPF_V3
  USE GeometryFieldsModule,     ONLY: &
    nGF, iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3, iGF_Gm_dd_11, &
    iGF_Gm_dd_22, iGF_Gm_dd_33
  USE TwoMoment_UtilitiesModule_Relativistic,  ONLY: &
    ComputePrimitive_TwoMoment
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels, nSpecies, nE
  USE MF_UtilitiesModule,                ONLY: &
    AMReX2thornado_Euler, &
    AMReX2thornado,       &
    thornado2AMReX,       &
    thornado2AMReX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeTimeStep_Fancy
  PUBLIC :: MF_ComputeTimeStep
  PUBLIC :: MF_ComputeFromConserved

CONTAINS




  SUBROUTINE MF_ComputeTimeStep_Fancy( MF_uGF, MF_uCR, CFL, TimeStepMin )

    TYPE(amrex_multifab),  INTENT(in)  :: MF_uGF(0:nLevels-1), MF_uCR(0:nLevels-1)
    REAL(amrex_real),     INTENT(in)  :: CFL
    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR(:,:,:,:)


    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:,:,:)

    REAL(amrex_real) :: TimeStep(0:nLevels-1)
    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
 
    TimeStepMin = HUGE( 1.0e0_amrex_real )    

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )


        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0(1:3) = BX % lo(1:3)
        iX_E0(1:3) = BX % hi(1:3)
        iX_B1(1:3) = BX % lo(1:3) - swX(1:3)
        iX_E1(1:3) = BX % hi(1:3) + swX(1:3)


        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U (1:nDOFZ,iE_B1:iE_E1,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCR,1:nSpecies) )

        CALL AMReX2thornado_Euler &
               ( nGF, iX_B1(1:3), iX_E1(1:3), &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iX_B1, iX_E1,                    &
                 uCR(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFZ*nCR*nSpecies*nE), &
                 U(1:nDOFZ,iE_B0:iE_E0, &
                           iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCR,1:nSpecies) )

!put compute timestep here

        TimeStep( iLevel ) = 0.01_amrex_real
 
        TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

        DEALLOCATE( U )
        DEALLOCATE( G )



      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

    END DO ! --- Loop over levels ---



  END SUBROUTINE MF_ComputeTimeStep_Fancy

  SUBROUTINE MF_ComputeTimeStep( nX, xR, xL, nNodes, CFL, TimeStepMin )

    INTEGER,              INTENT(in)  :: nX(:), nNodes
    REAL(amrex_real),     INTENT(in)  :: xR(:), xL(:), CFL
    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels-1)


    INTEGER          :: iLevel
    REAL(amrex_real) :: TimeStep(0:nLevels-1)


    TimeStepMin = HUGE( 1.0e0_amrex_real )    

    DO iLevel = 0, nLevels-1

      TimeStep( iLevel ) = CFL * MINVAL( (xR-xL) / DBLE(nX) ) &
                           / ( 2.0_amrex_real * DBLE(nNodes-1) + 1.0_amrex_real )

      TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

    END DO






  END SUBROUTINE MF_ComputeTimeStep

  SUBROUTINE MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

    TYPE(amrex_multifab), INTENT(in   ) :: &
      MF_uGF(0:nLevels-1), MF_uCR(0:nLevels-1), MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uPR(0:nLevels-1)


    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uPR(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: CF(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: PF(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: CR(:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: PR(:,:,:,:,:,:,:)

    INTEGER, ALLOCATABLE :: nIterations(:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), nPoints
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), i, iX1, iX2, iX3, iS, iE, iNodeX, iNodeZ, iPoint


    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        uPR  => MF_uPR (iLevel) % DataPtr( MFI )
        
        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX
 
        i=1          

        DO WHILE (i<=4)
          
          IF (i==1) THEN
          
            iZ_B0(i)=iE_B0 
            iZ_E0(i)=iE_E0
            iZ_B1(i)=iE_B1
            iZ_E1(i)=iE_E1

          ELSE 

            iZ_B0(i)=iX_B0(i-1) 
            iZ_E0(i)=iX_E0(i-1)
            iZ_B1(i)=iX_B1(i-1)
            iZ_E1(i)=iX_E1(i-1)

          END IF
          i = i + 1 
        END DO

        nPoints = nSpecies * PRODUCT( iX_E0 - iX_B0 + 1 ) &
              * ( iE_E0 - iE_B0 + 1 ) * nDOFZ

        ALLOCATE( nIterations(nPoints) )

        iPoint = 0
 
        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( CF (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( PF (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nPF) )

        ALLOCATE( CR (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        ALLOCATE( PR (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )


        CALL AMReX2thornado_Euler &
               ( nGF, iX_B1, iX_E1, &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado_Euler &
               ( nCF, iX_B1, iX_E1, &
                 uCF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nCF), &
                 CF(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF) )

        CALL AMReX2thornado &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iX_B1, iX_E1,                    &
                 uCR(      iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nDOFZ*nCR*nSpecies*nE), &
                 CR(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                           iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        CALL AMReX2thornado &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iX_B1, iX_E1,                    &
                 uPR(      iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nDOFZ*nPR*nSpecies*nE), &
                 PR(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                           iZ_B1(2):iZ_E1(2), &
                           iZ_B1(3):iZ_E1(3), &
                           iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )



        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputePrimitive_Euler_Relativistic &
               ( CF(1:nDOFX,iX1,iX2,iX3,iCF_D),         &
                 CF(1:nDOFX,iX1,iX2,iX3,iCF_S1),        &
                 CF(1:nDOFX,iX1,iX2,iX3,iCF_S2),        &
                 CF(1:nDOFX,iX1,iX2,iX3,iCF_S3),        &
                 CF(1:nDOFX,iX1,iX2,iX3,iCF_E),         &
                 CF(1:nDOFX,iX1,iX2,iX3,iCF_Ne),        &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_D),         &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_V1),        &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_V2),        &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_V3),        &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_E),         &
                 PF(1:nDOFX,iX1,iX2,iX3,iPF_Ne),        &
                 G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 G(1:nDOFX,iX1,iX2,iX3,iGF_Gm_dd_33) )
        END DO
        END DO
        END DO


        DO iS  = 1, nSpecies
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

        DO iE = iE_B0, iE_E0
          
          DO iNodeZ = 1, nDOFZ

            iPoint = iPoint + 1

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1


            CALL ComputePrimitive_TwoMoment &
                 ( CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                   CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                   CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                   CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                   PR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                   PR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                   PR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                   PR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                   PF(iNodeX,iX1,iX2,iX3,iPF_V1),     &
                   PF(iNodeX,iX1,iX2,iX3,iPF_V2),     &
                   PF(iNodeX,iX1,iX2,iX3,iPF_V3),     &
                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &                   
                   0.0_amrex_real, 0.0_amrex_real, 0.0_amrex_real, &
                   G(iNodeX,iX1,iX2,iX3,iGF_Alpha), &
                   G(iNodeX,iX1,iX2,iX3,iGF_Beta_1), &
                   G(iNodeX,iX1,iX2,iX3,iGF_Beta_2), &
                   G(iNodeX,iX1,iX2,iX3,iGF_Beta_3), &
                   nIterations(iPoint) )
        


          END DO

        END DO

      END DO
      END DO
      END DO
      END DO



        CALL thornado2AMReX &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, iX_B0, iX_E0, &
                 uPR(     iZ_B0(2):iZ_E0(2), &
                           iZ_B0(3):iZ_E0(3), &
                           iZ_B0(4):iZ_E0(4),1:nDOFZ*nPR*nSpecies*nE), &
                 PR(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                            iZ_B0(2):iZ_E0(2), &
                            iZ_B0(3):iZ_E0(3), &
                            iZ_B0(4):iZ_E0(4),1:nPR,1:nSpecies) )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO






  END SUBROUTINE MF_ComputeFromConserved



END MODULE MF_TwoMoment_UtilitiesModule
