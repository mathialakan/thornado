MODULE RefluxModule_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF
  USE MeshModule, ONLY: &
    MeshX

  ! --- Local Modules ---

  USE MF_FieldsModule, ONLY: &
    FluxRegister
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE AverageDownModule, ONLY: &
    AverageDownTo
  USE MF_UtilitiesModule, ONLY: &
    MultiplyWithMetric
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Reflux_Euler_MF

  INTERFACE Reflux_Euler_MF
    MODULE PROCEDURE Reflux_Euler_MF_SingleLevel
    MODULE PROCEDURE Reflux_Euler_MF_MultipleLevels
  END INTERFACE Reflux_Euler_MF

CONTAINS


  SUBROUTINE Reflux_Euler_MF_MultipleLevels( MF_uGF, MF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      IF( iLevel .GT. 0 ) &
        CALL Reflux_Euler_MF_SingleLevel( iLevel, MF_uGF, MF )

    END DO

  END SUBROUTINE Reflux_Euler_MF_MultipleLevels


  SUBROUTINE Reflux_Euler_MF_SingleLevel( FineLevel, MF_uGF, MF )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    CALL amrex_multifab_build &
           ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                  MF_uGF(FineLevel-1) % DM, nDOFX, swX )
    CALL SqrtGm(FineLevel-1) % COPY &
           ( MF_uGF(FineLevel-1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL amrex_multifab_build &
           ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                  MF_uGF(FineLevel  ) % DM, nDOFX, swX )
    CALL SqrtGm(FineLevel  ) % COPY &
           ( MF_uGF(FineLevel  ), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

    CALL CreateMesh_MF( FineLevel-1, MeshX )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

! todo: replace MF_uGF(FineLevel-1) with SqrtGm(FineLevel-1)
!       and remove iGF_SqrtGm from inputs to amrex
    CALL FluxRegister( FineLevel ) &
           % reflux_dg( MF_uGF(FineLevel-1), MF(FineLevel-1), &
                        nCF, dX1, dX2, dX3 )

    END ASSOCIATE

    CALL DestroyMesh_MF( MeshX )

    CALL MultiplyWithMetric( SqrtGm(FineLevel), MF_uGF(FineLevel), nGF, +1 )
    CALL MultiplyWithMetric( SqrtGm(FineLevel), MF    (FineLevel), nCF, +1 )

    CALL AverageDownTo( FineLevel-1, MF_uGF )
    CALL AverageDownTo( FineLevel-1, MF     )

    CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ), nGF, -1 )
    CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF    (FineLevel  ), nCF, -1 )
    CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1), nGF, -1 )
    CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF    (FineLevel-1), nCF, -1 )

    CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
    CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

  END SUBROUTINE Reflux_Euler_MF_SingleLevel


END MODULE RefluxModule_Euler
