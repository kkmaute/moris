/*
 * cl_DLA_Linear_Solver_PETSc.cpp
 *
 *  Created on: Dez 12, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Preconditioner_PETSc.hpp"
#include "cl_DLA_Linear_Solver_PETSc.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "fn_inv.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

using namespace moris;
using namespace dla;

Preconditioner_PETSc::Preconditioner_PETSc( Linear_Solver_PETSc *aLinearSolverAlgoritm )
        : mLinearSolverAlgoritm( aLinearSolverAlgoritm )
{
}

//----------------------------------------------------------------------------------------

void
Preconditioner_PETSc::build_multigrid_preconditioner( Linear_Problem *aLinearSystem )
{
    //----------------------------------------------------------------
    //                 Set general multigrid stuff
    //----------------------------------------------------------------
    // Build multigrid operators
    mLinearSolverAlgoritm->mSolverInterface->build_multigrid_operators();

    // get multigrid operators
    moris::Cell< sol::Dist_Matrix * > tProlongationList = mLinearSolverAlgoritm->mSolverInterface
                                                                  ->get_multigrid_operator_pointer()
                                                                  ->get_prolongation_list();

    sint tLevels = mLinearSolverAlgoritm->mParameterList.get< moris::sint >( "MultigridLevels" );

    // set multigrid levels and type
    PCMGSetLevels( mLinearSolverAlgoritm->mpc, tLevels, NULL );
    PCMGSetType( mLinearSolverAlgoritm->mpc, PC_MG_MULTIPLICATIVE );
    PCMGSetGalerkin( mLinearSolverAlgoritm->mpc, PC_MG_GALERKIN_BOTH );

    // get restriction operators
    moris::Cell< Mat > tTransposeOperators( tProlongationList.size() );
    for ( moris::uint Ik = 0; Ik < tProlongationList.size(); Ik++ )
    {
        MatTranspose( tProlongationList( Ik )->get_petsc_matrix(), MAT_INITIAL_MATRIX, &tTransposeOperators( Ik ) );

        if ( Ik == 0 )
        {
            std::cout << "writing interpolation operator" << std::endl;
            std::string tString = "Interpolation_Operator";
            tProlongationList( Ik )->save_matrix_to_matlab_file( tString.c_str() );
        }
    }

    // set restriction operators
    moris::sint tCounter = 0;
    for ( moris::sint Ik = tLevels - 1; Ik > 0; Ik-- )
    {
        PCMGSetInterpolation( mLinearSolverAlgoritm->mpc, Ik, tTransposeOperators( tCounter ) );
        //         PCMGSetRestriction( mLinearSolverAlgoritm->mpc, Ik, tTransposeOperators( tCounter ) );

        tCounter++;
    }

    // destroy original operators
    for ( moris::uint Ik = 0; Ik < tTransposeOperators.size(); Ik++ )
    {
        MatDestroy( &tTransposeOperators( Ik ) );
    }

    //----------------------------------------------------------------
    //                 Set coarsest level settings
    //----------------------------------------------------------------
    KSP tPetscKSPCoarseSolve;
    PCMGGetCoarseSolve( mLinearSolverAlgoritm->mpc, &tPetscKSPCoarseSolve );
    KSPSetType( tPetscKSPCoarseSolve, KSPRICHARDSON );
    PetscInt maxitss = 1000;
    KSPSetTolerances( tPetscKSPCoarseSolve, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxitss );
    PC cpc;
    KSPGetPC( tPetscKSPCoarseSolve, &cpc );
    PCSetType( cpc, PCLU );

    //----------------------------------------------------------------
    //                 Set schwarz preconditioner for finest level
    //----------------------------------------------------------------

    if ( mLinearSolverAlgoritm->mParameterList.get< bool >( "MG_use_schwarz_smoother" ) )
    {
        // set schwarz preconditiiner domains based on criteria
        moris::real tVolumeFractionThreshold =    //
                mLinearSolverAlgoritm->mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );
        moris::sint tSchwarzSmoothingIters =    //
                mLinearSolverAlgoritm->mParameterList.get< moris::sint >( "MG_schwarz_smoothing_iters" );

        moris::Cell< moris::Matrix< IdMat > > tCriteriaIds;
        mLinearSolverAlgoritm->mSolverInterface->get_adof_ids_based_on_criteria( tCriteriaIds,
                tVolumeFractionThreshold );

        //---------------------------------------------------------------

        uint tNumSerialDofs = mLinearSolverAlgoritm->mSolverInterface->get_num_my_dofs();

        moris::Matrix< DDSMat > tMat( tNumSerialDofs, 1, 0 );

        uint tNumBlocksInitialBlocks = tCriteriaIds.size();

        for ( uint Ik = 0; Ik < tNumBlocksInitialBlocks; Ik++ )
        {
            for ( uint Ii = 0; Ii < tCriteriaIds( Ik ).numel(); Ii++ )
            {
                if ( tMat( tCriteriaIds( Ik )( Ii ) ) == 0 )
                {
                    tMat( tCriteriaIds( Ik )( Ii ) ) = 1;
                }
            }
        }

        moris::Matrix< DDSMat > tNotInBlock( tNumSerialDofs, 1, 0 );
        tCounter = 0;
        for ( uint Ik = 0; Ik < tMat.numel(); Ik++ )
        {
            if ( tMat( Ik ) == 0 )
            {
                tNotInBlock( tCounter ) = Ik;
                tCounter++;
            }
        }

        tNotInBlock.resize( tCounter, 1 );

        tCriteriaIds.resize( tCriteriaIds.size() + tCounter );

        for ( uint Ik = 0; Ik < tNotInBlock.numel(); Ik++ )
        {
            tCriteriaIds( Ik + tNumBlocksInitialBlocks ).set_size( 1, 1, tNotInBlock( Ik ) );
        }

        //---------------------------------------------------------------

        KSP dKSPFirstDown;
        PCMGGetSmootherDown( mLinearSolverAlgoritm->mpc, tLevels - 1, &dKSPFirstDown );
        PCMGGetSmootherUp( mLinearSolverAlgoritm->mpc, tLevels - 1, &dKSPFirstDown );
        PC dPCFirstDown;
        KSPGetPC( dKSPFirstDown, &dPCFirstDown );
        KSPSetType( dKSPFirstDown, KSPRICHARDSON );
        moris::sint restart = 1;
        KSPGMRESSetRestart( dKSPFirstDown, restart );
        KSPSetTolerances( dKSPFirstDown, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, tSchwarzSmoothingIters );
        PCSetType( dPCFirstDown, PCASM );

        PCASMSetLocalType( dPCFirstDown, PC_COMPOSITE_MULTIPLICATIVE );
        PCASMSetType( dPCFirstDown, PC_ASM_BASIC );
        PCASMSetOverlap( dPCFirstDown, 0 );

        //-------------------------------------------------------------------------

        KSPSetOperators( dKSPFirstDown, aLinearSystem->get_matrix()->get_petsc_matrix(), aLinearSystem->get_matrix()->get_petsc_matrix() );


        uint tNumBlocks = tCriteriaIds.size();

        moris::Cell< IS > tIs( tNumBlocks );

        for ( uint Ik = 0; Ik < tNumBlocks; Ik++ )
        {
            ISCreateGeneral( PETSC_COMM_WORLD, tCriteriaIds( Ik ).numel(), tCriteriaIds( Ik ).data(), PETSC_COPY_VALUES, &tIs( Ik ) );
            //             ISView( tIs( Ik ) , PETSC_VIEWER_STDOUT_SELF );
        }

        PCASMSetLocalSubdomains( dPCFirstDown, tNumBlocks, NULL, tIs.data().data() );
        std::cout << "build-----" << std::endl;
        KSPSetUp( dKSPFirstDown );

        KSP *tSubksp;

        sint nlocal;
        sint first;

        PCASMGetSubKSP( dPCFirstDown, &nlocal, &first, &tSubksp );

        //         tKSPBlock.resize( tNumBlocks );
        Cell< PC > tPCBlock( nlocal );
        for ( sint Ik = 0; Ik < nlocal; Ik++ )
        {
            KSPSetType( tSubksp[ Ik ], KSPRICHARDSON );
            KSPGetPC( tSubksp[ Ik ], &tPCBlock( Ik ) );
            PCSetType( tPCBlock( Ik ), PCLU );
            PCSetUp( tPCBlock( Ik ) );

            //             KSPView( tSubksp[ Ik ], PETSC_VIEWER_STDOUT_SELF );
            //             PCView( tPCBlock( Ik ), PETSC_VIEWER_STDOUT_SELF );
        }

        std::cout << "----- ksps build -----" << std::endl;

        //         PCView( dPCFirstDown, PETSC_VIEWER_STDOUT_SELF );
    }
    else
    {
        KSP dkspDown1;
        PCMGGetSmootherDown( mLinearSolverAlgoritm->mpc, tLevels - 1, &dkspDown1 );
        PC dpcDown;
        KSPGetPC( dkspDown1, &dpcDown );
        KSPSetType( dkspDown1, KSPRICHARDSON );    // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
        moris::sint restart = 1;
        moris::sint tMaxit  = 1;
        KSPGMRESSetRestart( dkspDown1, restart );
        KSPSetTolerances( dkspDown1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, tMaxit );    // NOTE maxitr=restart;
        PCSetType( dpcDown, PCSOR );
        PCSORSetOmega( dpcDown, 0.8 );
        PCSORSetIterations( dpcDown, 1, 1 );

        KSP dkspUp1;
        PCMGGetSmootherUp( mLinearSolverAlgoritm->mpc, tLevels - 1, &dkspUp1 );
        PC dpcUp;
        KSPGetPC( dkspUp1, &dpcUp );
        KSPSetType( dkspUp1, KSPRICHARDSON );    // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
        KSPGMRESSetRestart( dkspUp1, restart );
        KSPSetTolerances( dkspUp1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, tMaxit );    // NOTE maxitr=restart;
        PCSetType( dpcUp, PCSOR );
        PCSORSetOmega( dpcUp, 0.8 );
        PCSORSetIterations( dpcUp, 1, 1 );
    }

    //----------------------------------------------------------------
    //                 Set jacobi smoother for coarser levels
    //----------------------------------------------------------------
    for ( PetscInt Ik = 1; Ik < tLevels - 1; Ik++ )
    {
        KSP dkspDown;
        PCMGGetSmootherDown( mLinearSolverAlgoritm->mpc, Ik, &dkspDown );
        PC dpcDown;
        KSPGetPC( dkspDown, &dpcDown );
        KSPSetType( dkspDown, KSPRICHARDSON );    // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
        moris::sint restart = 2;
        KSPGMRESSetRestart( dkspDown, restart );
        KSPSetTolerances( dkspDown, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, restart );    // NOTE maxitr=restart;
        PCSetType( dpcDown, PCSOR );                                                           // PCJACOBI, PCSOR for KSPCHEBYSHEV very good... Use KSPRICHARDSON for weighted Jacobi
    }

    for ( PetscInt k = 1; k < tLevels - 1; k++ )
    {
        KSP dkspUp;

        PCMGGetSmootherUp( mLinearSolverAlgoritm->mpc, k, &dkspUp );
        PC dpcUp;
        KSPGetPC( dkspUp, &dpcUp );
        KSPSetType( dkspUp, KSPRICHARDSON );
        //            KSPSetType(dkspUp,KSPGMRES);                                                                 // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
        moris::sint restart = 2;
        KSPGMRESSetRestart( dkspUp, restart );
        KSPSetTolerances( dkspUp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, restart );    // NOTE maxitr=restart;
        PCSetType( dpcUp, PCSOR );
        //         PCSetType(dpcUp,PCJACOBI);
    }
}

void
Preconditioner_PETSc::build_schwarz_preconditioner_petsc()
{
    // set schwarz preconditiiner domains based on criteria
    moris::real tVolumeFractionThreshold = mLinearSolverAlgoritm->mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );

    moris::Cell< moris::Matrix< IdMat > > tCriteriaIds;
    mLinearSolverAlgoritm->mSolverInterface->get_adof_ids_based_on_criteria( tCriteriaIds,
            tVolumeFractionThreshold );


    uint tNumSerialDofs = mLinearSolverAlgoritm->mSolverInterface->get_num_my_dofs();

    moris::Matrix< DDSMat > tMat( tNumSerialDofs, 1, 0 );

    uint tNumBlocksInitialBlocks = tCriteriaIds.size();

    for ( uint Ik = 0; Ik < tNumBlocksInitialBlocks; Ik++ )
    {
        for ( uint Ii = 0; Ii < tCriteriaIds( Ik ).numel(); Ii++ )
        {
            if ( tMat( tCriteriaIds( Ik )( Ii ) ) == 0 )
            {
                tMat( tCriteriaIds( Ik )( Ii ) ) = 1;
            }
        }
    }

    moris::Matrix< DDSMat > tNotInBlock( tNumSerialDofs, 1, 0 );
    uint                    tCounter = 0;
    for ( uint Ik = 0; Ik < tMat.numel(); Ik++ )
    {
        if ( tMat( Ik ) == 0 )
        {
            tNotInBlock( tCounter ) = Ik;
            tCounter++;
        }
    }

    tNotInBlock.resize( tCounter, 1 );

    tCriteriaIds.resize( tCriteriaIds.size() + tCounter );

    for ( uint Ik = 0; Ik < tNotInBlock.numel(); Ik++ )
    {
        tCriteriaIds( Ik + tNumBlocksInitialBlocks ).set_size( 1, 1, tNotInBlock( Ik ) );
    }

    //---------------------------------------------------------------

    sint tNumBlocks = tCriteriaIds.size();

    PCSetType( mLinearSolverAlgoritm->mpc, PCASM );

    PCASMSetLocalType( mLinearSolverAlgoritm->mpc, PC_COMPOSITE_ADDITIVE );
    PCASMSetType( mLinearSolverAlgoritm->mpc, PC_ASM_BASIC );
    PCASMSetOverlap( mLinearSolverAlgoritm->mpc, 0 );

    moris::Cell< IS > tIs( tNumBlocks );

    for ( sint Ik = 0; Ik < tNumBlocks; Ik++ )
    {
        ISCreateGeneral( PETSC_COMM_WORLD, tCriteriaIds( Ik ).numel(), tCriteriaIds( Ik ).data(), PETSC_COPY_VALUES, &tIs( Ik ) );
        //         ISView( tIs( Ik ) , PETSC_VIEWER_STDOUT_SELF );
    }

    PCASMSetLocalSubdomains( mLinearSolverAlgoritm->mpc, tNumBlocks, NULL, tIs.data().data() );

    KSPSetUp( mLinearSolverAlgoritm->mPetscKSPProblem );

    KSP *tSubksp;

    PCASMGetSubKSP( mLinearSolverAlgoritm->mpc, NULL, NULL, &tSubksp );

    Cell< PC > tPCBlock( tNumBlocks );
    for ( sint Ik = 0; Ik < tNumBlocks; Ik++ )
    {
        KSPSetType( tSubksp[ Ik ], KSPRICHARDSON );
        KSPGetPC( tSubksp[ Ik ], &tPCBlock( Ik ) );
        PCSetType( tPCBlock( Ik ), PCLU );
        PCSetUp( tPCBlock( Ik ) );

        //         KSPView( tSubksp[ Ik ], PETSC_VIEWER_STDOUT_SELF );
        //         PCView( tPCBlock( Ik ), PETSC_VIEWER_STDOUT_SELF );
    }

    //     PCView( mpc, PETSC_VIEWER_STDOUT_SELF );
}

void
Preconditioner_PETSc::build_schwarz_preconditioner( Linear_Problem *aLinearSystem )
{
    // set schwarz preconditiiner domains based on criteria
    moris::real tVolumeFractionThreshold = mLinearSolverAlgoritm->mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );

    moris::Cell< moris::Matrix< IdMat > > tCriteriaIds;
    mLinearSolverAlgoritm->mSolverInterface->get_adof_ids_based_on_criteria( tCriteriaIds,
            tVolumeFractionThreshold );


    uint tNumSerialDofs = mLinearSolverAlgoritm->mSolverInterface->get_num_my_dofs();

    moris::Matrix< DDSMat > tMat( tNumSerialDofs, 1, 0 );

    uint tNumBlocksInitialBlocks = tCriteriaIds.size();

    for ( uint Ik = 0; Ik < tNumBlocksInitialBlocks; Ik++ )
    {
        for ( uint Ii = 0; Ii < tCriteriaIds( Ik ).numel(); Ii++ )
        {
            if ( tMat( tCriteriaIds( Ik )( Ii ) ) == 0 )
            {
                tMat( tCriteriaIds( Ik )( Ii ) ) = 1;
            }
        }
    }

    moris::Matrix< DDSMat > tNotInBlock( tNumSerialDofs, 1, 0 );

    uint tCounter = 0;
    for ( uint Ik = 0; Ik < tMat.numel(); Ik++ )
    {
        if ( tMat( Ik ) == 0 )
        {
            tNotInBlock( tCounter ) = Ik;
            tCounter++;
        }
    }

    tNotInBlock.resize( tCounter, 1 );

    tCriteriaIds.resize( tCriteriaIds.size() + tCounter );

    for ( uint Ik = 0; Ik < tNotInBlock.numel(); Ik++ )
    {
        tCriteriaIds( Ik + tNumBlocksInitialBlocks ).set_size( 1, 1, tNotInBlock( Ik ) );
    }

    //---------------------------------------------------------------

    sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );

    mMapFree = tMatFactory.create_map( aLinearSystem->get_solver_input()->get_my_local_global_map(),
            aLinearSystem->get_solver_input()->get_constrained_Ids() );

    // Build matrix
    mPreconMat = tMatFactory.create_matrix( aLinearSystem->get_solver_input(), mMapFree );

    for ( uint Ik = 0; Ik < tCriteriaIds.size(); Ik++ )
    {
        uint tNumCriterias = tCriteriaIds( Ik ).numel();

        moris::Matrix< DDRMat > tValues( tNumCriterias, tNumCriterias, 0.0 );
        aLinearSystem->get_matrix()->get_matrix_values( tCriteriaIds( Ik ), tValues );

        moris::Matrix< DDRMat > tValuesInv;

        if ( tValues.numel() == 1 )
        {
            tValuesInv = { { 1.0 / std::sqrt( std::abs( tValues( 0 ) ) ) } };
        }
        else
        {
            tValuesInv = inv( tValues );
        }

        mPreconMat->fill_matrix(
                tNumCriterias,
                tValuesInv,
                tCriteriaIds( Ik ) );
    }

    mPreconMat->matrix_global_assembly();
}
