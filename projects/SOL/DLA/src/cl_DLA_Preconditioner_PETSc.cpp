/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_PETSc.cpp
 *
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
#include "HDF5_Tools.hpp"

using namespace moris;
using namespace dla;

//----------------------------------------------------------------------------------------

Preconditioner_PETSc::Preconditioner_PETSc(
        const Parameter_List &aParameterList )
        : Preconditioner( aParameterList )
{
}

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_ilu_preconditioner( Linear_Problem *aLinearSystem )
{
    // check for serail run; this preconditioner does not work in parallel
    MORIS_ERROR( par_size() == 1,
            "Preconditioner_PETSc::build_ilu_preconditioner - PETSc ILU preconditioner can only be used in serial." );

    // write preconditioner to log file
    MORIS_LOG_INFO( "KSP Preconditioner: ilu(%d)",
            mParameterList.get< moris::sint >( "ILUFill" ) );

    // Set PC type
    PCSetType( mpc, "ilu" );

    // Set operators
    PCSetOperators(
            mpc,
            aLinearSystem->get_matrix()->get_petsc_matrix(),
            aLinearSystem->get_matrix()->get_petsc_matrix() );

    // Set levels of fill for ILU
    PCFactorSetLevels(
            mpc,
            mParameterList.get< moris::sint >( "ILUFill" ) );

    // Set drop tolerance for Ilu
    PCFactorSetDropTolerance(
            mpc,
            mParameterList.get< moris::real >( "ILUTol" ),
            PETSC_DEFAULT,
            PETSC_DEFAULT );

    // Set preconditioner options from the options database
    PCSetFromOptions( mpc );

    // Finalize setup of preconditioner
    PCSetUp( mpc );
}

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_multigrid_preconditioner( Linear_Problem *aLinearSystem )
{
    // get number of mg levels
    sint tLevels = mParameterList.get< moris::sint >( "MultigridLevels" );

    // write preconditioner to log file
    MORIS_LOG_INFO( "KSP Preconditioner: mg(%d)", tLevels );

    // Set PC type
    PCSetType( mpc, "mg" );

    //----------------------------------------------------------------
    //                 Set general multigrid stuff
    //----------------------------------------------------------------

    // Build multigrid operators
    mGeoMultigrid = new Geometric_Multigrid( aLinearSystem->get_solver_input() );

    // get multigrid operators
    moris::Vector< sol::Dist_Matrix * > tProlongationList =
            mGeoMultigrid->get_prolongation_list();

    // set multigrid levels and type
    PCMGSetLevels( mpc, tLevels, nullptr );

    PCMGSetType( mpc, PC_MG_MULTIPLICATIVE );
    PCMGSetGalerkin( mpc, PC_MG_GALERKIN_BOTH );

    // get restriction operators
    moris::Vector< Mat > tTransposeOperators( tProlongationList.size() );
    for ( moris::uint Ik = 0; Ik < tProlongationList.size(); Ik++ )
    {
        MatTranspose( tProlongationList( Ik )->get_petsc_matrix(), MAT_INITIAL_MATRIX, &tTransposeOperators( Ik ) );

        if ( Ik == 0 )
        {
            MORIS_LOG_INFO( "Writing interpolation operator to MATLAB file" );

            std::string tString = "Interpolation_Operator";
            tProlongationList( Ik )->save_matrix_to_matlab_file( tString.c_str() );
        }
    }

    // set restriction operators
    moris::sint tCounter = 0;
    for ( moris::sint Ik = tLevels - 1; Ik > 0; Ik-- )
    {
        PCMGSetInterpolation( mpc, Ik, tTransposeOperators( tCounter ) );
        //         PCMGSetRestriction( mpc, Ik, tTransposeOperators( tCounter ) );

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
    PCMGGetCoarseSolve( mpc, &tPetscKSPCoarseSolve );
    KSPSetType( tPetscKSPCoarseSolve, KSPRICHARDSON );
    PetscInt maxitss = 1000;
    KSPSetTolerances( tPetscKSPCoarseSolve, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxitss );
    PC cpc;
    KSPGetPC( tPetscKSPCoarseSolve, &cpc );
    PCSetType( cpc, PCLU );

    //----------------------------------------------------------------
    //                 Set schwarz preconditioner for finest level
    //----------------------------------------------------------------

    if ( mParameterList.get< bool >( "MG_use_schwarz_smoother" ) )
    {
        // set schwarz preconditiiner domains based on criteria
        moris::real tVolumeFractionThreshold =    //
                mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );
        moris::sint tSchwarzSmoothingIters =    //
                mParameterList.get< moris::sint >( "MG_schwarz_smoothing_iters" );

        moris::Vector< moris::Matrix< IdMat > > tCriteriaIds;
        aLinearSystem->get_solver_input()->get_adof_ids_based_on_criteria( tCriteriaIds,
                tVolumeFractionThreshold );

        //---------------------------------------------------------------

        uint tNumSerialDofs = aLinearSystem->get_solver_input()->get_num_my_dofs();

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

        // construct string from output file name
        std::string tStrOutputFile = mParameterList.get< std::string >( "ASM_blocks_output_filename" );

        // if not empty
        if ( tStrOutputFile != "" )
        {
            // initialize hdf5 file
            hid_t  tFileID = create_hdf5_file( tStrOutputFile );
            herr_t tStatus = 0;

            for ( uint I = 0; I < tCriteriaIds.size(); I++ )
            {
                std::string str = ios::stringify( I );
                // write to file
                save_matrix_to_hdf5_file( tFileID, str, tCriteriaIds( I ), tStatus );
            }

            // close file
            close_hdf5_file( tFileID );
        }

        //---------------------------------------------------------------

        KSP dKSPFirstDown;
        PCMGGetSmootherDown( mpc, tLevels - 1, &dKSPFirstDown );
        PCMGGetSmootherUp( mpc, tLevels - 1, &dKSPFirstDown );
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

        moris::Vector< IS > tIs( tNumBlocks );

        for ( uint Ik = 0; Ik < tNumBlocks; Ik++ )
        {
            ISCreateGeneral( PETSC_COMM_WORLD, tCriteriaIds( Ik ).numel(), tCriteriaIds( Ik ).data(), PETSC_COPY_VALUES, &tIs( Ik ) );
            //             ISView( tIs( Ik ) , PETSC_VIEWER_STDOUT_SELF );
        }

        PCASMSetLocalSubdomains( dPCFirstDown, tNumBlocks, nullptr, tIs.data().data() );
        std::cout << "build-----" << '\n';
        KSPSetUp( dKSPFirstDown );

        KSP *tSubksp;

        sint nlocal;
        sint first;

        PCASMGetSubKSP( dPCFirstDown, &nlocal, &first, &tSubksp );

        //         tKSPBlock.resize( tNumBlocks );
        Vector< PC > tPCBlock( nlocal );
        for ( sint Ik = 0; Ik < nlocal; Ik++ )
        {
            KSPSetType( tSubksp[ Ik ], KSPRICHARDSON );
            KSPGetPC( tSubksp[ Ik ], &tPCBlock( Ik ) );
            PCSetType( tPCBlock( Ik ), PCLU );
            PCSetUp( tPCBlock( Ik ) );

            //             KSPView( tSubksp[ Ik ], PETSC_VIEWER_STDOUT_SELF );
            //             PCView( tPCBlock( Ik ), PETSC_VIEWER_STDOUT_SELF );
        }

        std::cout << "----- ksps build -----" << '\n';

        //         PCView( dPCFirstDown, PETSC_VIEWER_STDOUT_SELF );
    }
    else
    {
        KSP dkspDown1;
        PCMGGetSmootherDown( mpc, tLevels - 1, &dkspDown1 );
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
        PCMGGetSmootherUp( mpc, tLevels - 1, &dkspUp1 );
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
        PCMGGetSmootherDown( mpc, Ik, &dkspDown );
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

        PCMGGetSmootherUp( mpc, k, &dkspUp );
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

    delete mGeoMultigrid;
}

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_schwarz_preconditioner_petsc( Linear_Problem *aLinearSystem, KSP aPetscKSPProblem )
{
    // write preconditioner to log file
    MORIS_LOG_INFO( "KSP Preconditioner: asm" );

    // Set PC type
    PCSetType( mpc, "asm" );

    // set schwarz preconditioner domains based on criteria
    moris::real tVolumeFractionThreshold = mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );

    moris::Vector< moris::Matrix< IdMat > > tCriteriaIds;
    aLinearSystem->get_solver_input()->get_adof_ids_based_on_criteria( tCriteriaIds,
            tVolumeFractionThreshold );

    uint tNumSerialDofs = aLinearSystem->get_solver_input()->get_num_my_dofs();

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

    PCSetType( mpc, PCASM );

    PCASMSetLocalType( mpc, PC_COMPOSITE_ADDITIVE );
    PCASMSetType( mpc, PC_ASM_BASIC );
    PCASMSetOverlap( mpc, 0 );

    moris::Vector< IS > tIs( tNumBlocks );

    for ( sint Ik = 0; Ik < tNumBlocks; Ik++ )
    {
        ISCreateGeneral( PETSC_COMM_WORLD, tCriteriaIds( Ik ).numel(), tCriteriaIds( Ik ).data(), PETSC_COPY_VALUES, &tIs( Ik ) );
        //         ISView( tIs( Ik ) , PETSC_VIEWER_STDOUT_SELF );
    }

    PCASMSetLocalSubdomains( mpc, tNumBlocks, nullptr, tIs.data().data() );

    KSPSetUp( aPetscKSPProblem );

    KSP *tSubksp;

    PCASMGetSubKSP( mpc, nullptr, nullptr, &tSubksp );

    Vector< PC > tPCBlock( tNumBlocks );
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

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_schwarz_preconditioner( Linear_Problem *aLinearSystem )
{
    // write preconditioner to log file
    MORIS_LOG_INFO( "KSP Preconditioner: asm - mat" );

    // Set PC type
    PCSetType( mpc, "mat" );

    // set schwarz preconditiiner domains based on criteria
    moris::real tVolumeFractionThreshold = mParameterList.get< moris::real >( "ASM_volume_fraction_threshold" );

    moris::Vector< moris::Matrix< IdMat > > tCriteriaIds;
    aLinearSystem->get_solver_input()->get_adof_ids_based_on_criteria( tCriteriaIds,
            tVolumeFractionThreshold );

    uint tNumSerialDofs = aLinearSystem->get_solver_input()->get_num_my_dofs();

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

    mMapFree = tMatFactory.create_map(
            aLinearSystem->get_solver_input()->get_my_local_global_map(),
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

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_algebaric_multigrid_preconditioner( Linear_Problem *aLinearSystem )
{
    // write preconditioner to log file
    MORIS_LOG_INFO( "KSP Preconditioner: amg" );

    // Set PC type
    PCSetType( mpc, PCGAMG );

    // get number of pde equaitons for each grid
    sint tNumPDEEquaitons = mParameterList.get< moris::sint >( "num_pde_equations" );
    MatSetBlockSize( aLinearSystem->get_matrix()->get_petsc_matrix(), tNumPDEEquaitons );

    // see if the defualt values should be used then skip
    if ( mParameterList.get< bool >( "use_gamg_defaults" ) ) return;

    // set the values from the parameter list
    PCGAMGSetType( mpc, mParameterList.get< std::string >( "amg_type" ).c_str() );
    // real tThreshold = mParameterList.get< real >( "amg_threshold" );
    // PCGAMGSetThreshold( mpc,&tThreshold );

    // loop over the levels and set the smoothers(KSP,PC) for each level
    // KSPChebyshevEstEigGetKSP( mpc, &tKSPBlock( 0 ) );
}

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::build_preconditioner( Linear_Problem *aLinearSystem, KSP aPetscKSPProblem )
{
    // get preconditioner
    KSPGetPC( aPetscKSPProblem, &mpc );

    // set direct solver: superlu-dist
    if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "superlu-dist" ) )
    {
        // write solver to log file
        MORIS_LOG_INFO( "PC Solver: superlu-dist" );

        // set LU preconditioner
        PCSetType( mpc, PCLU );

        // set factorization method in preconditioner
        PCFactorSetMatSolverType( mpc, MATSOLVERSUPERLU_DIST );

        // set up the package to call for the factorization
        PCFactorSetUpMatSolverType( mpc );
    }

    // set direct solver: mumps
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "mumps" ) )
    {
#ifdef MORIS_USE_MUMPS
        // write solver to log file
        MORIS_LOG_INFO( "PC Solver: mumps" );

        // set LU preconditioner assuming that system is non-symmetric
        PCSetType( mpc, PCLU );

        // set factorization method in preconditioner
        PCFactorSetMatSolverType( mpc, MATSOLVERMUMPS );

        // set up the package to call for the factorization
        PCFactorSetUpMatSolverType( mpc );

        // get the factored matrix F from the preconditioner context
        Mat F;
        PCFactorGetMatrix( mpc, &F );

        // set MUMPS integer control parameters ICNTL to be passed to
        // MUMPS.  Setting entry 7 of MUMPS ICNTL array (of size 40) to a value
        // of 2. This sets use of Approximate Minimum Fill (AMF)
        PetscInt ival = 2, icntl = 7;

        // pass control parameters to MUMPS

        MatMumpsSetIcntl( F, icntl, ival );
#else
        MORIS_ERROR( false,
                "Linear_Solver_PETSc::construct_solver_and_preconditioner - MORIS installed without support for MUMPS" );
#endif
    }

    // set iterative solver: kspgmres
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "ilu" ) )
    {
        this->build_ilu_preconditioner( aLinearSystem );
    }
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "mg" ) )
    {
        this->build_multigrid_preconditioner( aLinearSystem );
    }
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "asm" ) )
    {
        // build schwarz preconditioner
        this->build_schwarz_preconditioner_petsc( aLinearSystem, aPetscKSPProblem );
    }
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "mat" ) )
    {
        // build schwarz preconditioner
        this->build_schwarz_preconditioner( aLinearSystem );
    }
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "gamg" ) )
    {
        this->build_algebaric_multigrid_preconditioner( aLinearSystem );
    }
    else if ( !strcmp( mParameterList.get< std::string >( "PCType" ).c_str(), "none" ) )
    {
        // Set PC type to none
        PCSetType( mpc, "none" );
    }
    else
    {
        MORIS_ERROR( false,
                "Linear_Solver_PETSc::construct_solver_and_preconditioner - no valid preconditioner was found." );
    }
}
