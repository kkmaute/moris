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
#include "cl_SOL_Dist_Map_Custom.hpp"

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

    // Generic LU path: choose the factor package through the parameter list.
    // Unlike the legacy mumps branch, availability is determined by PETSc itself.
    else if ( mParameterList.get< std::string >( "PCType" ) == "lu" )
    {
        PetscCallAbort( PETSC_COMM_WORLD, PCSetType( mpc, PCLU ) );
        const std::string tSolver = mParameterList.exists( "factor_solver_type" )
                                          ? mParameterList.get< std::string >( "factor_solver_type" )
                                          : "";

        if ( !tSolver.empty() )
        {
            PetscCallAbort( PETSC_COMM_WORLD, PCFactorSetMatSolverType( mpc, tSolver.c_str() ) );
        }
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
    else if ( mParameterList.get< std::string >( "PCType" ) == "bddc" )
    {
        // BDDC requires the overlapping subdomain matrices retained by MATIS.
        Mat       tMatrix  = aLinearSystem->get_matrix()->get_petsc_matrix();
        PetscBool tIsMATIS = PETSC_FALSE;

        PetscCallAbort( PETSC_COMM_WORLD, PetscObjectTypeCompare( (PetscObject)tMatrix, MATIS, &tIsMATIS ) );
        MORIS_ERROR( tIsMATIS,
                "BDDC requires MATIS. Set SOL_PETSc_Mat_Type to is before assembly (or use PETSc -mat_type is)." );

        // Set the preconditioner type before supplying BDDC metadata.
        MORIS_LOG_INFO( "KSP Preconditioner: bddc" );
        PetscCallAbort( PETSC_COMM_WORLD, PCSetType( mpc, PCBDDC ) );

        // Overlap maps may include ghost DOFs with no locally integrated terms.
        // Compress a preconditioning copy so BDDC never factors these empty rows;
        // retain the original map and matrix for subsequent MORIS reassembly.
        Mat tLocal;
        IS  tZeroRows;
        PetscCallAbort( PETSC_COMM_WORLD, MatISGetLocalMat( tMatrix, &tLocal ) );
        PetscCallAbort( PETSC_COMM_WORLD, MatFindZeroRows( tLocal, &tZeroRows ) );

        PetscInt tNumZeroRows = 0;
        if ( tZeroRows )
        {
            PetscCallAbort( PETSC_COMM_WORLD, ISGetLocalSize( tZeroRows, &tNumZeroRows ) );
        }

        PetscCallAbort( PETSC_COMM_WORLD, ISDestroy( &tZeroRows ) );

        // All ranks must participate in compaction if any rank has empty rows.
        PetscInt tGlobalZeroRows;
        MPI_Allreduce( &tNumZeroRows, &tGlobalZeroRows, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD );

        Mat tPreconditioningMatrix = tMatrix;

        if ( tGlobalZeroRows )
        {
            // PETSc 3.23's MatISFixLocalEmpty takes different MPIU_Allreduce
            // call sites when only some ranks have empty rows. Build the compact
            // MATIS explicitly so every rank follows the same collective sequence.
            // Retain the union of numerical row/column support for nonsymmetric A.
            PetscInt tLocalSize;
            PetscCallAbort( PETSC_COMM_WORLD, MatGetSize( tLocal, &tLocalSize, nullptr ) );

            std::vector< bool > tActive( tLocalSize, false );

            for ( PetscInt i = 0; i < tLocalSize; ++i )
            {
                PetscInt           tCount;
                const PetscInt    *tColumns;
                const PetscScalar *tValues;
                PetscCallAbort( PETSC_COMM_WORLD, MatGetRow( tLocal, i, &tCount, &tColumns, &tValues ) );

                for ( PetscInt j = 0; j < tCount; ++j )
                {
                    if ( PetscAbsScalar( tValues[ j ] ) > 0.0 )
                    {
                        tActive[ i ]             = true;
                        tActive[ tColumns[ j ] ] = true;
                    }
                }

                PetscCallAbort( PETSC_COMM_WORLD, MatRestoreRow( tLocal, i, &tCount, &tColumns, &tValues ) );
            }

            // Preserve the global IDs associated with the retained local indices.
            ISLocalToGlobalMapping tOriginalMap;
            PetscCallAbort( PETSC_COMM_WORLD, MatISGetLocalToGlobalMapping( tMatrix, &tOriginalMap, nullptr ) );

            const PetscInt *tOriginalIds;
            PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingGetIndices( tOriginalMap, &tOriginalIds ) );

            std::vector< PetscInt > tKeep, tGlobalIds;

            for ( PetscInt i = 0; i < tLocalSize; ++i )
            {
                if ( tActive[ i ] )
                {
                    tKeep.push_back( i );
                    tGlobalIds.push_back( tOriginalIds[ i ] );
                }
            }

            PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingRestoreIndices( tOriginalMap, &tOriginalIds ) );

            // Extract the active rows and columns from the subdomain matrix.
            IS tKeepIS;
            PetscCallAbort( PETSC_COMM_WORLD, ISCreateGeneral( PETSC_COMM_SELF, tKeep.size(), tKeep.data(), PETSC_COPY_VALUES, &tKeepIS ) );

            Mat tCompactLocal;
            PetscCallAbort( PETSC_COMM_WORLD, MatCreateSubMatrix( tLocal, tKeepIS, tKeepIS, MAT_INITIAL_MATRIX, &tCompactLocal ) );

            // Change the overlap map while preserving the owned global layout.
            ISLocalToGlobalMapping tCompactMap;
            PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingCreate( PETSC_COMM_WORLD, 1, tGlobalIds.size(), tGlobalIds.data(), PETSC_COPY_VALUES, &tCompactMap ) );

            PetscInt tOwnedSize, tGlobalSize;
            PetscCallAbort( PETSC_COMM_WORLD, MatGetLocalSize( tMatrix, &tOwnedSize, nullptr ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatGetSize( tMatrix, &tGlobalSize, nullptr ) );

            PetscCallAbort( PETSC_COMM_WORLD, MatCreateIS( PETSC_COMM_WORLD, 1, tOwnedSize, tOwnedSize, tGlobalSize, tGlobalSize, tCompactMap, tCompactMap, &tPreconditioningMatrix ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatISSetLocalMat( tPreconditioningMatrix, tCompactLocal ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatAssemblyBegin( tPreconditioningMatrix, MAT_FINAL_ASSEMBLY ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatAssemblyEnd( tPreconditioningMatrix, MAT_FINAL_ASSEMBLY ) );

            // Carry the near-nullspace metadata to the preconditioning copy.
            MatNullSpace tNearNullSpace;
            PetscCallAbort( PETSC_COMM_WORLD, MatGetNearNullSpace( tMatrix, &tNearNullSpace ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatSetNearNullSpace( tPreconditioningMatrix, tNearNullSpace ) );

            // MATIS retains references to its local matrix and map.
            PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingDestroy( &tCompactMap ) );
            PetscCallAbort( PETSC_COMM_WORLD, ISDestroy( &tKeepIS ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatDestroy( &tCompactLocal ) );
        }

        // Keep the assembled MORIS matrix as the operator for the linear solve.
        PetscCallAbort( PETSC_COMM_WORLD, PCSetOperators( mpc, tMatrix, tPreconditioningMatrix ) );

        if ( tGlobalZeroRows )
        {
            // PCSetOperators retains the preconditioning matrix.
            PetscCallAbort( PETSC_COMM_WORLD, MatDestroy( &tPreconditioningMatrix ) );
        }

        // Obtain the constrained DOFs in the original subdomain ordering.
        auto *tMorisMatrix = dynamic_cast< Matrix_PETSc * >( aLinearSystem->get_matrix() );
        MORIS_ERROR( tMorisMatrix, "BDDC requires a MORIS PETSc matrix." );

        auto                   tBoundaryDofs = tMorisMatrix->get_constrained_subdomain_dofs();
        ISLocalToGlobalMapping tMap;
        PetscCallAbort( PETSC_COMM_WORLD, MatGetLocalToGlobalMapping( tMatrix, &tMap, nullptr ) );

        // Global boundary numbering remains valid if the preconditioning copy
        // has removed unused local DOFs and changed its local ordering.
        PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingApply( tMap, tBoundaryDofs.size(), tBoundaryDofs.data(), tBoundaryDofs.data() ) );

        IS tBoundary;
        PetscCallAbort( PETSC_COMM_WORLD, ISCreateGeneral( PETSC_COMM_WORLD, tBoundaryDofs.size(), tBoundaryDofs.data(), PETSC_COPY_VALUES, &tBoundary ) );
        PetscCallAbort( PETSC_COMM_WORLD, PCBDDCSetDirichletBoundaries( mpc, tBoundary ) );
        PetscCallAbort( PETSC_COMM_WORLD, ISDestroy( &tBoundary ) );

        // Optionally identify each requested DOF type as a separate BDDC field.
        if ( mParameterList.exists( "bddc_use_dof_splitting" ) && mParameterList.get< bool >( "bddc_use_dof_splitting" ) )
        {
            auto *tInput  = aLinearSystem->get_solver_input();
            auto *tDofMap = dynamic_cast< Dist_Map_Custom * >( aLinearSystem->get_matrix()->get_map() );
            MORIS_ERROR( tInput && tDofMap, "BDDC DOF splitting requires a solver interface and a PETSc DOF map." );

            const auto       &tDofTypes = tInput->get_requested_dof_types();
            std::vector< IS > tFields( tDofTypes.size() );

            // The solver interface supplies owned MORIS IDs; PETSc expects global IDs.
            for ( uint i = 0; i < tDofTypes.size(); ++i )
            {
                Matrix< DDSMat > tIds = tInput->get_my_local_global_map( { tDofTypes( i ) } );
                tDofMap->map_from_moris_ids_to_petsc_ids( tIds );
                PetscCallAbort( PETSC_COMM_WORLD, ISCreateGeneral( PETSC_COMM_WORLD, tIds.numel(), tIds.data(), PETSC_COPY_VALUES, &tFields[ i ] ) );
            }

            PetscCallAbort( PETSC_COMM_WORLD, PCBDDCSetDofsSplitting( mpc, tFields.size(), tFields.data() ) );

            for ( auto &tField : tFields )
            {
                PetscCallAbort( PETSC_COMM_WORLD, ISDestroy( &tField ) );
            }
        }

        // Optionally supplement graph vertices with geometric mesh corners.
        // This is needed for floating elastic subdomains whose physical-boundary
        // endpoints are not graph vertices. The original AIJ path is unaffected.
        if ( mParameterList.exists( "bddc_use_subdomain_corners" ) && mParameterList.get< bool >( "bddc_use_subdomain_corners" ) )
        {
            auto *tInput = aLinearSystem->get_solver_input();
            MORIS_ERROR( tInput, "BDDC geometric corners require a solver interface." );

            auto  tCorners = tInput->get_subdomain_corner_dof_ids();
            auto *tDofMap  = dynamic_cast< Dist_Map_Custom * >( aLinearSystem->get_matrix()->get_map() );
            MORIS_ERROR( tDofMap, "BDDC geometric corners require a PETSc DOF map." );
            tDofMap->map_from_moris_ids_to_petsc_ids( tCorners );

            IS tVertices;
            PetscCallAbort( PETSC_COMM_WORLD, ISCreateGeneral( PETSC_COMM_WORLD, tCorners.numel(), tCorners.data(), PETSC_COPY_VALUES, &tVertices ) );
            PetscCallAbort( PETSC_COMM_WORLD, PCBDDCSetPrimalVerticesIS( mpc, tVertices ) );
            PetscCallAbort( PETSC_COMM_WORLD, ISDestroy( &tVertices ) );
        }

        // KSPSetFromOptions/SetUp follows this call. Managed BDDC options are
        // supplied by set_parameter_options; native PETSc overrides remain available.
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

//----------------------------------------------------------------------------------------

void Preconditioner_PETSc::set_parameter_options( PetscOptions aOptions )
{
    if ( mParameterList.get< std::string >( "PCType" ) != "bddc" )
    {
        return;
    }

    // PETSc exposes these BDDC switches through its options parser rather than
    // public setters. Keep the translation here, next to the BDDC implementation,
    // and honor the PC prefix so independent solvers do not share option names.
    const char *tPrefix = nullptr;
    PetscCallAbort( PETSC_COMM_WORLD, PCGetOptionsPrefix( mpc, &tPrefix ) );
    const std::string tOptionPrefix = "-" + std::string( tPrefix ? tPrefix : "" );

    if ( mParameterList.exists( "use_bddc_defaults" ) && !mParameterList.get< bool >( "use_bddc_defaults" ) )
    {
        const struct
        {
            const char *mParameter;
            const char *mOption;
        } tBDDCOptions[] = {
            { "bddc_symmetric", "pc_bddc_symmetric" },
            { "bddc_detect_disconnected", "pc_bddc_detect_disconnected" },
            { "bddc_use_vertices", "pc_bddc_use_vertices" },
            { "bddc_use_edges", "pc_bddc_use_edges" },
            { "bddc_use_faces", "pc_bddc_use_faces" },
            { "bddc_use_stiffness_scaling", "pc_is_use_stiffness_scaling" },
            { "bddc_use_deluxe_scaling", "pc_bddc_use_deluxe_scaling" },
            { "bddc_sub_schurs_symmetric", "sub_schurs_symmetric" },
            { "bddc_sub_schurs_posdef", "sub_schurs_posdef" }
        };

        // The sub_schurs options use the PC prefix directly; PETSc does not
        // prepend pc_bddc_ to their names.
        for ( const auto &tOption : tBDDCOptions )
        {
            const char *tValue = mParameterList.get< bool >( tOption.mParameter ) ? "true" : "false";
            PetscCallAbort( PETSC_COMM_WORLD, PetscOptionsSetValue( aOptions, ( tOptionPrefix + tOption.mOption ).c_str(), tValue ) );
        }
    }

    const char *tParameters[] = { "bddc_dirichlet_solver", "bddc_neumann_solver", "bddc_coarse_solver" };
    const char *tOptions[]    = { "pc_bddc_dirichlet_pc_factor_mat_solver_type",
           "pc_bddc_neumann_pc_factor_mat_solver_type",
           "pc_bddc_coarse_pc_factor_mat_solver_type" };

    for ( uint i = 0; i < 3; ++i )
    {
        if ( !mParameterList.exists( tParameters[ i ] ) )
        {
            continue;
        }

        const auto &tSolver = mParameterList.get< std::string >( tParameters[ i ] );

        if ( !tSolver.empty() )
        {
            PetscCallAbort( PETSC_COMM_WORLD, PetscOptionsSetValue( aOptions, ( tOptionPrefix + tOptions[ i ] ).c_str(), tSolver.c_str() ) );

            if ( i == 2 )
            {
                // An explicitly selected coarse factor package requires LU.
                PetscCallAbort( PETSC_COMM_WORLD, PetscOptionsSetValue( aOptions, ( tOptionPrefix + "pc_bddc_coarse_pc_type" ).c_str(), "lu" ) );
            }
        }
    }
}
