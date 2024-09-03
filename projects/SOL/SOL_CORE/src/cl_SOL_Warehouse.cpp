/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Warehouse.cpp
 *
 */

#include <functional>
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Communication_Tools.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_DLA_Eigen_Solver.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#ifdef MORIS_HAVE_PETSC
#include <petsc.h>
#ifdef MORIS_HAVE_SLEPC
#include <slepceps.h>
#endif
#endif

#include "cl_Library_IO.hpp"

using namespace moris;
using namespace sol;

// User-defined pointer function
typedef bool ( *Pointer_Function )( void* aPointer );

SOL_Warehouse::~SOL_Warehouse()
{
    for ( auto tLinearSolver : mLinearSolvers )
    {
        delete tLinearSolver;
    }
    for ( auto tNonLinearSolver : mNonlinearSolvers )
    {
        delete tNonLinearSolver;
    }
    for ( auto tTimeSolver : mTimeSolvers )
    {
        delete tTimeSolver;
    }

    for ( uint Ik = 0; Ik < mLinearSolverAlgorithms.size(); Ik++ )
    {
        mLinearSolverAlgorithms( Ik ).reset();
    }
    for ( uint Ik = 0; Ik < mNonlinearSolverAlgorithms.size(); Ik++ )
    {
        mNonlinearSolverAlgorithms( Ik ).reset();
    }
    for ( uint Ik = 0; Ik < mTimeSolverAlgorithms.size(); Ik++ )
    {
        mTimeSolverAlgorithms( Ik ).reset();
    }

#ifdef MORIS_HAVE_PETSC
    if ( mTPLType == moris::sol::MapType::Petsc )
    {
        PetscFinalize();
#ifdef MORIS_HAVE_SLEPC
	SlepcFinalize();
#endif
    }
#endif

    for ( auto& iPreconditioner : mPreconditioners )
    {
        delete iPreconditioner;
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::initialize()
{
    // Load parameters from parameterlist
    mTPLType = static_cast< moris::sol::MapType >( mParameterlist( 6 )( 0 ).get< moris::uint >( "SOL_TPL_Type" ) );

    mOperatorToMatlab      = mParameterlist( 6 )( 0 ).get< std::string >( "SOL_save_operator_to_matlab" );
    mSaveFinalSolVecToFile = mParameterlist( 6 )( 0 ).get< std::string >( "SOL_save_final_sol_vec_to_file" );

    mLoadSolVecFromFile    = mParameterlist( 6 )( 0 ).get< std::string >( "SOL_load_sol_vec_from_file" );
    mSolVecDataGroup       = mParameterlist( 6 )( 0 ).get< std::string >( "SOL_load_sol_vec_data_group" );
    mSolVecNumberOfVectors = mParameterlist( 6 )( 0 ).get< sint >( "SOL_load_sol_vec_num_vec" );

    mFilenameInitialGuess = mParameterlist( 6 )( 0 ).get< std::string >( "TSA_Initial_Sol_Vec" );

    mSaveFinalAdjointVecToFile = mParameterlist( 6 )( 0 ).get< std::string >( "SOL_save_final_adjoint_vec_to_file" );

    sol::SensitivityAnalysisType tSensitivityAnalysisType =
            static_cast< sol::SensitivityAnalysisType >( mParameterlist( 6 )( 0 ).get< moris::uint >( "Sensitivity_Analysis_Type" ) );

    mIsAdjointSensitivityAnalysis = tSensitivityAnalysisType == sol::SensitivityAnalysisType::ADJOINT;

#ifdef MORIS_HAVE_PETSC
    if ( mTPLType == moris::sol::MapType::Petsc )
    {
        PetscInitializeNoArguments();
#ifdef MORIS_HAVE_SLEPC
	SlepcInitializeNoArguments();
#endif
    
    }
#endif

    // create the requested preconditioner
    this->create_preconditioner_algorithms();

    // build solvers and solver algorithms
    this->create_linear_solver_algorithms();

    this->create_linear_solvers();

    this->create_nonlinear_solver_algorithms();

    this->create_nonlinear_solvers();

    this->create_time_solver_algorithms();

    this->create_time_solvers();
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_preconditioner_algorithms()
{
    uint tNumLinAlgorithms = mParameterlist( 7 ).size();

    mPreconditioners.resize( tNumLinAlgorithms );

    moris::dla::Solver_Factory tSolFactory;

    for ( uint Ik = 0; Ik < tNumLinAlgorithms; Ik++ )
    {
        mPreconditioners( Ik ) = tSolFactory.create_preconditioner( mParameterlist( 7 )( Ik ) );
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_linear_solver_algorithms()
{
    uint tNumLinAlgorithms = mParameterlist( 0 ).size();

    mLinearSolverAlgorithms.resize( tNumLinAlgorithms );

    moris::dla::Solver_Factory tSolFactory;

    for ( uint Ik = 0; Ik < tNumLinAlgorithms; Ik++ )
    {
        mLinearSolverAlgorithms( Ik ) = tSolFactory.create_solver( mParameterlist( 0 )( Ik ) );

        // get and set nonlinear sub-solvers for staggered methods
        Vector< uint > tPreconditionerIndices;
        string_to_cell( mParameterlist( 0 )( Ik ).get< std::string >( "preconditioners" ),
                tPreconditionerIndices );

        // loop over the list of preconditioners and assign it to the solver
        for ( auto tPrecIndex : tPreconditionerIndices )
        {
            mLinearSolverAlgorithms( Ik )->set_preconditioner( mPreconditioners( tPrecIndex ) );
        }

        // get and set nonlinear sub-solvers for staggered methods
        tPreconditionerIndices.clear();
        string_to_cell( mParameterlist( 0 )( Ik ).get< std::string >( "preconditioners_linear_operator" ),
                tPreconditionerIndices );

        // loop over the list of preconditioner and assign it to the solver
        for ( auto tPrecIndex : tPreconditionerIndices )
        {
            mLinearSolverAlgorithms( Ik )->set_left_hand_side_preconditioner( mPreconditioners( tPrecIndex ) );
        }

        // get and set nonlinear sub-solvers for staggered methods
        tPreconditionerIndices.clear();
        string_to_cell( mParameterlist( 0 )( Ik ).get< std::string >( "sub_linear_solver" ),
                tPreconditionerIndices );

        // set the sub linear solver for the iegn problem
        for ( auto tPrecIndex : tPreconditionerIndices )
        {
            MORIS_ASSERT( tPrecIndex < Ik, "SOL_Warehouse::create_linear_solver_algorithms - sub_linear_solver index %d is not defined yet", tPrecIndex );

            std::string tSubsolverPrec = mParameterlist( 0 )( tPrecIndex ).get< std::string >( "preconditioners" );

            if ( tSubsolverPrec not_eq "" )
            {
                mLinearSolverAlgorithms( Ik )->set_sublinear_solver_options( &mParameterlist( 0 )( tPrecIndex ), &mParameterlist( 7 )( std::stoi( tSubsolverPrec ) ) );
                continue;
            }

            mLinearSolverAlgorithms( Ik )->set_sublinear_solver_options( &mParameterlist( 0 )( tPrecIndex ), nullptr );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_linear_solvers()
{
    uint tNumLinSolvers = mParameterlist( 1 ).size();

    mLinearSolvers.resize( tNumLinSolvers );

    for ( uint Ik = 0; Ik < tNumLinSolvers; Ik++ )
    {
        mLinearSolvers( Ik ) = new dla::Linear_Solver( mParameterlist( 1 )( Ik ) );

        Vector< uint > tLinearSolverAlgorithmIndices;
        string_to_cell( mParameterlist( 1 )( Ik ).get< std::string >( "DLA_Linear_solver_algorithms" ),
                tLinearSolverAlgorithmIndices );

        // set output for LHS
        mLinearSolvers( Ik )->set_LHS_output_filename( mParameterlist( 1 )( Ik ).get< std::string >( "DLA_LHS_output_filename" ) );

        for ( uint Ii = 0; Ii < tLinearSolverAlgorithmIndices.size(); Ii++ )
        {
            mLinearSolvers( Ik )->set_linear_algorithm( Ii, mLinearSolverAlgorithms( tLinearSolverAlgorithmIndices( Ii ) ) );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_nonlinear_solver_algorithms()
{
    uint tNumNonLinAlgorithms = mParameterlist( 2 ).size();

    mNonlinearSolverAlgorithms.resize( tNumNonLinAlgorithms );

    NLA::Nonlinear_Solver_Factory tNonlinFactory;

    for ( uint Ik = 0; Ik < tNumNonLinAlgorithms; Ik++ )
    {
        mNonlinearSolverAlgorithms( Ik ) =
                tNonlinFactory.create_nonlinear_solver(    //
                        mParameterlist( 2 )( Ik ) );

        mNonlinearSolverAlgorithms( Ik )->set_linear_solver( mLinearSolvers( mParameterlist( 2 )( Ik ).get< moris::sint >( "NLA_Linear_solver" ) ) );

        if ( mParameterlist( 2 )( Ik ).get< moris::sint >( "NLA_linear_solver_for_adjoint_solve" ) == -1 )
        {
            mNonlinearSolverAlgorithms( Ik )    //
                    ->set_linear_solver_for_adjoint_solve( mLinearSolvers( mParameterlist( 2 )( Ik ).get< moris::sint >( "NLA_Linear_solver" ) ) );
        }
        else
        {
            mNonlinearSolverAlgorithms( Ik )    //
                    ->set_linear_solver_for_adjoint_solve( mLinearSolvers( mParameterlist( 2 )( Ik ).get< moris::sint >( "NLA_linear_solver_for_adjoint_solve" ) ) );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_nonlinear_solvers()
{
    // get number of nonlinear solvers
    uint tNumNonLinSolvers = mParameterlist( 3 ).size();

    mNonlinearSolvers.resize( tNumNonLinSolvers );

    // loop over all nonlinear solvers
    for ( uint Ik = 0; Ik < tNumNonLinSolvers; Ik++ )
    {
        // create nonlinear solver
        mNonlinearSolvers( Ik ) =
                new NLA::Nonlinear_Solver(
                        mParameterlist( 3 )( Ik ) );

        // set nonlinear algorithms
        moris::Matrix< DDSMat > tMat;
        string_to_mat( mParameterlist( 3 )( Ik ).get< std::string >( "NLA_Nonlinear_solver_algorithms" ),
                tMat );

        for ( uint Ii = 0; Ii < tMat.numel(); Ii++ )
        {
            // check that algorithm exists
            MORIS_ERROR( mNonlinearSolverAlgorithms.size() > (uint)tMat( Ii ),
                    "SOL_Warehouse::create_nonlinear_solvers - NLA_Nonlinear_solver_algorithm %d does not exist",
                    tMat( Ii ) );

            MORIS_ERROR( mNonlinearSolverAlgorithms( tMat( Ii ) ),
                    "SOL_Warehouse::create_nonlinear_solvers - NLA_Nonlinear_solver_algorithms %d has not been created",
                    tMat( Ii ) );

            mNonlinearSolvers( Ik )->set_nonlinear_algorithm( mNonlinearSolverAlgorithms( tMat( Ii ) ), Ii );
        }

        // get and set list of dof types
        Vector< Vector< MSI::Dof_Type > > tCellOfCells;

        map< std::string, enum MSI::Dof_Type > tMap = MSI::get_msi_dof_type_map();

        string_to_cell_of_cell(
                mParameterlist( 3 )( Ik ).get< std::string >( "NLA_DofTypes" ),
                tCellOfCells,
                tMap );

        for ( uint Ii = 0; Ii < tCellOfCells.size(); Ii++ )
        {
            MORIS_ERROR( tCellOfCells( Ii )( 0 ) != MSI::Dof_Type::UNDEFINED, "Dof types for nonlinear solver %-5i not specified", Ik );

            mNonlinearSolvers( Ik )->set_dof_type_list( tCellOfCells( Ii ) );
        }

        // get and set secondary dof types
        Vector< Vector< MSI::Dof_Type > > tCellOfCellsSecDofTypes;

        string_to_cell_of_cell(
                mParameterlist( 3 )( Ik ).get< std::string >( "NLA_Secondary_DofTypes" ),
                tCellOfCellsSecDofTypes,
                tMap );

        // if no secondary dof types are explicitly defined, build default list of dof types
        if ( tCellOfCellsSecDofTypes.size() > 0 )
        {
            if ( tCellOfCellsSecDofTypes( 0 )( 0 ) == MSI::Dof_Type::UNDEFINED )
            {
                this->get_default_secondary_dof_types( tCellOfCellsSecDofTypes, tCellOfCells );
            }
        }

        for ( uint Ii = 0; Ii < tCellOfCellsSecDofTypes.size(); Ii++ )
        {
            mNonlinearSolvers( Ik )->set_secondary_dof_type_list( tCellOfCellsSecDofTypes( Ii ) );
        }

        // get and set nonlinear sub-solvers for staggered methods
        moris::Matrix< DDSMat > tNonlinearSubSolvers;
        string_to_mat( mParameterlist( 3 )( Ik ).get< std::string >( "NLA_Sub_Nonlinear_Solver" ),
                tNonlinearSubSolvers );

        for ( uint Ii = 0; Ii < tNonlinearSubSolvers.numel(); Ii++ )
        {
            MORIS_ERROR( tNonlinearSubSolvers( Ii ) <= (sint)Ik,
                    "Nonlinear sub-solvers need to be defined before they can be used as sub-solvers in staggered methods." );

            mNonlinearSolvers( Ik )->set_sub_nonlinear_solver( mNonlinearSolvers( tNonlinearSubSolvers( Ii ) ) );
        }

        mNonlinearSolvers( Ik )->set_solver_warehouse( this );
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_time_solver_algorithms()
{
    uint tNumTimeAlgorithms = mParameterlist( 4 ).size();

    mTimeSolverAlgorithms.resize( tNumTimeAlgorithms );

    tsa::Time_Solver_Factory tTimeSolverFactory;

    for ( uint Ik = 0; Ik < tNumTimeAlgorithms; Ik++ )
    {
        mTimeSolverAlgorithms( Ik ) =
                tTimeSolverFactory.create_time_solver( static_cast< moris::tsa::TimeSolverType >( mParameterlist( 4 )( Ik ).get< moris::uint >( "TSA_Solver_Implementation" ) ),
                        mParameterlist( 4 )( Ik ) );

        mTimeSolverAlgorithms( Ik )->set_nonlinear_solver( mNonlinearSolvers( mParameterlist( 4 )( Ik ).get< moris::sint >( "TSA_Nonlinear_Solver" ) ) );

        if ( mParameterlist( 4 )( Ik ).get< moris::sint >( "TSA_Nonlinear_Sensitivity_Solver" ) == -1 )
        {
            mTimeSolverAlgorithms( Ik )->set_nonlinear_solver_for_sensitivity_analysis( mNonlinearSolvers( mParameterlist( 4 )( Ik ).get< moris::sint >( "TSA_Nonlinear_Solver" ) ) );
        }
        else
        {
            mTimeSolverAlgorithms( Ik )->set_nonlinear_solver_for_sensitivity_analysis( mNonlinearSolvers( mParameterlist( 4 )( Ik ).get< moris::sint >( "TSA_Nonlinear_Sensitivity_Solver" ) ) );
        }

        // set output file names
        mTimeSolverAlgorithms( Ik )->set_output_filename( mParameterlist( 6 )( 0 ).get< std::string >( "TSA_Save_Sol_Vecs_to_file" ) );
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::create_time_solvers()
{
    // get number if time solvers
    uint tNumTimeSolvers = mParameterlist( 5 ).size();

    // Resize time solver list
    mTimeSolvers.resize( tNumTimeSolvers );

    // loop over requested number of time solvers and create them
    for ( uint Ik = 0; Ik < tNumTimeSolvers; Ik++ )
    {
        // Create time solver with user defined parameter list
        mTimeSolvers( Ik ) = new tsa::Time_Solver( mParameterlist( 5 )( Ik ), this );

        // get tie solver algorithm indices for this time solver
        moris::Matrix< DDSMat > tMat;
        string_to_mat( mParameterlist( 5 )( Ik ).get< std::string >( "TSA_Solver_algorithms" ),
                tMat );

        // add these time solver algorithms to time solver
        for ( uint Ii = 0; Ii < tMat.numel(); Ii++ )
        {
            mTimeSolvers( Ik )->set_time_solver_algorithm( mTimeSolverAlgorithms( tMat( Ii ) ), Ii );
        }

        // get requested dof types for this time solver
        Vector< Vector< MSI::Dof_Type > > tCellOfCells;

        map< std::string, enum MSI::Dof_Type > tMap = MSI::get_msi_dof_type_map();

        string_to_cell_of_cell( mParameterlist( 5 )( Ik ).get< std::string >( "TSA_DofTypes" ),
                tCellOfCells,
                tMap );

        // add requested dof types to time solver
        for ( uint Ii = 0; Ii < tCellOfCells.size(); Ii++ )
        {
            MORIS_ERROR( tCellOfCells( Ii )( 0 ) != MSI::Dof_Type::UNDEFINED, "Dof types for time solver %-5i not specified", Ik );
            mTimeSolvers( Ik )->set_dof_type_list( tCellOfCells( Ii ) );
        }

        // get strings of output indices and criteria
        std::string tStringOutputInd      = mParameterlist( 5 )( Ik ).get< std::string >( "TSA_Output_Indices" );
        std::string tStringOutputCriteria = mParameterlist( 5 )( Ik ).get< std::string >( "TSA_Output_Criteria" );

        if ( tStringOutputInd.size() > 0 )
        {
            moris::Matrix< DDSMat > tOutputIndices;
            string_to_mat( tStringOutputInd, tOutputIndices );

            Vector< std::string > tOutputCriteria;
            string_to_cell( tStringOutputCriteria, tOutputCriteria );

            MORIS_ERROR( tOutputIndices.numel() == tOutputCriteria.size(),
                    "SOL_Warehouse::create_time_solvers(), Number of output indices and criteria must be the same" );

            for ( uint Ia = 0; Ia < tOutputCriteria.size(); Ia++ )
            {
                // if a criterion for the output has been defined
                if ( tOutputCriteria( Ia ) != "Default_Output_Criterion" )
                {
                    Pointer_Function tCriteriaFunc = mLibrary->load_function< Pointer_Function >( tOutputCriteria( Ia ) );
                    mTimeSolvers( Ik )->set_output( tOutputIndices( Ia ),
                            reinterpret_cast< bool ( * )( moris::tsa::Time_Solver* ) >( tCriteriaFunc ) );
                }
                // use default criterion (output after every time step) for output
                else
                {
                    mTimeSolvers( Ik )->set_output( tOutputIndices( Ia ), Default_Output_Criterion );
                }
            }
        }

        // Set pause function
        std::string tPauseFunctionName = mParameterlist( 5 )( Ik ).get< std::string >( "TSA_Pause_Function" );
        if ( not tPauseFunctionName.empty() )
        {
            mTimeSolvers( Ik )->set_pause_function( mLibrary->load_function< Void_Function >( tPauseFunctionName ) );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void SOL_Warehouse::get_default_secondary_dof_types(
        Vector< Vector< MSI::Dof_Type > >&        aCellOfCellsSecDofTypes,
        Vector< Vector< MSI::Dof_Type > > const & aCellOfCellDofTypes )
{
    // reset list of secondary dof types
    for ( auto& tCell : aCellOfCellsSecDofTypes )
    {
        tCell.clear();
    }
    aCellOfCellsSecDofTypes.clear();

    Vector< MSI::Dof_Type > tListOfAllPossibleDofTypes;

    // get all possible used dof types from the time solver parameter list
    for ( uint Ik = 0; Ik < mParameterlist( 5 ).size(); Ik++ )
    {
        Vector< Vector< MSI::Dof_Type > > tCellOfCells;

        map< std::string, enum MSI::Dof_Type > tMap = MSI::get_msi_dof_type_map();

        string_to_cell_of_cell( mParameterlist( 5 )( Ik ).get< std::string >( "TSA_DofTypes" ),
                tCellOfCells,
                tMap );

        // append dof types into non unique list
        for ( uint Ii = 0; Ii < tCellOfCells.size(); Ii++ )
        {
            tListOfAllPossibleDofTypes.append( tCellOfCells( Ii ) );
        }
    }

    // make list of dof types unique
    std::sort( ( tListOfAllPossibleDofTypes.data() ).data(),
            ( tListOfAllPossibleDofTypes.data() ).data() + tListOfAllPossibleDofTypes.size() );

    auto last = std::unique( ( tListOfAllPossibleDofTypes.data() ).data(),
            ( tListOfAllPossibleDofTypes.data() ).data() + tListOfAllPossibleDofTypes.size() );

    auto pos = std::distance( ( tListOfAllPossibleDofTypes.data() ).data(), last );

    tListOfAllPossibleDofTypes.resize( pos );

    // make list of residual dof types
    Vector< MSI::Dof_Type > tResDofTypes;
    for ( uint Ik = 0; Ik < aCellOfCellDofTypes.size(); Ik++ )
    {
        tResDofTypes.append( aCellOfCellDofTypes( Ik ) );
    }

    // find dof types which are not part of the residual dof types amd put them in list
    aCellOfCellsSecDofTypes.resize( 1 );
    for ( uint Ik = 0; Ik < tListOfAllPossibleDofTypes.size(); Ik++ )
    {
        bool tIsResDofType = false;

        for ( uint Ii = 0; Ii < tResDofTypes.size(); Ii++ )
        {
            if ( tResDofTypes( Ii ) == tListOfAllPossibleDofTypes( Ik ) )
            {
                tIsResDofType = true;
            }
        }

        // if dof type is not residual dof type, add to secondary dif types
        if ( not tIsResDofType )
        {
            aCellOfCellsSecDofTypes( 0 ).push_back( tListOfAllPossibleDofTypes( Ik ) );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

// void SOL_Warehouse::create_solver_manager_dependencies()
//{
//     // Set size of nonlinear solver manager list to number of nonlinear solver managers
//     mListSolverManagerDependencies.resize( mListNonlinearSolverManagers.size() );
//
//     // Loop over all nonlinear solver manager
//     for ( uint Ik = 0 ; Ik < mListNonlinearSolverManagers.size(); Ik++ )
//     {
//         // Get the list of list of dof types in which this particular nonlinear solver manager is operating on
//         Vector< Vector< enum MSI::Dof_Type > > tNonlinearSolverManagerDofTypeList = mListNonlinearSolverManagers( Ik )->get_dof_type_list();
//
//         // Set the size of the dependency list for this nonlinear solver manager == number sub-solvers
//         mListSolverManagerDependencies( Ik ).set_size( tNonlinearSolverManagerDofTypeList.size(), 1, -1 );
//
//         // Loop over all sub system dof type lists
//         for ( uint Ii = 0 ; Ii < tNonlinearSolverManagerDofTypeList.size(); Ii++ )
//         {
//             // Loop over all nonlinear solver managers
//             for ( uint Ij = 0 ; Ij < mListNonlinearSolverManagers.size(); Ij++ )
//             {
//                 // Check that the to be compared nonlinear solver manager indices are not equal
//                 if( Ij != Ik )
//                 {
//                     // Get dof type union of the downward nonlinear solver manager
//                     Vector< enum MSI::Dof_Type > tUnionDofTypeList = mListNonlinearSolverManagers( Ij )->get_dof_type_union();
//
//                     // Check if the subsystem doftype list and the downward solver union list are the same
//                     if( tNonlinearSolverManagerDofTypeList( Ii ).data() == tUnionDofTypeList.data() )
//                     {
//                         // Check that a possible downward solver is only found once
//                         MORIS_ERROR( mListSolverManagerDependencies( Ik )( Ii, 0 ) == -1,
//                                 "SOL_Warehouse::create_solver_manager_dependencies(): Two solvers are operating on the same dof types" );
//
//                         // Add downward solver nonlinear solver manager index to list
//                         mListSolverManagerDependencies( Ik )( Ii, 0 ) = Ij;
//                     }
//                 }
//             }
//         }
//     }
// }

//---------------------------------------------------------------------------------------------------------------------

// moris::sint SOL_Warehouse::get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
//                                                                     const moris::sint aDofTypeListIndex )
//{
//     // check if index was set.
//     MORIS_ERROR( mListSolverManagerDependencies( aSolverManagerIndex )( aDofTypeListIndex, 0 ) != -1,
//             "SOL_Warehouse::get_nonlinear_solver_manager_index(): Return nonlinear solver manager -1. Check create_solver_manager_dependencies()" );
//
//     return mListSolverManagerDependencies( aSolverManagerIndex )( aDofTypeListIndex , 0 );
// }

//---------------------------------------------------------------------------------------------------------------------

// void SOL_Warehouse::create_maps()
//{
//     // Create matrix vector factory
//     sol::Matrix_Vector_Factory tMatFactory( MapType::Epetra );
//
//     // Get number nonlinear solver managers
//     moris::uint tNumberNonlinSolverManager = mListNonlinearSolverManagers.size();
//
//     // Set size of list of maps
//     mListOfFreeMaps.resize( tNumberNonlinSolverManager + 1, nullptr );
//
//     // create full map pointer on last entry on map list
//     mListOfFreeMaps( tNumberNonlinSolverManager ) = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );
//
//     // Loop over all nonlinear solver managers
//     for ( uint Ik = 0 ; Ik < tNumberNonlinSolverManager; Ik++ )
//     {
//         // Get dof type union of the downward nonlinear solver manager
//         Vector< enum MSI::Dof_Type > tUnionDofTypeList = mListNonlinearSolverManagers( Ik )->get_dof_type_union();
//
//         // FIXME Set dof type on solver interface
//         mSolverInterface->set_requested_dof_types( tUnionDofTypeList );
//
//         // Create free map for this dof type list
//         mListOfFreeMaps( Ik ) = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );
//     }
// }

//---------------------------------------------------------------------------------------------------------------------

// Dist_Map * SOL_Warehouse::get_list_of_maps( const moris::sint aSolverManagerIndex )
//{
//     // return map for this nonlinear solver manager index
//     return mListOfFreeMaps( aSolverManagerIndex );
// }

//---------------------------------------------------------------------------------------------------------------------

// Dist_Map * SOL_Warehouse::get_full_maps( )
//{
//     // return map for this nonlinear solver manager index
//     moris::uint tMapListSize = mListOfFreeMaps.size();
//     return mListOfFreeMaps( tMapListSize - 1 );
// }

//---------------------------------------------------------------------------------------------------------------------

// void SOL_Warehouse::finalize()
//{
//     // Call to calculate nonlinear solver manager dependencies
//     this->create_solver_manager_dependencies();
//
//     // Build matrix vector factory
//     sol::Matrix_Vector_Factory    tMatFactory( MapType::Epetra );
//
//     mMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );
//
//     // Create Full Vector
//     mFullVector = tMatFactory.create_vector( mSolverInterface, mMap, VectorType::FREE );
//
//     // Initialize full vector with zeros
//     mFullVector->vec_put_scalar( 0.0 );
//
//     // Loop over all nonlinear solver managers and set pointer to the database
//     for ( uint Ik = 0; Ik < mListNonlinearSolverManagers.size(); Ik++ )
//     {
//         mListNonlinearSolverManagers( Ik )->set_solver_warehouse( this );
//     }
// }

//---------------------------------------------------------------------------------------------------------------------

// void SOL_Warehouse::set_nonlinear_solver_managers( Nonlinear_Solver * aNonlinearSolverManager )
//{
//     if( mCallCounter == 0 )
//     {
//         mListNonlinearSolverManagers.clear();
//
//         mListNonlinearSolverManagers.push_back( aNonlinearSolverManager );
//
//         mListNonlinearSolverManagers( mCallCounter )->set_nonlinear_solver_manager_index( mCallCounter );
//     }
//     else
//     {
//         mListNonlinearSolverManagers.push_back( aNonlinearSolverManager );
//
//         mListNonlinearSolverManagers( mCallCounter )->set_nonlinear_solver_manager_index( mCallCounter );
//     }
//
//     mCallCounter = mCallCounter + 1;
// }

//---------------------------------------------------------------------------------------------------------------------

// void SOL_Warehouse::solve()
//{
//     mListNonlinearSolverManagers( 0 )->solve( );
// }

//---------------------------------------------------------------------------------------------------------------------
