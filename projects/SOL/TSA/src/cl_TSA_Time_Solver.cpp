/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver.cpp
 *
 */
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"

#include "fn_Parsing_Tools.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include "HDF5_Tools.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

using namespace moris;
using namespace tsa;

//-------------------------------------------------------------------------------

Time_Solver::Time_Solver( const enum TimeSolverType aTimeSolverType )
        : mTimeSolverType( aTimeSolverType )
{
    // create solver factory
    Time_Solver_Factory tTimeSolFactory;

    // create solver object
    std::shared_ptr< Time_Solver_Algorithm > tTimeSolver = tTimeSolFactory.create_time_solver( aTimeSolverType );

    mTimeSolverAlgorithmList.clear();

    mTimeSolverAlgorithmList.resize( 1, nullptr );

    mTimeSolverAlgorithmList( 0 ) = tTimeSolver;

    mDofTypeList.resize( 0 );

    this->set_time_solver_parameters();
}

//--------------------------------------------------------------------------------------------------

Time_Solver::Time_Solver(
        const ParameterList            aParameterlist,
        sol::SOL_Warehouse*            aSolverWarehouse,
        const enum tsa::TimeSolverType aTimeSolverType )
        : mParameterListTimeSolver( aParameterlist )
        , mSolverWarehouse( aSolverWarehouse )
        , mSolverInterface( mSolverWarehouse->get_solver_interface() )
        , mTimeSolverType( aTimeSolverType )
{
    mDofTypeList.resize( 0 );

    this->initialize_time_levels();
}

//--------------------------------------------------------------------------------------------------

Time_Solver::Time_Solver(
        moris::Cell< std::shared_ptr< Time_Solver_Algorithm > >& aTimeSolverList,
        const enum TimeSolverType                                aTimeSolverType )
        : mTimeSolverType( aTimeSolverType )
{
    mTimeSolverAlgorithmList = aTimeSolverList;

    this->set_time_solver_parameters();
}

//--------------------------------------------------------------------------------------------------

Time_Solver::~Time_Solver()
{
    this->delete_pointers();
}

//--------------------------------------------------------------------------------------------------

void
Time_Solver::delete_pointers()
{
    if ( mIsMasterTimeSolver )
    {
        for ( auto tFullSolVec : mFullVector )
        {
            delete tFullSolVec;
        }

        mFullVector.clear();

        for ( auto tFullSolVec : mFullEigenVector )
        {
            delete tFullSolVec;
        }

        mFullEigenVector.clear();

        for ( auto tFullSolVec : mFullVectorSensitivity )
        {
            delete tFullSolVec;
        }

        mFullVectorSensitivity.clear();

        mTimeFrames.clear();

        delete mFullMap;
    }
}

//--------------------------------------------------------------------------------------------------

void
Time_Solver::set_dof_type_list(
        const moris::Cell< enum MSI::Dof_Type > aDofTypeList,
        const moris::sint                       aLevel )
{
    mDofTypeList.push_back( aDofTypeList );
}

//--------------------------------------------------------------------------------------------------

void
Time_Solver::set_time_solver_algorithm( std::shared_ptr< Time_Solver_Algorithm > aTimeSolver )
{
    if ( mCallCounter == 0 )
    {
        // removes all elements from the Cell and destroy them
        mTimeSolverAlgorithmList.clear();

        // Resize the Cell to size = 1
        mTimeSolverAlgorithmList.resize( 1, nullptr );

        // Set nonlinear solver on first entry
        mTimeSolverAlgorithmList( 0 ) = aTimeSolver;
    }
    else
    {
        // set nonlinear solver on next entry
        mTimeSolverAlgorithmList.push_back( aTimeSolver );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::set_time_solver_algorithm(
        std::shared_ptr< Time_Solver_Algorithm > aTimeSolver,
        const moris::uint                        aListEntry )
{
    // Check if list is smaller than given entry
    if ( aListEntry >= mTimeSolverAlgorithmList.size() )
    {
        // Resize to new entry value and set nullptr on new entries
        mTimeSolverAlgorithmList.resize( aListEntry + 1, nullptr );
    }

    // Set nonlinear solver on entry
    mTimeSolverAlgorithmList( aListEntry ) = aTimeSolver;
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::set_sub_time_solver( Time_Solver* aTimeSolver )
{
    if ( mCallCounterTimeSolver == 0 )
    {
        // removes all elements from the Cell and destroy them
        mTimeSubSolverList.clear();

        // Resize the Cell to size = 1
        mTimeSubSolverList.resize( 1, nullptr );

        // Set nonlinear solver on first entry
        mTimeSubSolverList( 0 ) = aTimeSolver;
    }
    else
    {
        // set nonlinear solver on next entry
        mTimeSubSolverList.push_back( aTimeSolver );
    }

    mCallCounterTimeSolver = mCallCounterTimeSolver + 1;
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::set_sub_time_solver(
        Time_Solver*      aTimeSolver,
        const moris::uint aListEntry )
{
    // Check if list is smaller than given entry
    if ( aListEntry >= mTimeSubSolverList.size() )
    {
        // Resize to new entry value and set nullptr on new entries
        mTimeSubSolverList.resize( aListEntry + 1, nullptr );
    }

    // Set nonlinear solver on entry
    mTimeSubSolverList( aListEntry ) = aTimeSolver;
}

//-------------------------------------------------------------------------------------------------------

moris::Cell< enum MSI::Dof_Type >
Time_Solver::get_dof_type_union()
{
    moris::sint tCounter = 0;

    // Loop over all dof type lists to determine the total number of dof types
    for ( moris::uint Ik = 0; Ik < mDofTypeList.size(); ++Ik )
    {
        tCounter = tCounter + mDofTypeList( Ik ).size();
    }

    // Create list of dof types with earlier determines size
    moris::Cell< enum MSI::Dof_Type > tUnionEnumList( tCounter );
    tCounter = 0;

    // Loop over all dof types. Add them to union list
    for ( moris::uint Ik = 0; Ik < mDofTypeList.size(); ++Ik )
    {
        for ( moris::uint Ii = 0; Ii < mDofTypeList( Ik ).size(); ++Ii )
        {
            tUnionEnumList( tCounter++ ) = mDofTypeList( Ik )( Ii );
        }
    }

    return tUnionEnumList;
}

//--------------------------------------------------------------------------------------------------

void
Time_Solver::set_solver_warehouse( sol::SOL_Warehouse* aSolverWarehouse )
{
    mSolverWarehouse = aSolverWarehouse;

    mSolverInterface = mSolverWarehouse->get_solver_interface();
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::set_output(
        const uint      aOutputIndex,
        Output_Criteria aOutputCriteria )
{
    mOutputIndices.push_back( aOutputIndex );
    mOutputCriteriaPointer.push_back( aOutputCriteria );
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::check_for_outputs(
        const moris::real& aTime,
        const bool         aEndOfTimeIteration )
{
    if ( mIsForwardSolve )
    {
        uint tCounter = 0;

        // loop over all outputs and check if it is triggered
        for ( Output_Criteria tOutputCriterias : mOutputCriteriaPointer )
        {
            bool tIsOutput = false;

            if ( tOutputCriterias == nullptr )
            {
                tIsOutput = true;
            }
            else
            {
                tIsOutput = tOutputCriterias( this );
            }

            if ( tIsOutput )
            {
                mSolverInterface->initiate_output( mOutputIndices( tCounter ), aTime, aEndOfTimeIteration );
            }

            tCounter++;
        }
    }
}

//-------------------------------------------------------------------------------------------------------

void
Time_Solver::solve( moris::Cell< sol::Dist_Vector* >& aFullVector )
{
    moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

    mSolverInterface->set_requested_dof_types( tDofTypeUnion );

    mFullVector = aFullVector;

    mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

    mTimeSolverAlgorithmList( 0 )->solve( aFullVector );

    // this->check_for_outputs();
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::solve()
{
    Tracer tTracer( "TimeSolver", "Forward Analysis", "Solve" );

    // flags if thats the master time solver and if this is a forward solve
    mIsMasterTimeSolver = true;
    mIsForwardSolve     = true;

    // delete pointers
    this->delete_pointers();

    // get solver interface
    mSolverInterface = mSolverWarehouse->get_solver_interface();

    mSolverInterface->set_is_forward_analysis();

    // create map object
    sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

    // build map for full vector
    mFullMap = tMatFactory.create_full_map(
            mSolverInterface->get_my_local_global_map(),
            mSolverInterface->get_my_local_global_overlapping_map() );

    // get number of RHS
    uint tNumRHMS = mSolverInterface->get_num_rhs();

    // set size for full solution vector on time step 0 and previous solution vector on timestep -1.
    mFullVector.resize( 2, nullptr );

    // full vector and prev full vector
    mFullVector( 0 ) = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );
    mFullVector( 1 ) = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

    // set time level 0 sol vec to interface
    mSolverInterface->set_solution_vector( mFullVector( 1 ) );
    mSolverInterface->set_solution_vector_prev_time_step( mFullVector( 0 ) );

    // get number of eigen vectors
    uint tNumEigenVectors = mSolverInterface->get_num_eigen_vectors();

    // create full vector for eigen vectors
    if ( tNumEigenVectors > 0 )
    {
        mFullEigenVector.resize( 1, nullptr );
        mFullEigenVector( 0 ) = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumEigenVectors );

        // set eigen vector in interface
        mSolverInterface->set_eigen_solution_vector( mFullEigenVector( 0 ) );
    }

    // initialize solution vector and prev solution vector
    this->initialize_sol_vec();
    this->initialize_prev_sol_vec();

    moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

    mSolverInterface->set_requested_dof_types( tDofTypeUnion );

    mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

    if ( not mSolverWarehouse->get_load_sol_vec_from_file().empty() )
    {
        std::string tSolVecPath = mSolverWarehouse->get_load_sol_vec_from_file();

        // detect file type
        std::string tType = tSolVecPath.substr( tSolVecPath.find_last_of( "." ) + 1, tSolVecPath.length() );

        if ( tType == "hdf5" || tType == "h5" )
        {
            // read solution vector from file
            mFullVector( 1 )->read_vector_from_HDF5( tSolVecPath.c_str() );
        }
        else
        {
            MORIS_ERROR( false, "Time_Solver::solve(), Solution vector input type unknown. New types can be implemented here." );
        }

        Matrix< DDRMat > tTime_0 = { { 0.0 }, { 0.0 } };
        Matrix< DDRMat > tTime_1 = { { 0.0 }, { 1.0 } };

        mSolverInterface->set_previous_time( tTime_0 );
        mSolverInterface->set_time( tTime_1 );

        mSolverInterface->compute_IQI();

        // input second time slap value for output
        this->check_for_outputs( 1.0, true );
    }
    else
    {
        // SOLVE CALL
        mTimeSolverAlgorithmList( 0 )->solve( mFullVector );

        // output solution vector to file
        if ( not mSolverWarehouse->get_save_final_sol_vec_to_file().empty() )
        {
            std::string tSolVecPath = mSolverWarehouse->get_save_final_sol_vec_to_file();

            // detect file type
            std::string tType = tSolVecPath.substr( tSolVecPath.find_last_of( "." ) + 1, tSolVecPath.length() );

            if ( tType == "hdf5" || tType == "h5" )
            {
                mFullVector( 1 )->save_vector_to_HDF5( tSolVecPath.c_str() );
            }
            else
            {
                MORIS_ERROR( false, "Time_Solver::solve(), Solution vector output type unknown. New types can be implemented here." );
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::solve_sensitivity()
{
    Tracer tTracer( "TimeSolver", "Sensitivity Analysis", "Solve" );
    mIsForwardSolve = false;

    // create map object
    sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

    uint tNumRHMS = mSolverInterface->get_num_rhs();

    mSolverInterface->set_is_sensitivity_analysis();

    // full vector and previous full vector
    mFullVectorSensitivity.resize( 2, nullptr );
    mFullVectorSensitivity( 0 ) = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );
    mFullVectorSensitivity( 1 ) = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

    mFullVectorSensitivity( 0 )->vec_put_scalar( 0.0 );
    mFullVectorSensitivity( 1 )->vec_put_scalar( 0.0 );

    mSolverInterface->set_adjoint_solution_vector( mFullVectorSensitivity( 0 ) );
    mSolverInterface->set_previous_adjoint_solution_vector( mFullVectorSensitivity( 1 ) );

    moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

    mSolverInterface->set_requested_dof_types( tDofTypeUnion );

    mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

    mTimeSolverAlgorithmList( 0 )->solve( mFullVectorSensitivity );

    // output solution vector to file
    if ( not mSolverWarehouse->get_save_final_adjoint_vec_to_file().empty() )
    {
        std::string tSolVecPath = mSolverWarehouse->get_save_final_adjoint_vec_to_file();

        // detect file type
        std::string tType = tSolVecPath.substr( tSolVecPath.find_last_of( "." ) + 1, tSolVecPath.length() );

        if ( tType == "hdf5" || tType == "h5" )
        {
            mFullVectorSensitivity( 0 )->save_vector_to_HDF5( tSolVecPath.c_str() );
        }
        else
        {
            MORIS_ERROR( false, "Time_Solver::solve_sensitivity(), Solution vector output type unknown. New types can be implemented here." );
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::initialize_sol_vec()
{
    // extract initialization string from parameter list
    moris::Cell< moris::Cell< std::string > > tDofTypeAndValuePair;

    // get string for initial guess from parameter list
    std::string tStrInitialGuess = mParameterListTimeSolver.get< std::string >( "TSA_Initialize_Sol_Vec" );

    // use string to determine initial guess
    // if parameter is not given, initialize with all zeros
    if ( tStrInitialGuess.size() == 0 )
    {
        mFullVector( 1 )->vec_put_scalar( 0.0 );
    }

    // if user defines parameter
    else
    {
        // check if file name with type is specified
        std::string tType = tStrInitialGuess.substr( tStrInitialGuess.find_last_of( "." ) + 1, tStrInitialGuess.length() );

        // if .hdf5 is defined as initial guess
        if ( tType == "hdf5" || tType == "h5" )
        {
            // log/print that the initial guess is read from file
            MORIS_LOG_INFO( "Reading initial guess for solution vector from file: ", tStrInitialGuess.c_str() );

            // FIXME: this option doesn't work in parallel, only for serial debugging purposes
            MORIS_ERROR( par_size() == 1, "Time_Solver::initialize_sol_vec() - Restarting from hdf5 file only possible in serial." );

            // read HDF5 file to moris matrix
            hid_t            tFileID = open_hdf5_file( tStrInitialGuess );
            herr_t           tStatus = 0;
            Matrix< DDRMat > tInitialGuess;
            load_matrix_from_hdf5_file( tFileID, "SolVec", tInitialGuess, tStatus );

            // get length of solution vector
            uint tSolVecLength = tInitialGuess.length();

            // initialize trivial map
            Matrix< DDSMat > tPos( tSolVecLength, 1, 0 );

            // create trivial map
            for ( uint Ik = 0; Ik < tSolVecLength; Ik++ )
            {
                tPos( Ik ) = Ik;
            }

            // populate solution vector
            mFullVector( 1 )->replace_global_values( tPos, tInitialGuess );
        }

        // else, assume matrix is specified and convert string to cell/matrix
        else
        {
            // convert the user defined string to cell/matrix
            string_to_cell_of_cell( tStrInitialGuess, tDofTypeAndValuePair );

            // create map object
            sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

            sol::Dist_Map* tFeeMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );

            uint tNumRHMS = mSolverInterface->get_num_rhs();

            // full vector and previous full vector
            sol::Dist_Vector* tFreeVector = tMatFactory.create_vector( mSolverInterface, tFeeMap, tNumRHMS );

            // initialize solution vector
            tFreeVector->vec_put_scalar( 0.0 );

            // get string to dof type map
            map< std::string, enum MSI::Dof_Type > tDofTypeMap = MSI::get_msi_dof_type_map();

            // loop over input dof types
            for ( uint Ik = 0; Ik < tDofTypeAndValuePair.size(); Ik++ )
            {
                // First string is dof type
                moris::Cell< enum MSI::Dof_Type > tDofType = { tDofTypeMap.find( tDofTypeAndValuePair( Ik )( 0 ) ) };

                // get local global ids for this dof type
                moris::Matrix< IdMat > tAdofIds = mSolverInterface->get_my_local_global_map( tDofType );

                // if one dof type is set to random all dof entries are set to random
                if ( tDofTypeAndValuePair( Ik )( 1 ) == "random" )
                {
                    tFreeVector->random();
                    break;
                }

                // get value from input
                moris::real tValue = std::stod( tDofTypeAndValuePair( Ik )( 1 ) );

                // add value into Matrix
                moris::Matrix< DDRMat > tValues( tAdofIds.numel(), 1, tValue );

                if ( tDofTypeAndValuePair( Ik ).size() == 2 )
                {
                    // replace initial values in solution vector
                    tFreeVector->replace_global_values(
                            tAdofIds,
                            tValues );
                }
                else if ( tDofTypeAndValuePair( Ik ).size() == 3 )
                {
                    // get solution vector index
                    moris::sint tVectorIndex = std::stoi( tDofTypeAndValuePair( Ik )( 2 ) );

                    // replace initial values in solution vector
                    tFreeVector->replace_global_values(
                            tAdofIds,
                            tValues,
                            (uint)tVectorIndex );
                }
                else
                {
                    MORIS_ERROR( false, " Time_Solver::initialize_sol_vec(), TSA_Initialize_Sol_Vec input not correct" );
                }
            }

            mFullVector( 1 )->import_local_to_global( *tFreeVector );

            delete ( tFreeVector );
            delete tFeeMap;
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::initialize_prev_sol_vec()
{
    // initialize prev solution vector with zero( time level -1 )
    mFullVector( 0 )->vec_put_scalar( 0.0 );

    //----------------------------------------------
    // get file name for initial guess from parameter list
    std::string tStrInitSol = mSolverWarehouse->get_TSA_initial_guess_input_filename();

    // check if user defined previous solution through file
    // FIXME: this is useless right now, as the time continuity IWG grabs the initial condition from the properties anyways
    if ( tStrInitSol.size() > 0 )
    {
        // log/print that the initial guess is read from file
        MORIS_LOG_INFO( "Reading initial solution from file: ", tStrInitSol.c_str() );

        // FIXME: this option doesn't work in parallel, only for serial debugging purposes
        MORIS_ERROR( par_size() == 1, "Time_Solver::initialize_sol_vec() - Restarting from hdf5 file only possible in serial." );

        // read HDF5 file to moris matrix
        hid_t            tFileID = open_hdf5_file( tStrInitSol );
        herr_t           tStatus = 0;
        Matrix< DDRMat > tInitSol;
        load_matrix_from_hdf5_file( tFileID, "SolVec", tInitSol, tStatus );

        // get length of solution vector
        uint tSolVecLength = tInitSol.length();

        // initialize trivial map
        Matrix< DDSMat > tPos( tSolVecLength, 1, 0 );

        // create trivial map
        for ( uint Ik = 0; Ik < tSolVecLength; Ik++ )
        {
            tPos( Ik ) = Ik;
        }

        // populate solution vector
        mFullVector( 0 )->replace_global_values( tPos, tInitSol );
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::prepare_sol_vec_for_next_time_step()
{
    if ( mIsMasterTimeSolver )
    {
        // get num RHS
        uint tNumRHMS = mSolverInterface->get_num_rhs();

        // create map object
        sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

        // full vector and prev full vector
        sol::Dist_Vector* tFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

        mFullVector.push_back( tFullVector );

        uint tNumSolVec = mFullVector.size();

        mFullVector( tNumSolVec - 1 )->vec_plus_vec( 1.0, *( mFullVector( tNumSolVec - 2 ) ), 0.0 );
        //        mFullVector( tNumSolVec-1 )->vec_put_scalar( 0.0 );
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::initialize_time_levels()
{
    // extract initialization string from parameter list
    moris::Cell< moris::Cell< std::string > > tDofTypeAndTimeLevelPair;

    string_to_cell_of_cell( mParameterListTimeSolver.get< std::string >( "TSA_time_level_per_type" ),
            tDofTypeAndTimeLevelPair );

    // get string to dof type map
    map< std::string, enum MSI::Dof_Type > tDofTypeMap = MSI::get_msi_dof_type_map();

    // loop over input dof types
    for ( uint Ik = 0; Ik < tDofTypeAndTimeLevelPair.size(); Ik++ )
    {
        // First string is dof type
        enum MSI::Dof_Type tDofType = tDofTypeMap.find( tDofTypeAndTimeLevelPair( Ik )( 0 ) );

        // get value from input
        moris::uint tValue = (moris::uint)std::stoi( tDofTypeAndTimeLevelPair( Ik )( 1 ) );

        mSolverInterface->set_time_levels_for_type( tDofType, tValue );
    }
}

//--------------------------------------------------------------------------------------------------------------------------
void
Time_Solver::set_time_solver_parameters()
{
    // Maximal number of linear solver restarts on fail
    mParameterListTimeSolver = prm::create_time_solver_parameter_list();
}

//--------------------------------------------------------------------------------------------------------------------------

void
Time_Solver::get_full_solution( moris::Matrix< DDRMat >& LHSValues )
{
    // get index of last solution vector
    uint tNumSolVec = mFullVector.size() - 2;

    // extract solution values
    mFullVector( tNumSolVec )->extract_copy( LHSValues );
}
