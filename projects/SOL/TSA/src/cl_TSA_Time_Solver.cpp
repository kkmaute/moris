/*
 * cl_TSA_Time_Solver.cpp
 *
 *  Created on: Okt 10, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"

#include "fn_Parsing_Tools.hpp"
#include "cl_PRM_SOL_Parameters.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"

// Detailed Logging package
//#include "cl_Tracer.hpp"

using namespace moris;
using namespace tsa;

    Time_Solver::Time_Solver( const enum TimeSolverType aTimeSolverType ) : mTimeSolverType( aTimeSolverType )
    {
        // create solver factory
        Time_Solver_Factory  tTimeSolFactory;

        // create solver object
        std::shared_ptr< Time_Solver_Algorithm > tTimeSolver = tTimeSolFactory.create_time_solver( aTimeSolverType );

        mTimeSolverAlgorithmList.clear();

        mTimeSolverAlgorithmList.resize( 1, nullptr );

        mTimeSolverAlgorithmList( 0 ) = tTimeSolver;

        mDofTypeList.resize( 0 );

        this->set_time_solver_parameters();
    }

//--------------------------------------------------------------------------------------------------

    Time_Solver::Time_Solver( const ParameterList         aParameterlist,
                                    sol::SOL_Warehouse  * aSolverWarehouse,
                              const enum TimeSolverType   aTimeSolverType ) : mParameterListTimeSolver( aParameterlist ),
                                                                              mSolverWarehouse( aSolverWarehouse ),
                                                                              mSolverInterface( mSolverWarehouse->get_solver_interface() ),
                                                                              mTimeSolverType( aTimeSolverType )
    {
        mDofTypeList.resize( 0 );

        this->initialize_time_levels();
    }

//--------------------------------------------------------------------------------------------------

    Time_Solver::Time_Solver(       moris::Cell< std::shared_ptr< Time_Solver_Algorithm > > & aTimeSolverList,
                              const enum TimeSolverType                                       aTimeSolverType ) : mTimeSolverType( aTimeSolverType )
    {
        mTimeSolverAlgorithmList = aTimeSolverList;

        this->set_time_solver_parameters();
    }

//--------------------------------------------------------------------------------------------------

    Time_Solver::~Time_Solver()
    {
        if( mIsMasterTimeSolver )
        {
            delete( mFullMap );
            delete( mFullVector );

            if ( mFullVectorSensitivity != nullptr )
            {
                delete( mFullVectorSensitivity );
            }
        }
    }

//--------------------------------------------------------------------------------------------------

    void Time_Solver::set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aDofTypeList,
                                         const moris::sint                       aLevel )
    {
        mDofTypeList.push_back( aDofTypeList );
    }

//--------------------------------------------------------------------------------------------------

    void Time_Solver::set_time_solver_algorithm( std::shared_ptr< Time_Solver_Algorithm > aTimeSolver )
    {
        if( mCallCounter == 0 )
        {
            // removes all elements from the Cell and destroy them
            mTimeSolverAlgorithmList.clear();

            // Resize the Cell to size = 1
            mTimeSolverAlgorithmList.resize( 1, nullptr );

            // Set nonlinear solver on first entry
            mTimeSolverAlgorithmList( 0 ) =  aTimeSolver;
        }
        else
        {
            // set nonlinear solver on next entry
            mTimeSolverAlgorithmList.push_back( aTimeSolver );
        }

        mCallCounter = mCallCounter + 1;
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::set_time_solver_algorithm(       std::shared_ptr< Time_Solver_Algorithm > aTimeSolver,
                                                 const moris::uint                              aListEntry )
    {
        // Check if list is smaller than given entry
        if( aListEntry >= mTimeSolverAlgorithmList.size() )
        {
            // Resize to new entry value and set nullptr on new entries
            mTimeSolverAlgorithmList.resize( aListEntry + 1, nullptr );
        }
        // Set nonlinear solver on entry
        mTimeSolverAlgorithmList( aListEntry ) = aTimeSolver;
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::set_sub_time_solver( Time_Solver * aTimeSolver )
    {
        if( mCallCounterTimeSolver == 0 )
        {
            // removes all elements from the Cell and destroy them
            mTimeSubSolverList.clear();

            // Resize the Cell to size = 1
            mTimeSubSolverList.resize( 1, nullptr );

            // Set nonlinear solver on first entry
            mTimeSubSolverList( 0 ) =  aTimeSolver;
        }
        else
        {
            // set nonlinear solver on next entry
            mTimeSubSolverList.push_back( aTimeSolver );
        }

        mCallCounterTimeSolver = mCallCounterTimeSolver + 1;
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::set_sub_time_solver(       Time_Solver * aTimeSolver,
                                           const moris::uint   aListEntry )
    {
        // Check if list is smaller than given entry
        if( aListEntry >= mTimeSubSolverList.size() )
        {
            // Resize to new entry value and set nullptr on new entries
            mTimeSubSolverList.resize( aListEntry + 1, nullptr );
        }
        // Set nonlinear solver on entry
        mTimeSubSolverList( aListEntry ) = aTimeSolver;
    }

//-------------------------------------------------------------------------------------------------------

    moris::Cell< enum MSI::Dof_Type > Time_Solver::get_dof_type_union()
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

    void Time_Solver::set_solver_warehouse( sol::SOL_Warehouse * aSolverWarehouse )
    {
        mSolverWarehouse = aSolverWarehouse;

        mSolverInterface = mSolverWarehouse->get_solver_interface() ;
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::set_output( const uint            aOutputIndex,
                                        Output_Criteria aOutputCriteria)
    {
        mOutputIndices.push_back( aOutputIndex );
        mOutputCriteriaPointer.push_back( aOutputCriteria );
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::check_for_outputs( const moris::real & aTime,
                                         const bool          aEndOfTimeIteration )
    {
         uint tCounter = 0;

         // loop over all outputs and check if it is triggered
        for( Output_Criteria tOutputCriterias : mOutputCriteriaPointer )
        {
            bool tIsOutput = false;

            if( tOutputCriterias == nullptr )
            {
                tIsOutput = true;
            }
            else
            {
                tIsOutput = tOutputCriterias( this );
            }

            if( tIsOutput )
            {
                MORIS_LOG_INFO(" Initiate output for output index %-5i", mOutputIndices( tCounter ) );

                mSolverInterface->initiate_output( mOutputIndices( tCounter ), aTime, aEndOfTimeIteration );
            }

            tCounter++;
        }
    }

//-------------------------------------------------------------------------------------------------------

    void Time_Solver::solve( Dist_Vector * aFullVector )
    {
        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mFullVector = aFullVector;

        mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

        mTimeSolverAlgorithmList( 0 )->solve( mFullVector );

        //this->check_for_outputs();
    }

    //--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::solve()
    {
        mIsMasterTimeSolver = true;

        // get solver interface
        mSolverInterface = mSolverWarehouse->get_solver_interface();

        // create map object
        Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );
        mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );

        uint tNumRHMS = mSolverInterface->get_num_rhs();

        // full vector and prev full vector
        mFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

        mSolverInterface->set_solution_vector( mFullVector );

        this->initialize_sol_vec();

        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

        mTimeSolverAlgorithmList( 0 )->solve( mFullVector );
    }

    //--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::solve_sensitivity()
    {
        // create map object
        Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

        uint tNumRHMS = mSolverInterface->get_num_rhs();

        // full vector and prev full vector
        mFullVectorSensitivity = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

        mSolverInterface->set_sensitivity_solution_vector( mFullVectorSensitivity );

        mFullVectorSensitivity->vec_put_scalar( 0.0 );

        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

        mTimeSolverAlgorithmList( 0 )->solve( mFullVectorSensitivity );
    }

//--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::initialize_sol_vec()
    {
        // initialize solution vector with zero
        mFullVector->vec_put_scalar( 0.0 );

        // extract initialization string from parameterlist
        moris::Cell< moris::Cell< std::string > > tDofTypeAndValuePair;
        string_to_cell_of_cell( mParameterListTimeSolver.get< std::string >( "TSA_Initialize_Sol_Vec" ),
                                tDofTypeAndValuePair );

        // get string to dof type map
        map< std::string, enum MSI::Dof_Type > tDofTypeMap = MSI::get_msi_dof_type_map();

        // loop over input dof types
        for( uint Ik = 0; Ik < tDofTypeAndValuePair.size(); Ik++ )
        {
            // First string is dof type
            moris::Cell< enum MSI::Dof_Type > tDofType = { tDofTypeMap.find( tDofTypeAndValuePair( Ik )( 0 ) ) };

            // get local global ids for this dof type
            moris::Matrix< IdMat > tAdofIds = mSolverInterface->get_my_local_global_map( tDofType );

            // get value from input
            moris::real tValue = std::stod( tDofTypeAndValuePair( Ik )( 1 ) );

            // add value into Matrix
            moris::Matrix< DDRMat > tValues( tAdofIds.numel(), 1, tValue );

            if( tDofTypeAndValuePair( Ik ).size() == 2 )
            {
                // replace inital values in solution vector
                mFullVector->replace_global_values( tAdofIds,
                                                    tValues );
            }
            else if( tDofTypeAndValuePair( Ik ).size() == 3 )
            {
                // get solution vector index
                moris::sint tVectorIndex = std::stoi( tDofTypeAndValuePair( Ik )( 2 ) );

                // replace inital values in solution vector
                mFullVector->replace_global_values( tAdofIds,
                                                    tValues,
                                                    ( uint ) tVectorIndex );
            }
            else
            {
                MORIS_ERROR( false, " Time_Solver::initialize_sol_vec(), TSA_Initialize_Sol_Vec input not correct" );
            }
        }
    }

//--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::initialize_time_levels()
    {
        // extract initialization string from parameterlist
        moris::Cell< moris::Cell< std::string > > tDofTypeAndTimeLevelPair;
        string_to_cell_of_cell( mParameterListTimeSolver.get< std::string >( "TSA_time_level_per_type" ),
                tDofTypeAndTimeLevelPair );

        // get string to dof type map
        map< std::string, enum MSI::Dof_Type > tDofTypeMap = MSI::get_msi_dof_type_map();

        // loop over input dof types
        for( uint Ik = 0; Ik < tDofTypeAndTimeLevelPair.size(); Ik++ )
        {
            // First string is dof type
            enum MSI::Dof_Type tDofType = tDofTypeMap.find( tDofTypeAndTimeLevelPair( Ik )( 0 ) );

            // get value from input
            moris::uint tValue = (moris::uint)std::stoi( tDofTypeAndTimeLevelPair( Ik )( 1 ) );

            mSolverInterface->set_time_levels_for_type( tDofType, tValue );
        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Time_Solver::set_time_solver_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListTimeSolver = prm::create_time_solver_parameter_list();
    }

//--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
    {
        mFullVector->extract_copy( LHSValues );
    }

