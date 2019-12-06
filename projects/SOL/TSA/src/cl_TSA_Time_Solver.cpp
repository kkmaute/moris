/*
 * cl_TSA_Time_Solver.cpp
 *
 *  Created on: Okt 10, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Map_Class.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"

#include "cl_SOL_Warehouse.hpp"

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"

using namespace moris;
using namespace tsa;

    Time_Solver::Time_Solver( const enum TimeSolverType aTimeSolverType ) : mTimeSolverType( aTimeSolverType )
    {
        // create solver factory
        Time_Solver_Factory  tTimeSolFactory;

        // create solver object
        std::shared_ptr< Time_Solver_Algorithm > tTimeSolver = tTimeSolFactory.create_time_solver( aTimeSolverType );

//        for( auto tTimeSolverAlgorithmList : mTimeSolverAlgorithmList )
//        {
////            delete tTimeSolverAlgorithmList;
//            tTimeSolverAlgorithmList = nullptr;
//        }
        mTimeSolverAlgorithmList.clear();

        mTimeSolverAlgorithmList.resize( 1, nullptr );

        mTimeSolverAlgorithmList( 0 ) = tTimeSolver;

        mDofTypeList.resize( 0 );

        this->set_time_solver_parameters();
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

    void Time_Solver::set_solver_warehouse( NLA::SOL_Warehouse * aSolverWarehouse )
    {
        mSolverWarehouse = aSolverWarehouse;

        mSolverInterface = mSolverWarehouse->get_solver_interface() ;
    }

    //-------------------------------------------------------------------------------------------------------

    void Time_Solver::set_output( const uint aOutputIndex,
                                        Output_Criteria aOutputCriteria)
    {
        mOutputIndices.push_back( aOutputIndex );
        mOutputCriteriaPointer.push_back( aOutputCriteria );
    }

    //-------------------------------------------------------------------------------------------------------

    void Time_Solver::check_for_outputs()
    {
         uint tCounter = 0;

        for( Output_Criteria tOutputCriterias : mOutputCriteriaPointer )
        {
            bool tIsOutput = tOutputCriterias( this );

            if( tIsOutput )
            {
                mSolverInterface->initiate_output( mOutputIndices( tCounter ) );
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
    }

    //--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::solve()
    {
        mIsMasterTimeSolver = true;

        mSolverInterface = mSolverWarehouse->get_solver_interface();

        // create map object
        Matrix_Vector_Factory tMatFactory( MapType::Epetra );
        mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map());

        // full vector and prev full vector
        mFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, VectorType::FREE );
        mFullVector->vec_put_scalar( 0.0 );

        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mTimeSolverAlgorithmList( 0 )->set_time_solver( this );

        mTimeSolverAlgorithmList( 0 )->solve( mFullVector );
    }
//--------------------------------------------------------------------------------------------------------------------------
    void Time_Solver::set_time_solver_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListTimeSolver.insert( "NLA_max_non_lin_solver_restarts" , 0 );
    }

//--------------------------------------------------------------------------------------------------------------------------

    void Time_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
    {
        mFullVector->extract_copy( LHSValues );
    }

