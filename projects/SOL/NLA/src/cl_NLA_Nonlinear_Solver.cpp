/*
 * cl_NLA_Nonlinear_Solver.cpp
 *
 *  Created on: Okt 10, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_NLA_Nonlinear_Database.hpp"

#include "cl_Logger.hpp"

using namespace moris;
using namespace NLA;

    Nonlinear_Solver::Nonlinear_Solver( const enum NonlinearSolverType aNonLinSolverType ) : mNonLinSolverType( aNonLinSolverType )
    {
        // create solver factory
        Nonlinear_Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Nonlinear_Algorithm > tNonLinSolver = tSolFactory.create_nonlinear_solver( aNonLinSolverType );

        mNonLinearSolverList.clear();

        mNonLinearSolverList.resize( 1, nullptr );

        mNonLinearSolverList( 0 ) = tNonLinSolver;

        mStaggeredDofTypeList.resize( 0 );

        this->set_nonlinear_solver_manager_parameters();
    }

    //--------------------------------------------------------------------------------------------------
    Nonlinear_Solver::Nonlinear_Solver(       moris::Cell< std::shared_ptr<Nonlinear_Algorithm > > & aNonlinerSolverList,
                                                        const enum NonlinearSolverType                            aNonLinSolverType ) : mNonLinSolverType( aNonLinSolverType )
    {
        mNonLinearSolverList = aNonlinerSolverList;

        this->set_nonlinear_solver_manager_parameters();
    }

    //--------------------------------------------------------------------------------------------------
    Nonlinear_Solver::~Nonlinear_Solver()
    {}

    //--------------------------------------------------------------------------------------------------
    void Nonlinear_Solver::set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
                                                      const moris::sint                       aLevel )
    {
        mStaggeredDofTypeList.push_back( aStaggeredDofTypeList );
    }

    //--------------------------------------------------------------------------------------------------
    void Nonlinear_Solver::set_nonlinear_solver( std::shared_ptr< Nonlinear_Algorithm > aNonLinSolver )
    {
        if( mCallCounter == 0 )
        {
            // removes all elements from the Cell and destroy them
            mNonLinearSolverList.clear();

            // Resize the Cell to size = 1
            mNonLinearSolverList.resize( 1, nullptr );

            // Set nonlinear solver on first entry
            mNonLinearSolverList( 0 ) =  aNonLinSolver;
        }
        else
        {
            // set nonlinear solver on next entry
            mNonLinearSolverList.push_back( aNonLinSolver );
        }

        mCallCounter = mCallCounter + 1;
    }

    //-------------------------------------------------------------------------------------------------------

    void Nonlinear_Solver::set_nonlinear_solver(       std::shared_ptr< Nonlinear_Algorithm > aNonLinSolver,
                                                         const moris::uint                         aListEntry )
    {
        // Check if list is smaller than given entry
        if( aListEntry >= mNonLinearSolverList.size() )
        {
            // Resize to new entry value and set nullptr on new entries
            mNonLinearSolverList.resize( aListEntry + 1, nullptr );
        }
        // Set nonlinear solver on entry
        mNonLinearSolverList( aListEntry ) = aNonLinSolver;
    }

    //-------------------------------------------------------------------------------------------------------

    moris::Cell< enum MSI::Dof_Type > Nonlinear_Solver::get_dof_type_union()
    {
        moris::sint tCounter = 0;

        // Loop over all dof type lists to determine the total number of dof types
        for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
        {
            tCounter = tCounter + mStaggeredDofTypeList( Ik ).size();
        }

        // Create list of dof types with earlier determines size
        moris::Cell< enum MSI::Dof_Type > tUnionEnumList( tCounter );
        tCounter = 0;

        // Loop over all dof types. Add them to union list
        for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
        {
            for ( moris::uint Ii = 0; Ii < mStaggeredDofTypeList( Ik ).size(); ++Ii )
            {
                tUnionEnumList( tCounter++ ) = mStaggeredDofTypeList( Ik )( Ii );
            }
        }

        return tUnionEnumList;
    }

    //--------------------------------------------------------------------------------------------------

    void Nonlinear_Solver::set_sonlinear_solver_manager_index( const moris::sint aNonlinearSolverManagerIndex )
    {
        mNonlinearSolverManagerIndex = aNonlinearSolverManagerIndex;
    }

    //--------------------------------------------------------------------------------------------------

    moris::sint Nonlinear_Solver::get_sonlinear_solver_manager_index()
    {
        MORIS_ERROR( mNonlinearSolverManagerIndex != -1,
                "Nonlinear_Solver_Manager::get_sonlinear_solver_manager_index(): mNonlinearSolverManagerIndex = -1. Solver manager index not set." );
        return mNonlinearSolverManagerIndex;
    }

    //--------------------------------------------------------------------------------------------------

    void Nonlinear_Solver::set_nonlinear_manager( Nonlinear_Database * aNonlinearDatabase )
    {
        mNonlinearDatabase = aNonlinearDatabase;

        mSolverInput = mNonlinearDatabase->get_solver_interface() ;
    }

    //-------------------------------------------------------------------------------------------------------

    void Nonlinear_Solver::solve()
    {
        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInput->set_requested_dof_types( tDofTypeUnion );

        if ( mNonLinSolverType == NonlinearSolverType::NLBGS_SOLVER )
        {
            mNonlinerProblem = new Nonlinear_Problem( mSolverInput, mNonlinearDatabase->get_full_vector(), mNonlinearSolverManagerIndex,  false );
        }
        else
        {
            mNonlinerProblem = new Nonlinear_Problem( mSolverInput, mNonlinearDatabase->get_full_vector(), mNonlinearSolverManagerIndex );
        }

        mNonLinearSolverList( 0 )->set_nonlinear_solver_manager( this );

        mNonLinearSolverList( 0 )->solver_nonlinear_system( mNonlinerProblem );
    }

    void Nonlinear_Solver::solve( Nonlinear_Problem * aNonlinearProblem )
    {
        moris::sint tErrorStatus = 0;
        moris::sint tMaxNumLinRestarts  = mParameterListNonLinearSolver.get< moris::sint >( "NLA_max_non_lin_solver_restarts" );
        moris::sint tTryRestartOnFailIt = 1;

        MORIS_ERROR( mNonLinSolverType != NonlinearSolverType::NLBGS_SOLVER, "Nonlinear_Solver::solve(); Nonliner Solver is NLBGS" );

        mNonLinearSolverList( 0 )->set_nonlinear_solver_manager( this );

        mNonLinearSolverList( 0 )->solver_nonlinear_system( aNonlinearProblem );

        // Restart the nonlinear solver
        while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && ( moris::sint )mNonLinearSolverList.size() <= tMaxNumLinRestarts )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( " ... Previous nonlinear solve failed. Trying restart %i of %i\n", tTryRestartOnFailIt, tMaxNumLinRestarts);
            }

            // Re-solve scaled linear system with current solution as an initial guess
            //tErrorStatus = mNonLinearSolverList( tTryRestartOnFailIt )->solver_nonlinear_system( aNonlinearProblem );
            mNonLinearSolverList( tTryRestartOnFailIt )->solver_nonlinear_system( aNonlinearProblem );

            // Iterate TryRestartOnFailIt counter
            tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Nonlinear_Solver::set_nonlinear_solver_manager_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListNonLinearSolver.insert( "NLA_max_non_lin_solver_restarts" , 0 );
    }

