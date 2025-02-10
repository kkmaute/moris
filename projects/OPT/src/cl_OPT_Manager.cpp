/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Manager.cpp
 *
 */

#include "cl_OPT_Manager.hpp"

#include <utility>
#include "fn_OPT_create_problem.hpp"
#include "fn_OPT_create_interface.hpp"
#include "fn_OPT_create_algorithm.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::opt
{
    // -------------------------------------------------------------------------------------------------------------

    Manager::Manager(
            const Module_Parameter_Lists&                          aParameterLists,
            const Vector< std::shared_ptr< Criteria_Interface > >& aInterfaces )
            : Manager(
                      aParameterLists( 2 ),
                      create_problem(
                              aParameterLists( 0 )( 0 ),
                              create_interface(
                                      aParameterLists( 1 ),
                                      aInterfaces ) ) )
    {
    }

    // -------------------------------------------------------------------------------------------------------------

    Manager::Manager(
            const Submodule_Parameter_Lists& aAlgorithmParameterLists,
            std::shared_ptr< Problem >       aProblem )
            : mProblem( std::move( aProblem ) )
    {
        // Construct Algorithm cell
        uint tNumAlgorithms = aAlgorithmParameterLists.size();

        for ( uint tAlgorithmIndex = 0; tAlgorithmIndex < tNumAlgorithms; tAlgorithmIndex++ )
        {
            mAlgorithms.push_back( create_algorithm( aAlgorithmParameterLists( tAlgorithmIndex ) ) );
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    Manager::~Manager()
    {
    }

    // -------------------------------------------------------------------------------------------------------------

    void Manager::perform()
    {
        // Trace optimization
        Tracer tTracer( "OPT", "Manager", "Perform" );

        gLogger.set_action_data( "GlobalClock", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR, "Epsilon", 0.0 );    // BRENDAN remove this hack to set action data
        gLogger.set_action_data( "GlobalClock", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR, "ADVID", 0 );        // BRENDAN remove this hack to set action data

        // initialize the problem
        mProblem->initialize();

        for ( uint i = 0; i < mAlgorithms.size(); i++ )
        {
            // solve the optimization problem based on the algorithm cell
            uint tOptIteration = mAlgorithms( i )->solve( i, mProblem );

            this->reinitialize( i, tOptIteration );

            // scale the solution of the optimization problem
            mProblem->scale_solution();

            // update the optimization problem
            mProblem->update_problem();
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    void Manager::reinitialize(
            uint aI,
            uint aOptIteration )
    {
        if ( mProblem->restart_optimization() )
        {
            mAlgorithms( aI )->set_restart_index( aOptIteration - 1 );

            mProblem->initialize();

            uint tOptIteration = mAlgorithms( aI )->solve( aI, mProblem );

            this->reinitialize( aI, tOptIteration );
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    Vector< real > Manager::get_advs()
    {
        return mProblem->get_advs();
    }

    // -------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Manager::get_objectives()
    {
        return mProblem->get_objectives();
    }

    // -------------------------------------------------------------------------------------------------------------
}    // namespace moris::opt
