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
#include "fn_OPT_create_problem.hpp"
#include "fn_OPT_create_interface.hpp"
#include "fn_OPT_create_algorithm.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace opt
    {
        // -------------------------------------------------------------------------------------------------------------

        Manager::Manager(
                const Vector<Vector<ParameterList>>           & aParameterLists,
                Vector<std::shared_ptr<Criteria_Interface>>   aInterfaces)
                : Manager(
                        aParameterLists(2),
                        create_problem(
                                aParameterLists(0)(0),
                                create_interface(
                                        aParameterLists(1),
                                        aInterfaces)))
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        Manager::Manager(
                const Vector<ParameterList>& aAlgorithmParameterLists,
                std::shared_ptr<Problem>   aProblem)
                : mProblem(aProblem)
        {
            // Construct Algorithm cell
            uint tNumAlgorithms = aAlgorithmParameterLists.size();

            for (uint tAlgorithmIndex = 0; tAlgorithmIndex < tNumAlgorithms; tAlgorithmIndex++)
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

            // initialize the problem
            mProblem->initialize();

            for (uint i = 0; i < mAlgorithms.size(); i++)
            {
                // solve the optimization problem based on the algorithm cell
                uint tOptIteration = mAlgorithms(i)->solve(i, mProblem);

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
            if( mProblem->restart_optimization())
            {
                mAlgorithms( aI )->set_restart_index( aOptIteration - 1 );

                mProblem->initialize();

                uint tOptIteration = mAlgorithms( aI )->solve( aI, mProblem );

                this->reinitialize( aI, tOptIteration );
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Manager::get_advs()
        {
            return mProblem->get_advs();
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Manager::get_objectives()
        {
            return mProblem->get_objectives();
        }

        // -------------------------------------------------------------------------------------------------------------
    }
}

