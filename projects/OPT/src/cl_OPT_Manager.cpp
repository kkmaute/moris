#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_problem.hpp"
#include "fn_OPT_create_interface.hpp"
#include "fn_OPT_create_algorithm.hpp"

namespace moris
{
    namespace opt
    {
        // -------------------------------------------------------------------------------------------------------------

        Manager::Manager(
                const Cell<Cell<ParameterList>>           & aParameterLists,
                Cell<std::shared_ptr<Criteria_Interface>>   aInterfaces)
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
                const Cell<ParameterList>& aAlgorithmParameterLists,
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
            // initialize the problem
            mProblem->initialize();

            for (uint i = 0; i < mAlgorithms.size(); i++)
            {
                // solve the optimization problem based on the algorithm cell
                uint tOptIteration = mAlgorithms(i)->solve(i, mProblem);

                this->restart_with_remesh( i, tOptIteration );

                // scale the solution of the optimization problem
                mProblem->scale_solution();

                // update the optimization problem
                mProblem->update_problem();
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Manager::restart_with_remesh(
                uint aI,
                uint aOptIteration )
        {
            if( mProblem->restart_optimization())
            {
                mProblem->initialize();

                mAlgorithms( aI )->set_restart_index( aOptIteration - 1 );

                uint tOptIteration = mAlgorithms( aI )->solve( aI, mProblem );

                this->restart_with_remesh( aI, tOptIteration );
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
