// Project header files
#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_problem.hpp"
#include "fn_OPT_create_interface.hpp"

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Manager::Manager(moris::Cell<moris::Cell<ParameterList>>& tParameterLists)
        {
            // Problem
            mProblem = create_problem(tParameterLists(0)(0), create_interface(tParameterLists(1)));

            // Algorithm Cell
            uint tNumAlgorithms = tParameterLists(2).size();
            for (uint tAlgorithmIndex = 0; tAlgorithmIndex < tNumAlgorithms; tAlgorithmIndex++)
            {
                Algorithm_API tNewAlgorithm(tParameterLists(2)(tAlgorithmIndex));
                mAlgorithms.push_back(tNewAlgorithm);
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
                mAlgorithms(i).solve(mProblem);

                // scale the solution of the optimization problem
                mProblem->scale_solution();

                // update the optimization problem
                mProblem->update_problem();
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
