// Project header files
#include "cl_OPT_Manager.hpp"

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Manager::Manager(Cell< Algorithm_API > aAlgorithms, Problem* aProblem) : mOptAlgCell(aAlgorithms), mProblem(aProblem)
        {
        }

        // ---------------------------------------------------------------------

        Manager::~Manager()
        {
        }

        // ---------------------------------------------------------------------

        void Manager::solve_opt_system()
        {
            // initialize the problem
            mProblem->initialize();

            for (uint i = 0; i < mOptAlgCell.size(); i++)
            {
                // solve the optimization problem based on the algorithm cell
                mOptAlgCell(i).solve(mProblem);

                // scale the solution of the optimization problem
                mProblem->scale_solution();

                // update the optimization problem
                mProblem->update_problem();
            }
        }

    }
}
