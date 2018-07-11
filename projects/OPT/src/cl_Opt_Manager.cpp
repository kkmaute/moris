// Project header files
#include "cl_Opt_Manager.hpp"
#include "Opt_Input.hpp"

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        OptManager::OptManager( OptProb & aOptProb )
        {
            mOptProb = aOptProb; // set the member variable mOptProb to aOptProb

            // call to input file to define the solution strategy i.e. decide
            // the optimization algorithms to be used
            mOptAlgCell = opt_input::define_opt_sol_strategy( );
        }

        // ---------------------------------------------------------------------

        OptManager::~OptManager()
        {
        }

        // ---------------------------------------------------------------------

        void OptManager::solve_opt_system()
        {
            for ( uint i=0; i < mOptAlgCell.size(); ++i )
            {
                // solve the optimization problem based on the optimization
                // algorithm Cell
                mOptAlgCell(i).solve( mOptProb );

                // scale the solution of the optimization problem
                mOptProb.scale_opt_sol();

                // update the optimization problem
                mOptProb.update_opt_prob();
            }
        }

    }
}
