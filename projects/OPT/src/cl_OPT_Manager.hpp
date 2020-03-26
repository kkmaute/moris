#ifndef MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_
#define MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_

// C++ header files.
#include <utility>

// MORIS project header files.
#include "cl_Cell.hpp" // CON/src
#include "core.hpp"
#include "cl_OPT_Algorithm_API.hpp" // OPT/src
#include "cl_OPT_Problem.hpp" // OPT/src

namespace moris
{
    namespace opt
    {
        class Manager
        {
        private:
            Cell< Algorithm_API > mOptAlgCell; // Cell of pointers of optimization algorithms
            Problem* mProblem;

        public:
            /**
             * Constructor
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            Manager(Cell< Algorithm_API > aAlgorithms, Problem* aProblem);

            /**
             * Destructor
             */
            ~Manager();

            /**
             * @brief Solve the Cell of optimization algorithms
             */
            void solve_opt_system( );
        };
    }  // namespace opt
}      // namespace moris

#endif  /* MORIS_OPTIMIZATION_CL_OPTDATA_HPP_ */
