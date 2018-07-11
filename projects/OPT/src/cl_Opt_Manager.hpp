#ifndef MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_
#define MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_

// C++ header files.
#include <utility>

// MORIS project header files.
#include "cl_Cell.hpp" // CON/src
#include "core.hpp"
#include "cl_Opt_Alg_API.hpp" // OPT/src
#include "cl_Opt_Prob.hpp" // OPT/src

namespace moris
{
    namespace opt
    {
        class OptManager
        {
        public:

            Cell< OptAlgAPI > mOptAlgCell; // Cell of pointers of optimization algorithms

            OptProb mOptProb; // Object of type OptProb

            /**
             * Constructor
             *
             * @param[in] aOptProb Object of type OptProb containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            OptManager( OptProb & aOptProb );

            /**
             * Destructor
             */
            ~OptManager( );

            /**
             * @brief Solve the Cell of optimization algorithms
             */
            void solve_opt_system( );
        };
    }  // namespace opt
}      // namespace moris

#endif  /* MORIS_OPTIMIZATION_CL_OPTDATA_HPP_ */
