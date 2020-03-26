#ifndef MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_OPT_Algorithm.hpp" // Base class // OPT/src

namespace moris
{
    namespace opt
    {
        class Algorithm_Sweep : public Algorithm
        {
        public:

            /**
             * Constructor
             */
            Algorithm_Sweep();

            /**
             * Destructor
             */
            ~Algorithm_Sweep();

            /**
             * @brief copy constructor through cloning
             */
            Algorithm* clone() const
            {
                return new Algorithm_Sweep(*this );
            }

            /**
             * @brief MORIS interface for solving of optimization problem using
             *        GCMMA
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            void solve(Problem* aOptProb );

            /**
             * @brief Call for evaluation of objectives and constraints
             *
             * @param[in]  aIter Optimization iteration number
             */
            void func( sint aIter );

            /**
             * @brief Call for evaluation of derivative of objectives and
             *        constraints
             */
            void grad( );

            /**
             *@brief Prints result of the  algorithm based on mStopFlag
             */
            void print_results( );
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_ */

