#ifndef MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Opt_Alg.hpp" // Base class // OPT/src

namespace moris
{
    namespace opt
    {
        class OptAlgSweep : public OptAlg
        {
        private:

        public:

            /**
             * Constructor
             */
            OptAlgSweep();

            /**
             * Destructor
             */
            ~OptAlgSweep();

            /**
             * @brief copy constructor through cloning
             */
            OptAlg*
            clone() const
            {
                return new OptAlgSweep( *this );
            }

            /**
             * @brief MORIS interface for solving of optimization problem using
             *        GCMMA
             *
             * @param[in] aOptProb Object of type OptProb containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            void solve( OptProb & aOptProb );

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
             *@brief Prints result of the GCMMA algorithm based on mStopFlag
             */
            void printresult( );
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGSWEEP_HPP_ */

