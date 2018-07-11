#ifndef MORIS_OPTIMIZATION_CL_OPTALGLBFGS_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGLBFGS_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Opt_Alg.hpp" // Base class // OPT/src

namespace moris
{
    namespace opt
    {
        class OptAlgLBFGS : public OptAlg
        {
        private:

            uint mOptIter; // optimization iteration counter
            uint mResFlag; // Flag from L-BFGS describing result of optimization algorithm

        public:

            /**
             * Constructor
             */
            OptAlgLBFGS();

            /**
             * Destructor
             */
            ~OptAlgLBFGS();

            /**
             * @brief copy constructor through cloning
             */
            OptAlg*
            clone() const
            {
                return new OptAlgLBFGS( *this );
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
             * @brief MORIS-GCMMA interface for evaluation of objectives and
             *        constraints
             *
             * @param[in]  aIter Optimization iteration number
             * @param[in]  aAdv Pointer to the abstract design variables
             *
             * @param[out] aObjval Value of the objective function
             */
            void func(
                    int     aIter,
                    double* aAdv,
                    double& aObjval );

            /**
             * @brief MORIS-GCMMA interface for evaluation of derivative of
             *        objectives and constraints
             *
             * @param[in]  aAdv Pointer to the abstract design variables
             *
             * @param[out] aD_Obj Pointer to the derivative of the objective
             *             function
             */
            void grad(
                    double*  aAdv,
                    double*  aD_Obj );

            /**
             *@brief Prints result of the L-BFGS algorithm based on mStopFlag
             */
            void printresult();
        };
    } // namespace opt
}     // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGLBFGS_HPP_ */
