#ifndef MORIS_OPTIMIZATION_CL_OPTALGLBFGS_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGLBFGS_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_OPT_Algorithm.hpp" // Base class // OPT/src

namespace moris
{
    namespace opt
    {
        class Algorithm_LBFGS : public Algorithm
        {
        private:

            uint mOptIter; // optimization iteration counter
            uint mResFlag; // Flag from L-BFGS describing result of optimization algorithm

        public:

            /**
             * Constructor
             */
            Algorithm_LBFGS();

            /**
             * Destructor
             */
            ~Algorithm_LBFGS();

            /**
             * @brief copy constructor through cloning
             */
            Algorithm*
            clone() const
            {
                return new Algorithm_LBFGS(*this );
            }

            /**
             * @brief MORIS interface for solving of optimization problem using
             *        GCMMA
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            void solve(std::shared_ptr<Problem> aOptProb );

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
