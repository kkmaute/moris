#ifndef MORIS_CL_OPT_ALGORITHM_LBFGS_HPP_
#define MORIS_CL_OPT_ALGORITHM_LBFGS_HPP_

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"

namespace moris
{
    namespace opt
    {
        class Algorithm_LBFGS : public Algorithm
        {
        private:

            uint mOptIter = 0;         // optimization iteration counter
            uint mResFlag;             // Flag from L-BFGS describing result of optimization algorithm
            sint mMaxIt;               // Maximum number of optimization iterations
            int mNumCorrections;       // Number of limited memory corrections
            double mNormDrop;          // Convergence criteria based on norm
            double mGradTolerance;      // Convergence criteria based on projected gradients

        public:

            /**
             * Constructor
             */
            Algorithm_LBFGS(ParameterList aParameterList);

            /**
             * Destructor
             */
            ~Algorithm_LBFGS();

            /**
             * @brief MORIS interface for solving of optimization problem using
             *        GCMMA
             *
             * @param[in] aCurrentOptAlgInd index of optimization algorithm
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            uint solve(
                    uint aCurrentOptAlgInd,
                    std::shared_ptr<Problem> aOptProb );

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
    }
}

#endif /* MORIS_CL_OPT_ALGORITHM_LBFGS_HPP_ */
