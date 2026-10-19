/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_ADAM.hpp
 *
 */

#pragma once

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"

namespace moris::opt
{
    class Algorithm_ADAM : public Algorithm
    {
      private:
        uint mOptIter = 0;    // optimization iteration counter

        uint mMaxIt;    // Maximum number of optimization iterations

        real mNormDrop;         // Convergence criteria based on norm drop of the optimization variables
        real mGradTolerance;    // Convergence criteria based on norm of the gradient

        real mMomentumBeta1;    // Decay rate for moving average of gradient
        real mMomentumBeta2;    // Decary rate for moving average of squared gradient
        real mLearningRate;     // Learning rate for ADAM update
        real mEpsilon;          // Small value to prevent division by zero

      public:
        /**
         * Constructor
         */
        Algorithm_ADAM( const Parameter_List& aParameterList );

        /**
         * Destructor
         */
        ~Algorithm_ADAM() override;

        /**
         * @brief MORIS interface for solving of optimization problem using
         *        ADAM
         *
         * @param[in] aCurrentOptAlgInd index of optimization algorithm
         * @param[in] aOptProb Object of type Problem containing relevant
         *            data regarding ADVs, the objective and constraints
         */
        uint solve(
                uint                       aCurrentOptAlgInd,
                std::shared_ptr< Problem > aOptProb ) override;

        /**
         * @brief Core ADAM algorithm to solve optimization problem
         */
        void ADAM_solve();

        /**
         * @brief Prints result of the ADAM algorithm based on mStopFlag
         */
        void printresult();

        //--------------------------------------------------------------------------------------------------------------
    };
}    // namespace moris::opt