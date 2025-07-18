/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_Sweep.hpp
 *
 */

#ifndef MORIS_CL_OPT_ALGORITHM_SWEEP_HPP_
#define MORIS_CL_OPT_ALGORITHM_SWEEP_HPP_

#include "cl_OPT_Algorithm.hpp"
#include "HDF5_Tools.hpp"

// forward declarations
namespace moris
{
    class Parameter_List;

    namespace opt
    {
        class Problem;
    }
}

namespace moris::opt
{
    class Algorithm_Sweep : public Algorithm
    {
      public:
        /**
         * Constructor
         */
        Algorithm_Sweep( const Parameter_List& aParameterList );

        /**
         * Destructor
         */
        ~Algorithm_Sweep() override;

        /**
         * @brief MORIS interface for solving of optimization problem using
         *        GCMMA
         *
         * @param[in] aCurrentOptAlgInd index of optimization algorithm
         * @param[in] aOptProb Object of type Problem containing relevant
         *            data regarding ADVs, the objective and constraints
         */
        uint solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb ) override;

      private:
        bool  mIncludeBounds;                  // whether or not to include upper/lower bounds in the sweep
        bool  mEvaluateObjectives;             // whether or not to compute new objectives when requested
        bool  mEvaluateConstraints;            // whether or not to compute new constraints when requested
        bool  mEvaluateObjectiveGradients;     // "                 " objective gradients
        bool  mEvaluateConstraintGradients;    // "                 " constraint gradients
        bool  mSave;                           // If saving the results of the sweep to an hdf5 file
        bool  mPrint;                          // If printing the results of the sweep to the screen
        hid_t mFileID;                         // Fild id for hdf5 file

        std::string mFiniteDifferenceType;    // Finite difference type

        Matrix< DDUMat > mNumEvaluations;      // Number of evaluations per ADV
        Matrix< DDRMat > mEvaluationPoints;    // Final evaluation points

        /**
         * Runs sweep algorithm on processor 0
         */
        void sweep_solve();

        /**
         * Sets up evaluation points
         */
        void set_up_evaluation_points();

        /**
         * Outputs the optimization problem at the current ADVs (objective and constraints)
         *
         * @param aEvaluationName The name to be printed/saved after the optimization variable type
         */
        void output_objectives_constraints( const std::string& aEvaluationName );

        /**
         * Evaluates the objective gradients at the current ADVs and outputs them
         *
         * @param aEvaluationName The name to be printed/saved after the optimization variable type
         */
        Matrix< DDRMat > evaluate_objective_gradients( const std::string& aEvaluationName );

        /**
         * Evaluates the constraint gradients at the current ADVs and outputs them
         *
         * @param aEvaluationName The name to be printed/saved after the optimization variable type
         */
        Matrix< DDRMat > evaluate_constraint_gradients( const std::string& aEvaluationName );

        /**
         * Outputs the given optimization variables based on printing/saving options
         *
         * @param aVariables Matrix of optimization variables
         * @param aFullEvaluationName Full name to be output to the screen/hdf5
         */
        void output_variables( const Matrix< DDRMat >& aVariables, const std::string& aFullEvaluationName );
    };
    }

#endif /* MORIS_CL_OPT_ALGORITHM_SWEEP_HPP_ */

