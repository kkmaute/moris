/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_Sweep.hpp
 *
 */

#pragma once

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "HDF5_Tools.hpp"

namespace moris::opt
{
    class Algorithm_Sweep : public Algorithm
    {
      public:
        /**
         * Constructor
         */
        Algorithm_Sweep( const Parameter_List& aParameterList, std::shared_ptr< Problem > aProblem = nullptr );

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

        /**
         * Given a direction vector, sets up the evaluation points for all ADVs along this direction based on the number of evaluations and step size
         */
        void setup_evaluation_points_from_direction( const Vector< real >& aDirection, real aStepSize );

        /**
         * Gets the direction vector between the current ADVs and the ADV values from the HDF5 file, and then sets up mNumEvaluationPoints based on this direction vector
         */
        void load_evaluation_points_from_file( const std::string& aFileName, real aStepSize = 1.0 );

        /**
         * Generates a random direction vector, and then sets up mNumEvaluationPoints to be values of the ADVs along the vector
         */
        void load_random_evaluation_points( real aStepSize = 1.0 );
    };
}    // namespace moris::opt