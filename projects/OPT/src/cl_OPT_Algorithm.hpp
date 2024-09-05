/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm.hpp
 *
 */

#ifndef MORIS_CL_OPT_ALGORITHM_HPP_
#define MORIS_CL_OPT_ALGORITHM_HPP_

#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris::opt
{
    /**
     * basic tasks for coordinating evaluation of design criteria and their gradients
     * in parallel
     */
    enum class Task
    {
        exit,
        wait,
        compute_criteria_forward_analysis,
        compute_criteria_finite_difference_analysis,
        compute_criteria_gradients_analytically,
        undefined
    };

    enum class SA_Type
    {
        analytical,
        forward,
        backward,
        central
    };

    //-----------------------------------------------------------------------

    class Algorithm
    {
      protected:
        uint mCurrentOptAlgInd;    // stores index of current optimization solver

        moris::uint mRestartIndex         = 0;    // iteration index to be set when restarting
        moris::sint mMaxIterations        = 0;    // maximum number of iterations
        moris::sint mMaxIterationsInitial = 0;    // maximum number of iterations

        std::shared_ptr< moris::opt::Problem > mProblem;    // pointer to problem algorithm operates on

        Matrix< DDSMat > mActive;    // flag for active/inactive constraints

        bool mGradientsHaveBeenComputed = false;    // flag whether gradients have been computed

        Task mRunning = Task::wait;    // task for coordinating parallel runs

        uint mPerturbationSizeIndex = 0;    // index that determines which perturbation size is used

        SA_Type mSAType = SA_Type::analytical;    // sensitivity analysis type; default analytical

        Matrix< DDRMat > mFDObjectiveGradients;      // FD objective gradients
        Matrix< DDRMat > mFDObjectiveConstraints;    // FD constraint gradients

        Matrix< DDRMat > mFiniteDifferenceEpsilons;    // Epsilon for finite differencing

        Matrix< DDUMat > mFiniteDifferenceADVs;    // Indices of ADVs with respect to which sensitivities are
                                                   // computed by finite differencing

        // Function pointer for computing gradients of design criteria either analytically or by finite differencing
        void ( Algorithm::*compute_design_criteria_gradients_by_type )( const Vector< real >& aADVs ) = nullptr;

        // Function pointers for computing the objective and constraint either analytically or by finite differencing
        const Matrix< DDRMat >& ( Algorithm::*get_objective_gradients_by_type )()  = nullptr;
        const Matrix< DDRMat >& ( Algorithm::*get_constraint_gradients_by_type )() = nullptr;

      public:
        /**
         * Constructor
         */
        Algorithm();

        /**
         * Destructor
         */
        virtual ~Algorithm();

        /**
         * @brief Calls the derived optimization algorithm
         *
         * @param[in] aOptProb Object of type Problem containing relevant
         *            data regarding ADVs, the objective and constraints
         */
        virtual uint solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb ) = 0;

        /**
         * Computes design criteria for a given vector of ADVs
         *
         * @param aADVs ADVs, empty if not on proc 0
         */
        void compute_design_criteria( Vector< real >& aADVs );

        /**
         * Computes objective values
         *
         */
        const Matrix< DDRMat >& get_objectives();

        /**
         * Computes constraint values
         *
         */
        const Matrix< DDRMat >& get_constraints();

        /**
         * Computes gradients of design criteria for a given vector of ADVs
         *
         * @param[in] aADVs         - ADVs, empty if not on proc 0
         */
        void compute_design_criteria_gradients( const Vector< real >& aADVs );

        /**
         * Returns gradients of objective(s) based on previous evaluations of sensitivities
         */
        const Matrix< DDRMat >& get_objective_gradients();

        /**
         * Returns gradients of constraints based on previous evaluations of sensitivities
         */
        const Matrix< DDRMat >& get_constraint_gradients();

        /**
         * Sets sensitivity analysis type
         *
         * @param aType - enum for sensitivity analysis type
         */
        void set_sensitivity_analysis_type( SA_Type aType = SA_Type::analytical );

        /**
         *  Initialize finite difference schemes
         */
        void initialize_finite_difference_schemes();

        /**
         * Sets perturbation size for finite differencing
         *
         * @param aEpsilons - vector of values for perturbing the ADVs; if only one value is given
         *                    this perturbation size is applied to all ADVs
         */
        void set_finite_difference_perturbation_size_index( uint aPerturbationSizeIndex );

        /**
         * Sets indices of ADVs with respect to which sensitivities are to be computed by FD
         */
        void set_finite_difference_advs( const Matrix< DDUMat >& aFiniteDifferenceADVs );

        /**
         * Communicates proc 0's running status to other processors so they know when to end.
         */
        void communicate_running_status();

        /**
         * Dummy solve on processors not running optimization algorithm.
         */
        void dummy_solve();

        /**
         * @brief write restart file with advs as well as upper and lower bounds
         */
        void write_advs_to_file( const Vector< real >& aADVs );

        /**
         * @brief set restart index
         */
        void set_restart_index( uint aRestartIndex );

      private:
        void compute_design_criteria_gradients_analytically( const Vector< real >& aADVs );
        void compute_design_criteria_gradients_fd_fwbw( const Vector< real >& aADVs );
        void compute_design_criteria_gradients_fd_central( const Vector< real >& aADVs );

        const Matrix< DDRMat >& get_objective_gradients_analytically();
        const Matrix< DDRMat >& get_objective_gradients_by_fd();

        const Matrix< DDRMat >& get_constraint_gradients_analytically();
        const Matrix< DDRMat >& get_constraint_gradients_by_fd();
    };
    }

#endif /* MORIS_CL_OPT_ALGORITHM_HPP_ */

