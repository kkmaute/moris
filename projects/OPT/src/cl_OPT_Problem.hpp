/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem.hpp
 *
 */

#ifndef MORIS_CL_OPT_PROBLEM_HPP_
#define MORIS_CL_OPT_PROBLEM_HPP_

#include "core.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace opt
    {
        class Problem
        {

            private:

                std::shared_ptr<Criteria_Interface> mInterface;

                Matrix<DDRMat> mUpperBounds;               // upper bounds on ADV vector
                Matrix<DDRMat> mLowerBounds;               // lower bounds on ADV vector
                Matrix<DDSMat> mConstraintTypes;           // flags for types of constraints
                Matrix<DDRMat> mObjectives = {{0.0}};      // objectives (always 1)
                Matrix<DDRMat> mConstraints;               // constraints
                Matrix<DDRMat> mObjectiveGradient;         // full gradient of the objectives with respect to the ADVs
                Matrix<DDRMat> mConstraintGradient;        // full gradient of the constraints with respect to the ADVs

                std::string mRestartFile;                  // Restart file

                real mADVNormTolerance = 1e-12;            // Tolerance for determining whether ADV vector has changed

                uint mNumObjectives;                       // Number of objectives (should be one)
                uint mNumConstraints;                      // Number of constraints

            protected:

                Matrix<DDRMat> mADVs;       // Abstract Design Variable vector
                Matrix<DDRMat> mCriteria;   // vector of criteria values
                Matrix<IdMat > mIjklIds;       // ijklIds for restart

            public:

                /**
                 * Constructor
                 *
                 * @param aParameterList Parameter list for constructing the problem
                 * @param aInterface Interface class written for other module
                 */
                Problem(
                      Parameter_List& aParameterList,
                        std::shared_ptr<Criteria_Interface> & aInterface);

                /**
                 * Destructor
                 */
                ~Problem();

                /**
                 * Get the number of advs
                 */
                uint get_num_advs()
                {
                    return mADVs.numel();
                }

                /**
                 * Get the number of objectives (should almost always be 1)
                 */
                uint get_num_objectives()
                {
                    return mNumObjectives;
                }

                /**
                 * Get the number of constraints
                 */
                uint get_num_constraints()
                {
                    return mNumConstraints;
                }

                /**
                 * Get number of equality constraints
                 */
                uint get_num_equality_constraints();

                /**
                 * Get reference to the adv vector
                 */
                Matrix<DDRMat> & get_advs()
                {
                    return mADVs;
                }

                /**
                 * Computes design criteria for given ADV vector.
                 *
                 * @param vector of ADV values
                 */
                void compute_design_criteria( Matrix< DDRMat >& aNewADVs );

                /**
                 * Computes analytically the criteria gradients for current ADV vector of problem
                 */
                void compute_design_criteria_gradients(const Matrix<DDRMat> & aNewADVs);

                /**
                 * Get reference to the adv upper bounds
                 *
                 * @return vector of upper bounds
                 */
                Matrix<DDRMat> & get_upper_bounds()
                {
                    return mUpperBounds;
                }

                /**
                 * Get reference to the adv lower bounds
                 *
                 * @return vector of lower bounds
                 */
                Matrix<DDRMat> & get_lower_bounds()
                {
                    return mLowerBounds;
                }

                Matrix<IdMat> & get_ijklIDs()
                {
                    return mIjklIds;
                }

                /**
                 * Initializes the ADVs and the upper and lower bounds, plus gets constraint types
                 */
                void initialize();

                /**
                 * get restart optimization flag
                 */
                bool restart_optimization()
                {
                    return mInterface->get_restart_optimization();
                }

                /**
                 * Checks that the constraints are given in the correct order (equality constraints first)
                 *
                 * @return bool, if order is correct
                 */
                bool check_constraint_order();

                /**
                 * Gets the objective values
                 *
                 * @return vector of objectives
                 */
                const Matrix<DDRMat> & get_objectives();

                /**
                 * Gets the constraint values
                 *
                 * @return vector of constraints
                 */
                const Matrix<DDRMat> & get_constraints();

                /**
                 * Sets objective and  constraint values
                 *
                 * @param[in] aObjectives - vector of objectives
                 * @param[in] aConstraints - vector of constraints
                 */
                void set_objectives_and_constraints(
                        const Matrix<DDRMat>& aObjectives,
                        const Matrix<DDRMat>& aConstraints);

                /**
                 * Returns the objective gradient, and computes it if not already available.
                 *
                 * @return full derivative matrix d(objective)_i/d(adv)_j
                 */
                const Matrix<DDRMat> & get_objective_gradients();

                /**
                 * Returns the constraint gradient, and computes it if not already available.
                 *
                 * @return full derivative matrix d(constraint)_i/d(adv)_j
                 */
                const Matrix<DDRMat> & get_constraint_gradients();

                /**
                 * Modifies the optimization solution. Not currently implemented.
                 */
                void scale_solution();

                /**
                 * Update the optimization problem. Not currently implemented.
                 */
                void update_problem();

                /**
                 * Gets the constraint types
                 *
                 * @return vector of integers, 0 = equality constraint, 1 = inequality constraint
                 */
                virtual Matrix<DDSMat> get_constraint_types() = 0;

            private:
                /**
                 * Read restart file
                 */
                void read_restart_file();

                /**
                 * Override the ADV values given by the interface. Optional.
                 */
                virtual void override_advs()
                {
                    // ADVs remain unchanged in base implementation
                };

                /**
                 * Calculates the objective value
                 *
                 * @return vector of objectives
                 */
                virtual Matrix<DDRMat> compute_objectives() = 0;

                /**
                 * Calculates the constraint values
                 *
                 * @return vector of constraints
                 */
                virtual Matrix<DDRMat> compute_constraints() = 0;

                /**
                 * Calculates the derivative of the objectives with respect to the advs
                 *
                 * @return matrix d(objective)_i/d(adv)_j
                 */
                virtual Matrix<DDRMat> compute_dobjective_dadv() = 0;

                /**
                 * Calculates the derivative of the objective with respect to the criteria.
                 *
                 * @return matrix d(objective)_i/d(criteria)_j
                 */
                virtual Matrix<DDRMat> compute_dobjective_dcriteria() = 0;

                /**
                 * Calculates the derivative of the constraints with respect to the advs
                 *
                 * @return matrix d(constraints)_i/d(adv)_j
                 */
                virtual Matrix<DDRMat> compute_dconstraint_dadv() = 0;

                /**
                 * Calculates the derivative of the constraints with respect to the criteria.
                 *
                 * @return matrix d(constraint)_i/d(criteria)_j
                 */
                virtual Matrix<DDRMat> compute_dconstraint_dcriteria() = 0;
        };
    }
}

#endif /* MORIS_CL_OPT_PROBLEM_HPP_ */

