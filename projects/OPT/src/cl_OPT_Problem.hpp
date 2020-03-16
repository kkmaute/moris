#ifndef MORIS_OPTIMIZATION_CL_OPTPROB_HPP_
#define MORIS_OPTIMIZATION_CL_OPTPROB_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_OPT_Interface.hpp"

namespace moris
{
    namespace opt
    {
        class Problem
        {

        private:
            Interface* mInterface;

            Matrix<DDRMat> mUpperBounds;  // upper bounds on ADV vector
            Matrix<DDRMat> mLowerBounds; // lower bounds on ADV vector
            Matrix<DDSMat> mConstraintTypes; // flags for types of constraints
            Matrix<DDRMat> mObjectives; // objective values (typically just 1, supporting multi-objective for future implementation)
            Matrix<DDRMat> mConstraints; // constraint values
            Matrix<DDRMat> mObjectiveGradient; // Full gradient of the objectives with respect to the ADVs
            Matrix<DDRMat> mConstraintGradient; // Full gradient of the constraints with respect to the ADVs

        protected:
            Matrix<DDRMat> mADVs;    // Abstract Design Variable vector
            Matrix<DDRMat> mCriteria; // vector of criteria values

        public:
            bool mUpdateObjectives = true;
            bool mUpdateConstraints = true;
            bool mUpdateObjectiveGradient = true;
            bool mUpdateConstraintGradient = true;

            /**
             * Constructor
             *
             * @param aInterface Interface class written for other module (e.g. GEN)
             */
            Problem(Interface* aInterface);

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
                return mObjectives.numel();
            }

            /**
             * Get the number of constraints
             */
            uint get_num_constraints()
            {
                return mConstraints.numel();
            }

            /**
             * Get number of equality constraints
             */
            uint get_num_equality_constraints();

            /**
             * Get the adv vector
             */
            Matrix<DDRMat> get_advs()
            {
                return mADVs;
            }

            /**
             * Sets the ADV vector. Note that this also indicates to the interface that values need to be recomputed.
             *
             * @param vector of ADV values
             */
            void set_advs(Matrix<DDRMat> aNewADVs);

            /**
             * Get the adv upper bounds
             */
            Matrix<DDRMat> get_upper_bounds()
            {
                return mUpperBounds;
            }

            /**
             * Get the adv lower bounds
             */
            Matrix<DDRMat> get_lower_bounds()
            {
                return mLowerBounds;
            }

            /**
             * Initializes the ADVs and the upper and lower bounds, plus gets constraint types
             */
            void initialize();

            /**
             * Checks that the constraints are given in the correct order (equality constraints first)
             */
            bool check_constraint_order();

            /**
             * Returns the objective gradient, and computes it if not already available.
             *
             * @return full derivative matrix d(objective)_i/d(adv)_j
             */
            Matrix<DDRMat> get_objective_gradient(); // TODO rename these

            /**
             * Returns the constraint gradient, and computes it if not already available.
             *
             * @return full derivative matrix d(constraint)_i/d(adv)_j
             */
            Matrix<DDRMat> get_constraint_gradient();

            /**
             * Modifies the optimization solution. Not currently implemented.
             */
            void scale_solution();

            /**
             * Update the optimization problem. Not currently implemented.
             */
            void update_problem();

            /**
             * Override the ADV values given by the interface. Optional.
             */
            virtual void override_advs()
            {
                // ADVs remain unchanged in base implementation
            };

            /**
             * Gets the constraint types
             *
             * @return vector of integers, 0 = equality constraint, 1 = inequality constraint
             */
            virtual Matrix<DDSMat> get_constraint_types() = 0;

            /**
             * Gets the objective values
             *
             * @return vector of objectives
             */
            virtual Matrix<DDRMat> get_objectives() = 0;

            /**
             * Gets the constraint values
             *
             * @return vector of constraints
             */
            virtual Matrix<DDRMat> get_constraints() = 0;

            /**
             * Gets the derivative of the objectives with respect to the advs
             *
             * @return matrix d(objective)_i/d(adv)_j
             */
            virtual Matrix<DDRMat> get_dobjective_dadv() = 0;

            /**
             * Gets the derivative of the constraints with respect to the advs
             *
             * @return matrix d(constraints)_i/d(adv)_j
             */
            virtual Matrix<DDRMat> get_dconstraint_dadv() = 0;

            /**
             * Gets the derivative of the objective with respect to the criteria.
             *
             * @return matrix d(objective)_i/d(criteria)_j
             */
            virtual Matrix<DDRMat> get_dobjective_dcriteria() = 0;

            /**
             * Gets the derivative of the constraints with respect to the criteria.
             *
             * @return matrix d(constraint)_i/d(criteria)_j
             */
            virtual Matrix<DDRMat> get_dconstraint_dcriteria() = 0;

        private:
            /**
             * Assembles the objective gradient from the explicit and implicit gradient terms.
             */
            void compute_objective_gradient();

            /**
             * Assembles the constraint gradient from the explicit and implicit gradient terms.
             */
            void compute_constraint_gradient();
        };
    }
}

#endif /* MORIS_OPTIMIZATION_CL_OPTPROB_HPP_ */
