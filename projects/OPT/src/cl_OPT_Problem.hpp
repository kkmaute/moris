#ifndef MORIS_OPTIMIZATION_CL_OPTPROB_HPP_
#define MORIS_OPTIMIZATION_CL_OPTPROB_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Problem
        {

        private:
            std::shared_ptr<Criteria_Interface> mInterface;

            Matrix<DDRMat> mUpperBounds;  // upper bounds on ADV vector
            Matrix<DDRMat> mLowerBounds; // lower bounds on ADV vector
            Matrix<DDSMat> mConstraintTypes; // flags for types of constraints
            Matrix<DDRMat> mObjectives; // objectives (always 1)
            Matrix<DDRMat> mConstraints; // constraints
            Matrix<DDRMat> mObjectiveGradient; // full gradient of the objectives with respect to the ADVs
            Matrix<DDRMat> mConstraintGradient; // full gradient of the constraints with respect to the ADVs
            Matrix<DDRMat> mFiniteDifferenceEpsilons; // Epsilon for finite differencing

            std::string mFiniteDifferenceType;
            real mADVNormTolerance = 1E-12;

        protected:
            Matrix<DDRMat> mADVs;    // Abstract Design Variable vector
            Matrix<DDRMat> mCriteria; // vector of criteria values

        public:
            bool mUpdateObjectives = true; // whether or not to compute new objectives when requested
            bool mUpdateConstraints = true; // whether or not to compute new constraints when requested
            bool mUpdateObjectiveGradients = true; // "                 " objective gradients
            bool mUpdateConstraintGradients = true; // "                 " constraint gradients

            /**
             * Constructor
             *
             * @param aParameterList Parameter list for constructing the problem
             * @param aInterface Interface class written for other module
             */
            Problem(ParameterList aParameterList, std::shared_ptr<Criteria_Interface> aInterface);

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
                return this->get_objectives().numel();
            }

            /**
             * Get the number of constraints
             */
            uint get_num_constraints()
            {
                return this->get_constraints().numel();
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
             * Sets the ADV vector. Note that this also indicates to the interface that values need to be recomputed.
             *
             * @param vector of ADV values
             */
            void set_advs(Matrix<DDRMat> aNewADVs);

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
            
            /**
             * Initializes the ADVs and the upper and lower bounds, plus gets constraint types
             */
            void initialize();

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
            Matrix<DDRMat> get_objectives();

            /**
             * Gets the constraint values
             *
             * @return vector of constraints
             */
            Matrix<DDRMat> get_constraints();

            /**
             * Returns the objective gradient, and computes it if not already available.
             *
             * @return full derivative matrix d(objective)_i/d(adv)_j
             */
            Matrix<DDRMat> get_objective_gradients();

            /**
             * Returns the constraint gradient, and computes it if not already available.
             *
             * @return full derivative matrix d(constraint)_i/d(adv)_j
             */
            Matrix<DDRMat> get_constraint_gradients();

            /**
             * Sets the type of finite differencing used for the objectives and an epsilon
             *
             * @param aType std::string defining type; forward, backward, central, or none
             * @param (optional) aEpsilons vector of values for perturbing the ADVs
             */
            void set_finite_differencing(std::string aType, Matrix<DDRMat> aEpsilons);
            void set_finite_differencing(std::string aType);

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

            /**
             * Function pointer which computes the objective gradient in a user-specified way
             */
            void (Problem::*compute_objective_gradient)() = &Problem::compute_objective_gradient_analytical;

            /**
             * Function pointer which computes the constraint gradients in a user-specified way
             */
            void (Problem::*compute_constraint_gradient)() = &Problem::compute_constraint_gradient_analytical;

            /**
             * Assembles the objective gradient from the explicit and implicit gradient terms.
             */
            void compute_objective_gradient_analytical();

            /**
             * Assembles the constraint gradient from the explicit and implicit gradient terms.
             */
            void compute_constraint_gradient_analytical();

            /**
             * Computes the objective gradient using forward or backward finite differencing
             */
            void compute_objective_gradient_fd_bias();

            /**
             * Computes the constraint gradient using forward or backward finite differencing
             */
            void compute_constraint_gradient_fd_bias();

            /**
             * Computes the objective gradient using central finite differencing
             */
            void compute_objective_gradient_fd_central();

            /**
             * Computes the constraint gradient using central finite differencing
             */
            void compute_constraint_gradient_fd_central();

        };
    }
}

#endif /* MORIS_OPTIMIZATION_CL_OPTPROB_HPP_ */
