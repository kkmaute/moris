//
// Created by christopherson on 2/7/20.
//

#ifndef MORIS_CL_OPT_PROBLEM_ROSENBROCK_HPP
#define MORIS_CL_OPT_PROBLEM_ROSENBROCK_HPP

#include "cl_OPT_Problem.hpp"
#include "cl_OPT_Interface.hpp"
#include "cl_Param_List.hpp" // CON/src

namespace moris
{
    namespace opt
    {
        class Problem_Rosenbrock : public moris::opt::Problem
        {

        public:

            bool mConstrained = false;

            /**
             * Constructor
             *
             * @param aInterface Interface class written for other module (e.g. GEN)
             */
            Problem_Rosenbrock(ParameterList aParameterList) : Problem(aParameterList)
            {}

            /**
             * Gets the constraint types
             *
             * @return vector of integers, 0 = equality constraint, 1 = inequality constraint
             */
            Matrix<DDSMat> get_constraint_types();

        protected:

            /**
             * Gets the objective values
             *
             * @return vector of objectives
             */
            Matrix<DDRMat> compute_objectives();

            /**
             * Gets the constraint values
             *
             * @return vector of constraints
             */
            Matrix<DDRMat> compute_constraints();

            /**
             * Gets the derivative of the objectives with respect to the advs
             *
             * @return matrix d(objective)_i/d(adv)_j
             */
            Matrix<DDRMat> compute_dobjective_dadv();

            /**
             * Gets the derivative of the constraints with respect to the advs
             *
             * @return matrix d(constraints)_i/d(adv)_j
             */
            Matrix<DDRMat> compute_dconstraint_dadv();

            /**
             * Gets the derivative of the objective with respect to the criteria.
             *
             * @return matrix d(objective)_i/d(criteria)_j
             */
            Matrix<DDRMat> compute_dobjective_dcriteria();

            /**
             * Gets the derivative of the constraints with respect to the criteria.
             *
             * @return matrix d(constraint)_i/d(criteria)_j
             */
            Matrix<DDRMat> compute_dconstraint_dcriteria();
        };
    }
}

#endif //MORIS_CL_OPT_PROBLEM_ROSENBROCK_HPP
