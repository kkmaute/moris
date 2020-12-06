#ifndef MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP
#define MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP

#include "cl_OPT_Problem.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace opt
    {
        class Problem_User_Defined : public moris::opt::Problem
        {

        public:

            /**
             * Constructor
             *
             * @param aParameterList parameter list for this problem specifying the needed library
             * @param aInterface Interface class written for other module
             */
            Problem_User_Defined(
                    ParameterList                       & aParameterList,
                    std::shared_ptr<Criteria_Interface> & aInterface);

            /**
             * Alternate constructor where the user-defined functions are provided directly. Used for OPT tests.
             *
             * @param aParameterList Parameter list for the base Problem class
             * @param aInterface Interface class written for other module
             * @param aConstraintTypesFunction Function for getting constraint types
             * @param aObjectiveFunction Objective function
             * @param aConstraintFunction Constraint function
             * @param aObjectiveADVGradientFunction
             * @param aObjectiveCriteriaGradientFunction
             * @param aConstraintADVGradientFunction
             * @param aConstraintCriteriaGradientFunction
             */
            Problem_User_Defined(
                    ParameterList                       aParameterList,
                    std::shared_ptr<Criteria_Interface> aInterface,
                    MORIS_CONSTRAINT_TYPES_FUNCTION     aConstraintTypesFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveADVGradientFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveCriteriaGradientFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintADVGradientFunction,
                    MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintCriteriaGradientFunction);

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

        private:
            MORIS_CONSTRAINT_TYPES_FUNCTION get_constraint_types_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_objectives_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_constraints_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_dobjective_dadv_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_dobjective_dcriteria_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_dconstraint_dadv_user_defined;
            MORIS_OBJECTIVE_CONSTRAINT_FUNCTION compute_dconstraint_dcriteria_user_defined;
        };
    }
}

#endif // MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP
