//
// Created by christopherson on 3/18/20.
//

#ifndef MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP
#define MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP

#include "cl_OPT_Problem.hpp"
#include "cl_OPT_Interface.hpp"
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
             * @param aInterface Interface class written for other module (e.g. GEN)
             */
            Problem_User_Defined(ParameterList aParameterList, std::shared_ptr<Interface> aInterface);

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
            MORIS_DDSMAT0_FUNCTION get_constraint_types_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_objectives_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_constraints_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_dobjective_dadv_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_dobjective_dcriteria_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_dconstraint_dadv_user_defined;
            MORIS_DDRMAT2_FUNCTION compute_dconstraint_dcriteria_user_defined;
        };
    }
}

#endif //MORIS_CL_OPT_PROBLEM_USER_DEFINED_HPP
