/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem_User_Defined.cpp
 *
 */

#include "cl_OPT_Problem_User_Defined.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Problem_User_Defined::Problem_User_Defined(
                ParameterList                       & aParameterList,
                std::shared_ptr<Criteria_Interface> & aInterface)
        : Problem(aParameterList, aInterface)
        {
            // Load library
            moris::Library_IO tLibrary(aParameterList.get<std::string>("library"));

            // Set user-defined functions
            get_constraint_types_user_defined          = tLibrary.load_function<Constraint_Types_Function>    ("get_constraint_types");
            compute_objectives_user_defined            = tLibrary.load_function<Objective_Constraint_Function>("compute_objectives");
            compute_constraints_user_defined           = tLibrary.load_function<Objective_Constraint_Function>("compute_constraints");
            compute_dobjective_dadv_user_defined       = tLibrary.load_function<Objective_Constraint_Function>("compute_dobjective_dadv");
            compute_dobjective_dcriteria_user_defined  = tLibrary.load_function<Objective_Constraint_Function>("compute_dobjective_dcriteria");
            compute_dconstraint_dadv_user_defined      = tLibrary.load_function<Objective_Constraint_Function>("compute_dconstraint_dadv");
            compute_dconstraint_dcriteria_user_defined = tLibrary.load_function<Objective_Constraint_Function>("compute_dconstraint_dcriteria");
        }

        //--------------------------------------------------------------------------------------------------------------

        Problem_User_Defined::Problem_User_Defined(
                ParameterList                       aParameterList,
                std::shared_ptr<Criteria_Interface> aInterface,
                Constraint_Types_Function     aConstraintTypesFunction,
                Objective_Constraint_Function aObjectiveFunction,
                Objective_Constraint_Function aConstraintFunction,
                Objective_Constraint_Function aObjectiveADVGradientFunction,
                Objective_Constraint_Function aObjectiveCriteriaGradientFunction,
                Objective_Constraint_Function aConstraintADVGradientFunction,
                Objective_Constraint_Function aConstraintCriteriaGradientFunction)
        : Problem(aParameterList, aInterface),
          get_constraint_types_user_defined(aConstraintTypesFunction),
          compute_objectives_user_defined(aObjectiveFunction),
          compute_constraints_user_defined(aConstraintFunction),
          compute_dobjective_dadv_user_defined(aObjectiveADVGradientFunction),
          compute_dobjective_dcriteria_user_defined(aObjectiveCriteriaGradientFunction),
          compute_dconstraint_dadv_user_defined(aConstraintADVGradientFunction),
          compute_dconstraint_dcriteria_user_defined(aConstraintCriteriaGradientFunction)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Problem_User_Defined::get_constraint_types()
        {
            return this->get_constraint_types_user_defined();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_objectives()
        {
            return this->compute_objectives_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_constraints()
        {
            return this->compute_constraints_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_dobjective_dadv()
        {
            return this->compute_dobjective_dadv_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_dobjective_dcriteria()
        {
            return this->compute_dobjective_dcriteria_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_dconstraint_dadv()
        {
            return this->compute_dconstraint_dadv_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_User_Defined::compute_dconstraint_dcriteria()
        {
            return this->compute_dconstraint_dcriteria_user_defined(mADVs, mCriteria);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

