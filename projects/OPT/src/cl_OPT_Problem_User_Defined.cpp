#include "cl_OPT_Problem_User_Defined.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Problem_User_Defined::Problem_User_Defined(ParameterList aParameterList,
                                                   std::shared_ptr<Criteria_Interface> aInterface)
                : Problem(aParameterList, aInterface)
        {
            // Load library
            moris::Library_IO tLibrary(aParameterList.get<std::string>("library"));

            // Set user-defined functions
            get_constraint_types_user_defined = tLibrary.load_constraint_types_function("get_constraint_types");
            compute_objectives_user_defined = tLibrary.load_objective_constraint_function("compute_objectives");
            compute_constraints_user_defined = tLibrary.load_objective_constraint_function("compute_constraints");
            compute_dobjective_dadv_user_defined = tLibrary.load_objective_constraint_function("compute_dobjective_dadv");
            compute_dobjective_dcriteria_user_defined = tLibrary.load_objective_constraint_function("compute_dobjective_dcriteria");
            compute_dconstraint_dadv_user_defined = tLibrary.load_objective_constraint_function("compute_dconstraint_dadv");
            compute_dconstraint_dcriteria_user_defined = tLibrary.load_objective_constraint_function("compute_dconstraint_dcriteria");
        }

        //--------------------------------------------------------------------------------------------------------------

        Problem_User_Defined::Problem_User_Defined(ParameterList                       aParameterList,
                                                   std::shared_ptr<Criteria_Interface> aInterface,
                                                   MORIS_CONSTRAINT_TYPES_FUNCTION     aConstraintTypesFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveADVGradientFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aObjectiveCriteriaGradientFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintADVGradientFunction,
                                                   MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aConstraintCriteriaGradientFunction)
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

