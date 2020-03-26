//
// Created by christopherson on 3/18/20.
//

#include "cl_OPT_Problem_User_Defined.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------
        Problem_User_Defined::Problem_User_Defined(ParameterList aParameterList) : Problem(aParameterList)
        {
            // Load library
            moris::Library_IO tLibrary(aParameterList.get<std::string>("library"));

            // Set user-defined functions
            get_constraint_types_user_defined = tLibrary.load_ddsmat0_function("get_constraint_types");
            compute_objectives_user_defined = tLibrary.load_ddrmat2_function("compute_objectives");
            compute_constraints_user_defined = tLibrary.load_ddrmat2_function("compute_constraints");
            compute_dobjective_dadv_user_defined = tLibrary.load_ddrmat2_function("compute_dobjective_dadv");
            compute_dobjective_dcriteria_user_defined = tLibrary.load_ddrmat2_function("compute_dobjective_dcriteria");
            compute_dconstraint_dadv_user_defined = tLibrary.load_ddrmat2_function("compute_dconstraint_dadv");
            compute_dconstraint_dcriteria_user_defined = tLibrary.load_ddrmat2_function("compute_dconstraint_dcriteria");
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

