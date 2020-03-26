//
// Created by christopherson on 3/16/20.
//

#include "cl_OPT_Interface_User_Defined.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace opt
    {
        //--------------------------------------------------------------------------------------------------------------

        Interface_User_Defined::Interface_User_Defined(ParameterList aParameterList)
        {
            // Load library
            moris::Library_IO tLibrary(aParameterList.get<std::string>("library"));

            // Set user-defined functions
            initialize_advs_user_defined = tLibrary.load_ddrmat0_function("initialize_advs");
            get_lower_adv_bounds_user_defined = tLibrary.load_ddrmat0_function("get_lower_adv_bounds");
            get_upper_adv_bounds_user_defined = tLibrary.load_ddrmat0_function("get_upper_adv_bounds");
            get_criteria_user_defined = tLibrary.load_ddrmat1_function("get_criteria");
            get_dcriteria_dadv_user_defined = tLibrary.load_ddrmat1_function("get_dcriteria_dadv");
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::initialize_advs()
        {
            return this->initialize_advs_user_defined();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::get_lower_adv_bounds()
        {
            return this->get_lower_adv_bounds_user_defined();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::get_upper_adv_bounds()
        {
            return this->get_upper_adv_bounds_user_defined();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interface_User_Defined::begin_new_analysis(Matrix<DDRMat> aNewADVs)
        {
            mADVs = aNewADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::get_criteria()
        {
            return this->get_criteria_user_defined(mADVs);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::get_dcriteria_dadv()
        {
            return this->get_dcriteria_dadv_user_defined(mADVs);
        }

        //--------------------------------------------------------------------------------------------------------------

    }   // namespace opt
}   // namespace moris
