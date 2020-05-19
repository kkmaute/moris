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

        Interface_User_Defined::Interface_User_Defined(ParameterList aParameterList) : mLibrary(aParameterList.get<std::string>("library"))
        {
            // Set user-defined functions
            initialize_user_defined = mLibrary.load_criteria_initialize_function("initialize");
            get_criteria_user_defined = mLibrary.load_criteria_function("get_criteria");
            get_dcriteria_dadv_user_defined = mLibrary.load_criteria_function("get_dcriteria_dadv");
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interface_User_Defined::initialize(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            initialize_user_defined(aADVs, aLowerBounds, aUpperBounds);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::perform(Matrix<DDRMat> aNewADVs)
        {
            mADVs = aNewADVs;
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
