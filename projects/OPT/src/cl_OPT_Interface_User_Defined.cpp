#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace opt
    {
        //--------------------------------------------------------------------------------------------------------------

        Interface_User_Defined::Interface_User_Defined(
                ParameterList aParameterList)
        {
            // set library with user defined functions
            mLibrary = std::make_shared<Library_IO>(aParameterList.get<std::string>("library"));

            // Set user-defined functions
            initialize_user_defined             = mLibrary->load_function<Criteria_Initialize_Function>("initialize");
            get_criteria_user_defined           = mLibrary->load_function<Criteria_Function>("get_criteria");
            compute_dcriteria_dadv_user_defined = mLibrary->load_function<Criteria_Function>("get_dcriteria_dadv");
        }

        //--------------------------------------------------------------------------------------------------------------

        Interface_User_Defined::Interface_User_Defined(
                Criteria_Initialize_Function aInitializationFunction,
                Criteria_Function aCriteriaEvaluationFunction,
                Criteria_Function aCriteriaGradientFunction)
        {
            // Set user-defined functions
            initialize_user_defined             = aInitializationFunction;
            get_criteria_user_defined           = aCriteriaEvaluationFunction;
            compute_dcriteria_dadv_user_defined = aCriteriaGradientFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interface_User_Defined::initialize(
                Matrix<DDRMat>& aADVs,
                Matrix<DDRMat>& aLowerBounds,
                Matrix<DDRMat>& aUpperBounds)
        {
            initialize_user_defined(aADVs, aLowerBounds, aUpperBounds);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_User_Defined::perform( Matrix< DDRMat >& aNewADVs )
        {
            mADVs = aNewADVs;

            return this->get_criteria_user_defined(mADVs);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_User_Defined::compute_dcriteria_dadv()
        {
            return this->compute_dcriteria_dadv_user_defined(mADVs);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
