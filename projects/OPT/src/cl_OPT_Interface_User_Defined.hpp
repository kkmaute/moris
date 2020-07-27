#ifndef MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP
#define MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP

#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace opt
    {
        class Interface_User_Defined : public Criteria_Interface
        {
        private:
            Matrix<DDRMat> mADVs;
            std::shared_ptr<Library_IO> mLibrary;

        public:

            /**
             * Constructor
             *
             * @param aParameterList Parameter list containing parameters for a user-defined interface.
             */
            Interface_User_Defined(ParameterList aParameterList);

            /**
             * Alternate constructor where the user-defined functions are provided directly. Used in the OPT tests.
             *
             * @param aInitializationFunction Function for initializing ADVs and lower/upper bounds.
             * @param aCriteriaEvaluationFunction Function for evaluating the criteria vector.
             * @param aCriteriaGradientFunction Function for evaluating the gradient of the criteria vector wrt ADVs.
             */
            Interface_User_Defined(MORIS_CRITERIA_INITIALIZE_FUNCTION aInitializationFunction,
                                   MORIS_CRITERIA_FUNCTION aCriteriaEvaluationFunction,
                                   MORIS_CRITERIA_FUNCTION aCriteriaGradientFunction);

            /**
             * Initializes the vectors of ADV values, lower bounds, and upper bounds
             *
             * @param aADVs Initial ADVs to be filled.
             * @param aLowerBounds Lower ADV bounds to be filled.
             * @param aUpperBounds Upper ADV bounds to be filled.
             */
            void initialize(Matrix<DDRMat>& aADVs,
                            Matrix<DDRMat>& aLowerBounds,
                            Matrix<DDRMat>& aUpperBounds);

            /**
             * Gets the criteria values.
             *
             * @return vector of criteria
             */
            Matrix<DDRMat> perform(Matrix<DDRMat> aNewADVs);

            /**
             * Gets the derivative of the criteria with respect to the advs.
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            Matrix<DDRMat> compute_dcriteria_dadv();
            
        private:
            // Loaded user-defined functions
            MORIS_CRITERIA_INITIALIZE_FUNCTION initialize_user_defined;
            MORIS_CRITERIA_FUNCTION get_criteria_user_defined;
            MORIS_CRITERIA_FUNCTION compute_dcriteria_dadv_user_defined;

        };
    }
}

#endif //MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP
