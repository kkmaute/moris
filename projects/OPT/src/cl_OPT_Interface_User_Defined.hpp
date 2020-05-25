//
// Created by christopherson on 3/16/20.
//

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
            moris::Library_IO mLibrary;

        public:

            /**
             * Constructor
             */
            Interface_User_Defined(ParameterList aParameterList);

            /**
             * Initializes the vectors of ADV values, lower bounds, and upper bounds
             */
            void initialize(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds);

            /**
             * Gets the criteria values
             *
             * @return vector of criteria
             */
            Matrix<DDRMat> perform(Matrix<DDRMat> aNewADVs);

            /**
             * Gets the derivative of the criteria with respect to the advs
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
    }   // namespace opt
}   // namespace moris

#endif //MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP
