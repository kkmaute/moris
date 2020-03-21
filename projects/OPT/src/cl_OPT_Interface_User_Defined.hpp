//
// Created by christopherson on 3/16/20.
//

#ifndef MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP
#define MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP

#include "cl_OPT_Interface.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace opt
    {
        class Interface_User_Defined : public Interface
        {
        private:
            Matrix<DDRMat> mADVs;

        public:

            /**
             * Constructor
             */
            Interface_User_Defined(ParameterList aParameterList);

            /**
             * Initializes the vector of ADV values
             */
            Matrix<DDRMat> initialize_advs();

            /**
             * Gets the lower bound values for the advs
             *
             * @return vector of lower bounds
             */
            Matrix<DDRMat> get_lower_adv_bounds();

            /**
             * Gets the upper bound values for the advs
             *
             * @return vector of upper bounds
             */
            Matrix<DDRMat> get_upper_adv_bounds();

            /**
             * Updates the ADVs to start an optimization step with a new analysis
             *
             * @param aNewADVs the new matrix of ADV values
             */
            void begin_new_analysis(Matrix<DDRMat> aNewADVs);

            /**
             * Gets the criteria values
             *
             * @return vector of criteria
             */
            Matrix<DDRMat> get_criteria();

            /**
             * Gets the derivative of the criteria with respect to the advs
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            Matrix<DDRMat> get_dcriteria_dadv();
            
        private:
            // Loaded user-defined functions
            MORIS_DDRMAT0_FUNCTION initialize_advs_user_defined;
            MORIS_DDRMAT0_FUNCTION get_lower_adv_bounds_user_defined;
            MORIS_DDRMAT0_FUNCTION get_upper_adv_bounds_user_defined;
            MORIS_DDSMAT0_FUNCTION get_constraint_types_user_defined;
            MORIS_DDRMAT1_FUNCTION get_criteria_user_defined;
            MORIS_DDRMAT1_FUNCTION get_dcriteria_dadv_user_defined;

        };
    }   // namespace opt
}   // namespace moris

#endif //MORIS_CL_OPT_INTERFACE_USER_DEFINED_HPP
