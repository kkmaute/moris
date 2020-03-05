//
// Created by christopherson on 2/7/20.
//

#ifndef MORIS_CL_OPT_INTERFACE_ROSENBROCK_HPP
#define MORIS_CL_OPT_INTERFACE_ROSENBROCK_HPP

#include "cl_OPT_Interface.hpp"
#include "cl_Param_List.hpp" // CON/src

namespace moris
{
    namespace opt
    {
        class Interface_Rosenbrock : public Interface
        {
        private:
            Matrix<DDRMat> mADVs;

        public:

            /**
             * Constructor
             */
            Interface_Rosenbrock(ParameterList aParameterList)
            {
            }

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
             * Gets the constraint types
             *
             * @return vector of integers, 0 = equality constraint, 1 = inequality constraint
             */
            Matrix<DDSMat> get_constraint_types();

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

        };
    }   // namespace opt
}   // namespace moris

#endif //MORIS_CL_OPT_INTERFACE_ROSENBROCK_HPP
