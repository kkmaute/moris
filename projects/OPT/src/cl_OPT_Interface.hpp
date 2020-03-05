//
// Created by christopherson on 1/16/20.
//

#ifndef MORIS_CL_OPT_INTERFACE_HPP
#define MORIS_CL_OPT_INTERFACE_HPP

#include "cl_Matrix.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Interface
        {
        public:

            /**
             * Constructor
             */
            Interface()
            {
            }

            /**
             * Destructor
             */
            virtual ~Interface()
            {
            }

            /**
             * Initializes the vector of ADV values
             */
            virtual Matrix<DDRMat> initialize_advs() = 0;

            /**
             * Gets the lower bound values for the advs
             *
             * @return vector of lower bounds
             */
            virtual Matrix<DDRMat> get_lower_adv_bounds() = 0;

            /**
             * Gets the upper bound values for the advs
             *
             * @return vector of upper bounds
             */
            virtual Matrix<DDRMat> get_upper_adv_bounds() = 0;

            /**
             * Updates the ADVs to start an optimization step with a new analysis
             *
             * @param aNewADVs the new matrix of ADV values
             */
            virtual void begin_new_analysis(Matrix<DDRMat> aNewADVs) = 0;

            /**
             * Gets the criteria values
             *
             * @return vector of criteria
             */
            virtual Matrix<DDRMat> get_criteria() = 0;

            /**
             * Gets the derivative of the criteria with respect to the advs
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            virtual Matrix<DDRMat> get_dcriteria_dadv() = 0;

        };
    }   // namespace opt
}   // namespace moris

#endif //MORIS_CL_OPT_INTERFACE_HPP
