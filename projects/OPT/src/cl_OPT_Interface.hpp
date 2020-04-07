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
             * Initializes the vectors of ADV values, lower bounds, and upper bounds
             */
            virtual void initialize(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds) = 0;

            /**
             * Gets the criteria values given a new set of ADVs
             *
             * @return vector of criteria
             */
            virtual Matrix<DDRMat> get_criteria(Matrix<DDRMat> aNewADVs) = 0;

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
