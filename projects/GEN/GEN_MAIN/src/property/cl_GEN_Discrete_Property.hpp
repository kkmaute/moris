//
// Created by christopherson on 5/19/20.
//

#ifndef MORIS_CL_GEN_DISCRETE_PROPERTY_HPP
#define MORIS_CL_GEN_DISCRETE_PROPERTY_HPP

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Discrete.hpp"

namespace moris
{
    namespace ge
    {
        class Discrete_Property : public Property, public Field_Discrete
        {
        public:

            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Discrete_Property(Matrix<DDRMat>& aADVs,
                              Matrix<DDUMat> aPropertyVariableIndices,
                              Matrix<DDUMat> aADVIndices,
                              Matrix<DDRMat> aConstantParameters);

            /**
             * Evaluate the property field based on the given coordinates
             *
             * @param aNodeIndex Node index
             */
            real evaluate_field_value(uint aNodeIndex);

            /**
             * Evaluate the sensitivities of the property value with respect to all internal variables
             *
             * @param aNodeIndex Node Index
             * @param aSensitivities Sensitivity matrix to fill
             */
            void evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities);

        };
    }   // end ge namespace
}   // end moris namespace

#endif //MORIS_CL_GEN_DISCRETE_PROPERTY_HPP
