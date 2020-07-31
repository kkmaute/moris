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
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Discrete_Property(Matrix<DDRMat>& aADVs,
                              Matrix<DDUMat> aPropertyVariableIndices,
                              Matrix<DDUMat> aADVIndices,
                              Matrix<DDRMat> aConstantParameters,
                              sint            aNumRefinements = 0,
                              sint            aRefinementFunctionIndex = -1,
                              sint            aBSplineMeshIndex = -1,
                              real            aBSplineLowerBound = -1.0,
                              real            aBSplineUpperBound = 1.0);

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
    }
}

#endif //MORIS_CL_GEN_DISCRETE_PROPERTY_HPP
