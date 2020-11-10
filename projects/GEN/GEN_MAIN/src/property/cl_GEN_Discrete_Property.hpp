#ifndef MORIS_CL_GEN_DISCRETE_PROPERTY_HPP
#define MORIS_CL_GEN_DISCRETE_PROPERTY_HPP

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Discrete_Integration.hpp"

namespace moris
{
    namespace ge
    {
        class Discrete_Property : public Property, public Field_Discrete_Integration
        {
        public:

            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementMeshIndices Indices of meshes to perform refinement on
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Discrete_Property(
                    Matrix<DDRMat>& aADVs,
                    Matrix<DDUMat>  aPropertyVariableIndices,
                    Matrix<DDUMat>  aADVIndices,
                    Matrix<DDRMat>  aConstantParameters,
                    std::string     aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aRefMeshIndex = {{}},
                    sint            aRefinementFunctionIndex = -1,
                    sint            aBSplineMeshIndex = -2,
                    real            aBSplineLowerBound = -1.0,
                    real            aBSplineUpperBound = 1.0);

            /**
             * Constructor
             *
             * @param aOwnedADVs Distributed owned ADVs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementMeshIndices Indices of meshes to perform refinement on
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Discrete_Property(
                    sol::Dist_Vector* aOwnedADVs,
                    Matrix<DDUMat>    aPropertyVariableIndices,
                    Matrix<DDUMat>    aADVIndices,
                    Matrix<DDRMat>    aConstantParameters,
                    std::string       aName = "",
                    Matrix<DDSMat>    aNumRefinements = {{}},
                    Matrix<DDSMat>    aRefMeshIndex = {{}},
                    sint              aRefinementFunctionIndex = -1,
                    sint              aBSplineMeshIndex = -2,
                    real              aBSplineLowerBound = -1.0,
                    real              aBSplineUpperBound = 1.0);

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Property value
             */
            real get_field_value(uint aNodeIndex);

            /**
             * Given a node index, evaluates the sensitivity of the property field with respect to all of the
             * property variables.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(uint aNodeIndex);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

        };
    }
}

#endif //MORIS_CL_GEN_DISCRETE_PROPERTY_HPP
