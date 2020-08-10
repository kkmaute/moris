#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Child_Node::Child_Node(Matrix<DDUMat>             aParentNodeIndices,
                               Cell<Matrix<DDRMat>>       aParentNodeCoordinates,
                               const xtk::Basis_Function& aBasisFunction,
                               Matrix<DDRMat>             aLocalCoordinates)
                : mParentNodeIndices(aParentNodeIndices),
                  mParentNodeCoordinates(aParentNodeCoordinates),
                  mLocalCoordinates(aLocalCoordinates)
        {
            // Shift local coordinates
            real tEpsilon = 1E-12;
            for (uint tDimension = 0; tDimension < mLocalCoordinates.length(); tDimension++)
            {
                mLocalCoordinates(tDimension) = std::min(mLocalCoordinates(tDimension), 1.0 - tEpsilon);
                mLocalCoordinates(tDimension) = std::max(mLocalCoordinates(tDimension), tEpsilon - 1.0);
            }

            // Evaluate basis function
            aBasisFunction.evaluate_basis_function(mLocalCoordinates, mBasisValues);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Child_Node::get_local_coordinates()
        {
            return mLocalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Child_Node::interpolate_geometry_field_value(Geometry* aGeometry)
        {
            // Get field values from parents
            Matrix<DDRMat> tGeometryFieldValues(mParentNodeIndices.length(), 1);
            for (uint tParentNode = 0; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                tGeometryFieldValues(tParentNode) = aGeometry->evaluate_field_value(mParentNodeIndices(tParentNode), mParentNodeCoordinates(tParentNode));
            }

            // Return interpolated value
            return Matrix<DDRMat>(mBasisValues * tGeometryFieldValues)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Child_Node::interpolate_geometry_sensitivity(Geometry* aGeometry, Matrix<DDRMat>& aSensitivities)
        {
            // Initialize using first parent
            aGeometry->evaluate_sensitivity(mParentNodeIndices(0), mParentNodeCoordinates(0), aSensitivities);
            aSensitivities = aSensitivities * mBasisValues(0);

            // Get sensitivity values from other parents
            for (uint tParentNode = 1; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                Matrix<DDRMat> tParentSensitivity(0, 0);
                aGeometry->evaluate_sensitivity(mParentNodeIndices(tParentNode), mParentNodeCoordinates(tParentNode), tParentSensitivity);
                aSensitivities = aSensitivities + mBasisValues(tParentNode) * tParentSensitivity;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
