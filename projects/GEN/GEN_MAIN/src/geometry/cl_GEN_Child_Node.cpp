#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"

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
                  mParentNodeCoordinates(aParentNodeCoordinates)
        {
            aBasisFunction.evaluate_basis_function(aLocalCoordinates, mBasisValues);
        }

        //--------------------------------------------------------------------------------------------------------------

        Child_Node::Child_Node(Matrix<DDUMat>             aParentNodeIndices,
                               Cell<Matrix<DDRMat>>       aParentNodeCoordinates,
                               const xtk::Basis_Function& aBasisFunction,
                               Matrix<DDRMat>             aNodeFieldValues,
                               real                       aIsocontourThreshold)
                : Child_Node(aParentNodeIndices,
                             aParentNodeCoordinates,
                             aBasisFunction,
                             Interpolation::linear_interpolation_value(aNodeFieldValues,
                                                                       aIsocontourThreshold))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Child_Node::interpolate_geometry_field_value(Geometry* aGeometry)
        {
            // Get field values from parent
            Matrix<DDRMat> tGeometryFieldValues(mParentNodeIndices.length(), 1);
            for (uint tParentNode = 0; tParentNode < mParentNodeIndices.length(); tParentNode++)
            {
                tGeometryFieldValues(tParentNode) = aGeometry->evaluate_field_value(mParentNodeIndices(tParentNode), mParentNodeCoordinates(tParentNode));
            }

            // Return interpolated value
            return Matrix<DDRMat>(mBasisValues * tGeometryFieldValues)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
