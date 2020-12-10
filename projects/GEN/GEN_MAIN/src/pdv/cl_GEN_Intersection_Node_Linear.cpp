#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_XTK_Linear_Basis_Functions.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node_Linear::Intersection_Node_Linear(
                uint                      aFirstNodeIndex,
                uint                      aSecondNodeIndex,
                const Matrix<DDRMat>&     aFirstNodeCoordinates,
                const Matrix<DDRMat>&     aSecondNodeCoordinates,
                std::shared_ptr<Geometry> aInterfaceGeometry,
                real                      aIsocontourThreshold,
                real                      aIsocontourTolerance)
                : Intersection_Node(
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        aFirstNodeCoordinates,
                        aSecondNodeCoordinates,
                        {{aFirstNodeIndex, aSecondNodeIndex}},
                        {aFirstNodeCoordinates, aSecondNodeCoordinates},
                        xtk::Linear_Basis_Function(),
                        interpolate_local_coordinates(
                                aFirstNodeIndex,
                                aSecondNodeIndex,
                                aFirstNodeCoordinates,
                                aSecondNodeCoordinates,
                                aInterfaceGeometry,
                                aIsocontourThreshold),
                        aInterfaceGeometry,
                        aIsocontourThreshold,
                        aIsocontourTolerance)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Linear::get_ancestor_coordinate_sensitivities(uint aAncestorIndex)
        {
            // Get geometry field values
            real tPhi0 = mInterfaceGeometry->get_field_value( mAncestorNodeIndices(0), mAncestorNodeCoordinates(0) );
            real tPhi1 = mInterfaceGeometry->get_field_value( mAncestorNodeIndices(1), mAncestorNodeCoordinates(1) );

            // Get geometry field sensitivity with respect to ADVs
            const Matrix<DDRMat>& tFieldSensitivity = mInterfaceGeometry->get_field_sensitivities(
                    mAncestorNodeIndices(aAncestorIndex),
                    mAncestorNodeCoordinates(aAncestorIndex));

            // Compute sensitivity of the global coordinate with respect to the field value
            Matrix<DDRMat> tCoordinateSensitivity =
                    (tPhi0 * (aAncestorIndex == 1) - tPhi1 * (aAncestorIndex == 0)) / std::pow((tPhi0 - tPhi1), 2)
                    * (mAncestorNodeCoordinates(1) - mAncestorNodeCoordinates(0));

            // Compute full sensitivity of global coordinates with respect to ADVs
            return trans(tCoordinateSensitivity) * tFieldSensitivity;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Linear::interpolate_local_coordinates(
                uint                      aFirstNodeIndex,
                uint                      aSecondNodeIndex,
                const Matrix<DDRMat>&     aFirstNodeCoordinates,
                const Matrix<DDRMat>&     aSecondNodeCoordinates,
                std::shared_ptr<Geometry> aInterfaceGeometry,
                real                      aIsocontourThreshold)
        {
            // Interface geometry values
            Matrix<DDRMat> tInterfaceGeometryValues = 
                    {{aInterfaceGeometry->get_field_value( aFirstNodeIndex, aFirstNodeCoordinates )},
                    { aInterfaceGeometry->get_field_value( aSecondNodeIndex, aSecondNodeCoordinates )}};

            // Interpolate
            Matrix<DDRMat> tLocalCoordinates = Interpolation::linear_interpolation_value(tInterfaceGeometryValues, aIsocontourThreshold);

            // Must store local coordinate in parent edge
            mLocalCoordinate = tLocalCoordinates(0);
            
            return tLocalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------
    
    }
}
