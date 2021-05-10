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
                std::shared_ptr<Intersection_Node> aFirstNode,
                std::shared_ptr<Intersection_Node> aSecondNode,
                uint                               aFirstNodeIndex,
                uint                               aSecondNodeIndex,
                const Matrix<DDRMat>&              aFirstNodeCoordinates,
                const Matrix<DDRMat>&              aSecondNodeCoordinates,
                std::shared_ptr<Geometry>          aInterfaceGeometry,
                real                               aIsocontourThreshold,
                real                               aIsocontourTolerance,
                real                               aIntersectionTolerance)
                : Intersection_Node(
                        get_local_coordinate(
                                aFirstNodeIndex,
                                aSecondNodeIndex,
                                aFirstNodeCoordinates,
                                aSecondNodeCoordinates,
                                aInterfaceGeometry,
                                aIsocontourThreshold),
                        aFirstNode,
                        aSecondNode,
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        {{-1}},
                        {{1}},
                        {{aFirstNodeIndex, aSecondNodeIndex}},
                        {aFirstNodeCoordinates, aSecondNodeCoordinates},
                        xtk::Linear_Basis_Function(),
                        aInterfaceGeometry,
                        aIsocontourThreshold,
                        aIsocontourTolerance,
                        aIntersectionTolerance)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Intersection_Node_Linear::get_dxi_dfield_from_ancestor(uint aAncestorIndex)
        {
            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get geometry field values
            real tPhi0 = tLockedInterfaceGeometry->get_field_value( mAncestorNodeIndices(0), mAncestorNodeCoordinates(0) );
            real tPhi1 = tLockedInterfaceGeometry->get_field_value( mAncestorNodeIndices(1), mAncestorNodeCoordinates(1) );

            // Compute sensitivity of the local coordinate with respect to the field value
            return 2 * ((tPhi0 - mIsocontourThreshold) * (aAncestorIndex == 1)
                      - (tPhi1 - mIsocontourThreshold) * (aAncestorIndex == 0)) / std::pow((tPhi1 - tPhi0), 2);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Linear::get_dxi_dcoordinate_first_parent()
        {
            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
            Matrix<DDRMat> tCoordinateSensitivities(1, 2);
            tLockedInterfaceGeometry->get_dfield_dcoordinates(
                    mAncestorNodeIndices(0),
                    mAncestorNodeCoordinates(0),
                    tCoordinateSensitivities);
            return this->get_dxi_dfield_from_ancestor(0) * tCoordinateSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Linear::get_dxi_dcoordinate_second_parent()
        {
            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
            Matrix<DDRMat> tCoordinateSensitivities(1, 2);
            tLockedInterfaceGeometry->get_dfield_dcoordinates(
                    mAncestorNodeIndices(1),
                    mAncestorNodeCoordinates(1),
                    tCoordinateSensitivities);
            return this->get_dxi_dfield_from_ancestor(1) * tCoordinateSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Intersection_Node_Linear::get_local_coordinate(
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

            return tLocalCoordinates(0);
        }

        //--------------------------------------------------------------------------------------------------------------
    
    }
}
