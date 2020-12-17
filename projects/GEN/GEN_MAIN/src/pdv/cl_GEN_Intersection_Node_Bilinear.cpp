#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_XTK_Quad_4_Basis_Function.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node_Bilinear::Intersection_Node_Bilinear(
                uint                        aFirstParentNodeIndex,
                uint                        aSecondParentNodeIndex,
                const Matrix<DDRMat>&       aFirstParentNodeLocalCoordinates,
                const Matrix<DDRMat>&       aSecondParentNodeLocalCoordinates,
                const Matrix<DDRMat>&       aFirstParentNodeGlobalCoordinates,
                const Matrix<DDRMat>&       aSecondParentNodeGlobalCoordinates,
                const Matrix<DDUMat>&       aAncestorNodeIndices,
                const Cell<Matrix<DDRMat>>& aAncestorNodeCoordinates,
                std::shared_ptr<Geometry>   aInterfaceGeometry,
                real                        aIsocontourThreshold,
                real                        aIsocontourTolerance)
                : Intersection_Node(
                        aFirstParentNodeIndex,
                        aSecondParentNodeIndex,
                        aFirstParentNodeGlobalCoordinates,
                        aSecondParentNodeGlobalCoordinates,
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        xtk::Quad_4_Basis_Function(),
                        interpolate_local_coordinates(
                                aFirstParentNodeLocalCoordinates,
                                aSecondParentNodeLocalCoordinates,
                                aAncestorNodeIndices,
                                aAncestorNodeCoordinates,
                                aInterfaceGeometry,
                                aIsocontourThreshold),
                        aInterfaceGeometry,
                        aIsocontourThreshold,
                        aIsocontourTolerance)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Bilinear::get_ancestor_coordinate_sensitivities(uint aAncestorIndex)
        {
            MORIS_ERROR(false, "Sensitivities not implemented yet for bilinear intersections.");
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node_Bilinear::interpolate_local_coordinates(
                const Matrix<DDRMat>&       aFirstParentNodeLocalCoordinates,
                const Matrix<DDRMat>&       aSecondParentNodeLocalCoordinates,
                const Matrix<DDUMat>&       aAncestorNodeIndices,
                const Cell<Matrix<DDRMat>>& aAncestorNodeCoordinates,
                std::shared_ptr<Geometry>   aInterfaceGeometry,
                real                        aIsocontourThreshold)
        {
            // Local coordinates
            real xi1 = aFirstParentNodeLocalCoordinates(0);
            real eta1 = aFirstParentNodeLocalCoordinates(1);
            real xi2 = aSecondParentNodeLocalCoordinates(0);
            real eta2 = aSecondParentNodeLocalCoordinates(1);

            // Geometry field values
            real phi1 = aInterfaceGeometry->get_field_value(aAncestorNodeIndices(0), aAncestorNodeCoordinates(0));
            real phi2 = aInterfaceGeometry->get_field_value(aAncestorNodeIndices(1), aAncestorNodeCoordinates(1));
            real phi3 = aInterfaceGeometry->get_field_value(aAncestorNodeIndices(2), aAncestorNodeCoordinates(2));
            real phi4 = aInterfaceGeometry->get_field_value(aAncestorNodeIndices(3), aAncestorNodeCoordinates(3));

            // Isocontour threshold
            real iso = aIsocontourThreshold;

            // Calculate local coordinates
            Matrix<DDRMat> tLocalCoordinates(2, 1);
            tLocalCoordinates(0) = -(4*sqrt((eta1*eta1*phi1*phi1*xi2*xi2)/16 - (eta1*eta1*phi1*phi1*xi2)/8 + (eta1*eta1*phi1*phi1)/16 - (eta1*eta1*phi1*phi2*xi2*xi2)/8 + (eta1*eta1*phi1*phi2)/8 + (eta1*eta1*phi1*phi3*xi2*xi2)/8 - (eta1*eta1*phi1*phi3)/8 - (eta1*eta1*phi1*phi4*xi2*xi2)/8 + (eta1*eta1*phi1*phi4*xi2)/4 - (eta1*eta1*phi1*phi4)/8 + (eta1*eta1*phi2*phi2*xi2*xi2)/16 + (eta1*eta1*phi2*phi2*xi2)/8 + (eta1*eta1*phi2*phi2)/16 - (eta1*eta1*phi2*phi3*xi2*xi2)/8 - (eta1*eta1*phi2*phi3*xi2)/4 - (eta1*eta1*phi2*phi3)/8 + (eta1*eta1*phi2*phi4*xi2*xi2)/8 - (eta1*eta1*phi2*phi4)/8 + (eta1*eta1*phi3*phi3*xi2*xi2)/16 + (eta1*eta1*phi3*phi3*xi2)/8 + (eta1*eta1*phi3*phi3)/16 - (eta1*eta1*phi3*phi4*xi2*xi2)/8 + (eta1*eta1*phi3*phi4)/8 + (eta1*eta1*phi4*phi4*xi2*xi2)/16 - (eta1*eta1*phi4*phi4*xi2)/8 + (eta1*eta1*phi4*phi4)/16 - (eta1*eta2*phi1*phi1*xi1*xi2)/8 + (eta1*eta2*phi1*phi1*xi1)/8 + (eta1*eta2*phi1*phi1*xi2)/8 - (eta1*eta2*phi1*phi1)/8 + (eta1*eta2*phi1*phi2*xi1*xi2)/4 - (eta1*eta2*phi1*phi2)/4 - (eta1*eta2*phi1*phi3*xi1*xi2)/4 + (eta1*eta2*phi1*phi3)/4 + (eta1*eta2*phi1*phi4*xi1*xi2)/4 - (eta1*eta2*phi1*phi4*xi1)/4 - (eta1*eta2*phi1*phi4*xi2)/4 + (eta1*eta2*phi1*phi4)/4 - (eta1*eta2*phi2*phi2*xi1*xi2)/8 - (eta1*eta2*phi2*phi2*xi1)/8 - (eta1*eta2*phi2*phi2*xi2)/8 - (eta1*eta2*phi2*phi2)/8 + (eta1*eta2*phi2*phi3*xi1*xi2)/4 + (eta1*eta2*phi2*phi3*xi1)/4 + (eta1*eta2*phi2*phi3*xi2)/4 + (eta1*eta2*phi2*phi3)/4 - (eta1*eta2*phi2*phi4*xi1*xi2)/4 + (eta1*eta2*phi2*phi4)/4 - (eta1*eta2*phi3*phi3*xi1*xi2)/8 - (eta1*eta2*phi3*phi3*xi1)/8 - (eta1*eta2*phi3*phi3*xi2)/8 - (eta1*eta2*phi3*phi3)/8 + (eta1*eta2*phi3*phi4*xi1*xi2)/4 - (eta1*eta2*phi3*phi4)/4 - (eta1*eta2*phi4*phi4*xi1*xi2)/8 + (eta1*eta2*phi4*phi4*xi1)/8 + (eta1*eta2*phi4*phi4*xi2)/8 - (eta1*eta2*phi4*phi4)/8 + (eta1*phi1*phi1*xi1*xi2)/8 - (eta1*phi1*phi1*xi1)/8 - (eta1*phi1*phi1*xi2*xi2)/8 + (eta1*phi1*phi1*xi2)/8 - (eta1*phi1*phi2*xi1*xi2)/4 + (eta1*phi1*phi2*xi2*xi2)/4 - (3*eta1*phi1*phi3*xi1)/4 + (3*eta1*phi1*phi3*xi2)/4 + iso*eta1*phi1*xi1 - iso*eta1*phi1*xi2 + (eta1*phi2*phi2*xi1*xi2)/8 + (eta1*phi2*phi2*xi1)/8 - (eta1*phi2*phi2*xi2*xi2)/8 - (eta1*phi2*phi2*xi2)/8 + (3*eta1*phi2*phi4*xi1)/4 - (3*eta1*phi2*phi4*xi2)/4 - iso*eta1*phi2*xi1 + iso*eta1*phi2*xi2 - (eta1*phi3*phi3*xi1*xi2)/8 - (eta1*phi3*phi3*xi1)/8 + (eta1*phi3*phi3*xi2*xi2)/8 + (eta1*phi3*phi3*xi2)/8 + (eta1*phi3*phi4*xi1*xi2)/4 - (eta1*phi3*phi4*xi2*xi2)/4 + iso*eta1*phi3*xi1 - iso*eta1*phi3*xi2 - (eta1*phi4*phi4*xi1*xi2)/8 + (eta1*phi4*phi4*xi1)/8 + (eta1*phi4*phi4*xi2*xi2)/8 - (eta1*phi4*phi4*xi2)/8 - iso*eta1*phi4*xi1 + iso*eta1*phi4*xi2 + (eta2*eta2*phi1*phi1*xi1*xi1)/16 - (eta2*eta2*phi1*phi1*xi1)/8 + (eta2*eta2*phi1*phi1)/16 - (eta2*eta2*phi1*phi2*xi1*xi1)/8 + (eta2*eta2*phi1*phi2)/8 + (eta2*eta2*phi1*phi3*xi1*xi1)/8 - (eta2*eta2*phi1*phi3)/8 - (eta2*eta2*phi1*phi4*xi1*xi1)/8 + (eta2*eta2*phi1*phi4*xi1)/4 - (eta2*eta2*phi1*phi4)/8 + (eta2*eta2*phi2*phi2*xi1*xi1)/16 + (eta2*eta2*phi2*phi2*xi1)/8 + (eta2*eta2*phi2*phi2)/16 - (eta2*eta2*phi2*phi3*xi1*xi1)/8 - (eta2*eta2*phi2*phi3*xi1)/4 - (eta2*eta2*phi2*phi3)/8 + (eta2*eta2*phi2*phi4*xi1*xi1)/8 - (eta2*eta2*phi2*phi4)/8 + (eta2*eta2*phi3*phi3*xi1*xi1)/16 + (eta2*eta2*phi3*phi3*xi1)/8 + (eta2*eta2*phi3*phi3)/16 - (eta2*eta2*phi3*phi4*xi1*xi1)/8 + (eta2*eta2*phi3*phi4)/8 + (eta2*eta2*phi4*phi4*xi1*xi1)/16 - (eta2*eta2*phi4*phi4*xi1)/8 + (eta2*eta2*phi4*phi4)/16 - (eta2*phi1*phi1*xi1*xi1)/8 + (eta2*phi1*phi1*xi1*xi2)/8 + (eta2*phi1*phi1*xi1)/8 - (eta2*phi1*phi1*xi2)/8 + (eta2*phi1*phi2*xi1*xi1)/4 - (eta2*phi1*phi2*xi1*xi2)/4 + (3*eta2*phi1*phi3*xi1)/4 - (3*eta2*phi1*phi3*xi2)/4 - iso*eta2*phi1*xi1 + iso*eta2*phi1*xi2 - (eta2*phi2*phi2*xi1*xi1)/8 + (eta2*phi2*phi2*xi1*xi2)/8 - (eta2*phi2*phi2*xi1)/8 + (eta2*phi2*phi2*xi2)/8 - (3*eta2*phi2*phi4*xi1)/4 + (3*eta2*phi2*phi4*xi2)/4 + iso*eta2*phi2*xi1 - iso*eta2*phi2*xi2 + (eta2*phi3*phi3*xi1*xi1)/8 - (eta2*phi3*phi3*xi1*xi2)/8 + (eta2*phi3*phi3*xi1)/8 - (eta2*phi3*phi3*xi2)/8 - (eta2*phi3*phi4*xi1*xi1)/4 + (eta2*phi3*phi4*xi1*xi2)/4 - iso*eta2*phi3*xi1 + iso*eta2*phi3*xi2 + (eta2*phi4*phi4*xi1*xi1)/8 - (eta2*phi4*phi4*xi1*xi2)/8 - (eta2*phi4*phi4*xi1)/8 + (eta2*phi4*phi4*xi2)/8 + iso*eta2*phi4*xi1 - iso*eta2*phi4*xi2 + (phi1*phi1*xi1*xi1)/16 - (phi1*phi1*xi1*xi2)/8 + (phi1*phi1*xi2*xi2)/16 - (phi1*phi2*xi1*xi1)/8 + (phi1*phi2*xi1*xi2)/4 - (phi1*phi2*xi2*xi2)/8 - (phi1*phi3*xi1*xi1)/8 + (phi1*phi3*xi1*xi2)/4 - (phi1*phi3*xi2*xi2)/8 + (phi1*phi4*xi1*xi1)/8 - (phi1*phi4*xi1*xi2)/4 + (phi1*phi4*xi2*xi2)/8 + (phi2*phi2*xi1*xi1)/16 - (phi2*phi2*xi1*xi2)/8 + (phi2*phi2*xi2*xi2)/16 + (phi2*phi3*xi1*xi1)/8 - (phi2*phi3*xi1*xi2)/4 + (phi2*phi3*xi2*xi2)/8 - (phi2*phi4*xi1*xi1)/8 + (phi2*phi4*xi1*xi2)/4 - (phi2*phi4*xi2*xi2)/8 + (phi3*phi3*xi1*xi1)/16 - (phi3*phi3*xi1*xi2)/8 + (phi3*phi3*xi2*xi2)/16 - (phi3*phi4*xi1*xi1)/8 + (phi3*phi4*xi1*xi2)/4 - (phi3*phi4*xi2*xi2)/8 + (phi4*phi4*xi1*xi1)/16 - (phi4*phi4*xi1*xi2)/8 + (phi4*phi4*xi2*xi2)/16) + eta1*phi1 + eta1*phi2 - eta2*phi1 - eta1*phi3 - eta2*phi2 - eta1*phi4 + eta2*phi3 + eta2*phi4 + phi1*xi1 - phi1*xi2 - phi2*xi1 + phi2*xi2 - phi3*xi1 + phi3*xi2 + phi4*xi1 - phi4*xi2 - 2*eta1*phi1*xi1 + eta1*phi1*xi2 + 2*eta1*phi2*xi1 + eta2*phi1*xi1 - eta1*phi2*xi2 - 2*eta1*phi3*xi1 - eta2*phi2*xi1 + eta1*phi3*xi2 + 2*eta1*phi4*xi1 + eta2*phi3*xi1 - eta1*phi4*xi2 - eta2*phi4*xi1)/(2*(eta1*phi1*xi1 - eta1*phi1*xi2 - eta1*phi2*xi1 - eta2*phi1*xi1 + eta1*phi2*xi2 + eta1*phi3*xi1 + eta2*phi1*xi2 + eta2*phi2*xi1 - eta1*phi3*xi2 - eta1*phi4*xi1 - eta2*phi2*xi2 - eta2*phi3*xi1 + eta1*phi4*xi2 + eta2*phi3*xi2 + eta2*phi4*xi1 - eta2*phi4*xi2));
            tLocalCoordinates(1) = (4*sqrt((eta1*eta1*phi1*phi1*xi2*xi2)/16 - (eta1*eta1*phi1*phi1*xi2)/8 + (eta1*eta1*phi1*phi1)/16 - (eta1*eta1*phi1*phi2*xi2*xi2)/8 + (eta1*eta1*phi1*phi2)/8 + (eta1*eta1*phi1*phi3*xi2*xi2)/8 - (eta1*eta1*phi1*phi3)/8 - (eta1*eta1*phi1*phi4*xi2*xi2)/8 + (eta1*eta1*phi1*phi4*xi2)/4 - (eta1*eta1*phi1*phi4)/8 + (eta1*eta1*phi2*phi2*xi2*xi2)/16 + (eta1*eta1*phi2*phi2*xi2)/8 + (eta1*eta1*phi2*phi2)/16 - (eta1*eta1*phi2*phi3*xi2*xi2)/8 - (eta1*eta1*phi2*phi3*xi2)/4 - (eta1*eta1*phi2*phi3)/8 + (eta1*eta1*phi2*phi4*xi2*xi2)/8 - (eta1*eta1*phi2*phi4)/8 + (eta1*eta1*phi3*phi3*xi2*xi2)/16 + (eta1*eta1*phi3*phi3*xi2)/8 + (eta1*eta1*phi3*phi3)/16 - (eta1*eta1*phi3*phi4*xi2*xi2)/8 + (eta1*eta1*phi3*phi4)/8 + (eta1*eta1*phi4*phi4*xi2*xi2)/16 - (eta1*eta1*phi4*phi4*xi2)/8 + (eta1*eta1*phi4*phi4)/16 - (eta1*eta2*phi1*phi1*xi1*xi2)/8 + (eta1*eta2*phi1*phi1*xi1)/8 + (eta1*eta2*phi1*phi1*xi2)/8 - (eta1*eta2*phi1*phi1)/8 + (eta1*eta2*phi1*phi2*xi1*xi2)/4 - (eta1*eta2*phi1*phi2)/4 - (eta1*eta2*phi1*phi3*xi1*xi2)/4 + (eta1*eta2*phi1*phi3)/4 + (eta1*eta2*phi1*phi4*xi1*xi2)/4 - (eta1*eta2*phi1*phi4*xi1)/4 - (eta1*eta2*phi1*phi4*xi2)/4 + (eta1*eta2*phi1*phi4)/4 - (eta1*eta2*phi2*phi2*xi1*xi2)/8 - (eta1*eta2*phi2*phi2*xi1)/8 - (eta1*eta2*phi2*phi2*xi2)/8 - (eta1*eta2*phi2*phi2)/8 + (eta1*eta2*phi2*phi3*xi1*xi2)/4 + (eta1*eta2*phi2*phi3*xi1)/4 + (eta1*eta2*phi2*phi3*xi2)/4 + (eta1*eta2*phi2*phi3)/4 - (eta1*eta2*phi2*phi4*xi1*xi2)/4 + (eta1*eta2*phi2*phi4)/4 - (eta1*eta2*phi3*phi3*xi1*xi2)/8 - (eta1*eta2*phi3*phi3*xi1)/8 - (eta1*eta2*phi3*phi3*xi2)/8 - (eta1*eta2*phi3*phi3)/8 + (eta1*eta2*phi3*phi4*xi1*xi2)/4 - (eta1*eta2*phi3*phi4)/4 - (eta1*eta2*phi4*phi4*xi1*xi2)/8 + (eta1*eta2*phi4*phi4*xi1)/8 + (eta1*eta2*phi4*phi4*xi2)/8 - (eta1*eta2*phi4*phi4)/8 + (eta1*phi1*phi1*xi1*xi2)/8 - (eta1*phi1*phi1*xi1)/8 - (eta1*phi1*phi1*xi2*xi2)/8 + (eta1*phi1*phi1*xi2)/8 - (eta1*phi1*phi2*xi1*xi2)/4 + (eta1*phi1*phi2*xi2*xi2)/4 - (3*eta1*phi1*phi3*xi1)/4 + (3*eta1*phi1*phi3*xi2)/4 + iso*eta1*phi1*xi1 - iso*eta1*phi1*xi2 + (eta1*phi2*phi2*xi1*xi2)/8 + (eta1*phi2*phi2*xi1)/8 - (eta1*phi2*phi2*xi2*xi2)/8 - (eta1*phi2*phi2*xi2)/8 + (3*eta1*phi2*phi4*xi1)/4 - (3*eta1*phi2*phi4*xi2)/4 - iso*eta1*phi2*xi1 + iso*eta1*phi2*xi2 - (eta1*phi3*phi3*xi1*xi2)/8 - (eta1*phi3*phi3*xi1)/8 + (eta1*phi3*phi3*xi2*xi2)/8 + (eta1*phi3*phi3*xi2)/8 + (eta1*phi3*phi4*xi1*xi2)/4 - (eta1*phi3*phi4*xi2*xi2)/4 + iso*eta1*phi3*xi1 - iso*eta1*phi3*xi2 - (eta1*phi4*phi4*xi1*xi2)/8 + (eta1*phi4*phi4*xi1)/8 + (eta1*phi4*phi4*xi2*xi2)/8 - (eta1*phi4*phi4*xi2)/8 - iso*eta1*phi4*xi1 + iso*eta1*phi4*xi2 + (eta2*eta2*phi1*phi1*xi1*xi1)/16 - (eta2*eta2*phi1*phi1*xi1)/8 + (eta2*eta2*phi1*phi1)/16 - (eta2*eta2*phi1*phi2*xi1*xi1)/8 + (eta2*eta2*phi1*phi2)/8 + (eta2*eta2*phi1*phi3*xi1*xi1)/8 - (eta2*eta2*phi1*phi3)/8 - (eta2*eta2*phi1*phi4*xi1*xi1)/8 + (eta2*eta2*phi1*phi4*xi1)/4 - (eta2*eta2*phi1*phi4)/8 + (eta2*eta2*phi2*phi2*xi1*xi1)/16 + (eta2*eta2*phi2*phi2*xi1)/8 + (eta2*eta2*phi2*phi2)/16 - (eta2*eta2*phi2*phi3*xi1*xi1)/8 - (eta2*eta2*phi2*phi3*xi1)/4 - (eta2*eta2*phi2*phi3)/8 + (eta2*eta2*phi2*phi4*xi1*xi1)/8 - (eta2*eta2*phi2*phi4)/8 + (eta2*eta2*phi3*phi3*xi1*xi1)/16 + (eta2*eta2*phi3*phi3*xi1)/8 + (eta2*eta2*phi3*phi3)/16 - (eta2*eta2*phi3*phi4*xi1*xi1)/8 + (eta2*eta2*phi3*phi4)/8 + (eta2*eta2*phi4*phi4*xi1*xi1)/16 - (eta2*eta2*phi4*phi4*xi1)/8 + (eta2*eta2*phi4*phi4)/16 - (eta2*phi1*phi1*xi1*xi1)/8 + (eta2*phi1*phi1*xi1*xi2)/8 + (eta2*phi1*phi1*xi1)/8 - (eta2*phi1*phi1*xi2)/8 + (eta2*phi1*phi2*xi1*xi1)/4 - (eta2*phi1*phi2*xi1*xi2)/4 + (3*eta2*phi1*phi3*xi1)/4 - (3*eta2*phi1*phi3*xi2)/4 - iso*eta2*phi1*xi1 + iso*eta2*phi1*xi2 - (eta2*phi2*phi2*xi1*xi1)/8 + (eta2*phi2*phi2*xi1*xi2)/8 - (eta2*phi2*phi2*xi1)/8 + (eta2*phi2*phi2*xi2)/8 - (3*eta2*phi2*phi4*xi1)/4 + (3*eta2*phi2*phi4*xi2)/4 + iso*eta2*phi2*xi1 - iso*eta2*phi2*xi2 + (eta2*phi3*phi3*xi1*xi1)/8 - (eta2*phi3*phi3*xi1*xi2)/8 + (eta2*phi3*phi3*xi1)/8 - (eta2*phi3*phi3*xi2)/8 - (eta2*phi3*phi4*xi1*xi1)/4 + (eta2*phi3*phi4*xi1*xi2)/4 - iso*eta2*phi3*xi1 + iso*eta2*phi3*xi2 + (eta2*phi4*phi4*xi1*xi1)/8 - (eta2*phi4*phi4*xi1*xi2)/8 - (eta2*phi4*phi4*xi1)/8 + (eta2*phi4*phi4*xi2)/8 + iso*eta2*phi4*xi1 - iso*eta2*phi4*xi2 + (phi1*phi1*xi1*xi1)/16 - (phi1*phi1*xi1*xi2)/8 + (phi1*phi1*xi2*xi2)/16 - (phi1*phi2*xi1*xi1)/8 + (phi1*phi2*xi1*xi2)/4 - (phi1*phi2*xi2*xi2)/8 - (phi1*phi3*xi1*xi1)/8 + (phi1*phi3*xi1*xi2)/4 - (phi1*phi3*xi2*xi2)/8 + (phi1*phi4*xi1*xi1)/8 - (phi1*phi4*xi1*xi2)/4 + (phi1*phi4*xi2*xi2)/8 + (phi2*phi2*xi1*xi1)/16 - (phi2*phi2*xi1*xi2)/8 + (phi2*phi2*xi2*xi2)/16 + (phi2*phi3*xi1*xi1)/8 - (phi2*phi3*xi1*xi2)/4 + (phi2*phi3*xi2*xi2)/8 - (phi2*phi4*xi1*xi1)/8 + (phi2*phi4*xi1*xi2)/4 - (phi2*phi4*xi2*xi2)/8 + (phi3*phi3*xi1*xi1)/16 - (phi3*phi3*xi1*xi2)/8 + (phi3*phi3*xi2*xi2)/16 - (phi3*phi4*xi1*xi1)/8 + (phi3*phi4*xi1*xi2)/4 - (phi3*phi4*xi2*xi2)/8 + (phi4*phi4*xi1*xi1)/16 - (phi4*phi4*xi1*xi2)/8 + (phi4*phi4*xi2*xi2)/16) - eta1*phi1 - eta1*phi2 + eta2*phi1 + eta1*phi3 + eta2*phi2 + eta1*phi4 - eta2*phi3 - eta2*phi4 - phi1*xi1 + phi1*xi2 + phi2*xi1 - phi2*xi2 + phi3*xi1 - phi3*xi2 - phi4*xi1 + phi4*xi2 + 2*eta1*phi1*xi1 - eta1*phi1*xi2 - 2*eta1*phi2*xi1 - eta2*phi1*xi1 + eta1*phi2*xi2 + 2*eta1*phi3*xi1 + eta2*phi2*xi1 - eta1*phi3*xi2 - 2*eta1*phi4*xi1 - eta2*phi3*xi1 + eta1*phi4*xi2 + eta2*phi4*xi1)/(2*(eta1*phi1*xi1 - eta1*phi1*xi2 - eta1*phi2*xi1 - eta2*phi1*xi1 + eta1*phi2*xi2 + eta1*phi3*xi1 + eta2*phi1*xi2 + eta2*phi2*xi1 - eta1*phi3*xi2 - eta1*phi4*xi1 - eta2*phi2*xi2 - eta2*phi3*xi1 + eta1*phi4*xi2 + eta2*phi3*xi2 + eta2*phi4*xi1 - eta2*phi4*xi2));

            // Must store local coordinate in parent edge
            mLocalCoordinate = 0.0; // FIXME, above equations are not working yet.
            
            return tLocalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------
    
    }
}
