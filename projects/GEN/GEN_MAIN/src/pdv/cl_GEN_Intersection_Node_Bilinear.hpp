#ifndef MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP

#include "cl_GEN_Intersection_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        class Intersection_Node_Bilinear : public Intersection_Node
        {

        public:

            /**
             * Constructor
             *
             * @param aFirstParentNodeIndex First parent node index
             * @param aSecondParentNodeIndex Second parent node index
             * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent node with respect to
             * the given ancestors
             * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent node with respect to
             * the given ancestors
             * @param aFirstParentNodeGlobalCoordinates Global coordinates of the first parent node
             * @param aSecondParentNodeGlobalCoordinates Global coordinates of the second parent node
             * @param aAncestorNodeIndices Ancestor node indices
             * @param aAncestorNodeCoordinates Ancestor node coordinates
             * @param aInterfaceGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining if parent nodes are on the interface or not
             */
            Intersection_Node_Bilinear(
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
                    real                        aIsocontourTolerance);

            /**
             * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
             * ancestor nodes.
             *
             * @return Sensitivities
             */
            Matrix<DDRMat> get_ancestor_coordinate_sensitivities(uint aAncestorIndex);

        private:

            /**
             * Interpolate and return the local coordinates of this intersection node. Used to clean up constructor.
             *
             * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent node with respect to
             * the given ancestors
             * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent node with respect to
             * the given ancestors
             * @param aAncestorNodeIndices Ancestor node indices
             * @param aAncestorNodeCoordinates Ancestor node coordinates
             * @param aInterfaceGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @return Local coordinates
             */
            Matrix<DDRMat> interpolate_local_coordinates(
                    const Matrix<DDRMat>&       aFirstParentNodeLocalCoordinates,
                    const Matrix<DDRMat>&       aSecondParentNodeLocalCoordinates,
                    const Matrix<DDUMat>&       aAncestorNodeIndices,
                    const Cell<Matrix<DDRMat>>& aAncestorNodeCoordinates,
                    std::shared_ptr<Geometry>   aInterfaceGeometry,
                    real                        aIsocontourThreshold);

        };
    }
}

#endif //MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP
