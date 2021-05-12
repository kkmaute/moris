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
        private:
            real xi1;
            real eta1;
            real xi2;
            real eta2;
            bool mLinear;
            bool mBilinearCase1;

        public:

            /**
             * Constructor
             *
             * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
             * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
             * @param aFirstParentNodeIndex Index of the first parent of this node
             * @param aSecondParentNodeIndex Index of the second parent of this node
             * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent node with respect to
             * the given ancestors
             * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent node with respect to
             * the given ancestors
             * @param aAncestorNodeIndices Ancestor node indices
             * @param aAncestorNodeCoordinates Ancestor node coordinates
             * @param aInterfaceGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining interface parent nodes based on geometry value
             * @param aIntersectionTolerance Tolerance for determining interface parent nodes with intersection distance
             */
            Intersection_Node_Bilinear(
                    std::shared_ptr<Intersection_Node> aFirstParentNode,
                    std::shared_ptr<Intersection_Node> aSecondParentNode,
                    uint                               aFirstParentNodeIndex,
                    uint                               aSecondParentNodeIndex,
                    const Matrix<DDRMat>&              aFirstParentNodeLocalCoordinates,
                    const Matrix<DDRMat>&              aSecondParentNodeLocalCoordinates,
                    const Matrix<DDUMat>&              aAncestorNodeIndices,
                    const Cell<Matrix<DDRMat>>&        aAncestorNodeCoordinates,
                    std::shared_ptr<Geometry>          aInterfaceGeometry,
                    real                               aIsocontourThreshold = 0.0,
                    real                               aIsocontourTolerance = 0.0,
                    real                               aIntersectionTolerance = 0.0);

        private:

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @return Local coordinate sensitivity
             */
            real get_dxi_dfield_from_ancestor(uint aAncestorIndex);

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix<DDRMat> get_dxi_dcoordinate_first_parent();

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix<DDRMat> get_dxi_dcoordinate_second_parent();

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
            real get_local_coordinate(
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
