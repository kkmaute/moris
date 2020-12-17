#ifndef MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP

#include "cl_GEN_Intersection_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        class Intersection_Node_Linear : public Intersection_Node
        {

        public:

            /**
             * Constructor
             *
             * @param aFirstNodeIndex Index of the first parent of this node
             * @param aSecondNodeIndex Index of the second parent of this node
             * @param aFirstNodeCoordinates Coordinates of the first parent of this node
             * @param aSecondNodeCoordinates Coordinates of the second parent of this node
             * @param aInterfaceGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining if parent nodes are on the interface or not
             */
            Intersection_Node_Linear(
                    uint                      aFirstNodeIndex,
                    uint                      aSecondNodeIndex,
                    const Matrix<DDRMat>&     aFirstNodeCoordinates,
                    const Matrix<DDRMat>&     aSecondNodeCoordinates,
                    std::shared_ptr<Geometry> aInterfaceGeometry,
                    real                      aIsocontourThreshold,
                    real                      aIsocontourTolerance);

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
             * @param aFirstNodeIndex Index of the first parent of this node
             * @param aSecondNodeIndex Index of the second parent of this node
             * @param aFirstNodeCoordinates Coordinates of the first parent of this node
             * @param aSecondNodeCoordinates Coordinates of the second parent of this node
             * @param aInterfaceGeometry Geometry that intersects the parent to create this node
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @return Local coordinates
             */
            Matrix<DDRMat> interpolate_local_coordinates(
                    uint                      aFirstNodeIndex,
                    uint                      aSecondNodeIndex,
                    const Matrix<DDRMat>&     aFirstNodeCoordinates,
                    const Matrix<DDRMat>&     aSecondNodeCoordinates,
                    std::shared_ptr<Geometry> aInterfaceGeometry,
                    real                      aIsocontourThreshold);

        };
    }
}

#endif //MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP
