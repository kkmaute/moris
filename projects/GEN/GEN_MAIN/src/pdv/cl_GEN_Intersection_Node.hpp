#ifndef MORIS_CL_GEN_INTERSECTION_NODE_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_HPP

#include "cl_GEN_Child_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        class Intersection_Node : public Child_Node
        {

        protected:
            real mLocalCoordinate;
            bool mIsIntersected;
            std::weak_ptr<Geometry> mInterfaceGeometry;
            real mIsocontourThreshold;

        private:
            std::shared_ptr<Intersection_Node> mFirstParentNode;
            std::shared_ptr<Intersection_Node> mSecondParentNode;
            bool mFirstParentOnInterface;
            bool mSecondParentOnInterface;
            Matrix<DDRMat> mGlobalCoordinates;
            Matrix<DDRMat> mParentVector;
            Matrix<DDRMat> mCoordinateSensitivities;
            Matrix<DDSMat> mCoordinateDeterminingADVIDs;

            moris_id mPDVStartingID;
            bool mPDVStartingIDSet = false;

            moris_id mNodeID = -1;
            moris_index mNodeOwner = -1;

        public:

            /**
             * Constructor
             *
             * @param aLocalCoordinate Local coordinate inside of parent edge
             * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
             * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
             * @param aFirstParentNodeIndex Index of the first parent of this node
             * @param aSecondParentNodeIndex Index of the second parent of this node
             * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent of this node
             * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent of this node
             * @param aAncestorNodeIndices Node indices of the ancestors of this intersection node
             * @param aAncestorNodeCoordinates Coordinates of the ancestors of this intersection node
             * @param aAncestorBasisFunction Basis function of the ancestor topology
             * @param aInterfaceGeometry Geometry that intersects the parent to create this node
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining interface parent nodes based on geometry value
             * @param aIntersectionTolerance Tolerance for determining interface parent nodes with intersection distance
             */
            Intersection_Node(
                    real                               aLocalCoordinate,
                    std::shared_ptr<Intersection_Node> aFirstParentNode,
                    std::shared_ptr<Intersection_Node> aSecondParentNode,
                    real                               aFirstParentNodeIndex,
                    real                               aSecondParentNodeIndex,
                    const Matrix<DDRMat>&              aFirstParentNodeLocalCoordinates,
                    const Matrix<DDRMat>&              aSecondParentNodeLocalCoordinates,
                    Matrix<DDUMat>                     aAncestorNodeIndices,
                    Cell<Matrix<DDRMat>>               aAncestorNodeCoordinates,
                    const xtk::Basis_Function&         aAncestorBasisFunction,
                    std::shared_ptr<Geometry>          aInterfaceGeometry,
                    real                               aIsocontourThreshold,
                    real                               aIsocontourTolerance,
                    real                               aIntersectionTolerance);

            /**
             * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
             * ancestor nodes.
             *
             * @return Sensitivities
             */
            Matrix<DDRMat> get_dcoordinate_dadv();

            /**
             * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
             *
             * @return ADV IDs
             */
            Matrix<DDSMat> get_coordinate_determining_adv_ids();

            /**
             * Returns if the parent edge is intersected (if the local coordinate of the intersection lies between
             * -1 and 1)
             *
             * @return If the edge is intersected
             */
            bool parent_edge_is_intersected();

            /**
             * Returns if the first parent used to create this node is on the geoemtry interface already.
             *
             * @return If the first parent is on the interface
             */
            bool first_parent_on_interface();

            /**
             * Returns if the second parent used to create this node is on the geometry interface already.
             *
             * @return If the second parent is on the interface
             */
            bool second_parent_on_interface();

            /**
             * Gets the local coordinate of this intersection node inside of the parent edge.
             *
             * @return Local coordinate
             */
            real get_local_coordinate();

            /**
             * Gets all global coordinate values for this intersection node.
             *
             * @return Global coordinates
             */
            const Matrix<DDRMat>& get_global_coordinates();

            /**
             * Get the value of a coordinate of this node
             *
             * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
             * @return Coordinate value
             */
            real get_coordinate_value(uint aCoordinateIndex);

            /**
             * Gets the number of PDVs on this intersection node.
             *
             * @return Number of PDVs
             */
            uint get_num_pdvs();

            /**
             * Sets the starting index to be able to use the intersection coordinates of this node as PDVs
             *
             * @param aPDVStartingID The global index of the first PDV on the host
             */
            void set_starting_pdv_id(moris_id aPDVStartingID);

            /**
             * Get the starting global index for the intersection coordinate PDVs
             *
             * @return The global index of the first PDV on the host
             */
            moris_id get_starting_pdv_id();

            /**
             * Set the node ID for this node.
             *
             * @param aNodeID Node ID
             */
            void set_id(moris_id aNodeID);

            /**
             * Set the owning processor for this node.
             *
             * @param aNodeOwner Owning processor
             */
            void set_owner(moris_index aNodeOwner);

            /**
             * Get the ID for this node.
             *
             * @return Node ID
             */
            moris_id get_id();

            /**
             * Get the owning processor for this node.
             *
             * @return Owning processor
             */
            moris_index get_owner();

        private:

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @param aAncestorIndex Ancestor index
             * @return Local coordinate sensitivity
             */
            virtual real get_dxi_dfield_from_ancestor(uint aAncestorIndex) = 0;

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix<DDRMat> get_dxi_dcoordinate_first_parent() = 0;

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix<DDRMat> get_dxi_dcoordinate_second_parent() = 0;

            /**
             * Function for appending to the coordinate sensitivities member variable, eliminating duplicate code
             *
             * @param aSensitivitiesToAdd Sensitivities to add
             */
            void join_coordinate_sensitivities(const Matrix<DDRMat>& aSensitivitiesToAdd);

            /**
             * Function for appending to the depending ADV IDs member variable, eliminating duplicate code
             *
             * @param aIDsToAdd IDs to add
             */
            void join_adv_ids(const Matrix<DDSMat>& aIDsToAdd);

        };
    }
}

#endif //MORIS_CL_GEN_INTERSECTION_NODE_HPP
