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
            std::weak_ptr<Geometry> mInterfaceGeometry;
            real mLocalCoordinate;

        private:
            bool mFirstParentOnInterface;
            bool mSecondParentOnInterface;
            Matrix<DDRMat> mGlobalCoordinates;
            Cell<std::shared_ptr<Intersection_Node>> mNodeDependencies; // TODO

            moris_id mPDVStartingID;
            bool mPDVStartingIDSet = false;

            moris_id mNodeID = -1;
            moris_index mNodeOwner = -1;

        public:

            /**
             * Constructor
             *
             * @param aFirstParentNodeIndex Index of the first parent of this node
             * @param aSecondParentNodeIndex Index of the second parent of this node
             * @param aFirstParentNodeCoordinates Coordinates of the first parent of this node
             * @param aSecondParentNodeCoordinates Coordinates of the second parent of this node
             * @param aAncestorNodeIndices Node indices of the ancestors of this intersection node
             * @param aAncestorNodeCoordinates Coordinates of the ancestors of this intersection node
             * @param aAncestorBasisFunction Basis function of the ancestor topology
             * @param aAncestorLocalCoordinates Local coordinate of this intersection inside of the ancestor element
             * @param aInterfaceGeometry Geometry that intersects the parent to create this node
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining if parent nodes are on the interface or not
             */
            Intersection_Node(
                    uint                       aFirstParentNodeIndex,
                    uint                       aSecondParentNodeIndex,
                    const Matrix<DDRMat>&      aFirstParentNodeCoordinates,
                    const Matrix<DDRMat>&      aSecondParentNodeCoordinates,
                    Matrix<DDUMat>             aAncestorNodeIndices,
                    Cell<Matrix<DDRMat>>       aAncestorNodeCoordinates,
                    const xtk::Basis_Function& aAncestorBasisFunction,
                    Matrix<DDRMat>             aAncestorLocalCoordinates,
                    std::shared_ptr<Geometry>  aInterfaceGeometry,
                    real                       aIsocontourThreshold,
                    real                       aIsocontourTolerance);

            /**
             * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
             * ancestor nodes.
             *
             * @return Sensitivities
             */
            virtual Matrix<DDRMat> get_ancestor_coordinate_sensitivities(uint aAncestorIndex) = 0;

            /**
             * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
             *
             * @return ADV IDs
             */
            Matrix<DDSMat> get_ancestor_coordinate_determining_adv_ids(uint aAncestorIndex);

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
            Matrix<DDRMat> get_global_coordinates();

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

        };
    }
}

#endif //MORIS_CL_GEN_INTERSECTION_NODE_HPP
