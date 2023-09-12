/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node.hpp
 *
 */

#ifndef MORIS_CL_GEN_INTERSECTION_NODE_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_HPP

#include "cl_GEN_Child_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        //------------------------------------------------------------------------------
        class Intersection_Node : public Child_Node
        {
          protected:
            real             mLocalCoordinate;
            bool             mIsIntersected;
            Matrix< DDRMat > mParentVector;
            Matrix< DDSMat > mCoordinateDeterminingADVIDs;

            std::shared_ptr< Intersection_Node > mFirstParentNode;
            std::shared_ptr< Intersection_Node > mSecondParentNode;
            moris_index                          mFirstParentNodeIndex;
            moris_index                          mSecondParentNodeIndex;
            bool                                 mFirstParentOnInterface;
            bool                                 mSecondParentOnInterface;
            Matrix< DDRMat >                     mGlobalCoordinates;


          private:
            moris_id mPDVStartingID;
            bool     mPDVStartingIDSet = false;

            moris_id    mNodeID    = -1;
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
             */
            Intersection_Node(
                    real                                 aLocalCoordinate,
                    std::shared_ptr< Intersection_Node > aFirstParentNode,
                    std::shared_ptr< Intersection_Node > aSecondParentNode,
                    moris_index                          aFirstParentNodeIndex,
                    moris_index                          aSecondParentNodeIndex,
                    const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                    Matrix< DDUMat >                     aAncestorNodeIndices,
                    Cell< Matrix< DDRMat > >             aAncestorNodeCoordinates,
                    const Element_Intersection_Type      aAncestorBasisFunction );
            /**
             * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
             * ancestor nodes.
             *
             * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
             * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
             * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
             */
            virtual void get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) = 0;

            /**
             * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
             *
             * @return ADV IDs
             */
            virtual Matrix< DDSMat > get_coordinate_determining_adv_ids() = 0;

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
            const Matrix< DDRMat >& get_global_coordinates();

            /**
             * Get the value of a coordinate of this node
             *
             * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
             * @return Coordinate value
             */
            real get_coordinate_value( uint aCoordinateIndex );

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
            void set_starting_pdv_id( moris_id aPDVStartingID );

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
            void set_id( moris_id aNodeID );

            /**
             * Set the owning processor for this node.
             *
             * @param aNodeOwner Owning processor
             */
            void set_owner( moris_index aNodeOwner );

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

            moris_index
            get_first_parent_node_index()
            {
                return mFirstParentNodeIndex;
            }

            moris_index
            get_second_parent_node_index()
            {
                return mSecondParentNodeIndex;
            }

          protected:
            /**
             * Computes basic member data for all intersection node derived classes.
             * Must be called by lowest level child class constructors.
             *
             */
            void initialize(
                    const Element_Intersection_Type aAncestorBasisFunction,
                    const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates );

            /**
             * Function for appending to the depending ADV IDs member variable, eliminating duplicate code
             *
             * @param aIDsToAdd IDs to add
             */
            void join_adv_ids( const Matrix< DDSMat >& aIDsToAdd );

          private:
            /**
             * Computes the global coordinates of the intersection and the parents.
             * Used by initialize() to set mGlobalCoordinates member data. Implementation provided by child class.
             *
             * @return Matrix< DDRMat > Global location of the intersection node and its parents
             */
            virtual Matrix< DDRMat > compute_global_coordinates() = 0;

            /**
             * Computes the vector from the first parent to the second parent
             * Used by initialize() to set mParentVector member data.
             *
             * @return Matrix< DDRMat > vector from the first parent to the second parent. Size determined by dimensionality of problem.
             */
            Matrix< DDRMat > compute_parent_vector(
                    const Element_Intersection_Type aAncestorBasisFunction,
                    const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates );

            /**
             * Determines if the parent nodes are intersected.
             * Used by initialize() to set mIsIntersected. Implementation provided by child class.
             *
             * @return if the parent nodes are intersected
             * @return false if there is no intersection detected
             */
            virtual bool determine_is_intersected(
                const Element_Intersection_Type aAncestorBasisFunction,
                const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates
            ) = 0;

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @param aAncestorIndex Ancestor index
             * @return Local coordinate sensitivity
             */
            virtual real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) = 0;

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix< DDRMat > get_dxi_dcoordinate_first_parent() = 0;

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix< DDRMat > get_dxi_dcoordinate_second_parent() = 0;

            /**
             * Function for appending to the coordinate sensitivities member variable, eliminating duplicate code
             *
             * @param aSensitivitiesToAdd Sensitivities to add
             */
            void join_coordinate_sensitivities( const Matrix< DDRMat >& aSensitivitiesToAdd );
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_HPP
