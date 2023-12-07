/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Discrete_Integration.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    class Field_Discrete_Integration : public Field
    {

      protected:
        mtk::Mesh_Pair mMeshPair;

      private:
        Cell< std::shared_ptr< Child_Node > > mChildNodes;

      public:

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        Field_Discrete_Integration(
                Matrix< DDRMat > aConstants,
                uint             aNumOriginalNodes );

        /**
         * Constructor that sets all field variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aSharedADVIds Shared ADV IDs needed for this field
         * @param aMeshPair Mesh pair for this field
         * @param aName Name of this field
         */
        Field_Discrete_Integration(
                const Matrix< DDSMat >& aSharedADVIds,
                mtk::Mesh_Pair          aMeshPair,
                std::string             aName );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) final;

        /**
         * Gets a field value of a derived node.
         *
         * @param aDerivedNode Derived node
         * @return Field value
         */
        real get_field_value( Derived_Node* aDerivedNode ) final;

        /**
         * Given a node index, returns the field value
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        virtual real get_field_value( uint aNodeIndex ) = 0;

        /**
         * Given a node index or coordinate, returns a matrix all sensitivities.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) final;

        /**
         * Gets a vector of the field derivatives with respect to ADVs of a derived node.
         *
         * @param aDerivedNode Derived node
         * @return d(field value)/d(ADV_j)
         */
        const Matrix< DDRMat >& get_dfield_dadvs( Derived_Node* aDerivedNode ) final;

        /**
         * Given a node index, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @return Vector of sensitivities
         */
        virtual const Matrix< DDRMat >& get_dfield_dadvs( uint aNodeIndex ) = 0;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) final;

        /**
         * Gets the IDs of ADVs that this field depends on for evaluations at a derived node.
         *
         * @param aDerivedNode Derived node
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids( Derived_Node* aDerivedNode ) final;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations for non-child nodes.
         *
         * @param aNodeIndex Node index
         * @return Determining ADV IDs at this node
         */
        virtual Matrix< DDSMat > get_determining_adv_ids( uint aNodeIndex );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;

        /**
         * Add a new child node for evaluation.
         *
         * @param aNodeIndex Index of the child node
         * @param aChildNode Contains information about how the child node was created
         */
        void add_child_node( uint aNodeIndex, std::shared_ptr< Child_Node > aChildNode ) override;

        /**
         * Resets all child nodes, called when a new XTK mesh is being created.
         */
        void reset_nodal_data() override;

        /**
         * Gets an MTK field, if this field needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;
    };
}
