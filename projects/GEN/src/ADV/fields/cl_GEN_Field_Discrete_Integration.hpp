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

namespace moris::gen
{
    class Field_Discrete_Integration : public Field
    {

      protected:
        mtk::Mesh_Pair mMeshPair;

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
         * @param aNodeManager Node manager
         * @return Field value
         */
        real get_field_value(
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) final;

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
         * @param aSensitivities Sensitivities to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_dfield_dadvs(
                Matrix< DDRMat >&   aSensitivities,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) final;

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
         * @param aDeterminingADVIDs Determining ADV IDs to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_determining_adv_ids(
                Matrix< DDSMat >&   aDeterminingADVIDs,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) final;

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
         * Gets an MTK field, if this field needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;
    };
}
