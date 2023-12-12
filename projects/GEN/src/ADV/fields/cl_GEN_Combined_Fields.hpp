/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Combined_Fields.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    class Combined_Fields : public Field
    {

      private:
        Cell< std::shared_ptr< Field > > mFields;
        real                             mScale;

      public:
        /**
         * Combined field constructor
         *
         * @param aFields Created fields
         * @param aUseMinimum Whether or not to use minimum or maximum value when combining fields
         */
        Combined_Fields(
                Cell< std::shared_ptr< Field > > aFields,
                bool                             aUseMinimum = true,
                std::string                      aName = "" );

        /**
         * Given a node coordinate, returns the minimum (or maximum) field value of all fields that have been combined.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates vector of coordinate values
         * @return Minimum field value, after scaling has been applied
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets a field value of a derived node.
         *
         * @param aDerivedNode Derived node
         * @return Field value
         */
        real get_field_value( Derived_Node* aDerivedNode ) override;

        /**
         * Given a node index or coordinate, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets a vector of the field derivatives with respect to ADVs of a derived node.
         *
         * @param aDerivedNode Derived node
         * @return d(field value)/d(ADV_j)
         */
        const Matrix< DDRMat >& get_dfield_dadvs( Derived_Node* aDerivedNode ) override;

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the IDs of ADVs that this field depends on for evaluations at a derived node.
         *
         * @param aDerivedNode Derived node
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids( Derived_Node* aDerivedNode ) override;

        /**
         * Gets an MTK field, if this field needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;
    };
}
