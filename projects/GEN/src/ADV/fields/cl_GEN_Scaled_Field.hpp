/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Scaled_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris::gen
{
    class Scaled_Field : public Field
    {

    private:
        std::shared_ptr< Field > mField;

    public:
        /**
         * Constructor
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         * @param aField Field that this field will scale
         */
        Scaled_Field(
                std::shared_ptr< Field > aField,
                ADV_ARG_TYPES )
                : Field( ADV_ARGS )
                , mField( aField )
        {
            VARIABLE_CHECK( 1 );
            MORIS_ERROR( mField, "A scaled field must be given an input field." );
        }

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets a field value of a derived node.
         *
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         * @return Field value
         */
        real get_field_value(
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) override;

        /**
         * Given a node index, evaluates the sensitivity of the field with respect to all of the
         * field variables.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

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
                const Node_Manager& aNodeManager ) override;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Vector< sint > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the IDs of ADVs that this field depends on for evaluations at a derived node.
         *
         * @param aDeterminingADVIDs Determining ADV IDs to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_determining_adv_ids(
                Vector< sint >&   aDeterminingADVIDs,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) override;

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
         * Sets the dependencies of this field after they have been found by update_dependencies().
         *
         * @param aDependencyFields Fields that this field depends on.
         */
        void set_dependencies( Vector< std::shared_ptr< Field > > aDependencyFields ) override;

        /**
         * Gets an MTK field, if this field needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;

    };
}
