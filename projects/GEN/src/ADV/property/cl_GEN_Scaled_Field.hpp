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

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Discrete_Interpolation.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris::ge
{
    // TODO a scaled field should not have this inheritance to begin with
    class Scaled_Field : public Field_Discrete_Interpolation
    {

    private:
        std::shared_ptr< Field > mField;

    public:
        /**
         * Constructor
         *
         * @param aADVs ADV vector
         * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the property variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         * @param aField Field that this property will scale
         */
        Scaled_Field(
                std::shared_ptr< Field > aField,
                ADV_ARG_TYPES )
                : Field_Discrete_Interpolation( mtk::Mesh_Pair( nullptr, nullptr ), ADV_ARGS )
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
         * @return Property value
         */
        real get_base_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index, evaluates the sensitivity of the property field with respect to all of the
         * property variables.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_base_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix<DDSMat> get_base_determining_adv_ids(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Sets the dependencies of this property after they have been found by update_dependencies().
         *
         * @param aDependencyFields Fields that this property depends on.
         */
        void set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields );

    };
}
