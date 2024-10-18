/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Signed_Distance_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"

namespace moris::gen
{
    class Signed_Distance_Field : public Field_Discrete_Integration
    {

      private:
        std::string    mObjectPath;
        Vector< real > mObjectOffset = { 0, 0, 0 };
        real           mShift        = 0;

        Matrix< DDRMat > mValues;

      public:
        /**
         * Constructor
         *
         * @param aMesh Mesh with the level set fields
         * @param aFieldNames Names of the fields
         */
        Signed_Distance_Field( const std::string& aObjectPath,
                Vector< real >                    aObjectOffset,
                real                              aSDFShift = 0 );

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        real get_field_value( uint aNodeIndex ) override;

        void evaluate_nodal_values();

        /**
         * Resets the nodal data for this signed distance field, evaluating the nodal values on the given mesh.
         *
         * @param aMesh New interpolation mesh for evaluating nodal values
         */
        void reset_nodal_data( mtk::Interpolation_Mesh* aMesh ) override;

        const Matrix< DDRMat >&
        get_dfield_dadvs( uint aNodeIndex ) override;

      private:
        /**
         * Given a node index, evaluates the sensitivity of the field with respect to all of the
         * field variables. This is currently not implemented for a signed distance field.
         *
         * @param aNodeIndex Node index
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >&
        get_base_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates )
        {
            MORIS_ERROR( false, "get_base_dfield_dadvs(), not implemented for a signed distance field" );
            return mValues;
        }
    };
}    // namespace moris::gen
