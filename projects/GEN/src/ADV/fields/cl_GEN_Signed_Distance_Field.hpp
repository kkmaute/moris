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
#include "cl_Cell.hpp"

namespace moris::ge
{
    class Signed_Distance_Field : public Field_Discrete_Integration
    {

      private:
        std::string      mObjectPath;
        Matrix< DDRMat > mObjectOffset = { { 0, 0, 0 } };
        real             mShift        = 0;

        Matrix< DDRMat > mValues;

      public:
        /**
         * Constructor
         *
         * @param aMesh Mesh with the level set fields
         * @param aFieldNames Names of the fields
         */
        Signed_Distance_Field( std::string         aObjectPath,
                Matrix< DDRMat >          aObjectOffset,
                real                      aSDFShift   = 0 );

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        real get_field_value( uint aNodeIndex );

        void evaluate_nodal_values();

        void reset_nodal_data();

        const Matrix< DDRMat >&
        get_dfield_dadvs( uint aNodeIndex );

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
            return mObjectOffset;
        }
    };
}
