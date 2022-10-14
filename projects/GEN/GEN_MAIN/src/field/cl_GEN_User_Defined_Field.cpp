/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Field.cpp
 *
 */

#include "cl_GEN_User_Defined_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Field::User_Defined_Field(
                Matrix< DDRMat > aConstants,
                Field_Function   aFieldFunction,
                Field_Parameters aParameters )
                : Field( aConstants, aParameters )
        {
            this->set_user_defined_functions( aFieldFunction, nullptr );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        User_Defined_Field::get_field_value( const Matrix< DDRMat >& aCoordinates )
        {
            return this->get_field_value_user_defined( aCoordinates, mFieldVariables );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        User_Defined_Field::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
        {
            this->get_dfield_dadvs_user_defined( aCoordinates, mFieldVariables, mSensitivities );
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        User_Defined_Field::get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities )
        {
            MORIS_ASSERT( aCoordinates.numel() == aSensitivities.numel(),
                    "User_Defined_Field::get_dfield_dcoordinates - inconsistent dimensions" );

            const real tPertubation = 1e-8;

            Matrix< DDRMat > tCoordinates( aCoordinates );

            for ( uint idim = 0; idim < aCoordinates.numel(); ++idim )
            {
                tCoordinates( idim ) += tPertubation;
                real tFieldValue = this->get_field_value_user_defined( tCoordinates, mFieldVariables );

                tCoordinates( idim ) -= 2.0 * tPertubation;
                tFieldValue -= this->get_field_value_user_defined( tCoordinates, mFieldVariables );

                aSensitivities( idim ) = tFieldValue / 2.0 / tPertubation;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        User_Defined_Field::set_user_defined_functions(
                Field_Function       aFieldFunction,
                Sensitivity_Function aSensitivityFunction )
        {
            // Set field evaluation function
            get_field_value_user_defined = aFieldFunction;

            // Check field evaluation function
            MORIS_ASSERT( std::isfinite( this->get_field_value_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables ) ),
                    "There is an error in a user-defined geometry field (field evaluates to nan/infinity)." );

            // Set sensitivity evaluation function
            if ( aSensitivityFunction == nullptr )
            {
                // Sensitivity function was not provided
                get_dfield_dadvs_user_defined = &( User_Defined_Field::no_sensitivities );
            }
            else
            {
                // Sensitivity function was provided
                get_dfield_dadvs_user_defined = aSensitivityFunction;

                // Check sensitivity function
                this->get_dfield_dadvs_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables, mSensitivities );

                // Check for row vector
                MORIS_ERROR( mSensitivities.n_rows() == 1,
                        "A user-defined geometry must provide a row vector for sensitivities." );

                // Check for size
                MORIS_ERROR( mSensitivities.n_cols() == mFieldVariables.size(),
                        "A user-defined geometry must have a sensitivity vector with a length equal to the total "
                        "number of geometry variables (ADVs + constants). sensitivities: %i   geom variables: %i \n",
                        mSensitivities.n_cols(),
                        mFieldVariables.size() );

                // Check for values not nan/infinity
                for ( uint tSensitivityIndex = 0; tSensitivityIndex < mSensitivities.n_cols(); tSensitivityIndex++ )
                {
                    MORIS_ASSERT( std::isfinite( mSensitivities( tSensitivityIndex ) ),
                            "There is an error in a user-defined geometry sensitivity (evaluates to nan/infinity)." );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        User_Defined_Field::no_sensitivities(
                const Matrix< DDRMat >& aCoordinates,
                const Cell< real* >&    aParameters,
                Matrix< DDRMat >&       aSensitivities )
        {
            MORIS_ERROR( false,
                    "A sensitivity evaluation function was not provided to a user-defined geometry. "
                    "Please make sure that you provide this function, or that sensitivities are not required." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        User_Defined_Field::compute_nodal_values()
        {
            // get interpolation mesh
            mtk::Mesh* tIPmesh = mMeshPair.get_interpolation_mesh();

            // make sure that nodal value matrix is properly sized
            mValues.resize( tIPmesh->get_num_nodes(), mNumberOfFields );

            // loop over all nodes
            for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
            {
                // loop over all fields
                for ( uint tFieldIndex = 0; tFieldIndex < mNumberOfFields; ++tFieldIndex )
                {
                    mValues( tNodeIndex, tFieldIndex ) = get_field_value_user_defined(
                            tIPmesh->get_node_coordinate( tNodeIndex ),
                            mFieldVariables );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
