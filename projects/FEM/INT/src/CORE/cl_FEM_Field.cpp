/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field.cpp
 *
 */

#include <iostream>
#include <cstdio>

#include "cl_FEM_Field.hpp"
#include "cl_Matrix.hpp"
// HD5 c-interface

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Field::Field( mtk::Mesh_Pair        aMeshPair,
                mtk::Field_Entity_Type aFieldEntityType,
                uint                        aDiscretizationMeshIndex )
                : mtk::Field( aMeshPair, 1, aFieldEntityType )
        {
            mUpdateNodalValues = false;

            mFieldImplementation = mtk::Field_Implementation::FEM;
        }

        //-----------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //-----------------------------------------------------------------------------

        void
        Field::set_field_type( const Vector< mtk::Field_Type >& aType )
        {
            mFieldType = aType;

            mValues.set_size( mValues.n_rows(), mFieldType.size(), MORIS_REAL_MIN );
        }

        //-----------------------------------------------------------------------------

        void
        Field::set_IQI_name( const std::string& aString )
        {
            mIQIName = aString;

            mPopulateFieldWithIQI = true;
        }

        //-----------------------------------------------------------------------------

        const std::string&
        Field::get_IQI_name()
        {
            return mIQIName;
        }

        //-----------------------------------------------------------------------------

        void
        Field::get_values(
                Matrix< IndexMat > const &      aIndex,
                Matrix< DDRMat >&               aValues,
                Vector< mtk::Field_Type > const & aFieldTypes )
        {
            // FIXME translate field types into index. implement map
            uint tNumFields = mFieldType.size();

            MORIS_ASSERT( tNumFields != 0, "Field::get_values(), No field types set for this field" );

            // FIXME this can be saved as a member variable or done in a smarter way
            moris::Matrix< IndexMat > tFieldIndex( tNumFields, 1 );
            for ( uint Ik = 0; Ik < tNumFields; Ik++ )
            {
                tFieldIndex( Ik ) = Ik;
            }

            this->get_value(
                    aIndex,
                    aValues,
                    tFieldIndex );
        }

        //-----------------------------------------------------------------------------

        void
        Field::set_field_from_file(
                const std::string& aString,
                const uint         aTimeIndex,
                const uint         aFieldIndex )
        {
            // detect file type
            std::string tType = aString.substr( aString.find_last_of( "." ) + 1, aString.length() );

            if ( tType == "hdf5" || tType == "h5" )
            {
                MORIS_ERROR( aTimeIndex == 0 && aFieldIndex == 0,
                        "Field::set_field_from_file - when loading from hdf5 time and vector indices need to be zero" );

                this->load_nodal_values_from_hdf5( aString );
            }
            else if ( tType == "exo" || tType == "e" )
            {
                load_field_from_exodus( aString, aTimeIndex, { { aFieldIndex } } );
            }
            else
            {
                MORIS_ERROR( false, "Field::set_field_from_file(), field type not known. New types can be implemented here." );
            }
        }

        //-----------------------------------------------------------------------------

        void
        Field::set_field_to_file( const std::string& aString )
        {
            mOutputFilePath = aString;
        }

        //------------------------------------------------------------------------------

        void
        Field::output_field_to_file()
        {
            if ( not mOutputFilePath.empty() )
            {
                // detect file type
                std::string tType = mOutputFilePath.substr( mOutputFilePath.find_last_of( "." ) + 1, mOutputFilePath.length() );

                if ( tType == "hdf5" || tType == "h5" )
                {
                    this->save_nodal_values_to_hdf5( mOutputFilePath );
                }
                else if ( tType == "exo" || tType == "e" )
                {
                    this->save_field_to_exodus( mOutputFilePath );
                }
                else
                {
                    MORIS_ERROR( false,
                            "Field::output_field_to_file(), field type not known. New types can be implemented here." );
                }
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
} /* namespace moris */
