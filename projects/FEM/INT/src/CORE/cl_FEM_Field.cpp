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

        Field::Field( mtk::Mesh_Pair aMeshPair,
                uint                 aDiscretizationMeshIndex )
        : mtk::Field( aMeshPair )
        {
            mtk::Interpolation_Mesh* tInterpolationMesh = aMeshPair.get_interpolation_mesh();

            // get number of nodes on block
            uint tNumberOfVertices = tInterpolationMesh->get_num_nodes();

            // set size of node values
            mNodalValues.set_size( tNumberOfVertices, 0, MORIS_REAL_MIN );

            mUpdateNodalValues = false;
        }

        //-----------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //-----------------------------------------------------------------------------

        void Field::set_field_type( const moris::Cell< mtk::Field_Type > & aType )
        {
            mFieldType = aType;

            mNodalValues.set_size( mNodalValues.n_rows(), mFieldType.size(), MORIS_REAL_MIN );
        }

        //-----------------------------------------------------------------------------

        void Field::set_IQI_name( const std::string & aString )
        {
            mIQIName = aString;

            mPopulateFieldWithIQI = true;
        }

        //-----------------------------------------------------------------------------

        const std::string & Field::get_IQI_name()
        {
            return mIQIName;
        }

        //-----------------------------------------------------------------------------

        void Field::get_nodal_values(
                Matrix< IndexMat > const      & aNodeIndex,
                Matrix< DDRMat >              & aNodalValues,
                Cell< mtk::Field_Type > const & aFieldTypes)
        {
            // FIXME translate field types into index. implement map
            uint tNumFields = mFieldType.size();

            MORIS_ASSERT( tNumFields != 0, "Field::get_nodal_values(), No field types set for this field" );

            //FIXME this can be saved as a member variable or done in a smarter way
            moris::Matrix< IndexMat > tFieldIndex( tNumFields, 1 );
            for( uint Ik = 0; Ik<tNumFields; Ik++ )
            {
                tFieldIndex( Ik ) = Ik;
            }

            this->get_nodal_value(
                    aNodeIndex,
                    aNodalValues,
                    tFieldIndex );
        }

        //-----------------------------------------------------------------------------

        void Field::set_field_from_file( const std::string & aString )
        {
            // detect file type
            std::string tType = aString.substr( aString.find_last_of(".")+1, aString.length() );

            if( tType == "hdf5" || tType == "h5" )
            {
                this->load_nodal_values_from_hdf5( aString );
            }
            else if ( tType == "exo" || tType == "e" )
            {
                load_field_from_exodus( aString );
            }
            else
            {
                MORIS_ERROR( false, "Field::set_field_from_file(), field type not known. New types can be implemented here.");
            }
        }

        //-----------------------------------------------------------------------------

        void Field::set_field_to_file( const std::string & aString )
        {
            mOutputFilePath = aString;
        }

        //------------------------------------------------------------------------------

        void Field::output_field_to_file()
        {
            if( not mOutputFilePath.empty() )
            {
                // detect file type
                std::string tType = mOutputFilePath.substr( mOutputFilePath.find_last_of(".")+1, mOutputFilePath.length() );

                if( tType == "hdf5" || tType == "h5" )
                {
                    this->save_nodal_values_to_hdf5( mOutputFilePath );
                }
                else if( tType == "exo" || tType == "e" )
                {
                    this->save_field_to_exodus( mOutputFilePath );
                }
                else
                {
                    MORIS_ERROR( false,
                            "Field::output_field_to_file(), field type not known. New types can be implemented here.");
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
