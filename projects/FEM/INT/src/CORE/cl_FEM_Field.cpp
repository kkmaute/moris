    #include <iostream>
#include <cstdio>

#include "cl_FEM_Field.hpp"
// HD5 c-interface


namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Field::Field( mtk::Mesh_Pair * aMeshPair,
                      uint const     & aDiscretizationMeshIndex )
        : mtk::Field( aMeshPair, aDiscretizationMeshIndex )
        {
            mtk::Interpolation_Mesh* tInterpolationMesh = aMeshPair->mInterpolationMesh;

            // get number of nodes on block
            uint tNumberOfVertices = tInterpolationMesh->get_num_nodes();

            // set size of node values
            mNodalValues.set_size( tNumberOfVertices, 1, MORIS_REAL_MIN );
        }

        Field::~Field()
        {

        }

        //-----------------------------------------------------------------------------

        void Field::set_field_type( const uint & aType )
        {
             //enum mtk::Field_Type tType = static_cast< mtk::Field_Type >( aType );
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
             else
             {
                 MORIS_ERROR( false, "Field::set_field_from_file(), field type not known. New types can be implemented here.");
             }
        }

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
