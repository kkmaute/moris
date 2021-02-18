    #include <iostream>
#include <cstdio>

#include "cl_FEM_Field.hpp"
// HD5 c-interface


namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Field::Field( std::shared_ptr<mtk::Mesh_Manager>   aMeshManager,
                      uint const                         & aMeshIndex,
                      uint const                         & aDiscretizationMeshIndex )
        : mtk::Field( aMeshManager, aMeshIndex, aDiscretizationMeshIndex )
        {
            mtk::Interpolation_Mesh* tInterpolationMesh =
                    mMeshManager->get_interpolation_mesh( mMeshIndex );

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

        }

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
