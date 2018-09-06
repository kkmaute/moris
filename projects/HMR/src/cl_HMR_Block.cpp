#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Block.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        std::string
        Block::get_label() const
        {
            return mLabel;
        }

//------------------------------------------------------------------------------

        void
        Block::set_label( const std::string & aLabel )
        {
            mLabel = aLabel;
        }

//------------------------------------------------------------------------------

        luint
        Block::get_number_of_vertices() const
        {
            return mMesh->get_number_of_nodes_on_proc();
        }


//------------------------------------------------------------------------------

        mtk::Vertex *
        Block::get_vertex_by_index( const luint & aIndex )
        {
            return mMesh->get_node_by_index( aIndex );
        }

//------------------------------------------------------------------------------

        luint
        Block::get_number_of_cells() const
        {
            return mMesh->get_number_of_elements();
        }

//------------------------------------------------------------------------------

        mtk::Cell *
        Block::get_cell_by_index( const luint & aIndex )
        {
            return mMesh->get_element( aIndex );
        }

//------------------------------------------------------------------------------

        luint
        Block::get_id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        sint
        Block::get_number_of_adofs_used_by_proc() const
        {
            return mMesh->get_number_of_bsplines_on_proc();
        }

//------------------------------------------------------------------------------

        void
        Block::get_adof_map( map< moris_id, moris_index > & aAdofMap ) const
        {
            aAdofMap.clear();

            moris_index tNumberOfBSplines = mMesh->get_number_of_bsplines_on_proc();

            for( moris_index k=0; k<tNumberOfBSplines; ++k )
            {
                aAdofMap[ mMesh->get_bspline( k )->get_id() ] = k;
            }
        }

//------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
