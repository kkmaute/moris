#include <string>
#include "cl_MTK_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR.hpp" //HMR/src
namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------

        Mesh::Mesh( HMR & aHMR,
                const uint & aActivationPattern ) : mHMR( aHMR )
        {

            // get number of meshes
            uint tNumberOfMeshes = mHMR.get_number_of_lagrange_meshes();

            // count number of meshes
            mNumberOfBlocks = 0;

            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = mHMR.get_lagrange_mesh_by_index( k );

                // test if mesh uses active pattern
                if ( tMesh->get_activation_pattern() == aActivationPattern )
                {
                    // increment counter
                    ++mNumberOfBlocks;
                }
            }

            // initialize block cell
            mBlocks.resize( mNumberOfBlocks, nullptr );

            // initialize counter
            uint tCount = 0;

            // create blocks
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = mHMR.get_lagrange_mesh_by_index( k );

                // test if mesh uses active pattern
                if ( tMesh->get_activation_pattern() == aActivationPattern )
                {
                    // create block counter
                    mBlocks( tCount ) = new hmr::Block( tMesh, k );

                    // create default name
                    std::string tLabel = "MESH_" +std::to_string( k );

                    // set default name
                    mBlocks( tCount++ )->set_label( tLabel );
                }
            }
        }

//-----------------------------------------------------------------------------

        Mesh::~Mesh()
        {
            // delete block pointers
            for( auto tBlock: mBlocks )
            {
                delete tBlock;
            }
        }

//-----------------------------------------------------------------------------

        uint
        Mesh::get_number_of_blocks() const
        {
            return mNumberOfBlocks;
        }

//-----------------------------------------------------------------------------

        mtk::Block *
        Mesh::get_block_by_index( const moris_index & aIndex )
        {
            return mBlocks( aIndex );
        }

//-----------------------------------------------------------------------------

        const mtk::Block *
        Mesh::get_block_by_index( const moris_index & aIndex ) const
        {
            return mBlocks( aIndex );
        }

//-----------------------------------------------------------------------------

        void
        Mesh::finalize()
        {
            mHMR.finalize();
        }

//-----------------------------------------------------------------------------

        Mat< uint >
        Mesh::get_communication_table() const
        {
            return mHMR.get_communication_table();
        }

//-----------------------------------------------------------------------------

        mtk::Field *
        Mesh::create_field( const std::string & aLabel )
        {
            // fixme: rethink the concept of multiple blocks on HMR

            // create field, by default, use block 0
            return new Field( aLabel, mBlocks( 0 ) );
        }

//-----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
