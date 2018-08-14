#include <string>
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Interface.hpp" //HMR/src

#include "cl_HMR.hpp" //HMR/src
namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------

        Interface::Interface( HMR & aHMR ) : mHMR( aHMR )
        {
            // get number of meshes
            uint tNumberOfMeshes = aHMR.get_number_of_lagrange_meshes();

            // count number of meshes that are non-zero
            mNumberOfBlocks = 0;
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = aHMR.get_lagrange_mesh_by_index( k );

                if( tMesh != NULL )
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
                auto tMesh = aHMR.get_lagrange_mesh_by_index( k );

                if( tMesh != NULL )
                {
                    // create block counter
                    mBlocks( tCount ) = new Block( tMesh, k );

                    // create default name
                    std::string tLabel = "ORDER_" +std::to_string( k );

                    // set default name
                    mBlocks( tCount++ )->set_label( tLabel );
                }
            }
        }

//-----------------------------------------------------------------------------

        Interface::~Interface()
        {
            // delete block pointers
            for( auto tBlock: mBlocks )
            {
                delete tBlock;
            }
        }

//-----------------------------------------------------------------------------

        uint
        Interface::get_number_of_blocks() const
        {
            return mNumberOfBlocks;
        }

//-----------------------------------------------------------------------------

        Block *
        Interface::get_block_by_index( const uint & aIndex )
        {
            return mBlocks( aIndex );
        }

//-----------------------------------------------------------------------------

        void
        Interface::finalize()
        {
            mHMR.finalize();
        }

//-----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
