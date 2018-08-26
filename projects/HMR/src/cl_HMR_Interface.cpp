#include <string>
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Interface.hpp" //HMR/src

#include "cl_HMR.hpp" //HMR/src
namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------

        Interface::Interface( HMR & aHMR, const uint & aActivationPattern ) : mHMR( aHMR )
        {

            // get number of meshes
            uint tNumberOfMeshes = mHMR.get_number_of_lagrange_meshes();

            // count number of meshes
            mNumberOfBlocks = 0;

            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = mHMR.get_lagrange_mesh_by_index( k );

                // test if mesh uses active pattern
                if ( tMesh->get_active_pattern() == aActivationPattern )
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
                if ( tMesh->get_active_pattern() == aActivationPattern )
                {
                    // create block counter
                    mBlocks( tCount ) = new Block( tMesh, k );

                    // create default name
                    std::string tLabel = "MESH_" +std::to_string( k );

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

        Mat< uint >
        Interface::get_communication_table() const
        {
            return mHMR.get_communication_table();
        }

//-----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
