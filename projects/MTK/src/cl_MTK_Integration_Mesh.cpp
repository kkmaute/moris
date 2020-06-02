#include <string>

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------

        Integration_Mesh::~Integration_Mesh()
        {
            for( auto tListofBlocks : mListofBlocks )
            {
                delete tListofBlocks;
            }
            mListofBlocks.clear();

            for( auto tListofSideSets : mListofSideSets )
            {
                delete tListofSideSets;
            }
            mListofSideSets.clear();
        }

        // ----------------------------------------------------------------------------

        moris::uint Integration_Mesh::get_num_sets() const
        {
            return mListOfAllSets.size();
        }

        // ----------------------------------------------------------------------------

        moris::mtk::Set * Integration_Mesh::get_set_by_name( std::string aSetLabel ) const
        {
            moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );

            return mListOfAllSets( tSetIndex );
        }

        // ----------------------------------------------------------------------------

        moris::mtk::Set * Integration_Mesh::get_set_by_index( moris_index aIndex ) const
        {
            return mListOfAllSets( aIndex );
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_set_index_by_name( std::string aSetLabel )
        {
            if ( ! mSetNameToIndexMap.key_exists( aSetLabel ) )
            {
                std::string tErrMsg = "Integration_Mesh::get_set_index_by_name - Set with name does not exists: " + aSetLabel;

                MORIS_ERROR( false, tErrMsg.c_str() );
            }

            return mSetNameToIndexMap.find( aSetLabel );
        }

        // ----------------------------------------------------------------------------

        std::string Integration_Mesh::get_block_set_label(moris_index aBlockSetOrdinal) const
        {
            MORIS_ERROR(0, "get_block_set_label has no default implementation");
            return "ERROR";
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_block_set_index(std::string aBlockSetLabel) const
        {
            MORIS_ERROR(0, "get_block_set_index has no default implementation");
            return MORIS_INDEX_MAX;
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_blocks() const
        {
            return mListofBlocks.size();
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_side_set() const
        {
            return mListofSideSets.size();
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_double_side_set() const
        {
            return mListofDoubleSideSets.size();
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_double_sided_set_index(std::string aDoubleSideSetLabel) const
        {
            MORIS_ERROR( false, "not implemented");
            return 0;
        }

        // ----------------------------------------------------------------------------

        void Integration_Mesh::collect_all_sets()
        {
            // reset
            mListOfAllSets.clear();
            mSetNameToIndexMap.clear();

            uint tCounter = 0;

            mListOfAllSets.append( mListofBlocks );

            for( uint Ik = 0; Ik < mListofBlocks.size(); Ik++ )
            {
                std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                mListOfAllSets( tCounter++ )->set_cell_topology( this->get_blockset_topology( tSetName ) );
            }

            mListOfAllSets.append( mListofSideSets );

            //FIXME implement cell tpopology for side set

            mListOfAllSets.append( mListofDoubleSideSets );

            for( uint Ik = 0; Ik < mListOfAllSets.size(); Ik++ )
            {
                mListOfAllSets( Ik )->set_set_index( Ik );

                mSetNameToIndexMap[ mListOfAllSets( Ik )->get_set_name() ] = Ik;
            }
        }

        // ----------------------------------------------------------------------------

    }
}
