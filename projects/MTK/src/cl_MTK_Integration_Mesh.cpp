#include <string>
#include <iomanip>
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

            for( auto tListofDoubleSideSets : mListofDoubleSideSets )
            {
                delete tListofDoubleSideSets;
            }
            mListofDoubleSideSets.clear();
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
 
        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_block_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mBlockSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------

         void
         Integration_Mesh::get_block_set_names_with_color(
                 moris_index const        & aColor,
                 moris::Cell<std::string> & aSetNames)
         {
             MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");

             uint tNumSets = mBlockSetToColor(aColor).size();
             aSetNames.resize( tNumSets );
             for(uint iS=0;iS<tNumSets;iS++)
             {
                 aSetNames(iS) = mBlockSetToColor(aColor)(iS)->get_set_name();
             }
         }

        // ----------------------------------------------------------------------------

         moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_side_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mSideSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------

         moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_double_side_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mDoubleSideSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------
 

        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_all_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mAllSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------
        void
        Integration_Mesh::print_sets_by_colors()
        {
            for(moris::uint i = 0; i < mBlockSetToColor.size(); i++)
            {
                std::cout<<"\n Color: "<<std::setw(8)<<i<<std::endl;
                std::cout<<"    Blocks: "<<std::endl;
                for(moris::uint iS = 0; iS <mBlockSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mBlockSetToColor(i)(iS)->get_set_name()<<std::endl;
                }

                std::cout<<"    Side Sets: "<<std::endl;
                for(moris::uint iS = 0; iS <mSideSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mSideSetToColor(i)(iS)->get_set_name()<<std::endl;
                }

                std::cout<<"    Double Side Sets: "<<std::endl;
                for(moris::uint iS = 0; iS <mDoubleSideSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mDoubleSideSetToColor(i)(iS)->get_set_name()<<std::endl;
                }
            }
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

        void
        Integration_Mesh::add_double_side_set(mtk::Set* aDblSideSet)
        {
            mListofDoubleSideSets.push_back(aDblSideSet); 
            this->collect_all_sets();  
        }


        void 
        Integration_Mesh::collect_all_sets( bool aSetShape )
        {
            // reset
            mListOfAllSets.clear();
            mSetNameToIndexMap.clear();

            uint tCounter = 0;

            // Append block sets to list of all sets
            mListOfAllSets.append( mListofBlocks );

            // Append side sets to list of all sets
            mListOfAllSets.append( mListofSideSets );

            //FIXME implement cell topology for side set

            // Append double side sets to list of all sets
            mListOfAllSets.append( mListofDoubleSideSets );

            // iterate through all sets and register their index in the map
            for( uint Ik = 0; Ik < mListOfAllSets.size(); Ik++ )
            {
                mListOfAllSets( Ik )->set_set_index( Ik );

                mSetNameToIndexMap[ mListOfAllSets( Ik )->get_set_name() ] = Ik;
            }

            // add the cell topology to the set
            for( uint Ik = 0; Ik < mListofBlocks.size(); Ik++ )
            {
                std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                // set the blockset cell topology and shape
                mListOfAllSets( tCounter )->set_cell_topology( this->get_blockset_topology( tSetName ) );
                if ( aSetShape == true )
                {
                    mListOfAllSets( tCounter )->set_IG_cell_shape( this->get_IG_blockset_shape( tSetName ) );
                    mListOfAllSets( tCounter )->set_IP_cell_shape( this->get_IP_blockset_shape( tSetName ) );
                }
                tCounter++;
            }

            if ( aSetShape == true )
            {
                // add the cell topology to the set
                for( uint Ik = 0; Ik < mListofSideSets.size(); Ik++ )
                {
                    std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                    // set the sideset cell topology and shape
                    mListOfAllSets( tCounter   )->set_IG_cell_shape( CellShape::STRAIGHT );
                    mListOfAllSets( tCounter++ )->set_IP_cell_shape( CellShape::STRAIGHT );
                }

                // add the cell topology to the set
                for( uint Ik = 0; Ik < mListofDoubleSideSets.size(); Ik++ )
                {
                    std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                    // set the sideset cell topology and shape
                    mListOfAllSets( tCounter   )->set_IG_cell_shape( CellShape::STRAIGHT );
                    mListOfAllSets( tCounter++ )->set_IP_cell_shape( CellShape::STRAIGHT );
                }
            }

            // setup color to set data
            this->setup_set_to_color();
        }

        // ----------------------------------------------------------------------------
        void
        Integration_Mesh::setup_set_to_color()
        {
            // determine maximum color for the sets (for data sizing)
            mMaxColor = 0;
            for(moris::uint i = 0; i < mListOfAllSets.size() ; i++)
            {
                Matrix<IndexMat> const & tSetColors = mListOfAllSets(i)->get_set_colors();

                // if the max set color is greater than current max
                // replace the current max
                if(tSetColors.numel() > 0)
                {
                    moris_index tMaxSetColor = tSetColors.max();
                    if(mMaxColor < tMaxSetColor)
                    {
                        mMaxColor = tMaxSetColor;
                    }
                }
            }

            // clear old data
            mBlockSetToColor.clear();
            mSideSetToColor.clear();
            mDoubleSideSetToColor.clear();
            mAllSetToColor.clear();

            // size outer cell size of member data
            mBlockSetToColor.resize(mMaxColor + 1);
            mSideSetToColor.resize(mMaxColor + 1);
            mDoubleSideSetToColor.resize(mMaxColor + 1);
            mAllSetToColor.resize(mMaxColor + 1);

            // iterate through block sets 
            for(moris::uint i = 0; i < mListofBlocks.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofBlocks(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mBlockSetToColor(tSetColors(j)).push_back(mListofBlocks(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofBlocks(i));
                }
            }

            // iterate through side sets 
            for(moris::uint i = 0; i < mListofSideSets.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofSideSets(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mSideSetToColor(tSetColors(j)).push_back(mListofSideSets(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofSideSets(i));
                }
            }

            // iterate through double side sets 
            for(moris::uint i = 0; i < mListofDoubleSideSets.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofDoubleSideSets(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mDoubleSideSetToColor(tSetColors(j)).push_back(mListofDoubleSideSets(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofDoubleSideSets(i));
                }
            }
        }
        // ----------------------------------------------------------------------------

        
        // ----------------------------------------------------------------------------

    }
}
