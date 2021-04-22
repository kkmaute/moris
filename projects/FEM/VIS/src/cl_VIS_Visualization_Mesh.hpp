/*
 * cl_VIS_Visualization_Mesh.hpp
 *
 *  Created on: Dec 02, 2019
 *      Author: Schmidt
 */

#ifndef PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_

#include "../../../MRS/CON/src/cl_Cell.hpp"
#include "../../../MTK/src/cl_MTK_Cell.hpp"
#include "../../../MTK/src/cl_MTK_Integration_Mesh.hpp"
#include "../../../MTK/src/cl_MTK_Mesh_Manager.hpp"       //MTK/src
#include "../../../MTK/src/cl_MTK_Set.hpp"
#include "../../../MTK/src/cl_MTK_Vertex.hpp"

#include "cl_VIS_Vertex_Visualization.hpp"
#include "cl_VIS_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"

namespace moris
{
namespace vis
{
class Visualization_Mesh : public mtk::Mesh
{
protected:
    moris::Cell< moris::mtk::Set * > mListofBlocks;
    moris::Cell< moris::mtk::Set * > mListofSideSets;

    moris::Cell< moris::mtk::Vertex const * >          mAllVertices;
    moris::Cell< mtk::Cell const * >                   mAllCells;
    moris::Cell< moris::Cell< mtk::Cell * > >          mCellsOnSet;
    moris::Cell< moris::Cell< mtk::Vertex * > >        mVerticesOnSet;
    moris::Cell< moris::Cell< const mtk::Cluster * > > mClustersOnBlock;

    const bool                                  mOnlyPrimary;

    moris::Cell< moris::mtk::Set * > mListOfAllSets;

    map< std::string, moris_index >  mSetNameToIndexMap;

public:

    Visualization_Mesh( moris::Cell< moris::mtk::Set * >                   aListofBlocks,
                        moris::Cell< moris::Cell< const mtk::Cluster * > > aClustersOnBlock,
                        moris::Cell< moris::Cell< mtk::Cell * > >          aCellsOnBlock,
                        moris::Cell< moris::Cell< mtk::Vertex * > >        aVerticesOnBlock,
                        const bool                                         aOnlyPrimary ) : mListofBlocks( aListofBlocks ),
                                                                                            mCellsOnSet( aCellsOnBlock ),
                                                                                            mVerticesOnSet( aVerticesOnBlock ),
                                                                                            mClustersOnBlock( aClustersOnBlock ),
                                                                                            mOnlyPrimary( aOnlyPrimary )
    {
        this->collect_all_sets();

        // All vertices/cells
        mAllCells = moris::Cell< mtk::Cell const * >( this->get_num_elems() );
        uint tNumBlocks = this->get_num_blocks();
        for(uint Ik=0; Ik < tNumBlocks; Ik ++)
        {
            for(uint Ii=0; Ii < mVerticesOnSet( Ik ).size(); Ii ++)
            {
                mAllVertices.push_back(mVerticesOnSet( Ik )( Ii ));
            }

            for(uint Ii=0; Ii < mCellsOnSet( Ik ).size(); Ii ++)
            {
                mAllCells( mCellsOnSet( Ik )( Ii )->get_index() ) = mCellsOnSet( Ik )( Ii );
            }
        }
    }

    // ----------------------------------------------------------------------------

    ~Visualization_Mesh()
    {
        for( auto tSet : mListOfAllSets )
        {
            delete tSet;
        }
        for( auto tCellOnSet : mCellsOnSet )
        {
            for( auto tCell : tCellOnSet )
            {
                delete tCell;
            }
        }
        for( auto tVerticesOnSet : mVerticesOnSet )
        {
            for( auto tVertex : tVerticesOnSet )
            {
                delete tVertex;
            }
        }
        for( auto tClustersOnSet : mClustersOnBlock )
        {
            for( auto tCluster : tClustersOnSet )
            {
                delete tCluster;
            }
        }
    };
    //##############################################
    // Cell Cluster Access
    //##############################################

    // ----------------------------------------------------------------------------

    moris::Cell<std::string> get_set_names(enum EntityRank aSetEntityRank ) const
    {
        moris::Cell<std::string> tSetNames;
        if (aSetEntityRank == EntityRank::ELEMENT)
        {
            uint tNumBlocks = this->get_num_sets();

            tSetNames.resize( tNumBlocks );

            for(uint Ik=0; Ik<tNumBlocks; Ik ++)
            {
                moris::mtk::Set * tSet = this->get_set_by_index( Ik);

                tSetNames( Ik ) = tSet->get_set_name();
            }
        }

        return tSetNames;
    }

    /**
     * Get the number of nodes on the mesh.
     *
     * @return Number of nodes
     */
    uint get_num_nodes() const
    {
        return mAllVertices.size();
    }

    // ----------------------------------------------------------------------------

    uint get_num_elems() const
    {
        uint tNumBlocks = this->get_num_blocks();

        uint tNumElems = 0;

        for(uint Ik=0; Ik<tNumBlocks; Ik ++)
        {
            moris::mtk::Set * tSet = this->get_set_by_index( Ik);

            tNumElems += tSet->get_num_cells_on_set( mOnlyPrimary );
        }
        return tNumElems;
    }

    // ----------------------------------------------------------------------------

    /*
     * Returns the mtk cells in a block set.
     */
    moris::Cell< mtk::Cell const * > get_set_cells( std::string aSetName ) const
    {
        uint tBlockIndex = mSetNameToIndexMap.find( aSetName );

        moris::mtk::Set * tSet = this->get_set_by_index( tBlockIndex );

        moris::Cell< const mtk::Cell * > tCellsOnSet( tSet->get_num_cells_on_set( mOnlyPrimary ) );

        for( uint Ik=0; Ik < mCellsOnSet( tBlockIndex ).size(); Ik ++ )
        {
            tCellsOnSet( Ik ) = mCellsOnSet( tBlockIndex )( Ik );
        }

        return tCellsOnSet;
    }

    // ----------------------------------------------------------------------------

    moris::Cell< moris::mtk::Vertex const * > get_all_vertices() const
    {
        return mAllVertices;
    }

    // ----------------------------------------------------------------------------
    /*
     * Get number of blocks
     */
    moris::uint get_num_blocks() const
    {
        return mListofBlocks.size();
    };

// ----------------------------------------------------------------------------

    moris::uint get_num_sets() const
    {
        return mListOfAllSets.size();
    }

// ----------------------------------------------------------------------------

    /*
     * Get set by index
     */
    moris::mtk::Set * get_set_by_index( moris::uint aSetIndex) const
    {
        MORIS_ASSERT( (uint)aSetIndex < mListOfAllSets.size(),"Set index out of bounds" );
        return mListOfAllSets( aSetIndex );
    };

    // ----------------------------------------------------------------------------

    moris::mtk::Set * get_set_by_name( std::string aSetLabel ) const
    {
        moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );

        MORIS_ASSERT( (uint)tSetIndex < mListOfAllSets.size(),"Set index out of bounds" );

        return mListOfAllSets( tSetIndex );
    }

// ----------------------------------------------------------------------------

    moris_index get_set_index_by_name( std::string aSetLabel ) const
    {
        return  mSetNameToIndexMap.find( aSetLabel );
    }

    // ----------------------------------------------------------------------------
    /*
     * Get number of side sets
     * Sometimes num side set
     */
    moris::uint get_num_side_set() const
    {
        return mListofSideSets.size();
    };

    /**
     * Gets element indices in a block set.
     *
     * @param aSetIndex Block set index
     * @return Element indices in the set
     */
    Matrix<IndexMat> get_element_indices_in_block_set(uint aSetIndex)
    {
        return mListofBlocks(aSetIndex)->get_cell_inds_on_block(false);
    }

    MeshType get_mesh_type() const
    {
        return MeshType::VIS;
    };

    Matrix< IdMat > get_communication_table() const
    {
        MORIS_ERROR( false, "get_communication_table(), not implemented for visualization mesh" );
        return Matrix<IdMat>( 0, 0 );
    }

    uint get_spatial_dim() const
    {
       return this->get_set_by_index( 0 )->get_spatial_dim();
    }

    uint get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const
    {
        MORIS_ERROR( false, "get_num_entities(), not implemented for visualization mesh" );
        return 0;
    }

    Matrix< IndexMat >
    get_nodes_connected_to_element_loc_inds(moris_index aElementIndex) const
    {
        return mAllCells(aElementIndex)->get_vertex_inds();
    }

    Matrix<IndexMat> get_entity_connected_to_entity_loc_inds (moris_index     aEntityIndex,
                                                             enum EntityRank aInputEntityRank,
                                                             enum EntityRank aOutputEntityRank,
															 const moris_index aIndex = 0) const
    {
        MORIS_ERROR( false, "get_entity_connected_to_entity_loc_inds(), not implemented for visualization mesh" );
        return Matrix<IndexMat>( 0, 0 );
    }

    Matrix< IndexMat > get_elements_connected_to_element_and_face_ord_loc_inds(moris_index aElementIndex) const
    {
        MORIS_ERROR( false, "get_entity_connected_to_entity_loc_inds(), not implemented for visualization mesh" );
        return Matrix<IndexMat>( 0, 0 );
    }

    Matrix< IndexMat > get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
    {
         MORIS_ERROR( false, "get_entity_connected_to_entity_loc_inds(), not implemented for visualization mesh" );
        return Matrix<IndexMat>( 0, 0 );
    }

    /**
     * Get the spatial coordinates of a node.
     *
     * @param aNodeIndex Node index
     * @return Node coordinates
     */
    Matrix< DDRMat >
    get_node_coordinate( moris_index aNodeIndex ) const
    {
        return mAllVertices(aNodeIndex)->get_coords();
    }


    /**
     * Get a global entity ID from an entity rank and local index.
     *
     * @param aEntityIndex Local entity index
     * @param aEntityRank Entity rank
     * @param aBSplineMeshIndex B-spline mesh Index
     * @return Global entity ID
     */
    moris_id
    get_glb_entity_id_from_entity_loc_index(
            moris_index        aEntityIndex,
            enum EntityRank    aEntityRank,
            const moris_index  aBSplineMeshIndex = 0) const
    {
        switch (aEntityRank)
        {
            case EntityRank::NODE:
                return mAllVertices(aEntityIndex)->get_id();
            case EntityRank::ELEMENT:
                return mAllCells(aEntityIndex)->get_id();
            default:
                MORIS_ERROR(false, "VIS mesh get_glb_entity_id_from_entity_loc_index() needs to be implemented.");
                return 0;
        }
    }

    /**
     * Gets the owner of a node.
     *
     * @param aNodeIndex Node index
     * @return Node owner
     */
    uint get_node_owner(moris_index aNodeIndex) const
    {
        return par_rank();
    }

    /**
     * Gets the owner of an element.
     *
     * @param aElementIndex Element index
     * @return Element owner
     */
    uint get_element_owner(moris_index aElementIndex) const
    {
        return par_rank();
    }

    enum CellTopology
    get_blockset_topology(const std::string & aSetName)
    {
        return mListofBlocks(mSetNameToIndexMap[aSetName])->get_cell_topology();
    }

    enum CellShape
    get_IG_blockset_shape(const std::string & aSetName)
    {
        return mListofBlocks(mSetNameToIndexMap[aSetName])->get_IG_cell_shape();
    }

    enum CellShape
    get_IP_blockset_shape(const std::string & aSetName)
    {
        return mListofBlocks(mSetNameToIndexMap[aSetName])->get_IP_cell_shape();
    }

protected:

    void collect_all_sets()
    {
        mListOfAllSets.append( mListofBlocks );
        mListOfAllSets.append( mListofSideSets );
//        mListOfAllSets.append( mListofDoubleSideSets );

        for( uint Ik = 0; Ik < mListOfAllSets.size(); Ik++ )
        {
           mListOfAllSets( Ik )->set_set_index( Ik );

           mSetNameToIndexMap[ mListOfAllSets( Ik )->get_set_name() ] = Ik;
        }
    }

};
}
}


#endif /* PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_ */
