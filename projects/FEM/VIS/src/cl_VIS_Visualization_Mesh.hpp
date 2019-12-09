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
#include "cl_VSI_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"

//#include "cl_Matrix.hpp"
//#include "cl_MTK_Block.hpp"

namespace moris
{
namespace vis
{
class Visualization_Mesh : public mtk::Mesh
{
protected:
    moris::Cell< moris::mtk::Set * > mListofBlocks;

    moris::Cell< moris::Cell< mtk::Cell * > >   mCellsOnBlock;
    moris::Cell< moris::Cell< mtk::Vertex * > >  mVerticesOnBlock;

    moris::map< std::string, moris::moris_index > mBlockNameToIndexMap;
//    moris::Cell< moris::mtk::Set * > mListofSideSets;
//    moris::Cell< moris::mtk::Set * > mListofDoubleSideSets;

    const bool                                  mOnlyPrimary;

public:

    Visualization_Mesh( moris::Cell< moris::mtk::Set * >            aListofBlocks,
                        moris::Cell< moris::Cell< mtk::Cell * > >   aCellsOnBlock,
                        moris::Cell< moris::Cell< mtk::Vertex * > > aVerticesOnBlock,
						const bool                                  aOnlyPrimary ) : mListofBlocks( aListofBlocks ),
                                                                                     mCellsOnBlock( aCellsOnBlock ),
                                                                                     mVerticesOnBlock( aVerticesOnBlock ),
                                                                                     mOnlyPrimary( aOnlyPrimary )
    {
        this->create_block_name_list();
    }
	
	// ----------------------------------------------------------------------------
		
    ~Visualization_Mesh()
    {
    };
    //##############################################
    // Cell Cluster Access
    //##############################################

    // ----------------------------------------------------------------------------
	
    moris::Cell<std::string> get_set_names(enum EntityRank aSetEntityRank) const
    {
        moris::Cell<std::string> tSetNames;
        if (aSetEntityRank == EntityRank::ELEMENT)
        {
            uint tNumBlocks = this->get_num_blocks();

            tSetNames.resize( tNumBlocks );

            for(uint Ik=0; Ik<tNumBlocks; Ik ++)
            {
                moris::mtk::Set * tSet = this->get_block_by_index( Ik);

                tSetNames( Ik ) = tSet->get_set_name();
            }
        }

        return tSetNames;
    }
	
    // ----------------------------------------------------------------------------
		
    void create_block_name_list()
    {
        uint tNumBlocks = this->get_num_blocks();

        for(uint Ik=0; Ik<tNumBlocks; Ik ++)
        {
            moris::mtk::Set * tSet = this->get_block_by_index( Ik);

            mBlockNameToIndexMap[ tSet->get_set_name() ] = Ik;
        }
    }
	
	// ----------------------------------------------------------------------------
//    /*
//     * Get a cell cluster related to an interpolation
//     * cell
//     */
//    virtual
//    Cell_Cluster const &
//    get_cell_cluster(Cell const & aInterpCell) const = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*
//     * Get block set names
//     */
//    virtual
//    moris::Cell<std::string>
//    get_block_set_names() const = 0;
//
//    /*!
//     * Returns the label
//     */
//    virtual
//    std::string
//    get_block_set_label(moris_index aBlockSetOrdinal) const
//    {
//        MORIS_ERROR(0, "get_block_set_label has no default implementation");
//        return "ERROR";
//    }

    //------------------------------------------------------------------------------
    /*
     * Get number of nodes
     */
    uint get_num_nodes() const
    {
        uint tNumBlocks = this->get_num_blocks();

        uint tNumNodes = 0;

        for(uint Ik=0; Ik<tNumBlocks; Ik ++)
        {
            moris::mtk::Set * tSet = this->get_block_by_index( Ik);

            tNumNodes += tSet->get_num_vertieces_on_set( mOnlyPrimary );
        }
        return tNumNodes;
    }

    // ----------------------------------------------------------------------------
	
    uint get_num_elems() const
    {
        uint tNumBlocks = this->get_num_blocks();

        uint tNumElems = 0;

        for(uint Ik=0; Ik<tNumBlocks; Ik ++)
        {
            moris::mtk::Set * tSet = this->get_block_by_index( Ik);

            tNumElems += tSet->get_num_cells_on_set( mOnlyPrimary );
        }
        return tNumElems;
    }

    // ----------------------------------------------------------------------------

    /*
     * Returns the mtk cells in a block set.
     */
    moris::Cell< mtk::Cell const * > get_block_set_cells( std::string aSetName) const
    {
        uint tBlockIndex = mBlockNameToIndexMap.find( aSetName );

        moris::mtk::Set * tSet = this->get_block_by_index( tBlockIndex );

        moris::Cell< const mtk::Cell * > tCellsOnSet( tSet->get_num_cells_on_set( mOnlyPrimary ) );

        for( uint Ik=0; Ik < mCellsOnBlock( tBlockIndex ).size(); Ik ++ )
        {
            tCellsOnSet( Ik ) = mCellsOnBlock( tBlockIndex )( Ik );
        }

        return tCellsOnSet;
    }

    // ----------------------------------------------------------------------------

    moris::Cell< moris::mtk::Vertex const * > get_all_vertices() const
    {
        uint tNumBlocks = this->get_num_blocks();

        moris::Cell< moris::mtk::Vertex const * > tVertices( this->get_num_nodes() );

        uint tCounter = 0;
        for(uint Ik=0; Ik < tNumBlocks; Ik ++)
        {

            for(uint Ii=0; Ii < mVerticesOnBlock( Ik ).size(); Ii ++)
            {
                tVertices( tCounter++ ) = mVerticesOnBlock( Ik )( Ii ) ;
            }
        }

        return tVertices;
    }

    // ----------------------------------------------------------------------------

    /*!
     * Returns the index given a label
     */
    virtual
    moris_index
    get_block_set_index(std::string aBlockSetLabel) const
    {
        MORIS_ERROR(0, "get_block_set_index has no default implementation");
        return MORIS_INDEX_MAX;
    }

    // ----------------------------------------------------------------------------
    /*
     * Get number of blocks
     */
    moris::uint
    get_num_blocks() const
    {
        return mListofBlocks.size();
    };

    // ----------------------------------------------------------------------------
    /*
     * Get block by index
     */
    moris::mtk::Set * get_block_by_index( moris::uint aBlockIndex) const
    {
        MORIS_ASSERT(aBlockIndex<mListofBlocks.size(),"Block index out of bounds");
        return mListofBlocks( aBlockIndex);
    };
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get number of blocks
//     * Sometimes num side set * 2. Ask Keenan
//     */
//    moris::uint
//    get_num_side_set() const
//    {
//        return mListofSideSets.size();
//    };
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get block by index
//     */
//    moris::mtk::Set *
//    get_side_set_by_index( moris::uint aSideSetIndex) const
//    {
//        MORIS_ASSERT(aSideSetIndex<mListofSideSets.size(),"Side set index out of bounds");
//        return mListofSideSets(aSideSetIndex);
//    };
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get number of blocks
//     * Sometimes num side set * 2. Ask Keenan
//     */
//    moris::uint
//    get_num_double_side_set() const
//    {
//        return mListofDoubleSideSets.size();
//    };
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get block by index
//     */
//    moris::mtk::Set *
//    get_double_side_set_by_index( moris::uint aSideSetIndex) const
//    {
//        MORIS_ASSERT(aSideSetIndex<mListofDoubleSideSets.size(),"Double side set index out of bounds");
//        return mListofDoubleSideSets(aSideSetIndex);
//    };
//
//    // ----------------------------------------------------------------------------
//
//    /*
//     * Get cell clusters within a block set
//     */
//    virtual
//    moris::Cell<Cluster const *>
//    get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const = 0;
//
//    //##############################################
//    // Side Set Cluster Access
//    //##############################################
//    /*!
//     * Get side clusters within a side set
//     */
//    virtual
//    moris::Cell<Cluster const *>
//    get_side_set_cluster(moris_index aSideSetOrdinal) const = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * get number of side sets
//     */
//    virtual
//    uint
//    get_num_side_sets() const = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the label
//     */
//    virtual
//    std::string
//    get_side_set_label(moris_index aSideSetOrdinal) const = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the index given a label
//     */
//    virtual
//    moris_index
//    get_side_set_index(std::string aSideSetLabel) const  = 0;
//
//    // ----------------------------------------------------------------------------
//
//    //##############################################
//    // Double Side Set Cluster Access
//    //##############################################
//
//    /*!
//     * Returns the number of double sided side sets in the mesh
//     */
//    virtual
//    uint
//    get_num_double_sided_sets() const  = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the label
//     */
//    virtual
//    std::string
//    get_double_sided_set_label(moris_index aSideSetOrdinal) const = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the index given a label
//     */
//    virtual
//    moris_index
//    get_double_sided_set_index(std::string aDoubleSideSetLabel) const  = 0;
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the double side clusters in the side set
//     */
//    virtual
//    moris::Cell<Cluster const*>
//    get_double_side_set_cluster(moris_index aSideSetOrdinal) const  = 0;
//
//    // ----------------------------------------------------------------------------


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
       return this->get_block_by_index( 0 )->get_spatial_dim();
    }

    uint get_num_entities( enum EntityRank aEntityRank ) const
    {
        MORIS_ERROR( false, "get_num_entities(), not implemented for visualization mesh" );
        return 0;
    }

    Matrix<IndexMat> get_entity_connected_to_entity_loc_inds (moris_index     aEntityIndex,
                                                             enum EntityRank aInputEntityRank,
                                                             enum EntityRank aOutputEntityRank) const
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
//

};
}
}


#endif /* PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_ */
