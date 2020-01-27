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
    moris::Cell< moris::mtk::Set * > mListofSideSets;

    moris::Cell< moris::Cell< mtk::Cell * > >   mCellsOnSet;
    moris::Cell< moris::Cell< mtk::Vertex * > >  mVerticesOnSet;

    const bool                                  mOnlyPrimary;

    moris::Cell< moris::mtk::Set * > mListOfAllSets;

    map< std::string, moris_index >  mSetNameToIndexMap;

public:

    Visualization_Mesh( moris::Cell< moris::mtk::Set * >            aListofBlocks,
                        moris::Cell< moris::Cell< mtk::Cell * > >   aCellsOnBlock,
                        moris::Cell< moris::Cell< mtk::Vertex * > > aVerticesOnBlock,
                        const bool                                  aOnlyPrimary ) : mListofBlocks( aListofBlocks ),
                                                                                     mCellsOnSet( aCellsOnBlock ),
                                                                                     mVerticesOnSet( aVerticesOnBlock ),
                                                                                     mOnlyPrimary( aOnlyPrimary )
    {
        this->collect_all_sets();
    }

    // ----------------------------------------------------------------------------

    ~Visualization_Mesh()
    {
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
            moris::mtk::Set * tSet = this->get_set_by_index( Ik);

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
        uint tNumBlocks = this->get_num_blocks();

        moris::Cell< moris::mtk::Vertex const * > tVertices( this->get_num_nodes() );

        uint tCounter = 0;
        for(uint Ik=0; Ik < tNumBlocks; Ik ++)
        {

            for(uint Ii=0; Ii < mVerticesOnSet( Ik ).size(); Ii ++)
            {
                tVertices( tCounter++ ) = mVerticesOnSet( Ik )( Ii ) ;
            }
        }

        return tVertices;
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

    moris::mtk::Set * get_set_by_name( std::string aSetLabel )
    {
        moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );

        MORIS_ASSERT( (uint)tSetIndex < mListOfAllSets.size(),"Set index out of bounds" );

        return mListOfAllSets( tSetIndex );
    }

// ----------------------------------------------------------------------------

    moris_index get_set_index_by_name( std::string aSetLabel )
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
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the label
//     */
//    virtual
//    std::string
//    get_side_set_label(moris_index aSideSetOrdinal) const = 0;

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
