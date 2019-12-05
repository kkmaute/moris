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
//    moris::Cell< moris::mtk::Set * > mListofSideSets;
//    moris::Cell< moris::mtk::Set * > mListofDoubleSideSets;

//    moris::Cell< moris::Cell< mtk::Cluster * > >        mClustersOnBlock;
//    moris::Cell< moris::Cell< const mtk::Cluster * > >  mClustersOnBlock_1;
//    moris::Cell< moris::Cell< mtk::Cell * > >   mCellsOnBlock;
//    moris::Cell< moris::Cell< mtk::Vertex * > > mVerticesOnBlock;
//
//    moris::Cell< Matrix< DDSMat > >             mVertexMapOnBlock;
//    moris::Cell< Matrix< DDSMat > >             mCellMapOnBlock;

public:
    Visualization_Mesh(){};
    // Functions only valid for visualization meshes

    Visualization_Mesh( moris::Cell< moris::mtk::Set * > aListofBlocks ) : mListofBlocks( aListofBlocks )
    {

    }

    Visualization_Mesh( mtk::Mesh_Manager* aMesh,
                        const uint         aMeshPairIndex )
    {
//        mListofBlocks.resize( 1, nullptr );
//        mClustersOnBlock.resize( 1 );
//        mCellsOnBlock.resize( 1 );
//        mVerticesOnBlock.resize( 1 );
//
//        mVertexMapOnBlock.resize( 1 );
//        mCellMapOnBlock.resize( 1 );
//
//        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
//        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
//        aMesh->get_mesh_pair( aMeshPairIndex, tInterpolationMesh, tIntegrationMesh );
//
//        moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( 0 );
//
//        uint tNumVerticesOnSet = tMeshSet->get_num_vertieces_on_set();
//
//        moris::Matrix< DDSMat > tVertexIndOnBlock = tMeshSet->get_vertieces_inds_on_block();
//        print(tVertexIndOnBlock,"tVertexIndOnBlock");
//
//        mVerticesOnBlock( 0 ).resize( tNumVerticesOnSet, nullptr );
//        mVertexMapOnBlock( 0 ).set_size( tVertexIndOnBlock.max() + 1, 1, -1 );
//
//        moris_index tVertexIndexCounter = 0;
//        moris_id tVertexIdCounter = 0;
//
//        for( uint Ik = 0; Ik < tNumVerticesOnSet; Ik++ )
//        {
//            mtk::Vertex tIntegrationVertex = tIntegrationMesh->get_mtk_vertex( tVertexIndOnBlock( Ik ) );
//
//            mVertexMapOnBlock( 0 )( tVertexIndOnBlock( Ik ) ) = Ik;
//
//            mVerticesOnBlock( 0 )( Ik ) = new Vertex_Visualization( tVertexIdCounter,
//                                                                    tVertexIndexCounter,
//                                                                    &tIntegrationVertex );
//
//            tVertexIndexCounter++;
//            tVertexIdCounter++;
//        }
//
//        std::cout<<"------------------------"<<std::endl;
////-----------------------------------------------------------------------------------------------------
//        uint tNumCellsOnSet = tMeshSet->get_num_cells_on_set();
//
//        moris::Matrix< DDSMat > tCellIndOnBlock = tMeshSet->get_cell_inds_on_block();
//
//        mCellsOnBlock( 0 ).resize( tNumCellsOnSet, nullptr );
//
//        mCellMapOnBlock( 0 ).set_size( tCellIndOnBlock.max() + 1, 1, -1 );
//
//        moris_index tCellIndexCounter = 0;
//        moris_id tCellIdCounter = 0;
//
//        std::cout<<tNumCellsOnSet<<" tNumCellsOnSet"<<std::endl;
//        std::cout<<mVertexMapOnBlock( 0 ).numel()<<" tIntegrationVertexInd"<<std::endl;
//
//        for( uint Ik = 0; Ik < tNumCellsOnSet; Ik++ )
//        {
////            mtk::Cell tIntegrationCell = tIntegrationMesh->get_mtk_cell( tCellIndOnBlock( Ik ) );
//
//            Matrix< IndexMat > tIntegrationVertexInd = tIntegrationMesh->get_mtk_cell( tCellIndOnBlock( Ik ) ).get_vertex_inds();
//
//            print(tIntegrationVertexInd , "tIntegrationVertexInd");
//            moris::Cell< mtk::Vertex* > tCellVerices( tIntegrationVertexInd.numel(), nullptr );
//
//
//
//            for( uint Ii = 0; Ii < tIntegrationVertexInd.numel(); Ii++ )
//            {
//                std::cout<<mVertexMapOnBlock( 0 )( tIntegrationVertexInd( Ii ) )<<std::endl;
//                tCellVerices( Ii ) = mVerticesOnBlock( 0 )( mVertexMapOnBlock( 0 )( tIntegrationVertexInd( Ii ) ) );
//            }
//
//            mCellMapOnBlock( 0 )( tCellIndOnBlock( Ik ) ) = Ik;
//
//            mCellsOnBlock( 0 )( Ik ) = new Cell_Visualization( tCellIdCounter,
//                                                               tCellIndexCounter,
//                                                               tCellVerices,
//                                                               &tIntegrationMesh->get_mtk_cell( tCellIndOnBlock( Ik ) ) );
//
//           tCellIndexCounter++;
//           tCellIdCounter++;
//        }
//
////-----------------------------------------------------------------------------------------------------
//        uint tNumClustersOnSet = tMeshSet->get_num_clusters_on_set();
//
//        moris::Cell< mtk::Cluster const * > tClustersOnBlock = tMeshSet->get_clusters_on_set();
//
//        mClustersOnBlock( 0 ).resize( tNumClustersOnSet, nullptr );
//
//
//        for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
//        {
//            moris::Cell<moris::mtk::Cell const *> tPrimaryCells = tClustersOnBlock( Ik )->get_primary_cells_in_cluster();
//
//            moris::Cell< mtk::Cell * > tClusterPrimaryCells( tPrimaryCells.size(), nullptr );
//
//            for( uint Ii = 0; Ii < tPrimaryCells.size(); Ii++ )
//            {
//                tClusterPrimaryCells( Ii ) = mCellsOnBlock( 0 )( mCellMapOnBlock( 0 )( tPrimaryCells( Ii )->get_index() ) );
//            }
//
//            moris::Cell< moris::mtk::Cell const * > tVoidCells = tClustersOnBlock( Ik )->get_void_cells_in_cluster();
//
//            moris::Cell< mtk::Cell * > tClusterVoidCells( tVoidCells.size(), nullptr );
//
//            for( uint Ii = 0; Ii < tVoidCells.size(); Ii++ )
//            {
//                tClusterVoidCells( Ii ) = mCellsOnBlock( 0 )( mCellMapOnBlock( 0 )( tVoidCells( Ii )->get_index() ) );
//            }
//
//            mClustersOnBlock( 0 )( Ik ) = new Cell_Cluster_Visualization;
//
//            if( tVoidCells.size() > 0 )
//            {
//                 mClustersOnBlock( 0 )( Ik )->mark_as_nontrivial();
//            }
//
//            mClustersOnBlock( 0 )( Ik )->add_primary_integration_cell( tPrimaryCells );
//            mClustersOnBlock( 0 )( Ik )->add_void_integration_cell( tVoidCells );
//            mClustersOnBlock( 0 )( Ik )->set_interpolation_cell( &tClustersOnBlock( Ik )->get_interpolation_cell() );
//        }
//
//        mClustersOnBlock_1.resize( 1 );
//        mClustersOnBlock_1( 0 ).resize( tNumClustersOnSet, nullptr );
//
//        for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
//        {
//            mClustersOnBlock_1( 0 )( Ik ) = mClustersOnBlock( 0 )( Ik );
//        }
//
//        mListofBlocks( 0 ) = new moris::mtk::Block( tMeshSet->get_set_name(),
//                                                    mClustersOnBlock_1( 0 ),
//                                                    tMeshSet->get_spatial_dim() );
    };

    ~Visualization_Mesh()
    {

    };
    //##############################################
    // Cell Cluster Access
    //##############################################

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
//
//    // ----------------------------------------------------------------------------
//
//    /*!
//     * Returns the index given a label
//     */
//    virtual
//    moris_index
//    get_block_set_index(std::string aBlockSetLabel) const
//    {
//        MORIS_ERROR(0, "get_block_set_index has no default implementation");
//        return MORIS_INDEX_MAX;
//    }
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get number of blocks
//     */
//    moris::uint
//    get_num_blocks() const
//    {
//        return mListofBlocks.size();
//    };
//
//    // ----------------------------------------------------------------------------
//    /*
//     * Get block by index
//     */
//    moris::mtk::Set *
//    get_block_by_index( moris::uint aBlockIndex) const
//    {
//        MORIS_ASSERT(aBlockIndex<mListofBlocks.size(),"Block index out of bounds");
//        return mListofBlocks(aBlockIndex);
//    };
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
        MORIS_ERROR( false, "get_mesh_type(), not implemented for visualization mesh" );
        return MeshType::END_ENUM;
    };

    Matrix< IdMat > get_communication_table() const
    {
        MORIS_ERROR( false, "get_communication_table(), not implemented for visualization mesh" );
        return Matrix<IdMat>( 0, 0 );
    }

    uint get_spatial_dim() const
    {
        MORIS_ERROR( false, "get_spatial_dim(), not implemented for visualization mesh" );
        return 0;
    }

    uint get_num_entities(enum EntityRank aEntityRank) const
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
