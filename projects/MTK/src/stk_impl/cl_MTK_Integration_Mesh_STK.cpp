/*
 * cl_MTK_Integration_Mesh_STK.cpp
 *
 *  Created on: Apr 29, 2019
 *      Author: doble
 */


#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Cell_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"
namespace moris
{
namespace mtk
{
// ----------------------------------------------------------------------------

Integration_Mesh_STK::Integration_Mesh_STK(std::shared_ptr<Mesh_Data_STK> aSTKMeshData):
    Mesh_Core_STK(aSTKMeshData)
{
}

// ----------------------------------------------------------------------------

Integration_Mesh_STK::Integration_Mesh_STK(
        std::string    aFileName,
        MtkMeshData*   aSuppMeshData,
        const bool     aCreateFacesAndEdges ):
            Mesh_Core_STK(aFileName,aSuppMeshData,aCreateFacesAndEdges)
{

}

// ----------------------------------------------------------------------------

Integration_Mesh_STK::Integration_Mesh_STK( MtkMeshData & aMeshData ):
                        Mesh_Core_STK(aMeshData){}

// ----------------------------------------------------------------------------

Integration_Mesh_STK::Integration_Mesh_STK( MtkMeshData &       aMeshData,
                                            Interpolation_Mesh* aInterpMesh,
                                            Cell_Cluster_Input* aCellClusterData,
                                            Side_Cluster_Input* aSideClusterData):
        Mesh_Core_STK(aMeshData)
{
    // setup cells and cell clusters
    this->setup_cell_clusters(*aInterpMesh,aCellClusterData);
    this->setup_blockset_with_cell_clusters();

    // setup side set clusters
    this->setup_side_set_clusters(*aInterpMesh,aSideClusterData);

}

// ----------------------------------------------------------------------------

Integration_Mesh_STK::Integration_Mesh_STK(Interpolation_Mesh & aInterpMesh,
                                           Cell_Cluster_Input * aCellClusterInput)
{
    MORIS_ERROR(aInterpMesh.get_mesh_type() == MeshType::STK,"operator= between an interpolation and integration mesh only valid between stk meshes");

    Interpolation_Mesh_STK* tInterpolationSTK = dynamic_cast<Interpolation_Mesh_STK*>(&aInterpMesh);

    // get the shared data from the stk interpolation mesh
    mSTKMeshData = tInterpolationSTK->get_stk_data_shared_pointer();

    this->setup_cell_clusters(aInterpMesh,aCellClusterInput);

    this->setup_blockset_with_cell_clusters();


}

// ----------------------------------------------------------------------------

Cell_Cluster const &
Integration_Mesh_STK::get_cell_cluster(Cell const & aInterpCell) const
{
    return mCellClusters(aInterpCell.get_index());
}

// ----------------------------------------------------------------------------

Cell_Cluster const &
Integration_Mesh_STK::get_cell_cluster(moris_index aInterpCellIndex) const
{
   MORIS_ASSERT(aInterpCellIndex<(moris_index)mCellClusters.size(),"Interpolation Cell index out of bounds");
   return mCellClusters(aInterpCellIndex);
}

// ----------------------------------------------------------------------------

moris::Cell<Cell_Cluster const *>
Integration_Mesh_STK::get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
{
    MORIS_ASSERT(aBlockSetOrdinal<(moris_index)mPrimaryBlockSetNames.size(),"Requested block set ordinal out of bounds.");

    moris::Cell<moris::moris_index> const & tClusterIndsInSet = mPrimaryBlockSetClusters(aBlockSetOrdinal);

    moris::Cell<Cell_Cluster const *> tClusterInSet(tClusterIndsInSet.size());

    for(moris::uint i = 0; i <tClusterIndsInSet.size(); i++)
    {
        tClusterInSet(i) = &this->get_cell_cluster(tClusterIndsInSet(i));
    }

    return tClusterInSet;

}

moris::Cell<Side_Cluster const *>
Integration_Mesh_STK::get_side_set_cluster(moris_index aSideSetOrdinal) const
{
    MORIS_ASSERT(aSideSetOrdinal < (moris_index)mSideSets.size(), "Side set ordinal out of bounds");

    moris::uint tNumSideClustersInSet = mSideSets(aSideSetOrdinal).size();

    moris::Cell<Side_Cluster const *> tSideClustersInSet(tNumSideClustersInSet);

    for(moris::uint i = 0; i <tNumSideClustersInSet; i++)
    {
        tSideClustersInSet(i) = & mSideSets(aSideSetOrdinal)(i);
    }

    return tSideClustersInSet;
}

// ----------------------------------------------------------------------------

void
Integration_Mesh_STK::setup_cell_clusters(Interpolation_Mesh & aInterpMesh,
                                          Cell_Cluster_Input * aCellClusterInput)
{
    if( aCellClusterInput != nullptr)
    {
        // Number of interpolation cells
        moris::uint tNumInterpCells = aInterpMesh.get_num_entities(EntityRank::ELEMENT);

        mCellClusters.resize(tNumInterpCells);


        //fixme: An assumption is made here that interp mesh is hex based
        moris::Matrix<moris::DDRMat> tHexLocalCoords = {{-1, -1, -1},{ 1, -1, -1},{ 1, 1, -1},{-1, 1, -1},{-1, -1, 1},{ 1, -1, 1},{ 1, 1, 1},{-1, 1, 1}};

        // iterate through cells
        for(moris::uint i = 0; i <tNumInterpCells; i++)
        {
            moris_id tCellId = aInterpMesh.get_glb_entity_id_from_entity_loc_index((moris_index)i,EntityRank::ELEMENT);

            moris_index tClusterIndex = aCellClusterInput->get_cluster_index(tCellId);

            if(tClusterIndex != MORIS_INDEX_MAX)
            {
                mCellClusters(i).set_interpolation_cell( aCellClusterInput->get_interp_cell(tClusterIndex) );

                moris::Matrix<IndexMat>* tPrimaryCellIds  = aCellClusterInput->get_primary_cell_ids(tClusterIndex);
                mCellClusters(i).add_primary_integration_cell(this->get_cell_pointers_from_ids(*tPrimaryCellIds));

                moris::Matrix<IndexMat>* tVoidCellIds  = aCellClusterInput->get_void_cell_ids(tClusterIndex);
                mCellClusters(i).add_void_integration_cell(this->get_cell_pointers_from_ids(*tVoidCellIds));

                moris::Matrix<IndexMat>* tVertexIds  = aCellClusterInput->get_vertex_in_cluster_ids(tClusterIndex);
                mCellClusters(i).add_vertex_to_cluster(this->get_vertex_pointers_from_ids(*tVertexIds));

                moris::Matrix<DDRMat>* tVertexLocalCoords = aCellClusterInput->get_vertex_local_coords_wrt_interpolation_cell(tClusterIndex);
                mCellClusters(i).add_vertex_local_coordinates_wrt_interp_cell(*tVertexLocalCoords);

            }
            else
            {
                // an assumption is made here that if the cell does not have a cluster then it is related 1 to 1
                // with the interpolation cell with the same id

                // interpolation cell
                mtk::Cell const * tInterpCell = &aInterpMesh.get_mtk_cell((moris_index) i);
                mCellClusters(i).set_interpolation_cell( tInterpCell );

                // integration cell (only primary cells here)
                moris_index tIntegCellIndex    = this->get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
                mtk::Cell const * tPrimaryCell = &this->get_mtk_cell(tIntegCellIndex);
                mCellClusters(i).add_primary_integration_cell(tPrimaryCell);

                moris::Matrix<IndexMat> tVertexIds  = this->get_entity_connected_to_entity_glob_ids(tCellId,EntityRank::ELEMENT,EntityRank::NODE);
                mCellClusters(i).add_vertex_to_cluster(this->get_vertex_pointers_from_ids(tVertexIds));

                mCellClusters(i).add_vertex_local_coordinates_wrt_interp_cell(tHexLocalCoords);

            }
        }
    }

}

void
Integration_Mesh_STK::setup_blockset_with_cell_clusters( )
{
   // construct integration to cell cluster relationship
    moris::Cell<moris::moris_index> tPrimaryIntegrationCellToClusterIndex(this->get_num_entities(EntityRank::ELEMENT),MORIS_INDEX_MAX);

    // iterate through clusters
    for(moris::uint  i = 0; i < mCellClusters.size(); i++)
    {
        Cell_Cluster_STK const & tCellCluster = mCellClusters(i);
        moris::Cell<moris::mtk::Cell const *> const & tPrimaryCells = tCellCluster.get_primary_cells_in_cluster();

        // iterate through primary cells
        for(moris::uint j = 0; j <tCellCluster.get_num_primary_cells(); j++)
        {
            moris::moris_index tCellIndex = tPrimaryCells(j)->get_index();

            MORIS_ASSERT(tPrimaryIntegrationCellToClusterIndex(tCellIndex) == MORIS_INDEX_MAX,"Integration cell can only appear as a primary cell in one cell cluster");
            tPrimaryIntegrationCellToClusterIndex(tCellIndex) = (moris_index) i;
        }
    }

    // get all block sets from the mesh
    moris::Cell<std::string> tBlockSetNames = this->get_set_names(EntityRank::ELEMENT);

    mPrimaryBlockSetClusters.resize(tBlockSetNames.size());
    mPrimaryBlockSetNames = tBlockSetNames;

    moris::Cell<uint> tSetsToRemove;

    for(moris::uint i = 0; i<tBlockSetNames.size(); i++)
    {

        moris::Matrix<moris::IndexMat> tCellsInSet = this->get_set_entity_loc_inds(EntityRank::ELEMENT,tBlockSetNames(i));

        bool tSetHasCluster = false;

        // add the cells in set's primary cluster to member data (make unique later)
        for(moris::uint j =0; j<tCellsInSet.numel(); j++)
        {
            moris::moris_index tCellIndex    = tCellsInSet(j);
            moris::moris_index tClusterIndex = tPrimaryIntegrationCellToClusterIndex(tCellIndex);

            if(tClusterIndex != MORIS_INDEX_MAX )
            {
                tSetHasCluster = true;
                mPrimaryBlockSetClusters(i).push_back(tClusterIndex);
            }
        }

        // if there are no elements which are a primary cell in one cluster for this block set, then we remove it.
        if( !tSetHasCluster )
        {
            tSetsToRemove.push_back(i);
        }

        // remove duplicates
        else
        {
            moris::unique( mPrimaryBlockSetClusters(i));
        }
    }


    // remove block sets which had not primary clusters
    for(moris::uint i = tSetsToRemove.size(); i>0; i--)
    {
        mPrimaryBlockSetClusters.erase(tSetsToRemove(i-1));
        mPrimaryBlockSetNames.erase(tSetsToRemove(i-1));
    }
}

void
Integration_Mesh_STK::setup_side_set_clusters(Interpolation_Mesh & aInterpMesh,
                                              Side_Cluster_Input * aSideClusterInput)
{
    moris::Cell<std::string> aSideSetNames = this->get_set_names(EntityRank::FACE);

    mSideSets.resize(aSideSetNames.size());

    // iterate through block sets
    for(moris::uint i = 0;  i < aSideSetNames.size(); i++)
    {
        // get the cells and side ordinals from the mesh for this side set
        moris::Cell< mtk::Cell const * > tCellsInSet(0);
        moris::Matrix<moris::IndexMat>   tSideOrdsInSet(0,0);
        this->get_sideset_cells_and_ords(aSideSetNames(i), tCellsInSet, tSideOrdsInSet);

        // integration cell to side set index

        // figure out which integration cells are in the side cluster input. these are assumed
        // the only non-trivial ones, all others will be marked as trivial
        if(aSideClusterInput !=nullptr)
        {
        moris::moris_index tSideClusterOrd = aSideClusterInput->get_side_label_ordinal(aSideSetNames(i));

        if(tSideClusterOrd == MORIS_INDEX_MAX)
        {
            // loop over cells in the side set and make sure they have all been included
            for(moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++)
            {
                // integration cell id
                moris_id tCellId = tCellsInSet(iIGCell)->get_id();

                // interpolation cell index
                moris_index tCellIndex = aInterpMesh.get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);

                // construct a trivial side cluster
                moris::mtk::Cell* tInterpCell = &aInterpMesh.get_mtk_cell(tCellIndex);
                mSideSets(i).push_back(Side_Cluster_STK(tInterpCell));
            }
        }
        else
        {
            // figure out which integration cells are in this mesh

            // access side set cluster data
            Side_Set_Cluster_Data const & tSideSetClusterData = aSideClusterInput->get_cluster_data(tSideClusterOrd);

            // mark integration cells which are not trivial
            std::unordered_map<moris_index,bool> tIntegrationCellsInSideSet;

            // collect integration cells in clusters for this side set
            uint tNumClusters = tSideSetClusterData.get_num_cell_clusters();

            for(moris::uint iC = 0; iC<tNumClusters; iC++)
            {
                // get data from the side set cluster data

                // cell ids and side ords
                moris::Matrix<moris::IdMat> const * tCellIdsAndOrds = tSideSetClusterData.get_integration_cell_ids_and_side_ords(iC);
                moris::Cell<moris::mtk::Cell const *> tCellPointers = this->get_cell_pointers_from_ids(tCellIdsAndOrds->get_column(0));
                moris::Matrix<moris::IdMat> tSideOrds = tCellIdsAndOrds->get_column(1);

                // vertex pointers
                moris::Matrix<moris::IdMat> const * tVertexInCluster = tSideSetClusterData.get_vertex_in_cluster_ids(iC);
                moris::Cell<moris::mtk::Vertex const *> tVertices = this->get_vertex_pointers_from_ids(*tVertexInCluster);

                // vertices in cluster
                // vertex parametric coordinate relative to
                mSideSets(i).push_back(Side_Cluster_STK(false,
                                       tSideSetClusterData.get_interp_cell(iC),
                                       tCellPointers,
                                       tSideOrds,
                                       tVertices,
                                       *tSideSetClusterData.get_vertex_local_coords_wrt_interpolation_cell(iC)));

                // mark all integration cells in this cluster
                for(moris::uint  iIGCell = 0; iIGCell<tCellIdsAndOrds->n_rows(); iIGCell++)
                {
                    tIntegrationCellsInSideSet[(*tCellIdsAndOrds)(iIGCell,0)] = true;
                }
            }

            // loop over cells in the side set and make sure they have all been included
            for(moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++)
            {
                moris_id tCellId = tCellsInSet(iIGCell)->get_id();

                if(tIntegrationCellsInSideSet.find(tCellId) == tIntegrationCellsInSideSet.end())
                {
                    // interpolation cell index
                    moris_index tCellIndex = aInterpMesh.get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);

                    // construct a trivial side cluster
                    moris::mtk::Cell* tInterpCell = &aInterpMesh.get_mtk_cell(tCellIndex);
                    mSideSets(i).push_back(Side_Cluster_STK(tInterpCell));
                }
            }


        }
        }

    }




}


moris::Cell<moris::mtk::Cell const *>
Integration_Mesh_STK::get_cell_pointers_from_ids(moris::Matrix<moris::IdMat> const & aCellIds) const
{
    moris::Cell<moris::mtk::Cell const *> tCellPtrs(aCellIds.numel());

    for(moris::uint i = 0; i < aCellIds.numel(); i++)
    {
        moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id(aCellIds(i),EntityRank::ELEMENT);
        tCellPtrs(i) = &this->get_mtk_cell(tCellIndex);
    }

    return tCellPtrs;
}

    // ----------------------------------------------------------------------------

moris::Cell<moris::mtk::Vertex const *>
Integration_Mesh_STK::get_vertex_pointers_from_ids(moris::Matrix<moris::IdMat> const & aVertexIds) const
{
    moris::Cell<moris::mtk::Vertex const *> tVertexPtrs(aVertexIds.numel());

    for(moris::uint i = 0; i < aVertexIds.numel(); i++)
    {
        moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id(aVertexIds(i),EntityRank::NODE);
        tVertexPtrs(i) = &this->get_mtk_vertex(tCellIndex);
    }

    return tVertexPtrs;
}    // ----------------------------------------------------------------------------

}
}
