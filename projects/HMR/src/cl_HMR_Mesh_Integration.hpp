/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh_Integration.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_

#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_HMR_Cell_Cluster.hpp"
#include "cl_HMR_Side_Cluster.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Cell_Info.hpp"
namespace moris
{
    namespace hmr
    {

        class Integration_Mesh_HMR : public Mesh, public mtk::Integration_Mesh
        {
            private:
                mtk::Cell_Cluster * mDummyCluster  = nullptr;

                // cell clusters
                moris::Cell<Cell_Cluster_HMR> mCellClusters;

                // Block sets containing Cell Clusters
                moris::Cell<std::string>                     mPrimaryBlockSetNames;
                moris::Cell<moris::Cell<moris::moris_index>> mPrimaryBlockSetClusters;

                // side sets
                std::unordered_map<std::string, moris_index> mSideSideSetLabelToOrd;
                moris::Cell<std::string>                     mSideSetLabels;
                moris::Cell<moris::Cell<Side_Cluster_HMR>>   mSideSets;

            public:

                Integration_Mesh_HMR(
                        const uint                  & aLagrangeOrder,
                        const uint                  & aLagrangePattern,
                        Interpolation_Mesh_HMR      * aInterpolationMesh )
            : Mesh( aInterpolationMesh->get_database(),
                    aLagrangeOrder,
                    aLagrangePattern )
            {
                    this->setup_cell_clusters( aInterpolationMesh );
                    this->setup_blockset_with_cell_clusters();
                    this->setup_side_set_clusters( aInterpolationMesh );
                    this->collect_all_sets();
            }

                Integration_Mesh_HMR(
                        const uint                               & aLagrangeMeshIndex,
                        Interpolation_Mesh_HMR                   * aInterpolationMesh  )
                : Mesh( aInterpolationMesh->get_database(),
                        aLagrangeMeshIndex )
                {
                    this->setup_cell_clusters( aInterpolationMesh );
                    this->setup_blockset_with_cell_clusters();
                    this->setup_side_set_clusters( aInterpolationMesh );
                    this->collect_all_sets();
                }

                mtk::Cell_Cluster const &
                get_cell_cluster(mtk::Cell const & aInterpCell) const
                {
                    return mCellClusters(aInterpCell.get_index());
                }

                /*
                 * Get a cell cluster related to an interpolation
                 * cell
                 */
                mtk::Cell_Cluster const &
                get_cell_cluster(moris_index aInterpCellIndex) const
                {
                    MORIS_ASSERT(aInterpCellIndex<(moris_index)mCellClusters.size(),"Interpolation Cell index out of bounds");
                    return mCellClusters(aInterpCellIndex);
                }

                /*
                 * Get block set names
                 */
                moris::Cell<std::string>
                get_block_set_names() const
                {
                    return mPrimaryBlockSetNames;
                }

                /*
                 * Get cell clusters within a block set
                 */
                moris::Cell<mtk::Cluster const *>
                get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
                {
                    MORIS_ASSERT(aBlockSetOrdinal < (moris_index)mPrimaryBlockSetNames.size(),"Requested block set ordinal out of bounds.");

                    moris::Cell<moris::moris_index> const & tClusterIndsInSet = mPrimaryBlockSetClusters(aBlockSetOrdinal);

                    moris::Cell<mtk::Cluster const *> tClusterInSet(tClusterIndsInSet.size());

                    for(moris::uint i = 0; i <tClusterIndsInSet.size(); i++)
                    {
                        tClusterInSet(i) = &this->get_cell_cluster(tClusterIndsInSet(i));
                    }

                    return tClusterInSet;
                }

                /*!
                 * get number of side sets
                 */
                uint get_num_side_sets() const
                {
                    return mSideSets.size();
                }

                moris::Cell<mtk::Cluster const *>
                get_side_set_cluster(moris_index aSideSetOrdinal) const
                {
                    MORIS_ASSERT(aSideSetOrdinal < (moris_index)mSideSets.size(), "Side set ordinal out of bounds");

                    moris::uint tNumSideClustersInSet = mSideSets(aSideSetOrdinal).size();

                    moris::Cell<mtk::Cluster const *> tSideClustersInSet(tNumSideClustersInSet);

                    for(moris::uint i = 0; i <tNumSideClustersInSet; i++)
                    {
                        tSideClustersInSet(i) = & mSideSets(aSideSetOrdinal)(i);
                    }

                    return tSideClustersInSet;
                }

                /*!
                 * Returns the label
                 */

                std::string
                get_side_set_label(moris_index aSideSetOrdinal) const
                {
                    MORIS_ASSERT(aSideSetOrdinal<(moris_index)mSideSetLabels.size(),
                            "Side set ordinal out of bounds");

                    return mSideSetLabels(aSideSetOrdinal);
                }

                /*!
                 * Returns the index given a label
                 */
                moris_index
                get_side_set_index(std::string aSideSetLabel) const
                {
                    auto tIter = mSideSideSetLabelToOrd.find(aSideSetLabel);

                    MORIS_ERROR(tIter != mSideSideSetLabelToOrd.end(),
                            "side side set label not found");

                    return tIter->second;
                }

                uint
                get_num_double_sided_sets() const
                {
                    return 0;
                }

                /*!
                 * Returns the label
                 */

                std::string
                get_double_sided_set_label(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR(0,"get_double_sided_set_label not implemented in HMR Integration mesh");
                    return "ERROR";
                }

                /*!
                 * Returns the double side clusters in the side set
                 */

                moris::Cell<moris::mtk::Cluster const *>
                get_double_side_set_cluster(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR(0,"get_double_side_set_cluster not implemented in HMR Integration mesh");
                    return moris::Cell<moris::mtk::Cluster const *>(0);
                }

                //-------------------------------------------------------------------------------
                /*
                 * Construct HMR Cell Clustering
                 */
                void
                setup_cell_clusters(Interpolation_Mesh_HMR * aInterpolationMesh)
                {
                    // check to see the meshes are the same (since all trivial)
                    MORIS_ASSERT(this->get_num_nodes() == aInterpolationMesh->get_num_nodes(),
                            "Mismatch nodes between integration and interpolation mesh nodes. %-5i | %-5i",
                            this->get_num_nodes(),
                            aInterpolationMesh->get_num_nodes());

                    MORIS_ASSERT(this->get_num_elemens_including_aura() == aInterpolationMesh->get_num_elemens_including_aura(),
                            "Mismatch elements between integration and interpolation mesh");

                    // get the cell rank
                    enum EntityRank tCellRank = this->get_cell_rank();

                    // number of interpolation cells
                    moris::uint tNumInterpCells = aInterpolationMesh->get_num_entities(tCellRank);
                    //        moris::uint tNumInterpCells = aInterpolationMesh.get_num_elemens_including_aura();

                    // size member data
                    mCellClusters.resize( tNumInterpCells );

                    moris::mtk::Cell_Info_Factory tCIFactory;

                    for(moris::uint i = 0; i < tNumInterpCells; i++)
                    {
                        // interpolation cell
                        mtk::Cell const * tInterpCell = &aInterpolationMesh->get_mtk_cell( (moris_index) i );
                        mCellClusters(i).set_interpolation_cell( tInterpCell );

                        // integration cell (only primary cells here)
                        //            moris_index tIntegCellIndex    = this->get_loc_entity_ind_from_entity_glb_id(tCellId,tCellRank);
                        mtk::Cell const * tPrimaryCell = &this->get_mtk_cell( i );
                        mCellClusters(i).add_primary_integration_cell( tPrimaryCell );

                        moris::Cell<moris::mtk::Vertex *> tVertexIds  = tPrimaryCell->get_vertex_pointers();
                        moris::Cell<moris::mtk::Vertex const *> tConstVertexPtrs(tVertexIds.size());
                        for(moris::uint iV = 0 ; iV < tVertexIds.size();iV ++)
                        {
                            tConstVertexPtrs(iV) = tVertexIds(iV);
                        }

                        // interpolation order of the interpolation cell
                        enum mtk::Interpolation_Order tInterpOrder = tInterpCell->get_interpolation_order();

                        // geometry
                        enum mtk::Geometry_Type tGeomType = tInterpCell->get_geometry_type();

                        // Cell info
                        moris::mtk::Cell_Info* tCellInfo = tCIFactory.create_cell_info(tGeomType, tInterpOrder);

                        Matrix<DDRMat> tXi;
                        tCellInfo->get_loc_coords_of_cell(tXi);

                        mCellClusters(i).add_vertex_to_cluster(tConstVertexPtrs);

                        mCellClusters(i).add_vertex_local_coordinates_wrt_interp_cell(tXi);

                        delete tCellInfo;

                    }
                }

                //-------------------------------------------------------------------------------

                void
                setup_blockset_with_cell_clusters()
                {
                    // construct integration to cell cluster relationship
                    moris::Cell<moris::moris_index> tPrimaryIntegrationCellToClusterIndex(
                            this->get_num_entities(EntityRank::ELEMENT),
                            MORIS_INDEX_MAX );

                    // iterate through clusters
                    for(moris::uint  i = 0; i < mCellClusters.size(); i++)
                    {
                        Cell_Cluster_HMR const & tCellCluster = mCellClusters(i);
                        moris::Cell<moris::mtk::Cell const *> const & tPrimaryCells = tCellCluster.get_primary_cells_in_cluster();

                        // iterate through primary cells
                        for(moris::uint j = 0; j <tCellCluster.get_num_primary_cells(); j++)
                        {
                            moris::moris_index tCellIndex = tPrimaryCells( j )->get_index();

                            MORIS_ASSERT( tPrimaryIntegrationCellToClusterIndex( tCellIndex ) == MORIS_INDEX_MAX,
                                    "Integration cell can only appear as a primary cell in one cell cluster" );

                            tPrimaryIntegrationCellToClusterIndex(tCellIndex) = (moris_index) i;
                        }
                    }

                    // get all block sets from the mesh
                    moris::Cell<std::string> tBlockSetNames = this->get_set_names(this->get_cell_rank());

                    mPrimaryBlockSetClusters.resize(tBlockSetNames.size());
                    mPrimaryBlockSetNames = tBlockSetNames;

                    moris::Cell<uint> tSetsToRemove;

                    for(moris::uint i = 0; i<tBlockSetNames.size(); i++)
                    {
                        moris::Matrix<moris::IndexMat> tCellsInSet =
                                this->get_set_entity_loc_inds(this->get_cell_rank(),tBlockSetNames(i));

                        bool tSetHasCluster = false;

                        // reserve memory
                        mPrimaryBlockSetClusters(i).reserve( tCellsInSet.numel() );

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

                    // trim outer and inner cells
                    shrink_to_fit_all( mPrimaryBlockSetClusters );

                    mListofBlocks.resize( mPrimaryBlockSetClusters.size(), nullptr );

                    for(moris::uint Ik = 0; Ik<mListofBlocks.size(); Ik++)
                    {
                        mListofBlocks( Ik ) = new moris::mtk::Block(
                                tBlockSetNames(Ik),
                                this->get_cell_clusters_in_set( Ik ),{{0}}, this->get_spatial_dim());
                    }
                }

                //-------------------------------------------------------------------------------

                /*
                 *  setup the side set cluster interface
                 */
                void setup_side_set_clusters( Interpolation_Mesh_HMR * aInterpMesh )
                {
                    moris::Cell<std::string> aSideSetNames = this->get_set_names( EntityRank::FACE );

                    mSideSets.resize( aSideSetNames.size() );

                    // copy strings labels
                    mSideSetLabels.append( aSideSetNames );

                    // cell info to use for setting up local coords
                    enum CellTopology tCellTopo = this->get_blockset_topology( "" );
                    mtk::Cell_Info_Factory tFactory;
                    std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = tFactory.create_cell_info_sp(tCellTopo);

                    // add to map
                    for(moris::uint i = 0; i < aSideSetNames.size(); i++)
                    {
                        MORIS_ASSERT(mSideSideSetLabelToOrd.find(mSideSetLabels(i)) == mSideSideSetLabelToOrd.end(),
                                "Duplicate side set label detected.");

                        mSideSideSetLabelToOrd[mSideSetLabels(i)] = i;
                    }

                    // iterate through block sets
                    for( moris::uint i = 0;  i < aSideSetNames.size(); i++ )
                    {
                        // get the cells and side ordinals from the mesh for this side set
                        moris::Cell< mtk::Cell const * > tCellsInSet( 0 );
                        moris::Matrix<moris::IndexMat>   tSideOrdsInSet( 0,0 );
                        this->get_sideset_cells_and_ords( aSideSetNames( i ), tCellsInSet, tSideOrdsInSet );

                        // reserve memory for side set information
                        mSideSets(i).reserve( tCellsInSet.size() );

                        // figure out which integration cells are in the side cluster input. these are assumed
                        // the only non-trivial ones, all others will be marked as trivial
                        // all trivial case
                        // loop over cells in the side set and make sure they have all been included
                        for( moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++ )
                        {
                            moris_id tCellId = tCellsInSet( iIGCell )->get_id();

                            // interpolation cell index
                            moris_index tCellIndex = aInterpMesh->get_loc_entity_ind_from_entity_glb_id( tCellId,
                                    EntityRank::ELEMENT );

                            // construct a trivial side cluster
                            moris::mtk::Cell* tInterpCell = &aInterpMesh->get_mtk_cell( tCellIndex );

                            // get local coordinates on the side
                            Matrix<DDRMat> tXi;
                            tCellInfo->get_loc_coord_on_side_ordinal(tSideOrdsInSet(iIGCell), tXi);

                            mSideSets(i).push_back(Side_Cluster_HMR(true,
                                                        tInterpCell,
                                                        {tCellsInSet(iIGCell)},
                                                        {{tSideOrdsInSet(iIGCell)}},
                                                        tCellsInSet(iIGCell)->get_vertices_on_side_ordinal(tSideOrdsInSet(iIGCell)),
                                                        tXi));
                        }
                    }

                    mListofSideSets.resize( mSideSets.size(), nullptr );

                    for(moris::uint Ik = 0; Ik<mListofSideSets.size(); Ik++)
                    {
                        mListofSideSets( Ik ) = new moris::mtk::Side_Set(aSideSetNames(Ik), this->get_side_set_cluster( Ik ),{{0}}, this->get_spatial_dim());
                    }
                }

                //-------------------------------------------------------------------------------

                enum EntityRank
                get_cell_rank()
                {
                    if(this->get_spatial_dim() == 1)
                    {
                        return EntityRank::FACE;
                    }
                    else if(this->get_spatial_dim() == 2)
                    {
                        return EntityRank::ELEMENT;
                    }
                    else if (this->get_spatial_dim() == 3)
                    {
                        return EntityRank::ELEMENT;
                    }
                    else
                    {
                        return EntityRank::INVALID;
                    }
                }
        };
    }
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_ */

