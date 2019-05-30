/*
 * cl_HMR_Mesh_Interpolation.hpp
 *
 *  Created on: Apr 19, 2019
 *      Author: doble
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_

#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_HMR_Cell_Cluster.hpp"
namespace moris
{
namespace hmr
{

class Integration_Mesh_HMR : public Mesh, public mtk::Integration_Mesh
{
public:
    Integration_Mesh_HMR(std::shared_ptr< Database > aDatabase,
                         const uint & aLagrangeOrder,
                         const uint & aLagrangePattern,
                         Interpolation_Mesh_HMR* aInterpolationMesh  ):
        Mesh( aDatabase, aLagrangeOrder, aLagrangePattern )
    {

    }




    mtk::Cell_Cluster const &
    get_cell_cluster(mtk::Cell const & aInterpCell) const
    {
        MORIS_ERROR(0,"Cell clusters not implemented in HMR");
        return *mDummyCluster;
    }

    /*
     * Get block set names
     */
    moris::Cell<std::string>
    get_block_set_names() const
    {
        MORIS_ERROR(0,"get_block_set_names not implemented in HMR");
        return moris::Cell<std::string>(0);
    }

    /*
     * Get cell clusters within a block set
     */
    moris::Cell<mtk::Cell_Cluster const *>
    get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
    {
        MORIS_ERROR(0,"get_block_set_names not implemented in HMR");
        moris::Cell<mtk::Cell_Cluster const *> tCellInCluster(0);
        return tCellInCluster;
    }

    /*!
     * get number of side sets
     */
    uint
    get_num_side_sets() const
    {
        MORIS_ERROR(0,"get_num_side_sets not implemented in HMR");
        return std::numeric_limits<moris_index>::max();
    }

    moris::Cell<mtk::Side_Cluster const *>
    get_side_set_cluster(moris_index aSideSetOrdinal) const
    {
        MORIS_ERROR(0,"get_side_set_cluster not implemented in HMR");
        moris::Cell<mtk::Side_Cluster const *> tSideClusters(0);
        return tSideClusters;
    }

    /*!
     * Returns the label
     */
    virtual
    std::string
    get_side_set_label(moris_index aSideSetOrdinal) const
    {
        MORIS_ERROR(0,"get_side_set_label not implemented in HMR Integration mesh");
        return "ERROR";
    }

    /*!
     * Returns the index given a label
     */
    virtual
    moris_index
    get_side_set_index(std::string aSideSetLabel) const
    {
        MORIS_ERROR(0,"get_side_set_index not implemented in HMR Integration mesh");
        return std::numeric_limits<moris_index>::max();
    }

    uint
    get_num_double_sided_sets() const
    {
        MORIS_ERROR(0,"get_num_double_sided_sets not implemented in HMR Integration mesh");
        return std::numeric_limits<uint>::max();
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
     * Returns the index given a label
     */

    moris_index
    get_double_sided_set_index(std::string aDoubleSideSetLabel) const
    {
        MORIS_ERROR(0,"get_double_sided_set_index not implemented in HMR Integration mesh");
        return MORIS_INDEX_MAX;
    }

    /*!
     * Returns the double side clusters in the side set
     */

    moris::Cell<moris::mtk::Double_Side_Cluster> const &
    get_double_side_set_cluster(moris_index aSideSetOrdinal) const
    {
        MORIS_ERROR(0,"get_double_side_set_cluster not implemented in HMR Integration mesh");
        return mDummyDoubleSide;
    }

private:
    mtk::Cell_Cluster * mDummyCluster     = nullptr;
    moris::Cell<moris::mtk::Double_Side_Cluster> mDummyDoubleSide;

    // cell clusters
    moris::Cell<Cell_Cluster_HMR> mCellClusters;


    // Block sets containing Cell Clusters
    moris::Cell<std::string>                     mPrimaryBlockSetNames;
    moris::Cell<moris::Cell<moris::moris_index>> mPrimaryBlockSetClusters;

    /*
     * Construct HMR Cell Clustering
     */
    void
    setup_cell_clusters(Interpolation_Mesh_HMR & aInterpolationMesh)
    {
        // check to see the meshes are the same (since all trivial)
        MORIS_ASSERT(this->get_num_nodes() == aInterpolationMesh.get_num_nodes(),"Mismatch nodes between integration and interpolation mesh");
        MORIS_ASSERT(this->get_num_elems() == aInterpolationMesh.get_num_elems(),"Mismatch elements between integration and interpolation mesh");

        // number of interpolation cells
        moris::uint tNumInterpCells = aInterpolationMesh.get_num_elems();

        for(moris::uint i = 0; i <tNumInterpCells; i++)
        {
                moris_id tCellId = aInterpolationMesh.get_glb_entity_id_from_entity_loc_index((moris_index)i,EntityRank::ELEMENT);

                // interpolation cell
                mtk::Cell const * tInterpCell = &aInterpolationMesh.get_mtk_cell((moris_index) i);
                mCellClusters(i).set_interpolation_cell( tInterpCell );

                // integration cell (only primary cells here)
                moris_index tIntegCellIndex    = this->get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
                mtk::Cell const * tPrimaryCell = &this->get_mtk_cell(tIntegCellIndex);
                mCellClusters(i).add_primary_integration_cell(tPrimaryCell);
        }
    }

    void
    setup_blockset_with_cell_clusters()
    {
        // construct integration to cell cluster relationship
        moris::Cell<moris::moris_index> tPrimaryIntegrationCellToClusterIndex(this->get_num_entities(EntityRank::ELEMENT),MORIS_INDEX_MAX);

        // iterate through clusters
        for(moris::uint  i = 0; i < mCellClusters.size(); i++)
        {
            Cell_Cluster_HMR const & tCellCluster = mCellClusters(i);
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

};
}
}


#endif /* PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_ */
