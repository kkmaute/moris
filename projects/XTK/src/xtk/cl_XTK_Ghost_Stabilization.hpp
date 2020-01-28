/*
 * cl_XTK_Ghost_Penalization.hpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_


#include "cl_XTK_Model.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
using namespace moris;

namespace moris
{
namespace mtk
{
    class Cell;
}
}
namespace xtk
{
/*
 * Data used while setting up the ghost stabilization (to not have to pass everything around in function inputs)
 */
struct Ghost_Setup_Data
{
    Cell<moris_index> mSubphaseIndexToInterpolationCellIndex;
    Cell<moris_index> mDblSideSetIndexInMesh;


    // interpolation cells
    Cell<Cell<moris_index>> mMasterSideIpCells;
    Cell<Cell<moris_index>> mSlaveSideIpCells;

    // Side ordinals
    Cell<Cell<moris_index>> mMasterSideIgCellSideOrds;
    Cell<Cell<moris_index>> mSlaveSideIgCellSideOrds;

    // trivial clusters are ones that are not on a hanging node interface
    Cell<Cell<moris_index>> mTrivialFlag;
    Cell<Cell<moris_index>> mTransitionLocation;
};

class Model;

class Ghost_Stabilization
{
public:
    Ghost_Stabilization():
        mXTKModel(nullptr)
    {}

    Ghost_Stabilization( Model* aXTKModel ):
        mXTKModel(aXTKModel)
    {}

    void
    setup_ghost_stabilization();

    void
    visualize_ghost_on_mesh(moris_index aBulkPhase);

    std::string
    get_ghost_dbl_side_set_name(moris_index aBulkPhase);

private:
    /*!
     * Construct interpolation cell for trivial cell cluster built on child mesh subphases.
     *
     * This is done so that we can still have 1 cluster per interpolation cell.
     */
    void
    construct_ip_ig_cells_for_ghost_side_clusters(Ghost_Setup_Data & aGhostSetupData);

    void
    prepare_ip_cell_id_answers( Cell<Matrix<IndexMat>>                 & aReceivedEnrCellIds,
                                Cell<moris_id>                         & aNewInterpCellIds,
                                Cell<Matrix<IndexMat>>                 & aEnrCellIds,
                                std::unordered_map<moris_id, moris_id> & aBaseEnrIdToIndexInNonTrivialOwned);

    void
    create_not_owned_ghost_ip_cells( Ghost_Setup_Data &                                      aGhostSetupData,
                                     Enriched_Interpolation_Mesh &                           aEnrInterpMesh,
                                     Cell<Cell<Interpolation_Cell_Unzipped  *>> const & aNonTrivialNotOwnedInterpCells,
                                     Cell<Matrix<IndexMat>>                          const & aReceivedEnrCellIds);

    void
    create_owned_ghost_ip_cells( Ghost_Setup_Data &                    aGhostSetupData,
                                 Enriched_Interpolation_Mesh &         aEnrInterpMesh,
                                 Cell<Interpolation_Cell_Unzipped *> & aNonTrivialOwnedInterpCells,
                                 Cell<moris_id>                      & aEnrCellIds);


    void
    declare_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData);

    void
    construct_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData);

    bool
    create_ghost(Ghost_Setup_Data &  aGhostSetupData,
                 moris_index const & aFirstSubphase,
                 moris_index const & aSecondSubphase,
                 moris_index &       aTrivialFlag);

    std::shared_ptr<Side_Cluster>
    create_slave_side_cluster(Ghost_Setup_Data &  aGhostSetupData,
                              Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                              uint const & aBulkIndex,
                              uint const & aCellIndex);

    std::shared_ptr<Side_Cluster>
    create_master_side_cluster(Ghost_Setup_Data &  aGhostSetupData,
                               Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                               uint const & aBulkIndex,
                               uint const & aCellIndex,
                               Side_Cluster* aSlaveSideCluster,
                               moris_index & aCurrentIndex,
                               moris_index & aCurrentId);

    moris::mtk::Cell*
    create_non_trivial_master_ig_cell(Ghost_Setup_Data &  aGhostSetupData,
                                      uint const & aBulkIndex,
                                      uint const & aCellIndex,
                                      Side_Cluster* aMasterSideCluster,
                                      Side_Cluster* aSlaveSideCluster,
                                      moris_index & aCurrentIndex,
                                      moris_index & aCurrentId);

    void
    permute_slave_vertices(moris::Cell<moris::mtk::Vertex const *> const & aSlaveVertices,
                           moris::Cell<moris::mtk::Vertex const *> const & aMasterVertices,
                           moris::Cell<moris::mtk::Vertex const *>  & aPermutedSlaveVertices,
                           moris::Cell<moris::mtk::Vertex const *>  & aPermutedAdjMastVertices);

    void
    get_local_coords_on_transition_side(moris_index const    & aMySideOrdinal,
                                        moris_index const    & aTransitionLoc,
                                        Cell<Matrix<DDRMat>> & aLocCoord);

    moris_index get_side_ordinals_for_non_trivial_master();

private:
    Model* mXTKModel; /*Pointer to the model*/

};

}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_ */
