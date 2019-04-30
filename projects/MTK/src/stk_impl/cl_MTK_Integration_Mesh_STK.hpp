/* cl_MTK_Integration_Mesh_STK.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Cell_Cluster_STK.hpp"
namespace moris
{
namespace mtk
{

class Interpolation_Mesh;

class MtkMeshData;

class Cell_Cluster_Input;

class Integration_Mesh_STK : public Mesh_Core_STK, public Integration_Mesh
{
    // Functions only valid for interpolation mIntegrationMeshes


public:
    /*!
     * Create STK integration mesh from an existing STK database
     */
    Integration_Mesh_STK(std::shared_ptr<Mesh_Data_STK> aSTKMeshData);

    /*!
     * Create a new integration mesh from file
     */
    Integration_Mesh_STK(
            std::string    aFileName,
            MtkMeshData*   aSuppMeshData,
            const bool     aCreateFacesAndEdges = true );

    /*!
     * Create a new integration mesh from data structure
     */
    Integration_Mesh_STK( MtkMeshData & aMeshData );

    /*!
     * Create a new integration mesh from data structure
     * with a link to an interpolation mesh
     */
    Integration_Mesh_STK( MtkMeshData & aMeshData,
                          Interpolation_Mesh* aInterpMesh,
                          Cell_Cluster_Input* aCellClusterData = nullptr);

    /*!
     * Create a integration mesh from an existing interpolation mesh
     */
    explicit
    Integration_Mesh_STK(Interpolation_Mesh & aInterpMesh,
                         Cell_Cluster_Input * aCellClusterInput);


    //##############################################
    // Cell Cluster Access
    //##############################################

    /*
     * Get a cell cluster related to an interpolation
     * cell
     */
    Cell_Cluster const &
    get_cell_cluster(Cell const & aInterpCell) const;

    /*
     * Get a cell cluster related to an interpolation
     * cell
     */
    Cell_Cluster const &
    get_cell_cluster(moris_index aInterpCellIndex) const;

    //##############################################
    // Block set with cluster access
    //##############################################
    moris::Cell<std::string>
    get_block_set_names()
    {
        return mPrimaryBlockSetNames;
    }

    moris::Cell<Cell_Cluster const *>
    get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const ;


private:
    // Cell Clusters
    moris::Cell<Cell_Cluster_STK> mCellClusters;

    // Block sets containing Cell Clusters
    moris::Cell<std::string>                     mPrimaryBlockSetNames;
    moris::Cell<moris::Cell<moris::moris_index>> mPrimaryBlockSetClusters;

    /*!
     * Setup the clustering interface
     */
    void
    setup_cell_clusters(Interpolation_Mesh & aInterpMesh,
                        Cell_Cluster_Input * aCellClusterInput);

    void
    setup_blockset_with_cell_clusters();

    moris::Cell<moris::mtk::Cell const *>
    get_cell_pointers_from_ids(moris::Matrix<moris::IdMat> const & aCellIds) const;

    moris::Cell<moris::mtk::Vertex const *>
    get_vertex_pointers_from_ids(moris::Matrix<moris::IdMat> const & aVertexIds) const;
};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_ */
