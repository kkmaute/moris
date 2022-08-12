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
#include "cl_MTK_Side_Cluster_STK.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

namespace moris
{
    namespace mtk
    {

        class Interpolation_Mesh;

        class MtkMeshData;

        class Cell_Cluster_Input;

        class Side_Cluster_Input;

        class Double_Side_Cluster_Input;

        class Integration_Mesh_STK : public Mesh_Core_STK
                , public Integration_Mesh
        {
            // Functions only valid for interpolation mIntegrationMeshes

          public:
            /*!
             * Create STK integration mesh from an existing STK database
             */
            Integration_Mesh_STK( std::shared_ptr< Mesh_Data_STK > aSTKMeshData );

            // ----------------------------------------------------------------------------

            /*!
             * Create a new integration mesh from file
             */
            Integration_Mesh_STK(
                    std::string  aFileName,
                    MtkMeshData *aSuppMeshData,
                    const bool   aCreateFacesAndEdges = true );

            // ----------------------------------------------------------------------------

            /*!
             * Create a new integration mesh from data structure
             */
            Integration_Mesh_STK( MtkMeshData &aMeshData );

            // ---------------------------------------------------------------------------

            /*!
             * Create a new integration mesh from data structure
             * with a link to an interpolation mesh
             */
            Integration_Mesh_STK(
                    MtkMeshData        &aMeshData,
                    Interpolation_Mesh *aInterpMesh );

            // ----------------------------------------------------------------------------

            /*!
             * Create a integration mesh from an existing interpolation mesh
             */
            explicit Integration_Mesh_STK(
                    Interpolation_Mesh &aInterpMesh,
                    Cell_Cluster_Input *aCellClusterInput );

            // ----------------------------------------------------------------------------

            ~Integration_Mesh_STK();

            //##############################################
            // Cell Cluster Access
            //##############################################

            /*
             * Get a cell cluster related to an interpolation
             * cell
             */
            Cell_Cluster const &
            get_cell_cluster( Cell const &aInterpCell ) const;

            // ----------------------------------------------------------------------------

            /*
             * Get a cell cluster related to an interpolation
             * cell
             */
            Cell_Cluster const &
            get_cell_cluster( moris_index aInterpCellIndex ) const;

            // ----------------------------------------------------------------------------

            //##############################################
            // Block set with cluster access
            //##############################################
            /*!
             * Returns the block set names
             */
            moris::Cell< std::string >
            get_block_set_names() const;

            // ----------------------------------------------------------------------------

            moris::Cell< Cluster const * >
            get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const;

            /*!
             * Returns the label
             */
            std::string
            get_block_set_label( moris_index aBlockSetOrdinal ) const;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the index given a label
             */
            moris_index
            get_block_set_index( std::string aBlockSetLabel ) const;

            // ----------------------------------------------------------------------------

            //##############################################
            // Side Set Cluster Access
            //##############################################

            /*!
             * Return a side set containing clusters
             */
            moris::Cell< Cluster const * >
            get_side_set_cluster( moris_index aSideSetOrdinal ) const;

            /*!
             * Retuns the number of side sets
             */
            uint
            get_num_side_sets() const;

            /*!
             * Returns the label
             */
            std::string
            get_side_set_label( moris_index aSideSetOrdinal ) const;

            /*!
             * Returns the index given a label
             */
            moris_index
            get_side_set_index( std::string aSideSetLabel ) const;

            //##############################################
            // Double Side Set Cluster Access
            //##############################################

            /*!
             * Returns the number of double sided side sets in the mesh
             */
            uint
            get_num_double_sided_sets() const;

            /*!
             * Returns the label
             */
            std::string
            get_double_sided_set_label( moris_index aSideSetOrdinal ) const;

            /*!
             * Returns the index given a label
             */
            moris_index
            get_double_sided_set_index( std::string aDoubleSideSetLabel ) const;

            /*!
             * Returns the double side clusters in the side set
             */
            moris::Cell< Cluster const * >
            get_double_side_set_cluster( moris_index aSideSetOrdinal ) const;

          private:
            // Cell Clusters
            moris::Cell< Cell_Cluster_STK > mCellClusters;

            // Block sets containing Cell Clusters
            std::unordered_map< std::string, moris_index >   mBlockSetLabelToOrd;
            moris::Cell< std::string >                       mPrimaryBlockSetNames;
            moris::Cell< moris::Cell< moris::moris_index > > mPrimaryBlockSetClusters;
            moris::Cell< moris::moris_index >                mIpCellToBlockSetOrd;

            // side sets
            std::unordered_map< std::string, moris_index > mSideSideSetLabelToOrd;
            moris::Cell< std::string >                     mSideSetLabels;
            moris::Cell< moris::Cell< Side_Cluster_STK > > mSideSets;

            // double side sets
            std::unordered_map< std::string, moris_index > mDoubleSideSetLabelToOrd;
            moris::Cell< std::string >                     mDoubleSideSetLabels;
            moris::Cell< moris::Cell< Cluster const * > >  mDoubleSideSets;
            moris::Cell< Side_Cluster_STK >                mDoubleSideSetSideClusters;

            /*!
             * Setup the clustering interface
             */
            void
            setup_cell_clusters( Interpolation_Mesh &aInterpMesh,
                    Cell_Cluster_Input              *aCellClusterInput );

            void setup_cell_clusters();

            /*!
             * Setup the blocksets which contain cell clusters
             */
            void
            setup_blockset_with_cell_clusters();

            void setup_blockset_with_cell_clusters_trivial();

            /*
             *  setup the side set cluster interface
             */
            void
            setup_side_set_clusters(
                    Interpolation_Mesh &aInterpMesh,
                    Side_Cluster_Input *aSideClusterInput );

            void setup_side_set_clusters_trivial();

            /*
             *  setup the double side set cluster interface
             */
            void
            setup_double_side_set_clusters(
                    Interpolation_Mesh        &aInterpMesh,
                    Double_Side_Cluster_Input *aSideClusterInput );

            void
            setup_double_side_set_clusters_all_trivial( Interpolation_Mesh &aInterpMesh );

            moris::Cell< moris::mtk::Cell const * >
            get_cell_pointers_from_ids( moris::Matrix< moris::IdMat > const &aCellIds ) const;

            moris::Cell< moris::mtk::Vertex const * >
            get_vertex_pointers_from_ids( moris::Matrix< moris::IdMat > const &aVertexIds ) const;
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_ */
