/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_Cluster.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_CLUSTER_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_CLUSTER_HPP_

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Function_Base.hpp"

// forward declare the mtk::Cluster_Group
namespace moris::mtk
{
    class Cluster_Group;
}

using namespace moris;

namespace moris::xtk
{
    class Interpolation_Cell_Unzipped;
    class Child_Mesh;

    struct IG_Cell_Group;
    struct IG_Vertex_Group;
    struct Facet_Based_Connectivity;

    class Cell_Cluster : public mtk::Cell_Cluster
    {
        friend class Enriched_Integration_Mesh;
        friend class Side_Cluster;

        //------------------------------------------------------------------------------

      protected:
        bool                                          mTrivial;
        bool                                          mVoid = false;
        bool                                          mInvalid;
        bool                                          mOnlyForVis = false;
        Interpolation_Cell_Unzipped const *           mInterpolationCell;
        Child_Mesh const *                            mChildMesh;    // FIXME: this doesn't seem to be used, should be removed
        Vector< moris::mtk::Cell const * >            mPrimaryIntegrationCells;
        Vector< moris::mtk::Cell const * >            mVoidIntegrationCells;
        Vector< moris::mtk::Vertex const * >          mVerticesInCluster;
        Matrix< DDRMat >                              mLocalCoords;
        Matrix< DDRMat >                              mQuadraturePoints; // Quadrature points.
        Matrix< DDRMat >                              mQuadratureWeights; // Quadrature weights
        std::shared_ptr< IG_Vertex_Group >            mVertexGroup;
        Vector< std::shared_ptr< IG_Cell_Group > >    mPrimaryIgCellGroup;
        Vector< std::shared_ptr< IG_Cell_Group > >    mVoidIgCellGroup;
        Vector< std::weak_ptr< mtk::Cluster_Group > > mClusterGroups;
        Vector< Vector< moris::mtk::Vertex* > >       mFacetVerticesOnSubphaseBoundary;
        Vector< Matrix< DDRMat > >                    mFacetNormals;
        Vector< Matrix< DDRMat > >                    mFacetVertexCoordinates;
        mtk::Interpolation_Function_Base*             mIPInterp; // Interpolation object to get the basis functions for moment fitting
        

        //------------------------------------------------------------------------------

      public:
        Cell_Cluster();
        ~Cell_Cluster() override;
        bool                                       is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        bool                                       is_full() const;
        bool                                       is_void() const;
        bool                                       is_invalid() const;
        Vector< moris::mtk::Cell const * > const & get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        Vector< moris::mtk::Cell const * > const & get_void_cells_in_cluster() const override;
        moris::mtk::Cell const &                   get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        Vector< moris::mtk::Vertex const * >       get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        Matrix< DDRMat >                           get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        Matrix< DDRMat >                           get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex, const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        moris_index                                get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;
        Matrix< DDRMat >                           get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const override;
        Matrix< DDRMat >                           get_quadrature_points() const override;
        Matrix< DDRMat >                           get_quadrature_weights() const override;

        //Vector< real >                             compute_quadrature_weights_moment_fitting( const uint, const uint ) ;

        // constructor for cell clusters that will only be used for visualization purposes
        // this results in cell-clusters that may miss some data and have reduced functionality
        Cell_Cluster( bool aOnlyForVisualization );

        //------------------------------------------------------------------------------

        // functions for internal XTK use
        Interpolation_Cell_Unzipped const *
        get_xtk_interpolation_cell() const;

        //------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_hanging_nodes() const override;

        //------------------------------------------------------------------------------

        // memory
        size_t
        capacity();

        //------------------------------------------------------------------------------

        void
        set_primary_integration_cell_group( const std::shared_ptr< IG_Cell_Group >& aPrimaryIgCells );

        //------------------------------------------------------------------------------

        void
        set_primary_integration_cell_groups( Vector< std::shared_ptr< IG_Cell_Group > > aPrimaryIgCells );

        //------------------------------------------------------------------------------

        void
        set_void_integration_cell_groups( Vector< std::shared_ptr< IG_Cell_Group > >& aVoidIgCells );

        //------------------------------------------------------------------------------

        void
        set_ig_vertex_group( const std::shared_ptr< IG_Vertex_Group >& aVertexGroup );

        //----------------------------------------------------------------

        bool
        has_cluster_group( const moris_index aDiscretizationMeshIndex ) const override;

        //----------------------------------------------------------------

        std::shared_ptr< mtk::Cluster_Group >
        get_cluster_group( const moris_index aDiscretizationMeshIndex ) const override;

        //------------------------------------------------------------------------------

        void
        set_cluster_group(
                const moris_index                            aDiscretizationMeshIndex,
                const std::shared_ptr< mtk::Cluster_Group >& aClusterGroupPtr ) override;
        
        //-------------------------------------------------------------------------------

        void 
        find_subphase_boundary_vertices(
                const std::shared_ptr< IG_Cell_Group >  aSubphaseIGCells,
                const std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
                const uint aDim
        );

        // ------------------------------------------------------------------------------

        void 
        set_quadrature_points( const uint aOrder , const uint aDim ) ;

        // -----------------------------------------------------------------------------

        void
        set_quadrature_weights( const uint aOrder , const uint aDim ) ;

        // -----------------------------------------------------------------------------

        void
        compute_quadrature_weights( const uint aOrder , const uint aDim ) ;

        //------------------------------------------------------------------------------

        std::shared_ptr< IG_Vertex_Group >
        get_ig_vertex_group();

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override;

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override;

        //------------------------------------------------------------------------------

    };    // class Cell_Cluster
}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_CLUSTER_HPP_ */
