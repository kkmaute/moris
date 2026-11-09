/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Cluster_Measure_Derivatives.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_MTK_Cell_ISC.hpp"
#include "cl_MTK_Cell_Info_Quad9.hpp"
#include "cl_MTK_Cell_Cluster_Proxy.hpp"
#include "cl_MTK_Vertex_Proxy.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    class Bulk_Measure_Test_Model : public FEM_Model
    {
      public:
        Matrix< DDSMat > mAssembly;
        uint             mRequests = 0;

        void get_integration_xyz_pdv_assembly_indices( const Matrix< IndexMat >& aNodes,
                const Vector< gen::PDV_Type >&                                   aTypes,
                Matrix< DDSMat >&                                                aAssembly ) override
        {
            ++mRequests;
            aAssembly = mAssembly.get_row( aNodes( 0 ) );
        }
    };

    class Bulk_Measure_Test_Mesh_Cluster : public moris::Cell_Cluster_Proxy
    {
      public:
        mutable uint                 mBulkCalls   = 0;
        mutable uint                 mScalarCalls = 0;
        mutable mtk::Leader_Follower mLastSide    = mtk::Leader_Follower::LEADER;
        mutable mtk::Primary_Void    mLastPhase   = mtk::Primary_Void::PRIMARY;

        Matrix< DDRMat > compute_cluster_cell_measure_derivatives( const Matrix< DDRMat >& aCoordinates,
                const Matrix< DDSMat >&                                                    aRequested,
                mtk::Primary_Void                                                          aPhase,
                mtk::Leader_Follower                                                       aSide ) const override
        {
            ++mBulkCalls;
            mLastSide  = aSide;
            mLastPhase = aPhase;

            return moris::Cell_Cluster_Proxy::compute_cluster_cell_measure_derivatives( aCoordinates, aRequested, aPhase, aSide );
        }

        real compute_cluster_cell_measure_derivative( const Matrix< DDRMat >& aCoordinates,
                uint                                                          aDirection,
                mtk::Primary_Void                                             aPhase,
                mtk::Leader_Follower                                          aSide ) const override
        {
            ++mScalarCalls;

            return moris::Cell_Cluster_Proxy::compute_cluster_cell_measure_derivative( aCoordinates, aDirection, aPhase, aSide );
        }
    };

    class Bulk_Measure_Test_Cluster : public Cluster
    {
      public:
        Bulk_Measure_Test_Cluster( mtk::Cluster const * aMeshCluster, Set* aSet )
        {
            mMeshCluster = aMeshCluster;
            mSet         = aSet;
        }
    };

    TEST_CASE( "FEM Bulk Cluster Measure Assembly", "[fem],[Bulk_Measure_Derivatives]" )
    {
        // A curved QUAD9 has active and inactive coordinate directions. Assembly
        // indices are deliberately permuted, with two nodes sharing one PDV.
        // FEM must make one bulk request and accumulate both contributions.
        auto             tInfo = std::make_shared< mtk::Cell_Info_Quad9 >();
        Matrix< DDRMat > tCoordinates;
        tInfo->get_loc_coords_of_cell( tCoordinates );

        Vector< mtk::Vertex_Proxy >    tVertices( 9 );
        Vector< mtk::Vertex* >         tPointers( 9 );
        Bulk_Measure_Test_Mesh_Cluster tMeshCluster;

        for ( uint i = 0; i < 9; ++i )
        {
            tCoordinates( i, 1 ) += 0.2 * ( 1.0 - tCoordinates( i, 0 ) * tCoordinates( i, 0 ) );
            tVertices( i ).mVertexCoord = tCoordinates.get_row( i );
            tVertices( i ).mVertexInd   = i;
            tPointers( i )              = &tVertices( i );
            tMeshCluster.mVerticesInCluster.push_back( &tVertices( i ) );
        }

        mtk::Cell_ISC tCell( 0, 0, 0, tInfo, tPointers );
        tMeshCluster.mPrimaryIntegrationCells = { &tCell };
        tMeshCluster.mVoidIntegrationCells    = { &tCell };

        Bulk_Measure_Test_Model tModel;
        tModel.mAssembly.set_size( 9, 2, -1 );
        tModel.mAssembly( 0, 0 ) = 3;
        tModel.mAssembly( 0, 1 ) = 0;
        tModel.mAssembly( 5, 0 ) = 1;
        tModel.mAssembly( 6, 1 ) = 2;
        tModel.mAssembly( 2, 0 ) = 3;

        Set tSet;
        tSet.set_equation_model( &tModel );
        tSet.get_geo_pdv_assembly_vector().set_size( 4, 1, 0 );

        Bulk_Measure_Test_Cluster tCluster( &tMeshCluster, &tSet );
        Matrix< DDRMat >          tExpected( 1, 4, 0.0 );

        for ( uint i = 0; i < 9; ++i )
        {
            for ( uint d = 0; d < 2; ++d )
            {
                const sint tIndex = tModel.mAssembly( i, d );

                if ( tIndex >= 0 )
                {
                    tExpected( tIndex ) += tCell.compute_cell_measure_deriv( i, d );
                }
            }
        }

        for ( auto tSide : { mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER } )
        {
            tMeshCluster.mBulkCalls   = 0;
            tMeshCluster.mScalarCalls = 0;

            const auto             tPhase  = tSide == mtk::Leader_Follower::LEADER ? mtk::Primary_Void::PRIMARY : mtk::Primary_Void::VOID;
            const Matrix< DDRMat > tActual = tCluster.compute_cluster_cell_measure_derivative( tPhase, tSide );

            CHECK( norm( tActual - tExpected ) < 1e-12 );
            CHECK( tMeshCluster.mBulkCalls == 1 );
            CHECK( tMeshCluster.mScalarCalls == 0 );
            CHECK( tMeshCluster.mLastSide == tSide );
            CHECK( tMeshCluster.mLastPhase == tPhase );
        }

        // Missing active directions and an empty PDV vector avoid bulk work.
        tModel.mAssembly.fill( -1 );
        tMeshCluster.mBulkCalls = 0;

        const Matrix< DDRMat > tInactive = tCluster.compute_cluster_cell_measure_derivative();

        CHECK( norm( tInactive ) == 0.0 );
        CHECK( tMeshCluster.mBulkCalls == 0 );
        tSet.get_geo_pdv_assembly_vector().set_size( 0, 1 );
        tModel.mRequests = 0;

        const Matrix< DDRMat > tEmpty = tCluster.compute_cluster_cell_measure_derivative();

        CHECK( tEmpty.numel() == 0 );
        CHECK( tModel.mRequests == 0 );
    }

}    // namespace moris::fem
