/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Checker.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_Vector.hpp"
#include "fn_isvector.hpp"
#include <unordered_map>

namespace moris::mtk
{

    class Serialized_Mesh_Data
    {
      public:
        // spatial dim
        moris_index mSpatialDim;

        // individual proc serial data

        // vertex related data
        Matrix< IdMat >            mVertexIds;                   //(x)
        Matrix< IdMat >            mVertexOwners;                //(x)
        Matrix< DDRMat >           mVertexCoordinates;           //(x)
        Vector< Matrix< DDRMat > > mVertexTMatrixWeights;        //(x)
        Vector< Matrix< IdMat > >  mVertexTMatrixBasisIds;       //(x)
        Vector< Matrix< IdMat > >  mVertexTMatrixBasisOwners;    //(x)

        // cell related data
        Matrix< IdMat >                                          mCellIds;       //(x)
        Matrix< IdMat >                                          mCellOwners;    //(x)
        Vector< std::unordered_map< moris_index, moris_index > > mCollectCellMaps;
        std::unordered_map< moris_index, moris_index >           mSerialCellMap;
        Vector< moris_index >                                    mCellSerialIndexToId;
        Vector< Matrix< IdMat > >                                mCollectCellIds;       //(x)
        Vector< Matrix< IdMat > >                                mCollectCellOwners;    //(x)

        Vector< Matrix< IdMat > >   mCellToVertex;    //(-) outer cell is the  cell topology
        Vector< enum CellTopology > mCellTopo;        //(-)

        // Block Sets
        Vector< std::string >     mCellSetNames;
        Vector< Matrix< IdMat > > mCellsInCellSet;

        // collected data on proc 0
        Vector< std::unordered_map< moris_index, moris_index > > mCollectVertexMaps;
        Vector< Matrix< IdMat > >                                mCollectVertexIds;                //(x)
        Vector< Matrix< IdMat > >                                mCollectVertexOwners;             //(x)
        Vector< Matrix< DDRMat > >                               mCollectVertexCoords;             //(x)
        Vector< Vector< Matrix< DDRMat > > >                     mCollectVertexTMatrixWeights;     //(x)
        Vector< Vector< Matrix< IdMat > > >                      mCollectVertexTMatrixBasisIds;    //(x)
        Vector< Vector< Matrix< IdMat > > >                      mCollectVertexTMatrixBasisOwners;

        // overall
        Vector< moris_index >                          mVertexSerialIndexToId;
        std::unordered_map< moris_index, moris_index > mSerialVertexMap;

        //
        moris_index
        get_vertex_index(
                const moris_index& aVertexId,
                const moris_index& aProcIndex ) const;

        //
        moris::Matrix< moris::DDRMat >
        get_vertex_coords(
                const moris_index& aVertexId,
                const moris_index& aProcIndex ) const;
    };

    class Mesh_Checker
    {
      public:
        Mesh_Checker();
        Mesh_Checker( moris_index   aMeshIndex,
                Interpolation_Mesh* aIpMesh,
                Integration_Mesh*   aIgMesh );
        ~Mesh_Checker();

        bool perform();
        void print_diagnostics();

        bool verify_double_side_sets( Integration_Mesh const * aIgMesh );
        bool verify_side_cluster( Cluster const * aCluster, enum Leader_Follower aLeaderFollower = Leader_Follower::LEADER );
        bool verify_basis_indexing( Interpolation_Mesh const * aIPMesh );

      private:
        moris_index         mMeshIndex;
        Interpolation_Mesh* mIpMesh;
        Integration_Mesh*   mIgMesh;

        Serialized_Mesh_Data mSerializedIpMesh;
        Serialized_Mesh_Data mSerializedIgMesh;

        // diagnostic flags
        // interp mesh
        bool mIpVertexDiag      = false;
        bool mIpVertexOwnerDiag = false;
        bool mIpVertexBasisDiag = false;
        bool mIpCellOwnerDiag   = false;

        // integ mesh
        bool mIgVertexDiag      = false;
        bool mIgVertexOwnerDiag = false;
        bool mIgCellOwnerDiag   = false;

        // Serialized mesh accessors
        moris_index
        get_vertex_proc_index_from_id( moris_id aId, moris_index aProc, Serialized_Mesh_Data* aSerializedMesh );

        // only proc 0
        moris_index
        get_vertex_serial_index_from_id( moris_id aId, Serialized_Mesh_Data* aSerializedMesh );

        moris_index
        get_cell_proc_index_from_id( moris_id aId, moris_index aProc, Serialized_Mesh_Data* aSerializedMesh );

        // only proc 0
        moris_index
        get_cell_serial_index_from_id( moris_id aId, Serialized_Mesh_Data* aSerializedMesh );

        void
        serialize_mesh();

        void
        serialize_mesh_core();

        void
        serialize_vertices(
                Mesh*                 aMesh,
                Serialized_Mesh_Data* aSerializedMesh );

        void
        serialize_vertex_t_matrices(
                Interpolation_Mesh*   aMesh,
                Serialized_Mesh_Data* aSerializedMesh );

        void
        serialize_cells(
                Mesh*                 aMesh,
                Serialized_Mesh_Data* aSerializedMesh );

        void
        gather_serialized_mesh( Serialized_Mesh_Data* aSerializedMesh );

        void
        gather_serialized_ip_mesh( Serialized_Mesh_Data* aSerializedMesh );

        void
        setup_vertex_maps( Serialized_Mesh_Data* aSerializedMesh );

        void
        setup_cell_maps( Serialized_Mesh_Data* aSerializedMesh );

        std::string
        bool_to_string( bool aBool );

        /*!
         * Verifies node coordinates in the mesh.
         * @param[in] aMesh - An integration or interpolation mesh
         * @param[in] aStackedVertexFlag - if true more than one vertex can be at the same spatial location
         */
        bool
        verify_vertex_coordinates(
                Serialized_Mesh_Data* aSerializedMesh,
                bool                  aStackedVertexFlag );

        /*!
         * Verifies node ownership in the mesh.
         * @param[in] aMesh - An integration or interpolation mesh
         * @param[in] aStackedVertexFlag - if true more than one vertex can be at the same spatial location
         */
        bool
        verify_vertex_ownership( Serialized_Mesh_Data* aSerializedMesh );

        /*!
         * Verifies cell ownership in the mesh.
         * @param[in] aMesh - An integration or interpolation mesh
         */
        bool
        verify_cell_ownership( Serialized_Mesh_Data* aSerializedMesh );

        template< typename MatrixType >
        Matrix< MatrixType >
        concatenate_cell_of_mats( Vector< Matrix< MatrixType > > aMat,
                moris_index                                      aFixedDim )
        {
            moris_index tFixedDimSize = 0;
            moris_index tTotalSize    = 0;

            Vector< moris_index > tMatSize( 2 );

            MORIS_ASSERT( aFixedDim <= 2, "Fixed dim can be 0 for row or 1 for col." );

            for ( moris::uint i = 0; i < aMat.size(); i++ )
            {

                tMatSize( 0 ) = aMat( i ).n_rows();
                tMatSize( 1 ) = aMat( i ).n_cols();

                if ( i == 0 )
                {
                    tFixedDimSize = tMatSize( aFixedDim );
                }

                else
                {
                    MORIS_ASSERT( tMatSize( aFixedDim ) == tFixedDimSize, "fixed dimension not consistent." );
                }

                tTotalSize = tTotalSize + aMat( i ).numel();
            }

            // size the concatenated matrix
            Matrix< MatrixType > tConcatenatedMat;

            moris_index tStart = 0;
            moris_index tEnd   = 0;

            if ( aFixedDim == 0 )
            {
                tConcatenatedMat.resize( tFixedDimSize, tTotalSize / tFixedDimSize );

                for ( moris::uint i = 0; i < aMat.size(); i++ )
                {
                    tEnd = tStart + aMat( i ).n_cols() - 1;

                    tConcatenatedMat( { 0, tFixedDimSize - 1 }, { tStart, tEnd } ) = aMat( i ).matrix_data();
                    tStart                                                         = tEnd + 1;
                }
            }
            else if ( aFixedDim == 1 )
            {
                tConcatenatedMat.resize( tTotalSize / tFixedDimSize, tFixedDimSize );
                for ( moris::uint i = 0; i < aMat.size(); i++ )
                {
                    tEnd = tStart + aMat( i ).n_rows() - 1;

                    tConcatenatedMat( { tStart, tEnd }, { 0, tFixedDimSize - 1 } ) = aMat( i ).matrix_data();
                    tStart                                                         = tEnd + 1;
                }
            }

            return tConcatenatedMat;
        }

        template< typename MatrixType >
        void
        cell_of_mats_to_flattened_mat(
                Vector< Matrix< MatrixType > > const & aCellOfMats,
                Matrix< MatrixType >&                  aData,
                Matrix< IndexMat >&                    aOffsets )
        {

            // setup offsets
            aOffsets.resize( 1, aCellOfMats.size() + 1 );
            aOffsets( 0 ) = 0;

            moris_index tCurrentIndex = 0;
            moris_index tDataSize     = 0;

            for ( moris::uint i = 0; i < aCellOfMats.size(); i++ )
            {

                MORIS_ASSERT( moris::isvector( aCellOfMats( i ) ), "Only implemented on vectors" );
                tDataSize = tDataSize + aCellOfMats( i ).numel();
            }

            aData.resize( 1, tDataSize );

            for ( moris::uint i = 0; i < aCellOfMats.size(); i++ )
            {

                for ( moris::uint j = 0; j < aCellOfMats( i ).numel(); j++ )
                {
                    aData( tCurrentIndex ) = aCellOfMats( i )( j );

                    tCurrentIndex++;
                }
                aOffsets( i + 1 ) = tCurrentIndex;
            }
        }

        template< typename MatrixType >
        void
        flattened_mat_to_cell_of_mats(
                Matrix< MatrixType > const &    aData,
                Matrix< IndexMat > const &      aOffsets,
                Vector< Matrix< MatrixType > >& aCellOfMats )
        {
            aCellOfMats.resize( aOffsets.numel() - 1, Matrix< MatrixType >( 0, 0 ) );

            moris_index tStart = 0;

            for ( moris::uint i = 0; i < aOffsets.numel() - 1; i++ )
            {

                // number of basis interpolating into the vertex
                moris::moris_index tNumel = aOffsets( i + 1 ) - tStart;

                aCellOfMats( i ).resize( 1, tNumel );

                // itere and grab  data
                for ( moris::moris_index j = 0; j < tNumel; j++ )
                {
                    aCellOfMats( i )( j ) = aData( tStart + j );
                }

                tStart = tStart + tNumel;
            }
        }
    };
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_ */
