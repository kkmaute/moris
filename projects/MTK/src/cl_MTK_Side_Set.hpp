/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Set.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_SIDE_SET_HPP_
#define SRC_MESH_CL_MTK_SIDE_SET_HPP_

#include <string>

#include "typedefs.hpp"               //MRS/COR/src
#include "fn_unique.hpp"              //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp"          //MTK/src
#include "cl_MTK_Cell.hpp"            //MTK/src

#include "cl_MTK_Cell_Cluster.hpp"    //MTK/src
#include "cl_MTK_Set.hpp"             //MTK/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------
        class Side_Set : public Set
        {
          private:
            uint                    mNumVerticesOnSet;
            moris::Matrix< DDSMat > mVerticesOnSet;

            uint                    mNumCellsOnBlock;
            moris::Matrix< DDSMat > mCellsOnBlock;

            //------------------------------------------------------------------------------

            void
            calculate_vertices_on_set()
            {
                uint tMaxNumVert = 0;

                for ( uint Ik = 0; Ik < mSetClusters.size(); Ik++ )
                {
                    const Matrix< IndexMat > tSideOrdinal = mSetClusters( Ik )->get_cell_side_ordinals();

                    const moris::Cell< moris::mtk::Cell const * > &tCellsInCluster =    //
                            mSetClusters( Ik )->get_primary_cells_in_cluster();

                    for ( uint Ij = 0; Ij < tCellsInCluster.size(); Ij++ )
                    {
                        tMaxNumVert = tMaxNumVert +    //
                                      tCellsInCluster( Ij )->get_vertices_ind_on_side_ordinal( tSideOrdinal( Ij ) ).numel();
                    }
                }

                moris::Matrix< DDSMat > tVerticesOnSet( 1, tMaxNumVert, -1 );

                uint tCounter = 0;

                for ( uint Ik = 0; Ik < mSetClusters.size(); Ik++ )
                {
                    Matrix< IndexMat > tSideOrdinal = mSetClusters( Ik )->get_cell_side_ordinals();

                    const moris::Cell< moris::mtk::Cell const * > &tCellsInCluster =    //
                            mSetClusters( Ik )->get_primary_cells_in_cluster();

                    for ( uint Ij = 0; Ij < tCellsInCluster.size(); Ij++ )
                    {
                        uint tNumVertices = tCellsInCluster( Ij )->get_vertices_ind_on_side_ordinal( tSideOrdinal( Ij ) ).numel();

                        tVerticesOnSet( { 0, 0 }, { tCounter, tCounter + tNumVertices - 1 } ) =
                                tCellsInCluster( Ij )->get_vertices_ind_on_side_ordinal( tSideOrdinal( Ij ) ).matrix_data();

                        tCounter = tCounter + tNumVertices;
                    }
                }

                unique( tVerticesOnSet, mVerticesOnSet );

                // FIXME delete this if statement the unique bud in ARMADILLO is fixed.
                if ( tVerticesOnSet.n_rows() != mVerticesOnSet.n_rows() )
                {
                    tVerticesOnSet = mVerticesOnSet;
                    mVerticesOnSet.set_size( 1, mVerticesOnSet.numel() );

                    for ( uint Ik = 0; Ik < tVerticesOnSet.numel(); Ik++ )
                    {
                        mVerticesOnSet( Ik ) = tVerticesOnSet( Ik );
                    }
                }

                mNumVerticesOnSet = mVerticesOnSet.numel();
            };

            //------------------------------------------------------------------------------

        void
        init_ig_geometry_type()
        {
            mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

                if ( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    mIGGeometryType = mSetClusters( 0 )->get_primary_cells_in_cluster()( 0 )->get_geometry_type();
                }

                mIGGeometryType = get_auto_side_geometry_type( mIGGeometryType );
            }

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Side_Set(
                    std::string const                    &aName,
                    moris::Cell< Cluster const * > const &aSideSetClusters,
                    Matrix< IndexMat > const             &aColors,
                    uint const                           &aSpatialDim )
                    : Set( aName, aSideSetClusters, aColors, aSpatialDim )
            {
                mSetType = mtk::SetType::SIDESET;

                this->calculate_vertices_on_set();

                this->init_ig_geometry_type();
            };

            //------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            ~Side_Set(){};

            //------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            // const moris::Matrix< DDUMat > &
            // get_list_of_block_cell_clusters() const
            // {
            //     return mMyBlockSetClusterInds;
            // }

            //------------------------------------------------------------------------------

            const Cluster *
            get_clusters_by_index( moris_index aClusterIndex ) const
            {
                return mSetClusters( aClusterIndex );
            }

            //------------------------------------------------------------------------------

            uint
            get_num_vertices_on_set( const bool aOnlyPrimary )
            {
                return mNumVerticesOnSet;
            }

            //------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_ig_vertices_inds_on_block( const bool aOnlyPrimary )
            {
                return mVerticesOnSet;
            }

            //------------------------------------------------------------------------------

            uint
            get_num_cells_on_set( const bool aOnlyPrimary )
            {
                MORIS_ERROR( false, "not implemented" );
                return mNumCellsOnBlock;
            }

            //------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_cell_inds_on_block( const bool aOnlyPrimary )
            {
                MORIS_ERROR( false, "not implemented" );
                return mCellsOnBlock;
            }

            //------------------------------------------------------------------------------

            moris::uint
            get_num_clusters_on_set() const
            {
                return mSetClusters.size();
            }

            //------------------------------------------------------------------------------

            moris::Cell< Cluster const * >
            get_clusters_on_set() const
            {
                return mSetClusters;
            }

            //------------------------------------------------------------------------------

            size_t
            capacity()
            {
                size_t tTotalSize = 0;

                // name of the set
                tTotalSize += sizeof( mNumVerticesOnSet );
                tTotalSize += sizeof( mNumCellsOnBlock );

                tTotalSize += sizeof( mVerticesOnSet ) + mVerticesOnSet.capacity();
                tTotalSize += sizeof( mCellsOnBlock ) + mCellsOnBlock.capacity();

                tTotalSize += Set::capacity();

                return tTotalSize;
            }

            //------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            // virtual std::string
            // get_label() const;

            //------------------------------------------------------------------------------

            /**
             * sets the name of a block
             */
            // virtual void
            // set_label( const std::string & aLabel );

            //------------------------------------------------------------------------------

            /**
             * returns the Id of this block
             */
            // virtual moris_id
            // get_id() const;

            //------------------------------------------------------------------------------

            /**
             * returns the number of vertices owned and shared on this proc
             */
            // virtual uint
            // get_number_of_vertices() const;

            //------------------------------------------------------------------------------

            /**
             * returns the number of element owned by this proc
             */
            // virtual uint
            // get_number_of_cells() const;

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to a vertex
             */
            // virtual Vertex *
            // get_vertex_by_index( const moris_index & aIndex );

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to a vertex ( const version )
             */
            // virtual const Vertex *
            // get_vertex_by_index( const moris_index & aIndex ) const;

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to a cell
             */
            // virtual Cell *
            // get_cell_by_index( const moris_index & aIndex );

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to a cell ( const version )
             */
            // virtual const Cell *
            // get_cell_by_index( const moris_index & aIndex ) const;

            //------------------------------------------------------------------------------

            // virtual sint
            // get_number_of_adofs_used_by_proc() const;

            //------------------------------------------------------------------------------

            // virtual void
            // get_adof_map( const uint aOrder, map< moris_id, moris_index > & aAdofMap  ) const;

            //------------------------------------------------------------------------------

            /**
             * return the interpolation order of the cells on the block
             */
            // virtual uint
            // get_interpolation_order() const;

            //------------------------------------------------------------------------------

        };    // end class: mtk::Side_Set

        //------------------------------------------------------------------------------

    }    // namespace mtk

}    // namespace moris

//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_SIDE_SET_HPP_ */
