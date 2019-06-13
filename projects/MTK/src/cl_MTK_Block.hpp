/*
 * cl_MTK_Element_Block.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_BLOCK_HPP_
#define SRC_MESH_CL_MTK_BLOCK_HPP_

#include <string>

#include "typedefs.hpp" //MRS/COR/src
#include "fn_unique.hpp" //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

#include "cl_MTK_Cell_Cluster.hpp" //MTK/src
#include "cl_MTK_Set.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
        class Block : public Set
        {
        private :
            uint                              mNumVerticesOnBlock;
            moris::Matrix< DDSMat >           mVerticesOnBlock;

//------------------------------------------------------------------------------

            void calculate_vertices_on_blocks()
            {
                uint tMaxNumVert = 0;

                for( uint Ik = 0; Ik < mBlockSetClusters.size(); Ik++)
                {
                    for( uint Ij = 0; Ij < mBlockSetClusters( Ik )->get_primary_cells_in_cluster().size(); Ij++)
                    {
                        tMaxNumVert = tMaxNumVert + mBlockSetClusters( Ik )->get_primary_cells_in_cluster()( Ij )
                                                                           ->get_vertex_inds().numel();
                    }
                }

                moris::Matrix< DDSMat > tVerticesOnBlock(1, tMaxNumVert, -1 );

                uint tCounter = 0;

                for( uint Ik = 0; Ik < mBlockSetClusters.size(); Ik++)
                {
                    for( uint Ij = 0; Ij < mBlockSetClusters( Ik )->get_primary_cells_in_cluster().size(); Ij++)
                    {
                        //FIXME rewrite for more readability
                        tVerticesOnBlock( { 0, 0 },{ tCounter, tCounter + mBlockSetClusters( Ik )->get_primary_cells_in_cluster()( Ij )
                                                                                        ->get_vertex_inds().numel() - 1 }) =
                                mBlockSetClusters( Ik )->get_primary_cells_in_cluster()( Ij )
                                                       ->get_vertex_inds().matrix_data();

                        tCounter =tCounter + mBlockSetClusters( Ik )->get_primary_cells_in_cluster()( Ij )
                                                                    ->get_vertex_inds().numel();
                    }
                }

                MORIS_ASSERT( tVerticesOnBlock.min() != -1, "calculate_vertices_on_blocks(): negative vertex index");

                unique( tVerticesOnBlock, mVerticesOnBlock);

//                print(mVerticesOnBlock,"mNumVerticesOnBlock");

                mNumVerticesOnBlock = mVerticesOnBlock.numel();
            };

//------------------------------------------------------------------------------

        protected :
            moris::Matrix< DDUMat >           mMyBlockSetClusterInds;
            moris::Cell<Cell_Cluster const *> mBlockSetClusters;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Block( moris::Cell<Cell_Cluster const *>  aBlockSetClusters ) : mBlockSetClusters(aBlockSetClusters)
            {
//                mMyBlockSetClusterInds.set_size( aBlockSetClusters.size(), 1 );
//
//                for( uint Ik = 0; Ik < aBlockSetClusters.size(); Ik++)
//                {
//                    mMyBlockSetClusterInds( Ik, 0 ) = aBlockSetClusterInd( Ik );
//                }

                this->calculate_vertices_on_blocks();
            };

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            ~Block(){};

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
//              const moris::Matrix< DDUMat > &
//              get_list_of_block_cell_clusters() const
//              {
//                  return mMyBlockSetClusterInds;
//              }

//------------------------------------------------------------------------------

              const Cell_Cluster  *
              get_cell_clusters_by_index( moris_index aCellClusterIndex ) const
              {
                  return mBlockSetClusters( aCellClusterIndex );
              }

//------------------------------------------------------------------------------

              const uint
              get_num_vertieces_on_set() const
              {
                  return mNumVerticesOnBlock;
              }

//------------------------------------------------------------------------------

              moris::Matrix< DDSMat >
              get_vertieces_inds_on_block() const
              {
                  return mVerticesOnBlock;
              }

//------------------------------------------------------------------------------

//              moris::Matrix< DDSMat >
//              get_vertieces_on_block() const
//              {
//                  return mVerticesOnBlock;
//              }

//------------------------------------------------------------------------------
            /**
             * return a label that describes the block
             */
//            virtual std::string
//            get_label() const;

//------------------------------------------------------------------------------

//            /**
//             * sets the name of a block
//             */
//            virtual void
//            set_label( const std::string & aLabel );
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns the Id of this block
//             */
//            virtual moris_id
//            get_id() const;
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns the number of vertices owned and shared on this proc
//             */
//            virtual uint
//            get_number_of_vertices() const;
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns the number of element owned by this proc
//             */
//            virtual uint
//            get_number_of_cells() const;
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns a pointer to a vertex
//             */
//            virtual Vertex *
//            get_vertex_by_index( const moris_index & aIndex );
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns a pointer to a vertex ( const version )
//             */
//            virtual const Vertex *
//            get_vertex_by_index( const moris_index & aIndex ) const;
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns a pointer to a cell
//             */
//            virtual Cell *
//            get_cell_by_index( const moris_index & aIndex );
//
////------------------------------------------------------------------------------
//
//            /**
//             * returns a pointer to a cell ( const version )
//             */
//            virtual const Cell *
//            get_cell_by_index( const moris_index & aIndex ) const;
//
////------------------------------------------------------------------------------
//
//            virtual sint
//            get_number_of_adofs_used_by_proc() const;
//
////------------------------------------------------------------------------------
//
//            virtual void
//            get_adof_map( const uint aOrder, map< moris_id, moris_index > & aAdofMap  ) const;
//
////------------------------------------------------------------------------------
//
//            /**
//             * return the interpolation order of the cells on the block
//             */
//            virtual uint
//            get_interpolation_order() const;

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_BLOCK_HPP_ */
