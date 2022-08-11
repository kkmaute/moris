/*
 * cl_MTK_Cell_STK.hpp
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Logger.hpp"
#include "cl_Cell.hpp"              //MRS/CON/src
#include "cl_MTK_Vertex_STK.hpp"    //MTK/src
#include "cl_MTK_Cell.hpp"          //MTK/src
#include "cl_MTK_Enums.hpp"         //MTK/src
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
        /**
         * \brief the mtk::Cell class provides the cell information that is
         * provided by the mesh.
         */

        class Cell_STK : public moris::mtk::Cell
        {
            moris::Cell< Vertex* > mCellVertices;
            Mesh*                  mSTKMeshData = nullptr;
            uint                   mBulkSetId   = MORIS_UINT_MAX;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell_STK(){};

            //------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Cell_STK(
                    std::shared_ptr< moris::mtk::Cell_Info > aCellConn,
                    moris_id                                 aCellId,
                    moris_index                              aCellInd,
                    const moris::Cell< Vertex* >&            aCellVertices,
                    Mesh*                                    aStkImplementation,
                    uint                                     aBulkSetId = MORIS_UINT_MAX )
                    : Cell( aCellId, aCellInd, aStkImplementation->get_entity_owner( aCellInd, EntityRank::ELEMENT ), aCellConn )
                    , mCellVertices( aCellVertices )
                    , mSTKMeshData( aStkImplementation )
                    , mBulkSetId( aBulkSetId ){};
            //------------------------------------------------------------------------------

            /**
             * Destructor. Must be virtual.
             */
            ~Cell_STK(){};


            //------------------------------------------------------------------------------

            /**
             * tells how many vertices are connected to this cell
             */
            uint
            get_number_of_vertices() const
            {
                return mCellVertices.size();
            };

            //------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            moris::Cell< Vertex* >
            get_vertex_pointers() const
            {
                return mCellVertices;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a Matrix of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const
            {
                size_t           tNumVertices = this->get_number_of_vertices();
                Matrix< DDRMat > tVertexCoords( tNumVertices, mSTKMeshData->get_spatial_dim() );
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    tVertexCoords.set_row( i, mCellVertices( i )->get_coords() );
                }
                return tVertexCoords;
            }

            //------------------------------------------------------------------------------

            /**
             * bulk set id
             *           */
            uint
            get_bulk_set_id() const
            {
                return mBulkSetId;
            }
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_ */
