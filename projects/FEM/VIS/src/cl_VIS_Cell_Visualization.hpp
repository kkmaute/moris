/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Cell_Visualization.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Logger.hpp"
#include "cl_Cell.hpp"        //MRS/CON/src
#include "cl_MTK_Cell.hpp"    //MRS/CON/src
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
//------------------------------------------------------------------------------
namespace moris
{
    namespace vis
    {

        /**
         * \brief the mtk::Cell class provides the cell information that is
         * provided by the mesh.
         */

        class Cell_Visualization : public moris::mtk::Cell
        {
            moris::Cell< mtk::Vertex* > mCellVertices;
            const moris::mtk::Cell*     mIntegrationcell = nullptr;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell_Visualization()    //: mCellInfo(nullptr)
                    {};

            //------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Cell_Visualization( moris_id              aCellId,
                    moris_index                       aCellInd,
                    const moris::Cell< mtk::Vertex* > aCellVertices,
                    const moris::mtk::Cell*           aIntegrationCell )
                    : Cell( aCellId,
                            aCellInd,
                            aIntegrationCell->get_owner(),
                            aIntegrationCell->get_cell_info_sp() )
                    , mCellVertices( aCellVertices )
                    , mIntegrationcell( aIntegrationCell ){};

            //------------------------------------------------------------------------------

            /**
             * Destructor. Must be virtual.
             */
            ~Cell_Visualization(){};

            //------------------------------------------------------------------------------

            /**
             * returns the id of the cell with which this cell was created
             *
             * @return luint ID
             */
            moris_id
            get_mesh_cell_id() const
            {
                return mIntegrationcell->get_id();
            };

            //------------------------------------------------------------------------------

            /**
             * returns the index of the cell with which this cell was created
             *
             * @return luint ID
             */
            moris_index
            get_mesh_cell_index() const
            {
                return mIntegrationcell->get_index();
            }

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
             * returns the proc id of the owner of this cell
             * ( this information is needed for STK )
             */
            moris_id
            get_owner() const
            {
                return mIntegrationcell->get_owner();
            }

            //------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            moris::Cell< mtk::Vertex* >
            get_vertex_pointers() const
            {
                return mCellVertices;
            }

            //------------------------------------------------------------------------------

            // TODO MESHCLEANUP
            void
            remove_vertex_pointer( moris_index aIndex )
            {
                // std::cout << "In Vis Cell Visu" << std::endl;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a Matrix with IDs of connected vertices
             */
            Matrix< IdMat >
            get_vertex_ids() const
            {
                size_t          tNumVertices = this->get_number_of_vertices();
                Matrix< IdMat > tVertexIds( 1, tNumVertices );
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    tVertexIds( i ) = mCellVertices( i )->get_id();
                }
                return tVertexIds;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a Mat with indices of connected vertices
             */
            Matrix< IndexMat >
            get_vertex_inds() const
            {
                size_t             tNumVertices = this->get_number_of_vertices();
                Matrix< IndexMat > tVertexInds( 1, tNumVertices );
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    tVertexInds( i ) = mCellVertices( i )->get_index();
                }
                return tVertexInds;
            }

            //------------------------------------------------------------------------------
            
            /**
             * returns a Matrix of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const
            {
                size_t tNumVertices = this->get_number_of_vertices();

                Matrix< DDRMat > tVertexCoords( tNumVertices, mCellVertices( 0 )->get_coords().numel() );
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    tVertexCoords.set_row( i, mCellVertices( i )->get_coords() );
                }
                return tVertexCoords;
            }

            //------------------------------------------------------------------------------

        };    // class Cell_Visualization

        //------------------------------------------------------------------------------

    } /* namespace vis */
} /* namespace moris */

//------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_ */
