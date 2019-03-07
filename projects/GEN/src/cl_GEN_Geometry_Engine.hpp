/*
 * cl_GEN_Geometry_Engine.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SRC_CL_GEN_GEOMETRY_ENGINE_HPP_
#define PROJECTS_GEN_SRC_CL_GEN_GEOMETRY_ENGINE_HPP_

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Block.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * This is a temporary class that provides the functionality of identifying
         * intersected and volume elements. This class will be merged into
         * something more general as the migration from LNA to LINALG continues
         */
        class Geometry_Engine
        {
//--------------------------------------------------------------------------------
        public:
//--------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Geometry_Engine(){};

//--------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Geometry_Engine(){};

//--------------------------------------------------------------------------------

//            void
//            find_cells_within_levelset(
//                          Cell< mtk::Cell * >      & aCells,
//                          Cell< mtk::Cell * >      & aCandidates,
//                          const  Matrix< DDRMat >  & aVertexValues,
//                   const              uint        aUpperBound = 0.0 )
//            {
//
//
//                // make sure that the field is a scalar field
//                MORIS_ASSERT( aVertexValues.n_cols() == 1,
//                        "find_cells_within_levelset() can only be performed on scalar fields" );
//
//                // make sure that node values are calculated
//                //MORIS_ASSERT( tVertexValues.length() == aScalarField->get_num_nodes(),
//                //        "number of field values does not match number of vertices on block" );
//
//                // initialize output cell
//                aCells.resize( aCandidates.size(), nullptr );
//
//                // initialize counter
//                uint tCount = 0;
//
//                // loop over all candidates
//                for( mtk::Cell * tCell : aCandidates )
//                {
//                    // get cell of vertex pointers
//                    Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();
//
//                    // get number of vertices on this element
//                    uint tNumberOfVertices = tVertices.size();
//
//                    // assign matrix with vertex values
//                    Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );
//
//                    // loop over all vertices and extract scalar field
//                    for( uint k=0; k<tNumberOfVertices; ++k )
//                    {
//                        // copy value from field into element local matrix
//                        tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
//                    }
//
//                    // test if cell is inside
//                    if(  tCellValues.max() <= aUpperBound )
//                    {
//                        // copy pointer to output
//                        aCells( tCount++ ) = tCell;
//                    }
//                }
//
//                // shrink output to fit
//                aCells.resize( tCount );
//            }

//--------------------------------------------------------------------------------

//            void
//            find_cells_intersected_by_levelset(
//                          Cell< mtk::Cell * > & aCells,
//                          Cell< mtk::Cell * > & aCandidates,
//                    const  Matrix< DDRMat >   & aVertexValues,
//                    const              real      aLowerBound = -0.0001,
//                    const              real      aUpperBound =  0.0001)
//            {
//                // make sure that input makes sense
//                MORIS_ASSERT( aLowerBound <= aUpperBound,
//                        "find_cells_intersected_by_levelset() : aLowerBound bound must be less or equal aUpperBound" );
//
//                // make sure that the field is a scalar field
//                MORIS_ASSERT( aVertexValues.n_cols() == 1,
//                        "find_cells_within_levelset() can only be performed on scalar fields" );
//
//                // initialize output cell
//                aCells.resize( aCandidates.size(), nullptr );
//
//                // initialize counter
//                uint tCount = 0;
//
//                // loop over all candidates
//                for( mtk::Cell * tCell : aCandidates )
//                {
//                    // get cell of vertex pointers
//                    Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();
//
//                    // get number of vertices on this element
//                    uint tNumberOfVertices = tVertices.size();
//
//                    // assign matrix with vertex values
//                    Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );
//
//                    // loop over all vertices and extract scalar field
//                    for( uint k=0; k<tNumberOfVertices; ++k )
//                    {
//                        // copy value from field into element local matrix
//                        tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
//                    }
//
//                    // test if cell is inside
//                    if ( tCellValues.min() <= aUpperBound && tCellValues.max() >= aLowerBound )
//                    {
//                        // copy pointer to output
//                        aCells( tCount++ ) = tCell;
//                    }
//                }
//
//                // shrink output to fit
//                aCells.resize( tCount );
//            }

//-------------------------------------------------------------------------------

        };
    }
}



#endif /* PROJECTS_GEN_SRC_CL_GEN_GEOMETRY_ENGINE_HPP_ */
