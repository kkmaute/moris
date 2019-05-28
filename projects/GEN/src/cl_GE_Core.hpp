/*
 * cl_GE_Core.hpp
 *
 *  Created on: May 17, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_CORE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_CORE_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"
#include "cl_GE_Nodal_Info.hpp"

// MTK includes
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Vertex.hpp"

// moris includes
#include "cl_Cell.hpp"

namespace moris
{
namespace ge
{
/*
 * @brief API to the geometry engine,
 *        * contains a member list of Nodal_Info objects who store all the relevant nodal information
 *        * each Nodal_Info object is associated with a {geom_rep,mesh} pair
 */
        class GE_Core
        {
        public:
            GE_Core(){};

            ~GE_Core(){};
            //------------------------------------------------------------------------------
            /**
             *  @brief set the current geometry and threshold, associate mesh with geometry, initialize
             *         nodal information container
             *
             *  @param[in] aGeometryPointer - pointer to geometry object
             *  @param[in] aMyManager       - pointer to mesh manager
             *  @param[in] aThreshold       - threshold value for current geometry LS function
             *
             */
            moris_index set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
                                      moris_index                   aMyMeshIndex = 0,
                                      real                          aThreshold = 0.0 )
            {
//                moris_index tGeomInd = mListOfGeoms.size();    //index of current geometry rep.
                mListOfGeoms.push_back( aGeomPointer );
                mThresholds.push_back( aThreshold );

                // initialize nodal information
                Nodal_Info tNodalInfo( aGeomPointer,aMyMeshIndex );

                mListOfNodalInfoObjects.push_back( tNodalInfo );
                return mListOfGeoms.size()-1;
            }
            //------------------------------------------------------------------------------
            std::shared_ptr< Geometry > get_geometry_pointer( moris_index aWhichGeometry )
            {
                return mListOfGeoms( aWhichGeometry );
            }
            //------------------------------------------------------------------------------
            /*
             * @brief return the nodal field vals from the specified NodalInfoObject
             *
             * @param[in] aWhichGeom - Nwhich geometry representation
             * @param[in] aWhichIndex - node/location to get information from
             *
             * @param[out] nodal information
             */
            Matrix< DDRMat > get_field_vals( moris_index aWhichGeom,
                                             moris_index aWhichIndex )
            {
                return mListOfNodalInfoObjects( aWhichGeom ).get_field_vals( aWhichIndex );
            };

            //------------------------------------------------------------------------------
            /*
             * @brief return the nodal sensitivity vals from the specified NodalInfoObject
             *
             * @param[in] aWhichGeom - which geometry representation
             * @param[in] aWhichIndex - node/location to get information from
             *
             * @param[out] nodal information
             */
            Matrix< DDRMat > get_sensitivity_vals( moris_index aWhichGeom,
                                                   moris_index aWhichIndex )
            {
                return mListOfNodalInfoObjects( aWhichGeom ).get_sensitivity_vals( aWhichIndex );
            };

            //------------------------------------------------------------------------------
            /*
             * @brief return the nodal normals from the specified NodalInfoObject
             *
             * @param[in] aWhichGeom - which geometry representation
             * @param[in] aWhichIndex - node/location to get information from
             *
             * @param[out] nodal information
             */
            Matrix< DDRMat > get_normal( moris_index aWhichGeom,
                                         moris_index aWhichIndex )
            {
                MORIS_ASSERT(false, "ge::GE_Core::get_nodal_normal(): currently not implemented yet");
                return mListOfNodalInfoObjects( aWhichGeom ).get_normal( aWhichIndex );
            };
            //------------------------------------------------------------------------------
            /*
             * @brief add a vertex to the nodal information object specified
             *
             * @param[in] aVertex - vertex to be added to nodal information tables
             * @param[in] aWhichGeom - which geometry representation
             */
            void add_vertex_and_value( mtk::Vertex &aVertex,
                                       moris_index  aWhichGeom )
            {
                mListOfNodalInfoObjects( aWhichGeom ).add_vertex_and_value( aVertex );
            };
            //------------------------------------------------------------------------------
            /*
             * @brief determine if there is an intersection
             *
             * @param[in] aWhichGeom - Which geometry representation are we checking for intersection with?
             * @param[in] aPrimitive - geometric primitive
             */
            //fixme need to loop over all the geometry representations to determine intersection with all of them
            void compute_intersection( moris_index aWhichGeom,
                                       Intersection_Object_Line*  aIntersectionObject )
            {
                return mListOfNodalInfoObjects( aWhichGeom ).compute_intersection(aIntersectionObject);
            };
            //------------------------------------------------------------------------------
            /*
             * @brief compute the LS value of a specified geometry representation at a specific coordinate
             *
             * @param[in] aCoordinate - coordinate where LS is to be evaluated at
             * @param[in] aWhichGeom - which geometry representation
             */
            real compute_value_at_point( Matrix< DDRMat > aCoordinate,
                                         moris_index  aWhichGeom )
            {
                MORIS_ASSERT( this->get_geometry_pointer(aWhichGeom)->get_geom_type() == GeomType::ANALYTIC, "ge::GE_Core::compute_value_at_point(): currently only implemented for analytic geometry representation" );
                return mListOfNodalInfoObjects( aWhichGeom ).compute_value_at_point( aCoordinate );
            };

//------------------------------------------------------------------------------
/*
 * ********************************************************************
 * BELOW FUNCTIONS ARE CURRENTLY BEING USED IN SOME HMR IMPLEMENTATIONS
 * ********************************************************************
 */
            /*
             * @brief determines elements (cells) intersected by the level set
             *
             * @param[in] aCells        - elements to be flagged for refinement
             * @param[in] aCandidates   - candidates for refinement
             * @param[in] aVertexValues - vertex values of scalar field
             * @param[in] aLowerBound   - lower bound of LS
             * @param[in] aUpperBound   - upper bound of LS
             */
            void
            find_cells_intersected_by_levelset(
                    Cell< mtk::Cell * >         & aCells,
                    Cell< mtk::Cell * >         & aCandidates,
                    const  Matrix< DDRMat >     & aVertexValues,
                    const  real                   aLowerBound = -0.0001,
                    const  real                   aUpperBound =  0.0001)
            {
                    // make sure that input makes sense
                    MORIS_ASSERT( aLowerBound <= aUpperBound,
                            "find_cells_intersected_by_levelset() : aLowerBound bound must be less or equal aUpperBound" );

                    // make sure that the field is a scalar field
                    MORIS_ASSERT( aVertexValues.n_cols() == 1,
                            "find_cells_within_levelset() can only be performed on scalar fields" );

                    // initialize output cell
                    aCells.resize( aCandidates.size(), nullptr );

                    // initialize counter
                    uint tCount = 0;

                    // loop over all candidates
                    for( mtk::Cell * tCell : aCandidates )
                    {
                        // get cell of vertex pointers
                        Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();

                        // get number of vertices on this element
                        uint tNumberOfVertices = tVertices.size();

                        // assign matrix with vertex values
                        Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );

                        // loop over all vertices and extract scalar field
                        for( uint k=0; k<tNumberOfVertices; ++k )
                        {
                            // copy value from field into element local matrix
                            tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
                        }

                        // test if cell is inside
                        if ( tCellValues.min() <= aUpperBound && tCellValues.max() >= aLowerBound )
                        {
                            // copy pointer to output
                            aCells( tCount++ ) = tCell;
                        }
                    }

                    // shrink output to fit
                    aCells.resize( tCount );
            }
            //------------------------------------------------------------------------------
                        /*
                         * @brief determines volume elements (cells)
                         *
                         * @param[in] aCells        - elements to be flagged for refinement
                         * @param[in] aCandidates   - candidates for refinement
                         * @param[in] aVertexValues - vertex values of scalar field
                         * @param[in] aUpperBound   - upper bound of LS
                         */
                        void
                        find_cells_within_levelset(
                                      Cell< mtk::Cell * >      & aCells,
                                      Cell< mtk::Cell * >      & aCandidates,
                                      const  Matrix< DDRMat >  & aVertexValues,
                                      const  uint                aUpperBound = 0.0 )
                        {


                            // make sure that the field is a scalar field
                            MORIS_ASSERT( aVertexValues.n_cols() == 1,
                                    "find_cells_within_levelset() can only be performed on scalar fields" );

                            // make sure that node values are calculated
                            //MORIS_ASSERT( tVertexValues.length() == aScalarField->get_num_nodes(),
                            //        "number of field values does not match number of vertices on block" );

                            // initialize output cell
                            aCells.resize( aCandidates.size(), nullptr );

                            // initialize counter
                            uint tCount = 0;

                            // loop over all candidates
                            for( mtk::Cell * tCell : aCandidates )
                            {
                                // get cell of vertex pointers
                                Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();

                                // get number of vertices on this element
                                uint tNumberOfVertices = tVertices.size();

                                // assign matrix with vertex values
                                Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );

                                // loop over all vertices and extract scalar field
                                for( uint k=0; k<tNumberOfVertices; ++k )
                                {
                                    // copy value from field into element local matrix
                                    tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
                                }

                                // test if cell is inside
                                if(  tCellValues.max() <= aUpperBound )
                                {
                                    // copy pointer to output
                                    aCells( tCount++ ) = tCell;
                                }
                            }

                            // shrink output to fit
                            aCells.resize( tCount );
                        }
//------------------------------------------------------------------------------
        private:
            // all corresponding information has the same index in each list (e.g. geom*(0) has corresponding threshold(0) and nodal information(0))
            moris::Cell< std::shared_ptr< Geometry > > mListOfGeoms;
            moris::Cell< real > mThresholds;
            moris::Cell< Nodal_Info > mListOfNodalInfoObjects;

//------------------------------------------------------------------------------
        protected:

        };

} // end namespace ge
} // end namespace moris




#endif /* PROJECTS_GEN_SRC_CL_GE_CORE_HPP_ */
