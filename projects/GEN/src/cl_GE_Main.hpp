/*
 * cl_GE_Main.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_HPP_

// moris includes
#include "cl_Cell.hpp"

// GE includes
#include "cl_GE_Geometry.hpp"

// MTK includes
#include "cl_MTK_Mesh.hpp"

// LinAlg includes
#include "op_minus.hpp"
#include "fn_norm.hpp"

namespace moris
{
	namespace ge
	{
		class GE
		{
//------------------------------------------------------------------------------

			// member variables
			moris::Cell< Geometry* > 		mListOfGeoms;
			moris::Cell< mtk::Mesh* > 		mListOfMeshs;
			moris::Cell< real > 			mThresholds;

//------------------------------------------------------------------------------
		private:
            real
            sphere_function( const Matrix< DDRMat > & aPoint, Cell< real > aInputs )
            {
                // aInputs(0) = x location of center
                // aInputs(1) = y location of center
                // aInputs(2) = z location of center
                // aInputs(3) = radius
                Matrix< DDRMat > tCenterVec(1,3);
                tCenterVec(0,0) = aInputs(0);
                tCenterVec(0,1) = aInputs(1);
                tCenterVec(0,2) = aInputs(2);
                return norm( aPoint - tCenterVec ) - aInputs(3);
            }

//------------------------------------------------------------------------------
		protected:
//------------------------------------------------------------------------------
		public:
			GE(){};

			~GE(){};

			/**
			 *  @brief  set the current geometry and threshold
			 *
			 *  @param[in]  aGeometryPointer        pointer to geometry object
			 *  @param[in]  aThreshold              threshold value for current geometry LS function
			 *
			 */
			void set_geometry( Geometry* & aGeomPointer,
			                   real        aThreshold   = 0.0)
			{
				mListOfGeoms.push_back( aGeomPointer );
				mThresholds.push_back(aThreshold);
			}

//------------------------------------------------------------------------------

			/**
			 *  @brief  set the current mesh
			 *
			 *  @param[in]  aMeshPointer           pointer to mesh object
			 *
			 */
			void set_mesh( mtk::Mesh* & aMeshPointer )
			{
				mListOfMeshs.push_back( aMeshPointer );
			}

//------------------------------------------------------------------------------
    	    /**
    	     * @brief function to create a list of flagged elements
    	     *
    	     * @param[in] aElementList          list of elements from mesh
    	     * @param[in] aConstant				cell of reals which are relevant to the geometry function (e.g. radius and center value for circle)
    	     * @param[in] aWhichGeometry		switch to determine which geometry function is being used from the list of geometries\
    	     *
    	     * @param[out] tRefineFlags			logical cell of 1s and 0s specifying which elements from the given list are flagged for refinement
    	     *
    	     */
			moris::Cell< uint >
			flag_element_list_for_refinement( moris::Cell< mtk::Cell* > & aElementList ,
											  moris::Cell< real >       & aConstant, // constants relative to LS function
											  uint                      aWhichGeometry )
			{
				uint mNumberOfElements = aElementList.size(); //number of elements
				moris::Cell< uint > tRefineFlags(mNumberOfElements); //flag list to be returned

				for(uint i=0; i<mNumberOfElements; i++ )
				{
					moris::uint j = aElementList(i)->get_number_of_vertices(); // number of nodes in the element

					Matrix< DDRMat > mFlag(j,1);
					for( uint k=0; k<j; k++ )
					{
						mFlag(k,0) = mListOfGeoms(aWhichGeometry)->get_val_at_vertex( aElementList(i)->get_vertex_pointers()(k)->get_coords() , aConstant );
					}

					if ( mFlag.max() > mThresholds(aWhichGeometry) && mFlag.min() < mThresholds(aWhichGeometry) )
					{
						tRefineFlags(i) = 1;
					}
					else
					{
						tRefineFlags(i) = 0;
					}
				}
				return tRefineFlags;
			};
//------------------------------------------------------------------------------
    	    /**
    	     * @brief similar to flag_element_list_for_refinement(), has the ability to check over multiple meshs
    	     *
    	     * @param[in] aConstant				cell of real consant relevant to the LS function (analytical)
    	     * @param[in] aWhichGeometry		switch to determine which geometry function is being used from the list of geometries
    	     * @param[in] aWhichMesh			switch to specify which mesh is being analyzed from the list of meshs
    	     *
    	     * @param[out] tRefFlags			logical list of elements which are flagged for refinement
    	     *
    	     */
			moris::Cell< uint >
			check_for_intersection( moris::Cell< real > &   aConstant,
										 uint               aWhichGeometry,
										 uint               aWhichMesh )
			{
				uint mNumElements = mListOfMeshs( 0 )->get_num_elems();
				moris::Cell< mtk::Cell* > aEleList( mNumElements );

				for(uint i=0; i<mNumElements; i++)
				{
					// loop through all elements and build the moris::Cell< mtk::Cell* >
					aEleList(i) = & mListOfMeshs(aWhichMesh)->get_mtk_cell(i);
				}
				moris::Cell< uint > tRefFlags(mNumElements);

				tRefFlags = this->flag_element_list_for_refinement( aEleList,
													                aConstant,
													                aWhichGeometry);
				return tRefFlags;
 		    };

//------------------------------------------------------------------------------
    	    /**
    	     * @brief function to get edge normal for 2D quad element
    	     *
    	     * @param[in] aElemGlobInd          global index of the element containing the edge
    	     * @param[in] aEdgeSideOrd			index of side ordinal local to the element, indices start at 0 on the bottom and continue CCW
    	     * @param[in] aMeshPointer			pointer to mesh object containing the element/edge
    	     * @param[out] tNormal			    edge normal
    	     *
    	     */
			Matrix< DDRMat >
			get_edge_normal_for_straight_edge_quad4( 	uint const & aElemGlobInd,
														uint const & aEdgeSideOrd,
														mtk::Mesh* & aMeshPointer 	)
			{
				Matrix< IdMat > tEdgesOnElem = aMeshPointer->get_edges_connected_to_element_glob_ids( aElemGlobInd );
				Matrix< IdMat > tNodesOnElem = aMeshPointer->get_nodes_connected_to_element_glob_ids( aElemGlobInd );
				MORIS_ASSERT( tEdgesOnElem.numel() == 4, "get_edge_normal_for_straight_edge_quad4() is only valid for a 2D Quad4 element");

				Matrix< DDRMat > tNodeCoord0( 1, 2 );
				Matrix< DDRMat > tNodeCoord1( 1, 2 );
				Matrix< DDRMat > tVec( 1, 2 );
				Matrix< DDRMat > tNormal( 1, 2 );

				switch ( aEdgeSideOrd )
				{
				case( 0 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
					tVec = tNodeCoord0 - tNodeCoord1;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
				case( 1 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
					tVec = tNodeCoord1 - tNodeCoord0;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
				case( 2 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
					tVec = tNodeCoord0 - tNodeCoord1;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
					break;
				}
				case( 3 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
					tVec = tNodeCoord1 - tNodeCoord0;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
                default		:
                {
                    MORIS_ERROR( false, "unknown side ordinal rank");
                    break;
                }
				}

				return tNormal;
			};

//------------------------------------------------------------------------------

			void
			seed_sphere_function( moris::Cell< mtk::Vertex * > & aPoints,
			                      moris::Matrix< DDRMat >      & aNumSphereList,
			                      real                         & aRadius)
			{
			    Matrix< DDRMat > tMinCoord;
			    tMinCoord = aPoints(0)->get_coords();

			    Matrix< DDRMat > tMaxCoord;
			    tMaxCoord = aPoints(1)->get_coords();

			    double tXDim = tMaxCoord(0,0) - tMinCoord(0,0);
			    double tYDim = tMaxCoord(0,1) - tMinCoord(0,1);
			    double tZDim = tMaxCoord(0,2) - tMinCoord(0,2);

			    std::cout<<"x dim = "<<tXDim<<std::endl;     std::cout<<"y dim = "<<tYDim<<std::endl;

			    double tXGap = tXDim/(aNumSphereList(0,0) - 1);
			    double tYGap = tYDim/(aNumSphereList(1,0) - 1);
			    double tZGap = tZDim/(aNumSphereList(2,0) - 1);

			    //-----generate the spheres-----
			    Matrix< DDRMat > tSphereVals(aNumSphereList(0), 1);
			    tSphereVals.fill( 0 );

			    // spheres in x-direction
			    for ( uint i=0; i<uint(aNumSphereList(1)); i++)
			    {
			        Cell< real > tSphereInputs(4);
			        tSphereInputs(0) = tXGap*i;
			        tSphereInputs(1) = tYGap*i;
			        tSphereInputs(2) = tZGap*i;
			        tSphereInputs(3) = aRadius;

			        tSphereVals(i)  = this->sphere_function( tMinCoord, tSphereInputs );
			        print(tSphereInputs, "inputs");
			    }

			    std::cout<<"xGap = "<<tXGap<<std::endl;
			    std::cout<<"yGap = "<<tYGap<<std::endl;

//			    print(tSphereVals, "vals");
			    print(aNumSphereList, "list");

			}

//------------------------------------------------------------------------------

            void
            find_cells_within_levelset(
                          Cell< mtk::Cell * >      & aCells,
                          Cell< mtk::Cell * >      & aCandidates,
                          const  Matrix< DDRMat >  & aVertexValues,
                   const              uint        aUpperBound = 0.0 )
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

            void
            find_cells_intersected_by_levelset(
                    Cell< mtk::Cell * > & aCells,
                    Cell< mtk::Cell * > & aCandidates,
                    const  Matrix< DDRMat >   & aVertexValues,
                    const              real      aLowerBound = -0.0001,
                    const              real      aUpperBound =  0.0001)
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

		};	/* ge class */
	}	/* ge namespace */
}	/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_HPP_ */
