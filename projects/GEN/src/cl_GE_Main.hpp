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

// FEM includes
#include "FEM/INT/src/cl_FEM_Integrator.hpp"
#include "FEM/INT/src/cl_FEM_Field_Interpolator.hpp"

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
			moris::Cell< real > 			mThresholds;

//------------------------------------------------------------------------------
		private:


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
						mFlag(k,0) = mListOfGeoms(aWhichGeometry)->get_field_val_at_coordinate( aElementList(i)->get_vertex_pointers()(k)->get_coords() , aConstant );
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
    	     * @brief similar to flag_element_list_for_refinement()
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
			                        mtk::Mesh*          &   aMeshPointer,
			                        uint                    aWhichGeometry)
			{
				uint mNumElements = aMeshPointer->get_num_elems();
				moris::Cell< mtk::Cell* > aEleList( mNumElements );

				for(uint i=0; i<mNumElements; i++)
				{
					// loop through all elements and build the moris::Cell< mtk::Cell* >
					aEleList(i) = & aMeshPointer->get_mtk_cell(i);
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
			/*
			 * @brief determines volume elements (cells)
			 *
			 *
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
            /*
             * @brief determines elements (cells) intersected by the level set
             *
             *
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
//fixme need to add default values for parts of this function (such as integration/interpolation rules, time evolution domain, etc.)

            /*@brief function to determine adof values (phi) at the nodal locations
             *
             * @param[in]    aTMatrix     - TMatrix to be used
             * @param[in]    aMeshPointer - pointer to mtk mesh
             * @param[in]    aGeomParams  - constant needed for specific geometry
             * @param[in]    aWhichGeom   - which geometry to analyze
             *
             * @param[out]   tPhiHat      - Bspline coefficient vector
             */
            Matrix< DDRMat >
            determine_phi_hat( Matrix< DDRMat >      & aTMatrix,
                                  mtk::Mesh*            & aMeshPointer,
                                  moris::Cell< real >   & aGeomParams,
                                  uint                  aWhichGeom )
            {
                // create a space time integration rule and interpolation rules for Geometry and LS field
                //------------------------------------------------------------------------------
                // set up the integration rule to be used
                fem::Integration_Rule tFieldIntegRule( mtk::Geometry_Type::QUAD,
                                                       fem::Integration_Type::GAUSS,
                                                       fem::Integration_Order::QUAD_3x3,
                                                       fem::Integration_Type::GAUSS,        // these last two are for the time integration
                                                       fem::Integration_Order::BAR_1);      // using BAR_1 as defailt since we are not evolving in time

                fem::Integrator tFieldIntegrator( tFieldIntegRule );                        // field integration object

                // interpolation rule for geometry
                fem::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                         fem::Interpolation_Type::LAGRANGE,
                                                         mtk::Interpolation_Order::QUADRATIC,
                                                         fem::Interpolation_Type::LAGRANGE,    // time integration
                                                         mtk::Interpolation_Order::LINEAR );   // time integration

                fem::Geometry_Interpolator* tGeomInterpolator = new fem::Geometry_Interpolator( tGeomInterpRule );

                // interpolation rule for the LS field
                fem::Interpolation_Rule tFieldRule( mtk::Geometry_Type::QUAD,
                                                    fem::Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::QUADRATIC,
                                                    fem::Interpolation_Type::CONSTANT,          // time integration (constant since we are not evolving in time)
                                                    mtk::Interpolation_Order::CONSTANT );       // time integration (constant since we are not evolving in time)

                uint tNumberOfFields = 1;
                // field interpolator object
                fem::Field_Interpolator tFieldInterpolator( tNumberOfFields,
                                                            tFieldRule,
                                                            tGeomInterpolator );

                // get GPs and weights from integrator
                Matrix< DDRMat > tSpaceIntegPoints   = tFieldIntegrator.get_points();
                Matrix< DDRMat > tSpaceIntegWeights  = tFieldIntegrator.get_weights();

                // evaluate LS at GP and Nodes
                //------------------------------------------------------------------------------
                uint tNumNodes = aMeshPointer->get_num_nodes();     // determine number of nodes in mesh
                Matrix< DDRMat > tTempCoord(1,2);                   // temporary coordinate used to fill coordinate matrix
                Matrix< DDRMat > tXHat(tNumNodes,2);                // global coordinate matrix

                for (uint j=0; j<tNumNodes; j++)                    // create global coordinate matrix
                {
                    tTempCoord = aMeshPointer->get_node_coordinate(j);
                    uint k = 0;

                    tXHat(j,k)   = tTempCoord(0,0);
                    tXHat(j,k+1) = tTempCoord(0,1);
                }

//                // create interpolation object
//                fem::Interpolation_Function_Base* tFunction = tGeomInterpRule.create_space_interpolation_function();

                //create time evolution domain (dummy in this case since we are not evolving in time)
                Matrix< DDRMat > tTHat( 2, 1 );
                tTHat( 0 ) = 0.0;
                tTHat( 1 ) = 1.0;

                tGeomInterpolator->set_coeff( tXHat, tTHat );       //set the coefficients xHat:(x,y), tHat:(t_0,t_f)

                Matrix< DDRMat > tPhiSAtNodes(1,tNumNodes);
                Matrix< DDRMat > tPoint;                            // temporary argument to be passed into LS function
                Matrix< DDRMat > tPhiS(1,tNumNodes);

                for (uint i=0; i<tNumNodes; i++)
                {
                    tPhiSAtNodes(i) = mListOfGeoms(aWhichGeom)->get_field_val_at_coordinate( aMeshPointer->get_node_coordinate(i) , aGeomParams );

                    //evaluate space and time at xi, tau
                    tPoint   = tGeomInterpolator->valx( tSpaceIntegPoints.get_column(i) );    // N * ( x,y ) (determines global coordinates of the GPs)
                    tPhiS(i) = mListOfGeoms(aWhichGeom)->get_field_val_at_coordinate( tPoint, aGeomParams );
                }
                //------------------------------------------------------------------------------

                // compute phi_hat
                //------------------------------------------------------------------------------

                Matrix< DDRMat > tPhiHat(tNumNodes,1);
                tPhiHat.fill( 1 );
                Matrix< DDRMat > tPhi(tNumNodes,1);
                tPhi.fill( 1 );

                real tTol = 1.0;

                for (uint k=0; k<10; k++)
                {
                    Matrix< DDRMat > tJacobian(tNumNodes,tNumNodes);
                    tJacobian.fill( 0 );
                    Matrix< DDRMat > tResidual(9,1);
                    tResidual.fill( 0 );

                    for(uint i=0; i<tNumNodes; i++)
                    {
                        tFieldInterpolator.set_space_time( tSpaceIntegPoints.get_column(i)); // sets up the integration points
                        Matrix< DDRMat > tN = tFieldInterpolator.N(); // create the shape function vector

                        real tMultFactor = tGeomInterpolator->det_J( tSpaceIntegPoints.get_column(i) ) * tSpaceIntegWeights(i); // multiplication factor for Jacobian and Residual
                        tJacobian        = tJacobian + tMultFactor*trans(tN)*tN;
                        tResidual        = tResidual + tMultFactor*trans(tN)*( tN*tPhi - tN*trans(tPhiSAtNodes));
                    }
                    Matrix< DDRMat > tJBar = trans(aTMatrix)*tJacobian*aTMatrix;  // in B-spline space (Tmatrix)
                    Matrix< DDRMat > tRBar = trans(aTMatrix)*tResidual;  // in B-spline space (Tmatrix)
                    if (k==1)
                    {
                        tTol = norm(tRBar*0.000001);
                    }

                    if (norm(tRBar)<tTol)
                    {
                        break;
                    }
                    tPhiHat = tPhiHat - inv(tJBar)*tRBar;  // solve Ax=B problem
                    tPhi    = aTMatrix*tPhiHat;            // determine tPhi for iteration loop
                }

                delete tGeomInterpolator;
                return tPhiHat;
            }

//------------------------------------------------------------------------------

		};	/* ge class */
	}	/* ge namespace */
}	/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_HPP_ */
