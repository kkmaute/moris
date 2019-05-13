/*
 * cl_GE_Core.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_CORE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_CORE_HPP_

// moris includes
#include "cl_Cell.hpp"
//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Geometry.hpp"
#include "cl_GE_Output_Object.hpp"
//------------------------------------------------------------------------------
// MTK includes
#include "cl_MTK_Mesh_Core.hpp"
//------------------------------------------------------------------------------
// FEM includes
#include "cl_FEM_Integrator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
//------------------------------------------------------------------------------
// LinAlg includes
#include "fn_linsolve.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"
#include "op_times.hpp"
#include "cl_Matrix.hpp"
//------------------------------------------------------------------------------

namespace moris
{
	namespace ge
	{
		class GE_Core : public Geometry_Engine_Interface
		{
//------------------------------------------------------------------------------
		public:
			GE_Core(){};

			/*
			 * initialize GE with a geometry pointer
			 */
			GE_Core( std::shared_ptr< Geometry > & aGeomPointer )
			{
			    this->set_geometry( aGeomPointer );
			};

			~GE_Core(){};

			/**
			 *  @brief  set the current geometry and threshold
			 *
			 *  @param[in]  aGeometryPointer        pointer to geometry object
			 *  @param[in]  aThreshold              threshold value for current geometry LS function
			 *
			 */
			void
			set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
			                                     real   aThreshold   = 0.0 );

//------------------------------------------------------------------------------
            /**
             *  @brief  get pointer to specified geometry object
             *
             *  @param[in]   aWhichGeometry                     which geometry pointer to get
             *
             *  @param[out]  mListOfGeoms( aWhichGeometry )     geometry class pointer
             */
			std::shared_ptr< Geometry >
			get_geometry_pointer( uint aWhichGeometry );

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
											  uint                      aWhichGeometry );

//------------------------------------------------------------------------------
    	    /**
    	     * @brief loops through a given mesh a returns a logical list of elements which are intersected with the specified geometry
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
			                        uint                    aWhichGeometry);

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
			get_edge_normal_for_straight_edge_quad4( uint const & aElemGlobInd,
													 uint const & aEdgeSideOrd,
													 mtk::Mesh* & aMeshPointer 	);

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
                          const  uint                aUpperBound = 0.0 );

//------------------------------------------------------------------------------
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
                    const  real                   aUpperBound =  0.0001);

//------------------------------------------------------------------------------

//fixme need to add default values for parts of this function (such as integration/interpolation rules, time evolution domain, etc.)
//fixme need to make geometry class robust enough so that this function knows how to get the field values at the nodes (for loop over nodes)

            /*@brief function to determine ADV's (phi) at the nodal locations
             *
             * @param[in]    aWhichGeom - which geometry representation from member list of geoms
             * @param[in]    aGeomParams - constant needed for specific geometry
             * @param[in]    aGeomType - elemental geometry type (e.g. QUAD)
             * @param[in]    aIntegType - integration type for the geom rep
             * @param[in]    aIntegOrder - integration order for the geom rep
             * @param[in]    aInterpType - interpolation type for the geom rep (e.g. LAGRANGE)
             * @param[in]    aInterpOrder - interpolation order for the geom rep (e.g. LINEAR)
             *
             * @param[out]   tADVs      - Bspline coefficient vector
             */
            Matrix< DDRMat >
            compute_nodal_advs( uint                          aWhichGeom,
                                moris::Cell< real >         & aGeomParams,
                                mtk::Geometry_Type            aGeomType,
                                fem::Integration_Type         aIntegType,
                                fem::Integration_Order        aIntegOrder,
                                fem::Interpolation_Type       aInterpType,
                                mtk::Interpolation_Order      aInterpOrder,
                                fem::Integration_Type         aTimeIntegType = fem::Integration_Type::GAUSS,
                                fem::Integration_Order        aTimeIntegOrder = fem::Integration_Order::BAR_1,
                                fem::Interpolation_Type       aTimeInterpTypeGeom = fem::Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order      aTimeInterpOrderGeom = mtk::Interpolation_Order::LINEAR,
                                fem::Interpolation_Type       aTimeInterpTypeTime = fem::Interpolation_Type::CONSTANT,
                                mtk::Interpolation_Order      aTImeInterpOrderTime = mtk::Interpolation_Order::CONSTANT,
                                uint                          aNumberOfFields = 1,
                                Matrix< DDRMat >      const & aTHat = Matrix< DDRMat >({{0},{1}}) )
            {
                // create a space time integration rule and interpolation rules
                //------------------------------------------------------------------------------
                // set up the integration rule to be used
                fem::Integration_Rule tFieldIntegRule( aGeomType,
                                                       aIntegType,
                                                       aIntegOrder,
                                                       aTimeIntegType,         // these last two are for the time integration
                                                       aTimeIntegOrder );      // using BAR_1 as default since we are not evolving in time

                fem::Integrator tFieldIntegrator( tFieldIntegRule );                        // field integration object

                // interpolation rule for geometry
                fem::Interpolation_Rule tGeomInterpRule( aGeomType,
                                                         aInterpType,
                                                         aInterpOrder,
                                                         aTimeInterpTypeGeom,      // time integration
                                                         aTimeInterpOrderGeom );   // time integration

                fem::Geometry_Interpolator* tGeomInterpolator = new fem::Geometry_Interpolator( tGeomInterpRule );

                // interpolation rule for the LS field
                fem::Interpolation_Rule tFieldRule( aGeomType,
                                                    aInterpType,
                                                    aInterpOrder,
                                                    aTimeInterpTypeTime,          // time integration (constant since we are not evolving in time)
                                                    aTImeInterpOrderTime );       // time integration (constant since we are not evolving in time)

                // field interpolator object
                fem::Field_Interpolator tFieldInterpolator( aNumberOfFields,
                                                            tFieldRule,
                                                            tGeomInterpolator );

                // get GPs and weights from integrator
                Matrix< DDRMat > tSpaceIntegPoints  = tFieldIntegrator.get_points();
                Matrix< DDRMat > tSpaceIntegWeights = tFieldIntegrator.get_weights();
                // evaluate LS at GP and Nodes
                //------------------------------------------------------------------------------
                uint tNumNodes = mListOfGeoms(aWhichGeom)->get_my_mesh()->get_num_nodes();    // determine number of nodes in mesh
                Matrix< DDRMat > tTempCoord(1,2);                       // temporary coordinate used to fill coordinate matrix
                Matrix< DDRMat > tXHat(tNumNodes,2);                    // global coordinate matrix

                for (uint j=0; j<tNumNodes; j++)                        // create global coordinate matrix
                {
                    tTempCoord = mListOfGeoms(aWhichGeom)->get_my_mesh()->get_node_coordinate(j);
                    uint k = 0;

                    tXHat(j,k)   = tTempCoord(0,0);
                    tXHat(j,k+1) = tTempCoord(0,1);
                }

//                // create interpolation object
//                fem::Interpolation_Function_Base* tFunction = tGeomInterpRule.create_space_interpolation_function();
                tGeomInterpolator->set_coeff( tXHat, aTHat );       //set the coefficients xHat:(x,y), tHat:(t_0,t_f)

                Matrix< DDRMat > tFieldValAtNodes(1,tNumNodes);
                Matrix< DDRMat > tPoint;                            // temporary argument to be passed into LS function
                Matrix< DDRMat > tPhiS(1,tNumNodes);                // field values at the Gauss Points

                for (uint i=0; i<tNumNodes; i++)
                {
                    if(mListOfGeoms(aWhichGeom)->is_analytic())     // right now only implemented for analytic and discrete geometry classes
                    {
                        tFieldValAtNodes(i) = mListOfGeoms(aWhichGeom)->get_field_val_at_coordinate( mListOfGeoms(aWhichGeom)->get_my_mesh()->get_node_coordinate(i) , aGeomParams );
                        //evaluate space and time at xi, tau
                        tPoint   = tGeomInterpolator->valx( tSpaceIntegPoints.get_column(i) );    // N * ( x,y ) (determines global coordinates of the GPs)

                        tPhiS(i) = mListOfGeoms(aWhichGeom)->get_field_val_at_coordinate( tPoint, aGeomParams );
                    }
                    else if(!mListOfGeoms(aWhichGeom)->is_analytic())
                    {
                        tFieldValAtNodes(i) = mListOfGeoms(aWhichGeom)->access_field_value_with_entity_index( i, EntityRank::NODE );

                        tPoint   = tGeomInterpolator->valx( tSpaceIntegPoints.get_column(i) );    // N * ( x,y ) (determines global coordinates of the GPs)

                        tPhiS(i) = mListOfGeoms(aWhichGeom)->access_field_value_with_entity_index( i, EntityRank::NODE );
                    }
                }
                //------------------------------------------------------------------------------
                // compute nodal advs from L2 projection
                //------------------------------------------------------------------------------
                Matrix< DDRMat > tADVs(tNumNodes,1, 1);

                Matrix< DDRMat > tPhi(tNumNodes,1, 1);

                real tTol = 1.0;

                for (uint k=0; k<10; k++)   // N-R loop
                {
                    Matrix< DDRMat > tJacobian(tNumNodes,tNumNodes);
                    tJacobian.fill( 0 );

                    Matrix< DDRMat > tResidual(mListOfGeoms(aWhichGeom)->get_my_t_matrix().n_cols(), 1, 0.0);

                    for(uint i=0; i<tNumNodes; i++)
                    {
                        tFieldInterpolator.set_space_time( tSpaceIntegPoints.get_column(i)); // sets up the integration points
                        Matrix< DDRMat > tN = tFieldInterpolator.N(); // create the shape function vector

                        real tMultFactor = tGeomInterpolator->det_J( tSpaceIntegPoints.get_column(i) ) * tSpaceIntegWeights(i); // multiplication factor for Jacobian and Residual
                        tJacobian        = tJacobian + tMultFactor*trans(tN)*tN;

                        tResidual        = tResidual + tMultFactor*trans(tN)*( tN*tPhi - tN*trans(tFieldValAtNodes));
                    }
                    Matrix< DDRMat > tJBar = trans(mListOfGeoms(aWhichGeom)->get_my_t_matrix())*tJacobian*mListOfGeoms(aWhichGeom)->get_my_t_matrix();  // in B-spline space (Tmatrix)
                    Matrix< DDRMat > tRBar = trans(mListOfGeoms(aWhichGeom)->get_my_t_matrix())*tResidual;                              // in B-spline space (Tmatrix)
                    if (k==0)
                    {
                        tTol = norm(tRBar*0.000001);
                    }

                    if (norm(tRBar)<tTol)
                    {
                        break;
                    }
                    tADVs = tADVs - solve(tJBar,tRBar);  // solve Ax=B problem
                    tPhi    = mListOfGeoms(aWhichGeom)->get_my_t_matrix()*tADVs;              // determine tPhi for iteration loop
                }
                delete tGeomInterpolator;
                return tADVs;
            }

        //------------------------------------------------------------------------------
		private:
            // member variables
            moris::Cell< std::shared_ptr< Geometry > > mListOfGeoms;
            moris::Cell< real >                        mThresholds;

            /*
             * should also store refinement list for certain geometry (if asked)
             */
        //------------------------------------------------------------------------------
        protected:

//------------------------------------------------------------------------------

		};	/* GE_Main class */
	}	/* ge namespace */
}	/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_CORE_HPP_ */
