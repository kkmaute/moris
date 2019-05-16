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

// HMR includes
#include "../projects/HMR/src/cl_HMR_Field.hpp"
//------------------------------------------------------------------------------
// MTK includes
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Mesh_Enums.hpp"
//------------------------------------------------------------------------------
// FEM includes
//#include "../projects/FEM/INT/src/cl_FEM_Integrator.hpp"
//#include "../projects/FEM/INT/src/cl_FEM_Field_Interpolator.hpp"
//------------------------------------------------------------------------------
// LinAlg includes
#include "cl_Matrix.hpp"
#include "fn_linsolve.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"
#include "op_times.hpp"
//------------------------------------------------------------------------------

namespace moris
{
	namespace ge
	{
		class GE_Core
		{
//------------------------------------------------------------------------------
		public:
			GE_Core(){};

			~GE_Core()
			{
			};
//------------------------------------------------------------------------------
			/**
			 *  @brief  set the current geometry and threshold
			 *
			 *  @param[in]  aGeometryPointer        pointer to geometry object
			 *  @param[in]  aThreshold              threshold value for current geometry LS function
			 *
			 */
			void
			set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
                          mtk::Mesh_Manager           & aMyManager,
			              real                          aThreshold = 0.0 )
			{
			    uint tGeomInd = mListOfGeoms.size();    //index of current geometry rep.
			    mListOfGeoms.push_back( aGeomPointer );
			    mThresholds.push_back( aThreshold );

                this->get_geometry_pointer( tGeomInd )->set_mesh( aMyManager );
			}

//------------------------------------------------------------------------------
            /**
             *  @brief  get pointer to specified geometry object
             *
             *  @param[in]   aWhichGeometry                     which geometry pointer to get
             *
             *  @param[out]  mListOfGeoms( aWhichGeometry )     geometry class pointer with index given
             */
			std::shared_ptr< Geometry >
			get_geometry_pointer( uint aWhichGeometry )
			{
			    return mListOfGeoms( aWhichGeometry );
			}

//------------------------------------------------------------------------------
			/*
			 * @brief initialize empty output object
			 *
			 *
			 */
			void
			initialize()
			{
                Output_Object tOutputObject;
                mListOfOutputObjects.push_back( tOutputObject );
			}
//------------------------------------------------------------------------------
			/*
			 * @brief   initialize output object
			 *
			 * @param[in]   aWhichGeometry - which geometry pointer the output object is associated with
			 * @param[in]   aMeshIndex     - associated mesh index for the geometry pointer of interest
			 * @param[in]   aField         - source field on mesh
			 * @param[in]   aTargetLabel   - target field label
			 * @param[in]   aTargetEntityRank - target entity rank
			 * @param[in]   aSourceEntityRank - source entity rank, defaulted to NODE
			 *
			 */
			void
			initialize( const uint        aWhichGeometry,
			            const moris_index aMeshIndex,
			                  std::shared_ptr<moris::hmr::Field>  & aField,
			                  std::string aTargetLabel,
			                  EntityRank  aTargetEntityRank,
			                  EntityRank  aSourceEntityRank = EntityRank::NODE)
			{
                MeshType tMeshType = this->get_geometry_pointer(aWhichGeometry)->get_my_mesh()->get_interpolation_mesh(aMeshIndex)->get_mesh_type();
                MORIS_ASSERT(tMeshType == MeshType::HMR, "GE_Core::initalize(): currently only set up to work with hmr mesh ");

                uint tBSplineOrder = aField->get_bspline_order();

			    enum type tType = this->get_geometry_pointer( aWhichGeometry )->get_geom_type();
			    if (tType == type::DISCRETE)
			    {

			        mapper::Mapper tMapper( this->get_geometry_pointer(aWhichGeometry)->get_my_mesh(),
			                                aMeshIndex,
			                                tBSplineOrder );

			        tMapper.perform_mapping( aField->get_label(),
			                                 aSourceEntityRank,
			                                 aTargetLabel,
			                                 aTargetEntityRank );
			        Output_Object tOutputObject( tMapper );
			        mListOfOutputObjects.push_back( tOutputObject );
			    }

			}

//------------------------------------------------------------------------------
			/**
			 *  @brief  get output object with all the queried information
			 *
			 *  @param[in]   aWhichGeometry                     output object associated with this geometry rep
			 *
			 *  @param[out]  mListOfOutputObjects( aWhichGeometry )     pointer to output object
			 */
			Output_Object
			get_output_object_pointer( uint aWhichGeometry )
			{
			    return mListOfOutputObjects( aWhichGeometry );
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
				}

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
			    }

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
			    case( 0 )   :
													         {
			        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
			        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
			        tVec = tNodeCoord0 - tNodeCoord1;
			        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
			        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
			        break;
													         }
			    case( 1 )   :
													         {
			        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
			        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
			        tVec = tNodeCoord1 - tNodeCoord0;
			        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
			        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
			        break;
													         }
			    case( 2 )   :
													         {
			        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
			        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
			        tVec = tNodeCoord0 - tNodeCoord1;
			        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
			        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
			        break;
													         }
			    case( 3 )   :
													         {
			        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
			        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
			        tVec = tNodeCoord1 - tNodeCoord0;
			        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
			        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
			        break;
													         }
			    default     :
			    {
			        MORIS_ERROR( false, "unknown side ordinal rank");
			        break;
			    }
			    }

			    return tNormal;
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


        //------------------------------------------------------------------------------
		private:
            // member variables
            moris::Cell< std::shared_ptr< Geometry > >      mListOfGeoms;
            moris::Cell< real >                             mThresholds;
            moris::Cell< Output_Object >                    mListOfOutputObjects;

        //------------------------------------------------------------------------------
        protected:

//------------------------------------------------------------------------------

		};	/* GE_Main class */
	}	/* ge namespace */
}	/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_CORE_HPP_ */
