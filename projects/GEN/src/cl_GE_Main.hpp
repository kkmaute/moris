/*
 * cl_GE_Main.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_HPP_

#include "cl_GE_Geometry.hpp"
#include "cl_Cell.hpp"


namespace moris
{
	namespace ge
	{
		class GE
		{
			//------------------------------------------------------------------------------

			// member variables
			moris::Cell< Geometry* > mListOfGeoms;
			moris::Cell< mtk::Mesh* > mListOfMeshs;
			moris::Cell< real > mThresholds;

			//------------------------------------------------------------------------------
		private:
			//------------------------------------------------------------------------------
		protected:
			//------------------------------------------------------------------------------
		public:
			GE(){};

			void set_geometry( Geometry* geomPointer )
			{
				mListOfGeoms.push_back( geomPointer );
			}

			//------------------------------------------------------------------------------
			void set_mesh( mtk::Mesh* meshPointer )
			{
				mListOfMeshs.push_back( meshPointer );
			}

			//------------------------------------------------------------------------------
			void
			set_threshold( moris::Cell< double > aThreshVals )
			{
				mThresholds.append(aThreshVals);
			}

			//------------------------------------------------------------------------------
    	    /**
    	     * @brief function to create a list of flagged elements
    	     *
    	     * @param[in] aElementList          list of elements from mesh
    	     * @param[in] aConstant				cell of reals which are relevant to the geometry function (e.g. radius and center value for circle)
    	     * @param[in] aWhichGeometry		switch to determine which geometry function is being used from the list of geometries
    	     * @param[in] aWhichThreshold		switch to determine which threshold to use for the LS function, from the list of thresholds
    	     * @param[out] tRefineFlags			logical cell of 1s and 0s specifying which elements from the given list are flagged for refinement
    	     *
    	     */
			moris::Cell< uint >
			flag_element_list_for_refinement( moris::Cell< mtk::Cell* > & aElementList ,
											  moris::Cell< real > & aConstant, // constants relative to LS function
											  uint aWhichGeometry,
											  uint aWhichThreshold )
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

					if ( mFlag.max() > mThresholds(aWhichThreshold) && mFlag.min() < mThresholds(aWhichThreshold) )
					{
//						std::cout<<"element "<<i+1<<" flagged for refinement"<<std::endl;
						tRefineFlags(i) = 1;
					}
					else
					{
//						std::cout<<"no need to refine element "<<i+1<<std::endl;
						tRefineFlags(i) = 0;
					}
				}
				return tRefineFlags;
			};
			//------------------------------------------------------------------------------
    	    /**
    	     * @brief similar to flag_element_list_for_refinement(), in fact uses said function, but has the ability to check over multiple meshs
    	     *
    	     * @param[in] aConstant				cell of real consant relevant to the LS function (analytical)
    	     * @param[in] aWhichGeometry		switch to determine which geometry function is being used from the list of geometries
    	     * @param[in] aWhichMesh			switch to specify which mesh is being analyzed from the list of meshs
    	     * @param[in] aWhichThreshold		switch to determine which threshold to use for the LS function, from the list of thresholds
    	     * @param[out] tRefFlags			logical list of elements which are flagged for refinement
    	     *
    	     */
			moris::Cell< uint >
			check_for_intersection( moris::Cell< real > & aConstant,
										 uint aWhichGeometry,
										 uint aWhichMesh,
										 uint aWhichThreshold)
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
														aWhichGeometry,
														aWhichThreshold);
				return tRefFlags;
 		    };

			//------------------------------------------------------------------------------
    	    /**
    	     * @brief function to get edge normal for 2D element
    	     *
    	     * @param[in] aElemGlobInd          global index of the element containing the edge
    	     * @param[in] aEdgeSideOrd			index of side ordinal local to the element, indices start at 0 on the bottom and continue CCW
    	     * @param[in] aMeshPointer			pointer to mesh object containing the element/edge
    	     * @param[out] tNormal			    edge normal
    	     *
    	     */
			Matrix< DDRMat >
			get_edge_normal_for_straight_edge( 	uint const & aElemGlobInd,
												uint const & aEdgeSideOrd,
												mtk::Mesh* & aMeshPointer )
			{
				Matrix< IndexMat > tEdgesOnElem = aMeshPointer->get_edges_connected_to_element_glob_ids( aElemGlobInd );
				print(tEdgesOnElem, "edges on element");
				Matrix< IdMat > tNodesOnElem = aMeshPointer->get_nodes_connected_to_element_glob_ids( aElemGlobInd );
				print(tNodesOnElem, "nodes on element");

				Matrix< DDRMat > tNodeCoord0( 1, 2 );
				Matrix< DDRMat > tNodeCoord1( 1, 2 );

				Matrix< DDRMat > tVec( 1, 2 );
				Matrix< DDRMat > tNormal( 1, 2 );

				switch ( aEdgeSideOrd )
				{
				case( 0 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
					print(tNodeCoord0, "coords of node 1");
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
					print(tNodeCoord1, "coords of node 2");
					tVec = tNodeCoord0 - tNodeCoord1;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
				case( 1 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
					print(tNodeCoord0, "coords of node 1");
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
					print(tNodeCoord1, "coords of node 2");
					tVec = tNodeCoord1 - tNodeCoord0;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
				case( 2 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
					print(tNodeCoord0, "coords of node 1");
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
					print(tNodeCoord1, "coords of node 2");
					tVec = tNodeCoord0 - tNodeCoord1;
					tVec(0,0) = tVec(0,0)/norm(tVec);	tVec(0,1) = tVec(0,1)/norm(tVec);
					tNormal(0,0) = tVec(0,1);	tNormal(0,1) = tVec(0,0);
                    break;
				}
				case( 3 )	:
				{
					tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
					print(tNodeCoord0, "coords of node 1");
					tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
					print(tNodeCoord1, "coords of node 2");
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

				print(tNormal, "normal");
				return tNormal;
			};

			//------------------------------------------------------------------------------

		};	/* ge class */
	}	/* ge namespace */
}	/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_HPP_ */
