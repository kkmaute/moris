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

			//create member variables
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
				mThresholds.append(aThreshVals);	// threshold for the level set field
			}

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

					Matrix< DDRMat > flag(j,1);
					for( uint k=0; k<j; k++ )
					{
						flag(k,0) = mListOfGeoms(aWhichGeometry)->get_val_at_vertex( aElementList(i)->get_vertex_pointers()(k)->get_coords() , aConstant );
					}

					if ( flag.max() > mThresholds(aWhichThreshold) && flag.min() < mThresholds(aWhichThreshold) )
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

		};
	}/* gen namespace */
}/* moris namespace */




#endif /* PROJECTS_GEN_SRC_CL_GE_HPP_ */
