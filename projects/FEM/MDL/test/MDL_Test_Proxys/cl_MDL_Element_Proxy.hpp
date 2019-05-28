/*
 * cl_GE_Element.hpp
 *
 *  Created on: Jan 2, 2019
 *      Author: sonne
 */
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MDL_Node_Proxy.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

#ifndef PROJECTS_GEN_SRC_CL_ELEMENT_HPP_
#define PROJECTS_GEN_SRC_CL_ELEMENT_HPP_

namespace moris{
	namespace ge{

	class Element : public mtk::Cell
	{
//		//------------------------------------------------------------------------------
//	private:
//
//		moris::Cell< mtk::Vertex* > mListOfNodes;
//
//		//------------------------------------------------------------------------------
//	public:
//		//------------------------------------------------------------------------------
//
//		Element( moris::Cell< mtk::Vertex* > aListOfNodes )
//		{
//			mListOfNodes = aListOfNodes;
//		};
//		~Element(){};
//
//		//------------------------------------------------------------------------------
//		// fills a moris::cell with pointers to connected nodes
//		moris::Cell< mtk::Vertex* >
//		get_vertex_pointers() const
//	    {
//			return mListOfNodes;
//		};
//
//		//------------------------------------------------------------------------------
//		// gives the number of nodes on the element
//		uint get_number_of_vertices() const
//		{
//			return mListOfNodes.size();
//		}
//
//		//------------------------------------------------------------------------------
//		// domain wide id of the cell
//        moris_id
//        get_id() const
//        {
//        	MORIS_ERROR( false,"get_id() not implemented " );
//        	moris_id tId = 0;
//            return tId;
//        };
//
//		//------------------------------------------------------------------------------
//        // local index of the cell
//        moris_index
//		get_index() const
//        {
//        	MORIS_ERROR( false, "get_index() not implemented " );
//        	moris_index mIndex = 0;
//        	return mIndex;
//        };
//
//		//------------------------------------------------------------------------------
//        // return the proc id of the owner of this element
//        moris_id
//		get_owner() const
//        {
//        	MORIS_ERROR( false, "get_owner() not available for this element " );
//        	moris_id tId = 0;
//        	return tId;
//        };
//
//		//------------------------------------------------------------------------------
//        // returns matrix with IDs of connected nodes
//        Matrix< IdMat >
//        get_vertex_ids() const
//		{
//        	MORIS_ERROR( false, "get_vertex_ids() not available for this element " );
//        	return Matrix< IdMat > (0,0);
//		};
//
//		//------------------------------------------------------------------------------
//        // returns a matrix with indices of connected vertices
//        Matrix< IndexMat >
//        get_vertex_inds() const
//		{
//        	MORIS_ERROR( false, "get_vertex_inds() not available for this element " );
//        	return Matrix< IndexMat > (0,0);
//		};
//
//		//------------------------------------------------------------------------------
//        // returns matrix of dimension: < # indices * # dimensions >
//        Matrix< DDRMat >
//        get_vertex_coords() const
//		{
//        	MORIS_ERROR( false, "get_vertex_coords() not available for this element " );
//        	return Matrix< DDRMat > (0,0);
//		};
//
//		//------------------------------------------------------------------------------
//        // return an enum which defines the geometry type of the element
//        mtk::Geometry_Type
//		get_geometry_type() const
//        {
//        	MORIS_ERROR( false, "get_geometry_type() not available for this element " );
//        	return mtk::Geometry_Type::UNDEFINED;
//        };
//
//		//------------------------------------------------------------------------------
//        // returns the order of the element
//        mtk::Interpolation_Order
//		get_interpolation_order() const
//        {
//        	MORIS_ERROR( false, "get_interpolation_order() not available for this element " );
//        	return mtk::Interpolation_Order::UNDEFINED;
//        };
//
//		//------------------------------------------------------------------------------
//        // special hmr function
//        luint
//		get_memory_index_of_background_element() const
//        {
//        	MORIS_ERROR( false, "get_memory_index_of_background_element() not available for this element " );
//        	return 0;
//        };

	protected:
	};

	}/* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_ELEMENT_HPP_ */
