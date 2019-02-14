/*
 * cl_FEM_Element_Proxy.hpp
 *
 *  Created on: Jan 23, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENT_PROXY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENT_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp" //MTK/src

namespace moris
{
namespace fem
{

	class ElementProxy : public mtk::Cell
	{
	private:
		moris::Cell< mtk::Vertex* > 	mNodeList;
		enum mtk::Geometry_Type 		mGeometryType = mtk::Geometry_Type::UNDEFINED;

	public:
		ElementProxy(moris::Cell< mtk::Vertex* > aNodeList, enum mtk::Geometry_Type aGeometryType)
		{
			mNodeList = aNodeList;
			mGeometryType = aGeometryType;
		}

		~ElementProxy(){};

		moris_id get_id() const
		{
			MORIS_ERROR( false,"get_id() not implemented " );
			moris_id tId = 0;
			return tId;
		};

		//------------------------------------------------------------------------------
		moris_index get_index() const
		{
			MORIS_ERROR( false,"get_index() not implemented " );
			moris_index mIndex = 0;
			return mIndex;
		};

		//------------------------------------------------------------------------------
		uint get_number_of_vertices() const
		{
			return mNodeList.size();
		};


		//------------------------------------------------------------------------------
		moris_id get_owner() const
		{
			MORIS_ERROR( false,"get_owner() not implemented " );
		    moris_id tId = 0;
		    return tId;
		};

		//------------------------------------------------------------------------------
		moris::Cell< mtk::Vertex* > get_vertex_pointers() const
		{
			return mNodeList;
		};

		//------------------------------------------------------------------------------
		Matrix< IdMat > get_vertex_ids() const
		{
			MORIS_ERROR( false, "get_vertex_ids() not available for this element " );
			return Matrix< IdMat > (0,0);
		};

		//------------------------------------------------------------------------------
		Matrix< IndexMat > get_vertex_inds() const
		{
			MORIS_ERROR( false, "get_vertex_ids() not available for this element " );
			return Matrix< IndexMat > (0,0);
		};

		//------------------------------------------------------------------------------
		Matrix< DDRMat > get_vertex_coords() const
		{
			MORIS_ERROR( false, "get_vertex_coords() not available for this element " );
			return Matrix< DDRMat > (0,0);
		};

		//------------------------------------------------------------------------------
		mtk::Geometry_Type get_geometry_type() const
		{
			return mGeometryType;
		};

		//------------------------------------------------------------------------------
		mtk::Interpolation_Order get_interpolation_order() const
		{
			MORIS_ERROR( false, "get_interpolation_order() not available for this element " );
			return mtk::Interpolation_Order::UNDEFINED;
		};

	};

}/* namespace fem */
}/* namespace moris */



#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENT_PROXY_HPP_ */
