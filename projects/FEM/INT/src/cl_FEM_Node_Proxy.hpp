/*
 * cl_FEM_Node_Proxy.hpp
 *
 *  Created on: Jan 23, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_NODE_PROXY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_NODE_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
namespace fem
{

	class NodeProxy : public mtk::Vertex
	{
	private:
		moris::Matrix< DDRMat > coord;
		sint mindex = -1;

	public:
		NodeProxy( double xval, double yval , uint index){
			coord.set_size(1,2,1.0);
			coord(0,0) = xval; coord(0,1) = yval;
			mindex = index;

		}

		//------------------------------------------------------------------------------
		~NodeProxy(){};

		//------------------------------------------------------------------------------
		Matrix< DDRMat > get_coords() const
		{
			return coord;
		};

	    moris_id
		get_id() const
		{
	    	MORIS_ERROR( false,"Function not implemented in base vertex");
		    return gNoID;
		};

		//------------------------------------------------------------------------------

		moris_index
		get_index() const
		{
			return mindex;
		};

		//------------------------------------------------------------------------------
		moris_index
		get_owner() const
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		    return 0;
		};

		//------------------------------------------------------------------------------
	    mtk::Vertex_Interpolation *
		get_interpolation( const uint aOrder )
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		    return nullptr;
		}

		//------------------------------------------------------------------------------
		mtk::Vertex_Interpolation *
		get_interpolation( const uint aOrder ) const
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		    return nullptr;
		 }

		//------------------------------------------------------------------------------
		uint
		get_level() const
		{
			return 0;
		};

		//------------------------------------------------------------------------------
		void
		flag()
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		};

		//------------------------------------------------------------------------------
		void
		unflag()
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		};

		//------------------------------------------------------------------------------
		bool
		is_flagged() const
		{
			MORIS_ERROR( false,"Function not implemented in base vertex" );
		    return false;
		};
	};

}/* namespace fem */
}/* namespace moris */



#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_NODE_PROXY_HPP_ */
