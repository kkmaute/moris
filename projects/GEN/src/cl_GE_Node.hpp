/*
 * cl_GE_Node.hpp
 *
 *  Created on: Jan 2, 2019
 *      Author: sonne
 */
#ifndef PROJECTS_GEN_SRC_CL_NODE_HPP_
#define PROJECTS_GEN_SRC_CL_NODE_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
	namespace ge
	{
		class Node : public mtk::Vertex
		{
			//------------------------------------------------------------------------------
		private:
			moris::Matrix< DDRMat > coord;

			//------------------------------------------------------------------------------

		public:
			//------------------------------------------------------------------------------
			Node(double xval, double yval){
				coord.set_size(1,2,1.0);
				coord(0,0) = xval; coord(0,1) = yval;
			}

            Node(double xval, double yval, double zval){
                coord.set_size(1,3,1.0);
                coord(0,0) = xval; coord(0,1) = yval; coord(0,2) = zval;
            }
			//------------------------------------------------------------------------------
			~Node(){};

			//------------------------------------------------------------------------------
			// gives the coordinates of node
			Matrix< DDRMat > get_coords() const
			{
				return coord;
			};

			//------------------------------------------------------------------------------
			// gives domain wide ID of node
			moris_id
			get_id() const
			{
				MORIS_ERROR( false, "get_id() not implemented " );
				return gNoID;
			};

			//------------------------------------------------------------------------------
			// gives index of node
			moris_index
			get_index() const
			{
				MORIS_ERROR( false, "get_index() not implemented " );
				return 0;
			};

			//------------------------------------------------------------------------------
			// gives owner of node
			moris_index
			get_owner() const
			{
				MORIS_ERROR( false, "get_index() not implemented " );
				return 0;
			};

			//------------------------------------------------------------------------------
			mtk::Vertex_Interpolation*
			get_interpolation( const uint aOrder )
			{
				MORIS_ERROR( false, "get_interpolation() not implemented " );
				return nullptr;
			};

			//------------------------------------------------------------------------------
			const mtk::Vertex_Interpolation*
			get_interpolation( const uint aOrder ) const
			{
				MORIS_ERROR( false, "get_interpolation() not implemented " );
				return nullptr;
			};

			//------------------------------------------------------------------------------
			uint
			get_level() const
			{
				return( 0 );
			};

			//------------------------------------------------------------------------------
			void
			flag()
			{
				MORIS_ERROR( false, "flag() not implemented " );
			};

			//------------------------------------------------------------------------------
			void
			unflag() const
			{
				MORIS_ERROR( false, "unflag() not implemented " );
			};

			//------------------------------------------------------------------------------
			bool
			is_flagged() const
			{
				MORIS_ERROR( false, "is_flagged() not implemented " );
				return false;
			};


		};

	} /* namespace ge */
} /* namespace moris */




#endif /* PROJECTS_GEN_SRC_CL_NODE_HPP_ */
