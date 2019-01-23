/*
 * cl_GE_Geometry.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_

#include "cl_GE_Element.hpp"

namespace moris
{
namespace ge
{
	class Geometry{
	private:

	protected:
		uint mNumEle = 0;
	public:
		Geometry(){};

		~Geometry(){};
        //------------------------------------------------------------------------------

		virtual real get_val_at_vertex( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConstant )
		{
			std::cout<<"yikes, you shouldn't be seeing this....."<<std::endl;
			return 0;
		};

        //------------------------------------------------------------------------------

		virtual void flag_intersection_locations()
		{
			std::cout<<"yikes, you shouldn't be seeing this....."<<std::endl;
		};

        //------------------------------------------------------------------------------

		virtual void set_analytical_function( real ( *mFunctionAnalytical )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ) )
		{
			std::cout<<"yikes, you shouldn't be seeing this....."<<std::endl;
		};

        //------------------------------------------------------------------------------

	};
} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_GEN_SRC_CL_GE_BASE_HPP_ */
