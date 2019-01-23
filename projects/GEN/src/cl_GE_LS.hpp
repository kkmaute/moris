/*
 * cl_GE_LS.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_LS_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_LS_HPP_

#include "fn_norm.hpp"
#include "cl_GE_Geometry.hpp"

namespace moris
{
namespace ge
{
	class LS : public Geometry
	{
	private:

	protected:

	public:
		LS()
	{
			std::cout<<"LS constructor"<<std::endl;
	};
		~LS(){};

	};
} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_MTK_GE_SRC_CL_GE_LS_HPP_ */
