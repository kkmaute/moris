/*
 * cl_GE_Enums.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_ENUMS_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_ENUMS_HPP_

namespace moris
{
	namespace ge
	{
		enum class type{
// geometry pointer types
			ANALYTIC,
			DISCRETE,
// analytic geometry types
			SPHERE,
			PLANE,
			SPIRAL,
			GYROID,
			COMPOSITE_FIBER,
			COMPOSITE_FIBER_WAVE_1,
			COMPOSITE_FIBER_WAVE_2,
			COMPOSITE_FIBER_WAVE_3,
			COMPOSITE_FIBER_STRAIGHT_1,
			COMPOSITE_FIBER_STRAIGHT_2,
			COMPOSITE_FIBER_STRAIGHT_3,
			MULTI_CYLINDER,
// end of enums
			END_ENUM
		};
	} /* namespace gen */
} /* namespace moris */



#endif /* PROJECTS_MTK_GE_SRC_CL_GE_ENUMS_HPP_ */
