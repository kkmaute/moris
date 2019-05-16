/*
 * cl_GE_SDF.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_

#include <iostream>
#include <string>
#include <memory>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_GE_Geometry.hpp"

#include "MTK/src/cl_MTK_Mesh.hpp"
#include "GEN/SDF/src/cl_SDF_Object.hpp"
#include "GEN/SDF/src/cl_SDF_Mesh.hpp"

namespace moris
{
namespace ge
{
	class SDF : public Geometry
	{
	public:

		SDF(){};
		~SDF(){};

        //------------------------------------------------------------------------------
        /*
         * @brief function to report if geometry representation type
         */
        enum type get_geom_type() const
        {
            return type::SDF;
        }

    private:

    protected:

	};
} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_ */
