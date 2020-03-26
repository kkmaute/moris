/*
 * cl_PRM_GEN_Parameters.hpp
 *
 *  Created on: Feb 18, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
//#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Param_List.hpp"


namespace moris
{
    namespace prm
    {
    //------------------------------------------------------------------------------
    moris::ParameterList create_gen_parameter_list()
    {
        ParameterList tParameterList;

        tParameterList.insert( "geometries", std::string("") );

        return tParameterList;
    }

    //------------------------------------------------------------------------------
    }   // end prm namespace
}       // end moris namespace



#endif /* PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_ */
