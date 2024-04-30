/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Parameters.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_PARAMETERS_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_PARAMETERS_HPP_

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_XML_Parser.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

//        typedef Param_List< boost::variant< sint, real, std::string  > > ParameterList;

Parameter_List
create_sdf_parameter_list();

//-------------------------------------------------------------------------------

Parameter_List
create_sdf_object_parameter_list();

//-------------------------------------------------------------------------------

        void
        load_sdf_parameter_list_from_xml(
                const std::string            & aFilePath,
        Parameter_List           & aGlobalParameters,
                Vector< Parameter_List > & aObjectParameters );

//-------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_PARAMETERS_HPP_ */

