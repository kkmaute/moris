/*
 * cl_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_

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

    // creates a parameter list with default inputs
    moris::ParameterList create_msi_parameter_list()
    {
        ParameterList mMSIParameterList;

        mMSIParameterList.insert( "UX"         , 0 );
        mMSIParameterList.insert( "UY"         , 0 );
        mMSIParameterList.insert( "UZ"         , 0 );
        mMSIParameterList.insert( "TEMP"       , 0 );
        mMSIParameterList.insert( "L2"         , 0 );
        mMSIParameterList.insert( "MAPPING_DOF", 0 );
        mMSIParameterList.insert( "LS1"        , 0 );
        mMSIParameterList.insert( "LS2"        , 0 );
        mMSIParameterList.insert( "NLSX"       , 0 );
        mMSIParameterList.insert( "NLSY"       , 0 );
        mMSIParameterList.insert( "NLSZ"       , 0 );
        mMSIParameterList.insert( "VX"         , 0 );
        mMSIParameterList.insert( "VY"         , 0 );
        mMSIParameterList.insert( "VZ"         , 0 );
        mMSIParameterList.insert( "P"          , 0 );

        return mMSIParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_ */
