/*
 * fn_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_

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
    inline
    moris::ParameterList create_msi_parameter_list()
    {
        ParameterList mMSIParameterList;

        // Adof type interpolation index
        mMSIParameterList.insert( "TEMP"       , 0 );
        mMSIParameterList.insert( "P"          , 0 );
        mMSIParameterList.insert( "RHO"        , 0 );
        mMSIParameterList.insert( "E"          , 0 );
        mMSIParameterList.insert( "EVP"        , 0 );
        mMSIParameterList.insert( "EVT"        , 0 );

        mMSIParameterList.insert( "UX"         , 0 );
        mMSIParameterList.insert( "UY"         , 0 );
        mMSIParameterList.insert( "UZ"         , 0 );
        mMSIParameterList.insert( "VX"         , 0 );
        mMSIParameterList.insert( "VY"         , 0 );
        mMSIParameterList.insert( "VZ"         , 0 );
        mMSIParameterList.insert( "MX"         , 0 );
        mMSIParameterList.insert( "MY"         , 0 );
        mMSIParameterList.insert( "MZ"         , 0 );
        mMSIParameterList.insert( "EVX"        , 0 );
        mMSIParameterList.insert( "EVY"        , 0 );
        mMSIParameterList.insert( "EVZ"        , 0 );
        mMSIParameterList.insert( "NLSX"       , 0 );
        mMSIParameterList.insert( "NLSY"       , 0 );
        mMSIParameterList.insert( "NLSZ"       , 0 );

        mMSIParameterList.insert( "L2"         , 0 );
        mMSIParameterList.insert( "MAPPING_DOF", 0 );
        mMSIParameterList.insert( "LS1"        , 0 );
        mMSIParameterList.insert( "LS2"        , 0 );

        mMSIParameterList.insert( "THETA"      , 0 );
        mMSIParameterList.insert( "PHID"       , 0 );
        mMSIParameterList.insert( "PHISD"      , 0 );

        mMSIParameterList.insert( "VISCOSITY"  , 0 );
        mMSIParameterList.insert( "STRESS_DOF" , 0 );

        // General MSI parameters
        mMSIParameterList.insert( "order_adofs_by_host", false );

        // Geometric multigrid parameters
        mMSIParameterList.insert( "multigrid"  , false );
        mMSIParameterList.insert( "level"      , 2 );

        return mMSIParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_ */
