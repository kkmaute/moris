/*
 * cl_MSI_Parameters.hpp
 *
 *  Created on: May 5, 2018
 *      Author: schmidt
 */

#ifndef SRC_HMR_CL_MSI_PARAMETERS_HPP_
#define SRC_HMR_CL_MSI_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"

#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Param_List.hpp"       //CON/src

#include "cl_XML_Parser.hpp"

namespace moris
{
    namespace MSI
    {

// -----------------------------------------------------------------------------

        // fixme: to be deleted soon
        // creates a parameter list with default inputs
        void load_hmr_parameter_list_from_xml( const std::string   & aFilePath,
                                                     ParameterList & aParameterList );

// -----------------------------------------------------------------------------

        // creates a parameter list with default options
        ParameterList create_hmr_parameter_list();

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_MSI_PARAMETERS_HPP_ */
