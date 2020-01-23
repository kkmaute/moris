/*
 * cl_MSI_Parameters.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: schmidt
 */

#include "cl_MSI_Parameters.hpp" //HMR/src

#include "assert.hpp"
#include "fn_Parsing_Tools.hpp"

#include "fn_unique.hpp"

namespace moris
{
    namespace MSI
    {

// -----------------------------------------------------------------------------

    // creates a parameter list with default inputs
    ParameterList create_hmr_parameter_list()
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

//--------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    void load_hmr_parameter_list_from_xml( const std::string   & aFilePath,
                                                 ParameterList & aParameterList )
    {
        MORIS_ERROR( false, "not implemented");
    }

//--------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
