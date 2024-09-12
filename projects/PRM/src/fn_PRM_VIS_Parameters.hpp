/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_VIS_Parameters.hpp
 *
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_VIS_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_VIS_PARAMETERS_HPP_

#include "cl_Parameter_List.hpp"

#include "cl_VIS_Output_Enums.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline moris::Parameter_List create_vis_parameter_list()
    {
        Parameter_List mVISParameterList;

        mVISParameterList.insert( "Output_Index", 0 );
        mVISParameterList.insert( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        mVISParameterList.insert( "File_Name", std::pair< std::string, std::string >( "", "" ) );
        mVISParameterList.insert( "Temp_Name", std::pair< std::string, std::string >( "./", "temp.exo" ) );
        mVISParameterList.insert( "Save_Frequency", MORIS_SINT_MAX );
        mVISParameterList.insert( "Time_Offset", 0.0 );
        mVISParameterList.insert( "Set_Names", "" );
        mVISParameterList.insert( "Field_Names", "" );
        mVISParameterList.insert( "Field_Type", "" );
        mVISParameterList.insert( "IQI_Names", "" );

        return mVISParameterList;
    }
    //------------------------------------------------------------------------------

}    // namespace moris::prm

#endif /* PROJECTS_PRM_SRC_FN_PRM_VIS_PARAMETERS_HPP_ */

