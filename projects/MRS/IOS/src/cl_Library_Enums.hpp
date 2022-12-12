/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_Library_Enums.hpp
 *
 */
#ifndef SRC_cl_Library_Enums
#define SRC_cl_Library_Enums

#include <string>
#include "typedefs.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    // Static constants
    const std::string PARAM_LIST_FUNC_NAME_ENDING = "ParameterList";
    const std::string XML_PARAMETER_FILE_ROOT = "ParameterLists";
    const std::string PRINT_ERROR = "print_error";

    // -----------------------------------------------------------------------------

    // define enums that define the state of an xml parser
    enum class XML_Mode
    {
        UNINITIALIZED,
        READ,
        WRITE,
        END_ENUM
    };

    // -----------------------------------------------------------------------------

    // define enums for types of library
    enum class Library_Type
    {
        UNDEFINED,
        STANDARD,
        MESHGEN,
        END_ENUM
    };

    // -----------------------------------------------------------------------------

    // define enums for input file types
    enum class File_Type
    {
        UNDEFINED,
        SO_FILE,
        XML_FILE,
        END_ENUM
    };

    // -----------------------------------------------------------------------------

    // define enums for input file types
    enum class Parameter_List_Type
    {
        //! only have types here for which we have a standard parameter list (with the exception of "END_ENUM")
        //! Make sure to include all enums defined here in the convert to string function below (with the exception of "END_ENUM")
        OPT,
        HMR,
        STK,
        XTK,
        GEN,
        FEM,
        SOL,
        MSI,
        VIS,
        MIG,
        WRK,
        MORISGENERAL,
        END_ENUM
    };

    // -----------------------------------------------------------------------------

    /**
     * @brief converts an enum of type Parameter_List_Type to a string spelled equivalent to the enum
     *
     * @param aParameterListType enum to convert to string
     * @return std::string enum spelled out as string
     */
    std::string
    convert_parameter_list_enum_to_string( Parameter_List_Type aParameterListType );

    // -----------------------------------------------------------------------------

    /**
     * @brief Get the name for parameter list function expected for a given module
     *
     * @param aParameterListType Parameter_List_Type enum naming the module for which parameters should be parsed
     * @return std::string name of the parameter list funciton expected in an .so input file
     */
    std::string
    get_name_for_parameter_list_type( Parameter_List_Type aParameterListType );

    // -----------------------------------------------------------------------------

    std::string
    get_outer_sub_parameter_list_name( 
            Parameter_List_Type aModule,
            uint                aParamListIndex );

    // -----------------------------------------------------------------------------
    
    std::string
    get_inner_sub_parameter_list_name( 
            Parameter_List_Type aModule,
            uint                aParamListIndex );

    // -----------------------------------------------------------------------------

}    // namespace moris

#endif /* cl_Library_Enums.hpp */