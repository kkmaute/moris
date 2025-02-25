/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_Library_Enums.hpp
 *
 */
#pragma once

#include <string>
#include "moris_typedefs.hpp"
#include "fn_enum_macros.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    // Static constants
    const std::string PARAM_LIST_FUNC_NAME_ENDING = "ParameterList";
    const std::string XML_PARAMETER_FILE_ROOT     = "ParameterLists";
    const std::string PRINT_ERROR                 = "print_error";

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
    enum class Module_Type
    {
        //! only have types here for which we have a standard parameter list (with the exception of "UNDEFINED")
        //! Make sure to include all enums defined here in the convert to string function below (with the exception of "UNDEFINED")
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

    ENUM_MACRO( OPT_Submodule,
        OPTIMIZATION_PROBLEMS,
        INTERFACE,
        ALGORITHMS )

    ENUM_MACRO( HMR_Submodule,
        GENERAL,
        LAGRANGE_MESHES,
        BSPLINE_MESHES )

    enum class STK_Submodule
    {
        GENERAL
    };

    enum class XTK_Submodule
    {
        GENERAL
    };

    ENUM_MACRO( GEN_Submodule,
        GENERAL,
        GEOMETRIES,
        PROPERTIES )

    ENUM_MACRO( FEM_Submodule,
        PROPERTIES,
        CONSTITUTIVE_MODELS,
        STABILIZATION,
        IWG,
        IQI,
        COMPUTATION,
        FIELDS,
        PHASES,
        MATERIAL_MODELS )

    ENUM_MACRO( SOL_Submodule,
        LINEAR_ALGORITHMS,
        LINEAR_SOLVERS,
        NONLINEAR_ALGORITHMS,
        NONLINEAR_SOLVERS,
        TIME_SOLVER_ALGORITHMS,
        TIME_SOLVERS,
        SOLVER_WAREHOUSE,
        PRECONDITIONERS )

    enum class MSI_Submodule
    {
        GENERAL
    };

    enum class VIS_Submodule
    {
        OUTPUT_MESHES
    };

    enum class MIG_Submodule
    {
        GENERAL
    };

    enum class WRK_Submodule
    {
        GENERAL
    };

    enum class MORISGENERAL_Submodule
    {
        REMESHING,
        REFINEMENT,
        MAPPING
    };

    // -----------------------------------------------------------------------------

    /**
     * @brief converts an enum of type Parameter_List_Type to a string spelled equivalent to the enum
     *
     * @param aModuleType enum to convert to string
     * @return std::string enum spelled out as string
     */
    std::string
    convert_parameter_list_enum_to_string( Module_Type aModuleType );

    // -----------------------------------------------------------------------------

    /**
     * @brief Get the name for parameter list function expected for a given module
     *
     * @param aModuleType Parameter_List_Type enum naming the module for which parameters should be parsed
     * @return std::string name of the parameter list function expected in an .so input file
     */
    std::string
    get_name_for_parameter_list_type( Module_Type aModuleType );

    /**
     * Gets the submodule names for a particular MORIS module.
     * 
     * @param aModuleType Module parameter list type
     * @return Submodule names
     */
    Vector< std::string > get_submodule_names( Module_Type aModuleType );

    /**
     * @brief Get the number of sub parameter lists that can be defined in a given module
     *
     * @param aModule the module
     * @return uint number of sub parameter lists in Module
     */
    uint
    get_number_of_sub_parameter_lists_in_module( Module_Type aModule );

    // -----------------------------------------------------------------------------

    std::string
    get_inner_sub_parameter_list_name(
            Module_Type aModule,
            uint                aParamListIndex );

    // -----------------------------------------------------------------------------

    // define enums for geometry types
    enum class Geometry_Type
    {
        UNDEFINED,
        PRE_DEFINED,
        USER_DEFINED,
        IMAGE_FILE,
        END_ENUM
    };

    // -----------------------------------------------------------------------------

}    // namespace moris
