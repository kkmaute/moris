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
    enum class Parameter_List_Type
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

    enum class OPT_SubModule
    {
        OPTIMIZATION_PROBLEMS,
        INTERFACE,
        ALGORITHMS
    };

    enum class HMR_SubModule
    {
        GENERAL
    };

    enum class STK_SubModule
    {
        GENERAL
    };

    enum class XTK_SubModule
    {
        GENERAL
    };

    enum class GEN_SubModule
    {
        GENERAL,
        GEOMETRIES,
        PROPERTIES
    };

    enum class FEM_SubModule
    {
        PROPERTIES,
        CONSTITUTIVE_MODELS,
        STABILIZATION_PARAMETER,
        IWG,
        IQI,
        COMPUTATION_PARAMETERS,
        FIELDS,
        MATERIAL_MODELS
    };

    enum class SOL_SubModule
    {
        LINEAR_ALGORITHMS,
        LINEAR_SOLVERS,
        NONLINEAR_ALGORITHMS,
        NONLINEAR_SOLVERS,
        TIME_SOLVER_ALGORITHMS,
        TIME_SOLVERS,
        SOLVER_WAREHOUSE,
        PRECONDITIONERS
    };

    enum class MSI_SubModule
    {
        GENERAL
    };

    enum class VIS_SubModule
    {
        OUTPUT_MESHES
    };

    enum class MIG_SubModule
    {
        GENERAL
    };

    enum class WRK_SubModule
    {
        GENERAL
    };

    enum class MORISGENERAL_SubModule
    {
        REMESHING,
        REFINEMENT,
        MAPPING
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
     * @return std::string name of the parameter list function expected in an .so input file
     */
    std::string
    get_name_for_parameter_list_type( Parameter_List_Type aParameterListType );

    // -----------------------------------------------------------------------------

    /**
     * @brief Get the number of sub parameter lists that can be defined in a given module
     *
     * @param aModule the module
     * @return uint number of sub parameter lists in Module
     */
    uint
    get_number_of_sub_parameter_lists_in_module( Parameter_List_Type aModule );

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
