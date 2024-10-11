/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_Library_Enums.cpp
 *
 */

#include "cl_Library_Enums.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    std::string
    convert_parameter_list_enum_to_string( Module_Type aParameterListType )
    {
        switch ( aParameterListType )
        {
            case Module_Type::OPT:
                return "OPT";
            case Module_Type::HMR:
                return "HMR";
            case Module_Type::STK:
                return "STK";
            case Module_Type::XTK:
                return "XTK";
            case Module_Type::GEN:
                return "GEN";
            case Module_Type::FEM:
                return "FEM";
            case Module_Type::SOL:
                return "SOL";
            case Module_Type::MSI:
                return "MSI";
            case Module_Type::VIS:
                return "VIS";
            case Module_Type::MIG:
                return "MIG";
            case Module_Type::WRK:
                return "WRK";
            case Module_Type::MORISGENERAL:
                return "MORISGENERAL";

            default:
                MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
                return "";
        }
    }

    // -----------------------------------------------------------------------------

    std::string
    get_name_for_parameter_list_type( Module_Type aParameterListType )
    {
        return convert_parameter_list_enum_to_string( aParameterListType ) + "ParameterList";
    }
    
    // -----------------------------------------------------------------------------
    
    Vector< std::string > get_submodule_names( Module_Type aModuleType )
    {
        // get the names of the sub-parameter lists for each of the modules
        switch ( aModuleType )
        {
            case Module_Type::OPT:
                return OPT_Submodule_String::values;

            case Module_Type::GEN:
                return GEN_Submodule_String::values;

            case Module_Type::FEM:
                return FEM_Submodule_String::values;

            case Module_Type::SOL:
                return SOL_Submodule_String::values;

            case Module_Type::MORISGENERAL:
                // needs reworking anyways
                return { "Remeshing", "Refinement", "Mapping" };
                
            case Module_Type::END_ENUM:
                return {};

            default:
                return { "General" };
        }
    }

    // -----------------------------------------------------------------------------

    uint
    get_number_of_sub_parameter_lists_in_module( Module_Type aModuleType )
    {
        return get_submodule_names( aModuleType ).size();
    }

    // -----------------------------------------------------------------------------

    std::string
    get_inner_sub_parameter_list_name(
            Module_Type aModuleType,
            uint                aParamListIndex )
    {
        // initialize the names with the standard
        Vector< std::string > tNames = { "" };

        // get the names of the sub-parameter lists for each of the modules
        switch ( aModuleType )
        {
            case Module_Type::OPT:
                tNames = { "OptimizationProblem", "Interface", "Algorithms" };
                break;

            case Module_Type::HMR:
                break;    // standard name

            case Module_Type::STK:
                break;    // standard name

            case Module_Type::XTK:
                break;    // standard name

            case Module_Type::GEN:
                tNames = { "General", "Geometry", "Field" };
                break;

            case Module_Type::FEM:
                tNames = {
                    "Property",         // 0
                    "CM",               // 1
                    "SP",               // 2
                    "IWG",              // 3
                    "IQI",              // 4
                    "",                 // 5
                    "Field",            // 6
                    "MaterialPhase",    // 7
                    "MM"                // 8
                };
                break;

            case Module_Type::SOL:
                tNames = {
                    "LinearAlgorithm",        // 0
                    "LinearSolver",           // 1
                    "NonLinearAlgorithm",     // 2
                    "NonLinearSolver",        // 3
                    "TimeSolverAlgorithm",    // 4
                    "TimeSolver",             // 5
                    "SolverWarehouse",        // 6
                    "Preconditioner"          // 7
                };
                break;

            case Module_Type::MSI:
                break;    // standard name

            case Module_Type::VIS:
                tNames = { "OutputMesh" };
                break;

            case Module_Type::MIG:
                break;    // standard name

            case Module_Type::WRK:
                break;    // standard name

            case Module_Type::MORISGENERAL:
                tNames = { "", "", "" };
                break;

            default:
                MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
                break;
        }

        // check validity of the input
        uint tNumSubParamLists = tNames.size();
        MORIS_ERROR( aParamListIndex < tNumSubParamLists,
                "Library_Enums::get_inner_sub_parameter_list_name() - "
                "Trying to retrieve names of parameters in sub-parameter list %i for module %s, but it has only %i sub-parameter lists.",
                aParamListIndex,
                convert_parameter_list_enum_to_string( aModuleType ).c_str(),
                tNumSubParamLists );

        // retrieve the name for the specific sub-parameter list requested
        return tNames( aParamListIndex );
    }

    // -----------------------------------------------------------------------------

}    // namespace moris
