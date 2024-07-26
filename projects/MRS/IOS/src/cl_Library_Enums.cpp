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
    convert_parameter_list_enum_to_string( Parameter_List_Type aParameterListType )
    {
        switch ( aParameterListType )
        {
            case Parameter_List_Type::OPT:
                return "OPT";
            case Parameter_List_Type::HMR:
                return "HMR";
            case Parameter_List_Type::STK:
                return "STK";
            case Parameter_List_Type::XTK:
                return "XTK";
            case Parameter_List_Type::GEN:
                return "GEN";
            case Parameter_List_Type::FEM:
                return "FEM";
            case Parameter_List_Type::SOL:
                return "SOL";
            case Parameter_List_Type::MSI:
                return "MSI";
            case Parameter_List_Type::VIS:
                return "VIS";
            case Parameter_List_Type::MIG:
                return "MIG";
            case Parameter_List_Type::WRK:
                return "WRK";
            case Parameter_List_Type::MORISGENERAL:
                return "MORISGENERAL";

            default:
                MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
                return "";
        }
    }

    // -----------------------------------------------------------------------------

    std::string
    get_name_for_parameter_list_type( Parameter_List_Type aParameterListType )
    {
        return convert_parameter_list_enum_to_string( aParameterListType ) + "ParameterList";
    }

    // -----------------------------------------------------------------------------

    uint
    get_number_of_sub_parameter_lists_in_module( Parameter_List_Type aModule )
    {
        // get the names of the sub-parameter lists for each of the modules
        switch ( aModule )
        {
            case Parameter_List_Type::OPT:
                return 1;

            case Parameter_List_Type::HMR:
                return 1;

            case Parameter_List_Type::STK:
                return 1;

            case Parameter_List_Type::XTK:
                return 1;

            case Parameter_List_Type::GEN:
                return 3;

            case Parameter_List_Type::FEM:
                return 9;

            case Parameter_List_Type::SOL:
                return 8;

            case Parameter_List_Type::MSI:
                return 1;

            case Parameter_List_Type::VIS:
                return 1;

            case Parameter_List_Type::MIG:
                return 1;

            case Parameter_List_Type::WRK:
                return 1;

            case Parameter_List_Type::MORISGENERAL:
                return 3;

            default:
                MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
                return 0;
        }
    }

    // -----------------------------------------------------------------------------

    std::string
    get_outer_sub_parameter_list_name(
            Parameter_List_Type aModule,
            uint                aParamListIndex )
    {
        // initialize the names with the standard
        Vector< std::string > tNames = { "General" };

        // get the names of the sub-parameter lists for each of the modules
        switch ( aModule )
        {
            case Parameter_List_Type::OPT:
                tNames = { "OptimizationProblems", "Interface", "Algorithms" };
                break;

            case Parameter_List_Type::HMR:
                break;    // standard name

            case Parameter_List_Type::STK:
                break;    // standard name

            case Parameter_List_Type::XTK:
                break;    // standard name

            case Parameter_List_Type::GEN:
                tNames = { "General", "Geometries", "Properties" };
                break;

            case Parameter_List_Type::FEM:
                tNames = {
                    "Properties",                 // 0
                    "ConstitutiveModels",         // 1
                    "StabilizationParameters",    // 2
                    "WeakForms",                  // 3
                    "QuantitiesOfInterest",       // 4
                    "ComputationParameters",      // 5
                    "Fields",                     // 6
                    "Materials",                  // 7
                    "MaterialModels"              // 8
                };
                break;

            case Parameter_List_Type::SOL:
                tNames = {
                    "LinearAlgorithms",        // 0
                    "LinearSolvers",           // 1
                    "NonLinearAlgorithms",     // 2
                    "NonLinearSolvers",        // 3
                    "TimeSolverAlgorithms",    // 4
                    "TimeSolvers",             // 5
                    "SolverWarehouse",         // 6
                    "Preconditioners"          // 7
                };
                break;

            case Parameter_List_Type::MSI:
                break;    // standard name

            case Parameter_List_Type::VIS:
                tNames = { "OutputMeshes" };
                break;

            case Parameter_List_Type::MIG:
                break;    // standard name

            case Parameter_List_Type::WRK:
                break;    // standard name

            case Parameter_List_Type::MORISGENERAL:
                tNames = { "Remeshing", "Refinement", "Mapping" };
                break;

            default:
                MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
                break;
        }

        // check validity of the input
        uint tNumSubParamLists = tNames.size();
        MORIS_ERROR( aParamListIndex < tNumSubParamLists,
                "Library_Enums::get_outer_sub_parameter_list_name() - "
                "Trying to retrieve name for sub-parameter list %i for module %s, but it has only %i sub-parameter lists.",
                aParamListIndex,
                convert_parameter_list_enum_to_string( aModule ).c_str(),
                tNumSubParamLists );

        // retrieve the name for the specific sub-parameter list requested
        return tNames( aParamListIndex );
    }

    // -----------------------------------------------------------------------------

    std::string
    get_inner_sub_parameter_list_name(
            Parameter_List_Type aModule,
            uint                aParamListIndex )
    {
        // initialize the names with the standard
        Vector< std::string > tNames = { "" };

        // get the names of the sub-parameter lists for each of the modules
        switch ( aModule )
        {
            case Parameter_List_Type::OPT:
                tNames = { "OptimizationProblem", "Interface", "Algorithms" };
                break;

            case Parameter_List_Type::HMR:
                break;    // standard name

            case Parameter_List_Type::STK:
                break;    // standard name

            case Parameter_List_Type::XTK:
                break;    // standard name

            case Parameter_List_Type::GEN:
                tNames = { "", "Geometry", "Field" };
                break;

            case Parameter_List_Type::FEM:
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

            case Parameter_List_Type::SOL:
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

            case Parameter_List_Type::MSI:
                break;    // standard name

            case Parameter_List_Type::VIS:
                tNames = { "OutputMesh" };
                break;

            case Parameter_List_Type::MIG:
                break;    // standard name

            case Parameter_List_Type::WRK:
                break;    // standard name

            case Parameter_List_Type::MORISGENERAL:
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
                convert_parameter_list_enum_to_string( aModule ).c_str(),
                tNumSubParamLists );

        // retrieve the name for the specific sub-parameter list requested
        return tNames( aParamListIndex );
    }

    // -----------------------------------------------------------------------------

}    // namespace moris
