/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.cpp
 *
 */

#include "cl_Library_IO_Standard.hpp"
#include "enums.hpp"
#include "parameters.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::Library_IO_Standard()
            : Library_IO()    // initialize base class data as usual
    {
        // list of supported parameter list types
        mSupportedParamListTypes = {
            Module_Type::OPT,
            Module_Type::HMR,
            Module_Type::STK,
            Module_Type::XTK,
            Module_Type::GEN,
            Module_Type::FEM,
            Module_Type::SOL,
            Module_Type::MSI,
            Module_Type::VIS,
            Module_Type::MIG,
            Module_Type::WRK,
            Module_Type::MORISGENERAL
        };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::~Library_IO_Standard()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::load_all_standard_parameters()
    {
        // go over all modules and check that
        for ( uint iModule = 0; iModule < (uint)( Module_Type::END_ENUM ); iModule++ )
        {
            // get the current module
            Module_Type tParamListType = (Module_Type)( iModule );

            // get access to this module's parameter list
            Module_Parameter_Lists tModuleParamList = mParameterLists( iModule );

            // for each parameter list type, initialize it with a default
            switch ( tParamListType )
            {
                case Module_Type::OPT:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
                    break;

                case Module_Type::HMR:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
                    break;

                case Module_Type::STK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_stk_parameter_list() );
                    break;

                case Module_Type::XTK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
                    break;

                case Module_Type::GEN:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
                    break;

                case Module_Type::FEM:
                    tModuleParamList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
                    break;

                case Module_Type::SOL:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
                    tModuleParamList( 1 ).add_parameter_list( prm::create_linear_solver_parameter_list() );
                    tModuleParamList( 2 ).add_parameter_list( prm::create_nonlinear_algorithm_parameter_list() );
                    tModuleParamList( 3 ).add_parameter_list( prm::create_nonlinear_solver_parameter_list() );
                    tModuleParamList( 4 ).add_parameter_list( prm::create_time_solver_algorithm_parameter_list() );
                    tModuleParamList( 5 ).add_parameter_list( prm::create_time_solver_parameter_list() );
                    tModuleParamList( 6 ).add_parameter_list( prm::create_solver_warehouse_parameterlist() );
                    tModuleParamList( 7 ).add_parameter_list( prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
                    break;

                case Module_Type::MSI:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
                    break;

                case Module_Type::VIS:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
                    tModuleParamList( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
                    tModuleParamList( 0 )( 0 ).set( "Save_Frequency", 1 );
                    break;

                case Module_Type::MIG:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_mig_parameter_list() );
                    break;

                case Module_Type::WRK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_wrk_parameter_list() );
                    break;

                case Module_Type::MORISGENERAL:
                    break;

                case Module_Type::END_ENUM:
                    MORIS_ERROR( false,
                            "Library_IO_Standard::load_all_standard_parameters() - "
                            "Module_Type is UNDEFINED. This loop shouldn't get here." );
                    break;

                default:
                    MORIS_ERROR( false,
                            "Library_IO_Standard::load_all_standard_parameters() - "
                            "No standard library defined for module parameter list type #%i. "
                            "Most likely a new enum has been defined that isn't used in this function.",
                            (uint)( tParamListType ) );
                    break;
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
