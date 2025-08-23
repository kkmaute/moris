/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Module_Parameter_Lists.cpp
 *
 */

#include "cl_Module_Parameter_Lists.hpp"
#include "parameters.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Module_Parameter_Lists::Module_Parameter_Lists( Module_Type aModuleType )
            : mParameterListType( aModuleType )
    {
        // Get submodule names
        Vector< std::string > tSubmoduleNames = get_submodule_names( aModuleType );

        // Add new submodules
        for ( const std::string& iSubmoduleName : tSubmoduleNames )
        {
            mSubmoduleParameterLists.emplace_back( iSubmoduleName );
        }

        // Create default parameter lists and set uniqueness flag
        switch ( mParameterListType )
        {
            case Module_Type::OPT:
                mSubmoduleParameterLists( index( OPT::OPTIMIZATION_PROBLEMS, std::true_type() ) ).add_parameter_list( prm::create_opt_problem_parameter_list(), true );
                break;
            case Module_Type::HMR:
                mSubmoduleParameterLists( index( HMR::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_hmr_parameter_list(), true );
                break;
            case Module_Type::STK:
                mSubmoduleParameterLists( index( STK_Submodule::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_stk_parameter_list(), true );
                break;
            case Module_Type::XTK:
                mSubmoduleParameterLists( index( XTK_Submodule::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_xtk_parameter_list(), true );
                break;
            case Module_Type::GEN:
                mSubmoduleParameterLists( index( GEN::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_gen_parameter_list(), true );
                break;
            case Module_Type::FEM:
                mSubmoduleParameterLists( index( FEM_Submodule::COMPUTATION, std::true_type() ) ).add_parameter_list( prm::create_computation_parameter_list(), true );
                break;
            case Module_Type::SOL:
                mSubmoduleParameterLists( index( SOL::SOLVER_WAREHOUSE, std::true_type() ) ).add_parameter_list( prm::create_solver_warehouse_parameterlist(), true );
                break;
            case Module_Type::MSI:
                mSubmoduleParameterLists( index( MSI_Submodule::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_msi_parameter_list(), true );
                break;
            case Module_Type::VIS:
                mSubmoduleParameterLists( index( VIS_Submodule::OUTPUT_MESHES, std::true_type() ) ).add_parameter_list( prm::create_vis_parameter_list(), false);
                break;
            case Module_Type::MIG:
                mSubmoduleParameterLists( index( MIG_Submodule::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_mig_parameter_list(), true );
                break;
            case Module_Type::WRK:
                mSubmoduleParameterLists( index( WRK_Submodule::GENERAL, std::true_type() ) ).add_parameter_list( prm::create_wrk_parameter_list(), true );
                break;
            default:    // Do nothing
                break;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Module_Parameter_Lists::size() const
    {
        return mSubmoduleParameterLists.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Module_Parameter_Lists::clear()
    {
        mSubmoduleParameterLists.clear();
        mParameterListType = Module_Type::END_ENUM;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Module_Parameter_Lists::insert(
            const std::string&     aName,
            const Design_Variable& aDesignVariable )
    {
        mSubmoduleParameterLists( mCurrentSubmoduleIndex ).insert( aName, aDesignVariable );
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Module_Parameter_Lists::begin() -> decltype( mSubmoduleParameterLists.begin() )
    {
        return mSubmoduleParameterLists.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Module_Parameter_Lists::end() -> decltype( mSubmoduleParameterLists.end() )
    {
        return mSubmoduleParameterLists.end();
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris
