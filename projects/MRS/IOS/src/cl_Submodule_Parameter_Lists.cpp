/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Submodule_Parameter_Lists.cpp
 *
 */

#include "cl_Submodule_Parameter_Lists.hpp"
#include "parameters.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Submodule_Parameter_Lists::Submodule_Parameter_Lists( std::string aType )
            : mType( std::move( aType ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::string& Submodule_Parameter_Lists::get_type()
    {
        return mType;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Submodule_Parameter_Lists::size() const
    {
        return mParameterLists.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( const Parameter_List& aParameterList )
    {
        mParameterLists.push_back( aParameterList );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list()
    {
        if ( mType == HMR_Submodule_String::values( static_cast< uint >( HMR::LAGRANGE_MESHES ) ) )
        {
            mParameterLists.push_back( prm::create_lagrange_mesh_parameter_list() );
        }
        else if ( mType == HMR_Submodule_String::values( static_cast< uint >( HMR::BSPLINE_MESHES ) ) )
        {
            mParameterLists.push_back( prm::create_bspline_mesh_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::PROPERTIES ) ) )
        {
            mParameterLists.push_back( prm::create_property_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::CONSTITUTIVE_MODELS ) ) )
        {
            mParameterLists.push_back( prm::create_constitutive_model_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::STABILIZATION ) ) )
        {
            mParameterLists.push_back( prm::create_stabilization_parameter_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::IWG ) ) )
        {
            mParameterLists.push_back( prm::create_IWG_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::IQI ) ) )
        {
            mParameterLists.push_back( prm::create_IQI_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::FIELDS ) ) )
        {
            mParameterLists.push_back( prm::create_fem_field_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::PHASES ) ) )
        {
            mParameterLists.push_back( prm::create_phase_parameter_list() );
        }
        else if ( mType == FEM_Submodule_String::values( static_cast< uint >( FEM::MATERIAL_MODELS ) ) )
        {
            mParameterLists.push_back( prm::create_material_model_parameter_list() );
        }
        else if ( mType == SOL_Submodule_String::values( static_cast< uint >( SOL::LINEAR_SOLVERS ) ) )
        {
            mParameterLists.push_back( prm::create_linear_solver_parameter_list() );
        }
        else if ( mType == SOL_Submodule_String::values( static_cast< uint >( SOL::NONLINEAR_ALGORITHMS ) ) )
        {
            mParameterLists.push_back( prm::create_nonlinear_algorithm_parameter_list() );
        }
        else if ( mType == SOL_Submodule_String::values( static_cast< uint >( SOL::NONLINEAR_SOLVERS ) ) )
        {
            mParameterLists.push_back( prm::create_nonlinear_solver_parameter_list() );
        }
        else if ( mType == SOL_Submodule_String::values( static_cast< uint >( SOL::TIME_SOLVER_ALGORITHMS ) ) )
        {
            mParameterLists.push_back( prm::create_time_solver_algorithm_parameter_list() );
        }
        else if ( mType == SOL_Submodule_String::values( static_cast< uint >( SOL::TIME_SOLVERS ) ) )
        {
            mParameterLists.push_back( prm::create_time_solver_parameter_list() );
        }
        else
        {
            MORIS_ERROR( false, "The %s submodule does not support adding additional parameter lists.", mType.c_str() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( opt::Optimization_Algorithm_Type aOptimizationAlgorithmType )
    {
        // Check for correct submodule type
        this->check_submodule_type( OPT_Submodule_String::values( static_cast< uint >( OPT::ALGORITHMS ) ) );

        // Add new parameter list
        mParameterLists.push_back( prm::create_optimization_algorithm_parameter_list( aOptimizationAlgorithmType ) );
    }

    //--------------------------------------------------------------------------------------------------------------

//    void Submodule_Parameter_Lists::add_parameter_list( uint aGeometryType, gen::Field_Type aFieldType )
//    {
//        // Check for correct submodule type
//        this->check_submodule_type( GEN_Submodule_String::values( static_cast< uint >( GEN::GEOMETRIES ) ) );
//
//        // Add new parameter list TODO BRENDAN: add overall GEN function for this
//        MORIS_ERROR( false, "GEN geometries cannot be created automatically yet" );
//    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( gen::Field_Type aFieldType )
    {
        // Check for correct submodule type
        this->check_submodule_type( GEN_Submodule_String::values( static_cast< uint >( GEN::PROPERTIES ) ) );

        // Add new parameter list
        mParameterLists.push_back( prm::create_gen_property_parameter_list( aFieldType ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( sol::SolverType aSolverType )
    {
        // Check for correct submodule type
        this->check_submodule_type( SOL_Submodule_String::values( static_cast< uint >( SOL::LINEAR_ALGORITHMS ) ) );

        // Add new parameter list
        mParameterLists.push_back( prm::create_linear_algorithm_parameter_list( aSolverType ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( sol::PreconditionerType aPreconditionerType )
    {
        // Check for correct submodule type
        this->check_submodule_type( SOL_Submodule_String::values( static_cast< uint >( SOL::PRECONDITIONERS ) ) );

        // Add new parameter list
        mParameterLists.push_back( prm::create_preconditioner_parameter_list( aPreconditionerType ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::append( const Submodule_Parameter_Lists& aSubmoduleToAppend )
    {
        mParameterLists.append( aSubmoduleToAppend.mParameterLists );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Submodule_Parameter_Lists::empty()
    {
        return mParameterLists.empty();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::insert(
        const std::string&     aName,
        const Design_Variable& aDesignVariable )
    {
        // Insert into parameter list
        mParameterLists.back().insert( aName, aDesignVariable );
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List& Submodule_Parameter_Lists::operator()( uint aParameterListIndex )
    {
        return mParameterLists( aParameterListIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Parameter_List& Submodule_Parameter_Lists::operator()( uint aParameterListIndex ) const
    {
        return mParameterLists( aParameterListIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Submodule_Parameter_Lists::begin()->decltype( mParameterLists.begin() )
    {
        return mParameterLists.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Submodule_Parameter_Lists::end()->decltype( mParameterLists.end() )
    {
        return mParameterLists.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::check_submodule_type( const std::string& aExpectedType )
    {
        MORIS_ERROR( mType == aExpectedType, "A %s submodule cannot create a %s parameter list.", mType.c_str(), aExpectedType.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------
}
