/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Constitutive_Model.cpp
 *
 */

#include "cl_FEM_Constitutive_Model.hpp"

#include <utility>
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    Constitutive_Model::Constitutive_Model()
    {
        // set storage for evaluation
        mdStraindx.resize( mMaxSpaceDerOrder );

        // set flag for evaluation
        mdStraindxEval.set_size( mMaxSpaceDerOrder, 1, true );
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::print_names()
    {
        std::cout << "----------" << '\n';
        std::cout << "CM: " << mName << '\n';

        // properties
        for ( uint iProp = 0; iProp < mProperties.size(); iProp++ )
        {
            if ( mProperties( iProp ) != nullptr )
            {
                std::cout << "Property: " << mProperties( iProp )->get_name() << '\n';
            }
        }
        std::cout << "----------" << '\n';
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::reset_eval_flags()
    {
        // reset the flux value and derivative flags
        mFluxEval = true;
        mdFluxdDofEval.fill( true );

        // reset the divergence of the flux value and derivative flags
        mDivFluxEval = true;
        mddivfluxduEval.fill( true );

        // reset the gradient of the divergence of the flux value and derivative flags
        mGradDivFluxEval = true;
        mGradDivFluxDofEval.fill( true );

        // reset the traction value and derivative flags
        mTractionEval = true;
        mdTractiondDofEval.fill( true );

        // reset the test traction value and derivative flags
        mTestTractionEval.fill( true );
        mTestTractionTransEval.fill( true );
        mdTestTractiondDofEval.fill( true );

        // reset the stress value and derivative flags
        mStressEval = true;
        mdStressdDofEval.fill( true );

        // reset the strain value and derivative flags
        mStrainEval = true;
        mdStraindDofEval.fill( true );
        mdStraindxEval.fill( true );

        // reset the divergence of the strain value and derivative flags
        mDivStrainEval = true;
        mddivstrainduEval.fill( true );

        // reset the test strain value and derivative flags
        mTestStrainEval      = true;
        mTestStrainTransEval = true;
        mdTestStraindDofEval.fill( true );

        // reset the constitutive matrix value and derivative flags
        mConstEval = true;
        mdConstdDofEval.fill( true );

        // reset the energy value and derivative flags
        mEnergyEval = true;
        mEnergyDofEval.fill( true );

        // reset the rate of the energy value and derivative flags
        mEnergyDotEval = true;
        mEnergyDotDofEval.fill( true );

        // reset the gradient of the energy value and derivative flags
        mGradEnergyEval = true;
        mGradEnergyDofEval.fill( true );

        // reset the rate of the gradient of the energy value and derivative flags
        mGradEnergyDotEval = true;
        mGradEnergyDotDofEval.fill( true );

        // reset underlying properties
        for ( const std::shared_ptr< Property >& tProp : mProperties )
        {
            if ( tProp != nullptr )
            {
                tProp->reset_eval_flags();
            }
        }

        // reset underlying material model
        for ( const std::shared_ptr< Material_Model >& tMaterialModel : mMaterialModels )
        {
            if ( tMaterialModel != nullptr )
            {
                tMaterialModel->reset_eval_flags();
            }
        }

        // reset evaluation flags for specific constitutive model
        this->reset_specific_eval_flags();
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes )
    {
        // set the dof types
        mDofTypes = aDofTypes;

        // build a map for the dof types
        this->build_dof_type_map();
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // set parameters
        mParameters = aParameters;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_property(
            std::shared_ptr< fem::Property > aProperty,
            const std::string&               aPropertyString )
    {
        // check that aPropertyString makes sense
        MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                "Constitutive_Model::set_property - CM %s - Unknown aPropertyString : %s \n",
                mName.c_str(),
                aPropertyString.c_str() );

        // set the property in the property cell
        mProperties( mPropertyMap[ aPropertyString ] ) = std::move( aProperty );
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< fem::Property >& Constitutive_Model::get_property(
            const std::string& aPropertyString )
    {
        // check that aPropertyString makes sense
        MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                "Constitutive_Model::get_property - CM %s - Unknown aPropertyString : %s \n",
                mName.c_str(),
                aPropertyString.c_str() );

        // get the property in the property cell
        return mProperties( mPropertyMap[ aPropertyString ] );
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_material_model(
            std::shared_ptr< fem::Material_Model > aMaterialModel,
            const std::string&                     aMaterialModelString )
    {
        // check that aMaterialModelString makes sense
        MORIS_ERROR( mMaterialModelMap.find( aMaterialModelString ) != mMaterialModelMap.end(),
                "Constitutive_Model::set_material_model - CM %s - Unknown aMaterialModelString : %s \n",
                mName.c_str(),
                aMaterialModelString.c_str() );

        // set the MM in the MM cell
        mMaterialModels( mMaterialModelMap[ aMaterialModelString ] ) = std::move( aMaterialModel );
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< fem::Material_Model >& Constitutive_Model::get_material_model(
            const std::string& aMaterialModelString )
    {
        // check that aPropertyString makes sense
        MORIS_ERROR( mMaterialModelMap.find( aMaterialModelString ) != mPropertyMap.end(),
                "Constitutive_Model::get_material_model - CM %s - Unknown aMaterialModelString : %s \n",
                mName.c_str(),
                aMaterialModelString.c_str() );

        // get the material model in the material model cell
        return mMaterialModels( mMaterialModelMap[ aMaterialModelString ] );
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_dof_type_map()
    {
        // get number of dof types
        uint tNumDofTypes = mDofTypes.size();

        // determine the max Dof_Type enum
        sint tMaxEnum = 0;
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDOF )( 0 ) ) );
        }
        tMaxEnum++;

        // set map size
        mDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over the dof types
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            // fill the dof type map
            mDofTypeMap( static_cast< int >( mDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_dof_type_list()
    {
        // get number of dof types
        uint tNumDofTypes = mDofTypes.size();

        // set the size of the dof type list
        uint tCounterMax = tNumDofTypes;

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dof_type_list().size();
            }
        }
        mGlobalDofTypes.resize( tCounterMax );
        Vector< sint > tCheckList( tCounterMax, -1 );

        // initialize total dof counter
        uint tCounter = 0;

        // get active dof type for constitutive model
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            tCheckList( tCounter )      = static_cast< uint >( mDofTypes( iDOF )( 0 ) );
            mGlobalDofTypes( tCounter ) = mDofTypes( iDOF );
            tCounter++;
        }

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get active dof types
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tProperty->get_dof_type_list();

                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        tCheckList( tCounter )      = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
                        mGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dof type groups, i.e. the number of interpolators
        mGlobalDofTypes.resize( tCounter );

        // number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();
        uint tNumDirectDofTypes = mDofTypes.size();

        // set flag for evaluation
        mEnergyDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mEnergyDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mGradEnergyDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mGradEnergyDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mGradDivFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mTestTractionEval.set_size( tNumDirectDofTypes, 1, true );
        mTestTractionTransEval.set_size( tNumDirectDofTypes, 1, true );
        mdFluxdDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mddivfluxduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdTractiondDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mdTestTractiondDofEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mdStressdDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mdStraindDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mdTestStraindDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mddivstrainduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdConstdDofEval.set_size( tNumGlobalDofTypes, 1, true );

        // set storage for evaluation
        mEnergyDof.resize( tNumGlobalDofTypes );
        mEnergyDotDof.resize( tNumGlobalDofTypes );
        mGradEnergyDotDof.resize( tNumGlobalDofTypes );
        mGradEnergyDof.resize( tNumGlobalDofTypes );
        mGradDivFluxDof.resize( tNumGlobalDofTypes );
        mTestTraction.resize( tNumDirectDofTypes );
        mTestTractionTrans.resize( tNumDirectDofTypes );
        mdFluxdDof.resize( tNumGlobalDofTypes );
        mddivfluxdu.resize( tNumGlobalDofTypes );
        mdTractiondDof.resize( tNumGlobalDofTypes );
        mdTestTractiondDof.resize( tNumDirectDofTypes );
        for ( uint iDirectDof = 0; iDirectDof < tNumDirectDofTypes; iDirectDof++ )
        {
            mdTestTractiondDof( iDirectDof ).resize( tNumGlobalDofTypes );
        }
        mdStressdDof.resize( tNumGlobalDofTypes );
        mdStraindDof.resize( tNumGlobalDofTypes );
        mdTestStraindDof.resize( tNumGlobalDofTypes );
        mddivstraindu.resize( tNumGlobalDofTypes );
        mdConstdDof.resize( tNumGlobalDofTypes );

        // initialize storage variables specific to child CMs
        this->initialize_spec_storage_vars_and_eval_flags();
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< MSI::Dof_Type > >& Constitutive_Model::get_global_dof_type_list()
    {
        if ( mGlobalDofBuild )
        {
            // build the stabilization parameter global dof type list
            this->build_global_dof_type_list();

            // update build flag
            mGlobalDofBuild = false;
        }

        if ( mGlobalDofMapBuild )
        {
            // build the stabilization parameter global dof type map
            this->build_global_dof_type_map();

            // update build flag
            mGlobalDofMapBuild = false;
        }

        return mGlobalDofTypes;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_dof_type_map()
    {
        if ( mGlobalDofBuild )
        {
            // build the stabilization parameter global dof type list
            this->build_global_dof_type_list();

            // update build flag
            mGlobalDofBuild = false;
        }

        // get number of global dof types
        uint tNumDofTypes = mGlobalDofTypes.size();

        // determine the max Dof_Type enum
        sint tMaxEnum = 0;
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Dof_Type map size
        mGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Dof_Type map
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            // fill the DoF map
            mGlobalDofTypeMap( static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDSMat >& Constitutive_Model::get_global_dof_type_map()
    {
        if ( mGlobalDofMapBuild )
        {
            // build the stabilization parameter global dof type map
            this->build_global_dof_type_map();

            // update build flag
            mGlobalDofMapBuild = false;
        }

        return mGlobalDofTypeMap;
    }

    //------------------------------------------------------------------------------

    bool Constitutive_Model::check_dof_dependency(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // set bool for dependency
        bool tDofDependency = false;

        // get dof type index
        uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

        // if aDofType is an active dv type for the constitutive model
        if ( tDofIndex < this->get_global_dof_type_map().numel() && this->get_global_dof_type_map()( tDofIndex ) != -1 )
        {
            // bool is set to true
            tDofDependency = true;
        }

        // return bool for dependency
        return tDofDependency;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_dv_type_list(
            const Vector< Vector< gen::PDV_Type > >& aDvTypes )
    {
        // set the dv types
        mDvTypes = aDvTypes;

        // build a map for the dv types
        this->build_dv_type_map();
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_dv_type_map()
    {
        // get number of dv types
        uint tNumDvTypes = mDvTypes.size();

        // determine the max Dv_Type enum
        sint tMaxEnum = 0;
        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDvTypes( iDV )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Dv_Type map size
        mDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over the dv types
        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            // fill the dv type map
            mDvTypeMap( static_cast< int >( mDvTypes( iDV )( 0 ) ), 0 ) = iDV;
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_dv_type_list()
    {
        // get number of dv types
        uint tNumDvTypes = mDvTypes.size();

        // set the size of the dv type list
        uint tCounterMax = tNumDvTypes;

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dv_type_list().size();
            }
        }
        mGlobalDvTypes.resize( tCounterMax );
        Vector< sint > tCheckList( tCounterMax, -1 );

        // initialize total dv counter
        uint tCounter = 0;

        // get active dv type for constitutive model
        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            tCheckList( tCounter )     = static_cast< uint >( mDvTypes( iDV )( 0 ) );
            mGlobalDvTypes( tCounter ) = mDvTypes( iDV );
            tCounter++;
        }

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get active dv types
                Vector< Vector< gen::PDV_Type > > tActiveDvType = tProperty->get_dv_type_list();

                for ( uint iDV = 0; iDV < tActiveDvType.size(); iDV++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDV )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        tCheckList( tCounter )     = static_cast< uint >( tActiveDvType( iDV )( 0 ) );
                        mGlobalDvTypes( tCounter ) = tActiveDvType( iDV );
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dv type groups, i.e. the number of interpolators
        mGlobalDvTypes.resize( tCounter );

        // build global dv type map
        this->build_global_dv_type_map();

        // number of dof types
        uint tNumGlobalDvTypes = mGlobalDvTypes.size();
        uint tNumDirectDvTypes = mDvTypes.size();

        // set flag for evaluation
        mdFluxdDvEval.set_size( tNumGlobalDvTypes, 1, true );
        mdTractiondDvEval.set_size( tNumGlobalDvTypes, true );
        mdTestTractiondDvEval.set_size( tNumDirectDvTypes, tNumGlobalDvTypes, true );
        mdStraindDvEval.set_size( tNumGlobalDvTypes, true );
        mdConstdDvEval.set_size( tNumGlobalDvTypes, true );

        // set storage for evaluation
        mdFluxdDv.resize( tNumGlobalDvTypes );
        mdTractiondDv.resize( tNumGlobalDvTypes );
        mdTestTractiondDv.resize( mDvTypes.size() );
        for ( uint iDirectDv = 0; iDirectDv < mDvTypes.size(); iDirectDv++ )
        {
            mdTestTractiondDv( iDirectDv ).resize( tNumGlobalDvTypes );
        }
        mdStraindDv.resize( tNumGlobalDvTypes );
        mdConstdDv.resize( tNumGlobalDvTypes );
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< gen::PDV_Type > >& Constitutive_Model::get_global_dv_type_list()
    {
        if ( mGlobalDvBuild )
        {
            // build the stabilization parameter global dv type list
            this->build_global_dv_type_list();

            // build the stabilization parameter global dv type map
            this->build_global_dv_type_map();

            // update build flag
            mGlobalDvBuild = false;
        }

        return mGlobalDvTypes;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_dv_type_map()
    {
        // get number of global dv types
        uint tNumDvTypes = mGlobalDvTypes.size();

        // determine the max Dv_Type enum
        sint tMaxEnum = 0;
        for ( uint iDOF = 0; iDOF < tNumDvTypes; iDOF++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDvTypes( iDOF )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Dv_Type map size
        mGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over the dv types
        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            // fill the dv type map
            mGlobalDvTypeMap( static_cast< int >( mGlobalDvTypes( iDV )( 0 ) ), 0 ) = iDV;
        }
    }

    //------------------------------------------------------------------------------

    bool Constitutive_Model::check_dv_dependency(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // set bool for dependency
        bool tDvDependency = false;

        // get dv type index
        uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

        // if aDvType is an active dv type for the constitutive model
        if ( tDvIndex < mGlobalDvTypeMap.numel() && mGlobalDvTypeMap( tDvIndex ) != -1 )
        {
            // bool is set to true
            tDvDependency = true;
        }
        // return bool for dependency
        return tDvDependency;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_field_type_list(
            const Vector< Vector< mtk::Field_Type > >& aFieldTypes )
    {
        // set the field types
        mFieldTypes = aFieldTypes;

        // build a map for the field types
        this->build_field_type_map();
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_field_type_map()
    {
        // get number of field types
        uint tNumFieldTypes = mFieldTypes.size();

        // determine the max Dof_Type enum
        sint tMaxEnum = 0;
        for ( uint iField = 0; iField < tNumFieldTypes; iField++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFieldTypes( iField )( 0 ) ) );
        }
        tMaxEnum++;

        // set map size
        mFieldTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over the field types
        for ( uint iField = 0; iField < tNumFieldTypes; iField++ )
        {
            // fill the field type map
            mFieldTypeMap( static_cast< int >( mFieldTypes( iField )( 0 ) ), 0 ) = iField;
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_field_type_list()
    {
        // get number of field types
        uint tNumFieldTypes = mFieldTypes.size();

        // set the size of the field type list
        uint tCounterMax = tNumFieldTypes;

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_field_type_list().size();
            }
        }
        mGlobalFieldTypes.resize( tCounterMax );
        Vector< sint > tCheckList( tCounterMax, -1 );

        // initialize total dv counter
        uint tCounter = 0;

        // get active Field type for constitutive model
        for ( uint iField = 0; iField < tNumFieldTypes; iField++ )
        {
            tCheckList( tCounter )        = static_cast< uint >( mFieldTypes( iField )( 0 ) );
            mGlobalFieldTypes( tCounter ) = mFieldTypes( iField );
            tCounter++;
        }

        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get active Field types
                Vector< Vector< mtk::Field_Type > > tActiveFieldType = tProperty->get_field_type_list();

                for ( uint iField = 0; iField < tActiveFieldType.size(); iField++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveFieldType( iField )( 0 ) ) );
                    }

                    // if Field enum not in the list
                    if ( !tCheck )
                    {
                        tCheckList( tCounter )        = static_cast< uint >( tActiveFieldType( iField )( 0 ) );
                        mGlobalFieldTypes( tCounter ) = tActiveFieldType( iField );
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique Field type groups, i.e. the number of interpolators
        mGlobalFieldTypes.resize( tCounter );

        //            // build global Field type map
        //            this->build_global_field_type_map();
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< mtk::Field_Type > >& Constitutive_Model::get_global_field_type_list()
    {
        if ( mGlobalFieldBuild )
        {
            // build the stabilization parameter global field type list
            this->build_global_field_type_list();

            // update build flag
            mGlobalFieldBuild = false;
        }

        if ( mGlobalFieldMapBuild )
        {
            // build the stabilization parameter global dv type map
            this->build_global_field_type_map();

            // update build flag
            mGlobalFieldMapBuild = false;
        }

        return mGlobalFieldTypes;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::build_global_field_type_map()
    {
        if ( mGlobalFieldBuild )
        {
            // build the stabilization parameter global dof type list
            this->build_global_field_type_list();

            // update build flag
            mGlobalFieldBuild = false;
        }

        // get number of global dof types
        uint tNumFieldTypes = mGlobalFieldTypes.size();

        // determine the max Field_Type enum
        sint tMaxEnum = 0;
        for ( uint iField = 0; iField < tNumFieldTypes; iField++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalFieldTypes( iField )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Field_Type map size
        mGlobalFieldTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Field_Type map
        for ( uint iField = 0; iField < tNumFieldTypes; iField++ )
        {
            // fill the Field map
            mGlobalFieldTypeMap( static_cast< int >( mGlobalFieldTypes( iField )( 0 ) ), 0 ) = iField;
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDSMat >& Constitutive_Model::get_global_field_type_map()
    {
        if ( mGlobalFieldMapBuild )
        {
            // build the global field type map
            this->build_global_field_type_map();

            // update build flag
            mGlobalFieldMapBuild = false;
        }

        return mGlobalFieldTypeMap;
    }

    //------------------------------------------------------------------------------

    bool Constitutive_Model::check_field_dependency(
            const Vector< mtk::Field_Type >& aFieldType )
    {
        // set bool for dependency
        bool tFieldDependency = false;

        // get dof type index
        uint tFieldIndex = static_cast< uint >( aFieldType( 0 ) );

        // if aDofType is an active dv type for the constitutive model
        if ( tFieldIndex < this->get_global_field_type_map().numel() && this->get_global_field_type_map()( tFieldIndex ) != -1 )
        {
            // bool is set to true
            tFieldDependency = true;
        }

        // return bool for dependency
        return tFieldDependency;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::set_field_interpolator_manager(
            Field_Interpolator_Manager* aFieldInterpolatorManager )
    {
        // set the field interpolator manager for the constitutive model
        mFIManager = aFieldInterpolatorManager;

        // loop over the underlying properties
        for ( const std::shared_ptr< Property >& tProp : this->get_properties() )
        {
            if ( tProp != nullptr )
            {
                // set the field interpolator manager for the property
                tProp->set_field_interpolator_manager( mFIManager );
            }
        }

        // loop over the underlying material models
        for ( const std::shared_ptr< Material_Model >& tMM : this->get_material_models() )
        {
            if ( tMM != nullptr )
            {
                // set the field interpolator manager for the property
                tMM->set_field_interpolator_manager( mFIManager );
            }
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::get_non_unique_dof_types(
            Vector< MSI::Dof_Type >& aDofTypes )
    {
        // initialize dof counter
        uint tCounter = 0;

        // loop over direct dof dependencies
        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
        {
            // update counter
            tCounter += mDofTypes( iDOF ).size();
        }

        // loop over properties
        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get property dof type list
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tProperty->get_dof_type_list();

                // loop over property dof types
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // update counter
                    tCounter += tActiveDofType( iDOF ).size();
                }
            }
        }

        // loop over material models
        for ( const std::shared_ptr< Material_Model >& tMM : mMaterialModels )
        {
            if ( tMM != nullptr )
            {
                // get MM dof type list
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tMM->get_dof_type_list();

                // loop over MM dof types
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // update counter
                    tCounter += tActiveDofType( iDOF ).size();
                }
            }
        }

        // reserve memory for the non unique dof type list
        aDofTypes.reserve( tCounter );

        // loop over direct dof dependencies
        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
        {
            // populate the dof type list
            aDofTypes.append( mDofTypes( iDOF ) );
        }

        // loop over the properties
        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get property dof type list
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tProperty->get_dof_type_list();

                // loop over property dof types
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // populate dof type list
                    aDofTypes.append( tActiveDofType( iDOF ) );
                }
            }
        }

        // loop over the material models
        for ( const std::shared_ptr< Material_Model >& tMM : mMaterialModels )
        {
            if ( tMM != nullptr )
            {
                // get MM dof type list
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tMM->get_dof_type_list();

                // loop over MM dof types
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // populate dof type list
                    aDofTypes.append( tActiveDofType( iDOF ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::get_non_unique_dof_dv_and_field_types(
            Vector< MSI::Dof_Type >&   aDofTypes,
            Vector< gen::PDV_Type >&   aDvTypes,
            Vector< mtk::Field_Type >& aFieldTypes )
    {
        // initialize dof counter
        uint tDofCounter   = 0;
        uint tDvCounter    = 0;
        uint tFieldCounter = 0;

        // loop over direct dof dependencies
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            // update counter
            tDofCounter += mDofTypes( iDof ).size();
        }

        // loop over direct dv dependencies
        for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
        {
            // update counter
            tDvCounter += mDvTypes( iDv ).size();
        }

        // loop over direct field dependencies
        for ( uint iField = 0; iField < mFieldTypes.size(); iField++ )
        {
            // update counter
            tFieldCounter += mFieldTypes( iField ).size();
        }

        // loop over properties
        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get property dof type list
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;
                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update counter
                tDofCounter += tActiveDofTypes.size();
                tDvCounter += tActiveDvTypes.size();
                tFieldCounter += tActiveDvTypes.size();
            }
        }

        // loop over material models
        for ( const std::shared_ptr< Material_Model >& tMM : mMaterialModels )
        {
            if ( tMM != nullptr )
            {
                // get MM dof type list
                Vector< MSI::Dof_Type > tActiveDofTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // update counter
                tDofCounter += tActiveDofTypes.size();
            }
        }

        // reserve memory for the non unique dof and dv types
        aDofTypes.reserve( tDofCounter );
        aDvTypes.reserve( tDvCounter );
        aFieldTypes.reserve( tFieldCounter );

        // loop over direct dof dependencies
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            // populate the dof type list
            aDofTypes.append( mDofTypes( iDof ) );
        }

        // loop over direct dv dependencies
        for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
        {
            // populate the dv type list
            aDvTypes.append( mDvTypes( iDv ) );
        }

        // loop over direct field dependencies
        for ( uint iField = 0; iField < mFieldTypes.size(); iField++ )
        {
            // populate the dv type list
            aFieldTypes.append( mFieldTypes( iField ) );
        }

        // loop over the properties
        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != nullptr )
            {
                // get property dof and dv type list
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // populate the dof and dv type lists
                aDofTypes.append( tActiveDofTypes );
                aDvTypes.append( tActiveDvTypes );
                aFieldTypes.append( tActiveFieldTypes );
            }
        }

        // loop over the MMs
        for ( const std::shared_ptr< Material_Model >& tMM : mMaterialModels )
        {
            if ( tMM != nullptr )
            {
                // get MM dof and dv type list
                Vector< MSI::Dof_Type > tActiveDofTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // populate the dof type lists
                aDofTypes.append( tActiveDofTypes );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Constitutive_Model::eval_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            Matrix< DDRMat >&              aDerivativeFD,
            const Vector< MSI::Dof_Type >& aDofTypes,
            real                           aPerturbation,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check requested derivative
        MORIS_ERROR( aCMRequestType != CM_Request_Type::END_CM_REQUEST_TYPE,
                "Constitutive_Model::eval_derivative_FD - aCMRequestType needs to be defined." );

        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // reset properties
        this->reset_eval_flags();

        // evaluate unperturbed value
        const Matrix< DDRMat >& tUnperturbed = this->select_derivative_FD(
                aCMRequestType,
                aTestDofTypes,
                aNormal,
                aJump,
                aCMFunctionType );

        // set size for derivative
        aDerivativeFD.set_size( tUnperturbed.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        const Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed contribution to derivative
                    aDerivativeFD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tUnperturbed / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // assemble the derivative
                    aDerivativeFD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->select_derivative_FD( aCMRequestType, aTestDofTypes, aNormal, aJump, aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }

        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        this->set_derivative_FD(
                aCMRequestType,
                aDerivativeFD,
                aDofTypes,
                aTestDofTypes,
                aCMFunctionType );
    }

    /**
     * select derivative wrt to a dof type
     * @param[ out ] aCMRequestType
     */
    const Matrix< DDRMat >&
    Constitutive_Model::select_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMRequestType )
        {
            case CM_Request_Type::STRAIN:
            {
                return this->strain( aCMFunctionType );
                break;
            }
            case CM_Request_Type::FLUX:
            {
                return this->flux( aCMFunctionType );
                break;
            }
            case CM_Request_Type::TRACTION:
            {
                return this->traction( aNormal, aCMFunctionType );
                break;
            }
            case CM_Request_Type::TEST_TRACTION:
            {
                mTraction = this->testTraction_trans( aNormal, aTestDofTypes, aCMFunctionType ) * aJump;
                return mTraction;
                break;
            }
            default:
                MORIS_ERROR( false, "Constitutive_Model::select_derivative_FD: aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    void
    Constitutive_Model::set_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            Matrix< DDRMat >&              aDerivativeFD,
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        switch ( aCMRequestType )
        {
            case CM_Request_Type::STRAIN:
            {
                // set value to storage
                mdStraindDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::FLUX:
            {
                // set value to storage
                mdFluxdDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TRACTION:
            {
                // set value to storage
                mdTractiondDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TEST_TRACTION:
            {
                // get the test dof index
                uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

                // set value to storage
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) = aDerivativeFD;
                break;
            }
            default:
                MORIS_ERROR( false, "Constitutive_Model::set_derivative_FD: aCMRequestType undefined" );
        }
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dFluxdDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adFluxdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // evaluate unperturbed flux
        Matrix< DDRMat > tFlux = this->flux( aCMFunctionType );

        // set size for derivative
        adFluxdDOF_FD.set_size( tFlux.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed flux contribution to dfluxdu
                    adFluxdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tFlux / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // assemble the jacobian
                    adFluxdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->flux( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mdFluxdDof( tDofIndex ) = adFluxdDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dtractiondu_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adtractiondu_FD,
            real                           aPerturbation,
            Matrix< DDRMat >&              aNormal,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFIDerivative =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
        uint tDerNumFields = tFIDerivative->get_number_of_fields();

        // compute unperturbed traction
        Matrix< DDRMat > tTraction = this->traction( aNormal, aCMFunctionType );

        // set size for derivative
        adtractiondu_FD.set_size( tTraction.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed traction contribution to dtractiondu
                    adtractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tTraction / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFIDerivative->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add contribution of point to dtractiondu
                    adtractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->traction( aNormal, aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFIDerivative->set_coeff( tCoeff );

        // set value for storage
        mdTractiondDof( tDofIndex ) = adtractiondu_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dtesttractiondu_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            Matrix< DDRMat >&              adtesttractiondu_FD,
            real                           aPerturbation,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the derivative dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFIDerivative =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
        uint tDerNumFields = tFIDerivative->get_number_of_fields();

        // compute unperturbed test traction
        Matrix< DDRMat > tUnperturbedTestTraction =
                trans( this->testTraction( aNormal, aTestDofTypes, aCMFunctionType ) ) * aJump;
        adtesttractiondu_FD.set_size(
                tUnperturbedTestTraction.n_rows(),
                tDerNumDof,
                0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed test traction contribution to dtesttractiondu
                    adtesttractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tUnperturbedTestTraction / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFIDerivative->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add unperturbed test traction contribution to dtesttractiondu
                    adtesttractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * trans( this->testTraction( aNormal, aTestDofTypes, aCMFunctionType ) ) * aJump / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFIDerivative->set_coeff( tCoeff );

        // set value for storage
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ) = adtesttractiondu_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dtesttractiondu_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            Matrix< DDRMat >&              adtesttractiondu_FD,
            real                           aPerturbation,
            const Matrix< DDRMat >&        aNormal,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the derivative dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFIDerivative =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
        uint tDerNumFields = tFIDerivative->get_number_of_fields();

        // compute unperturbed test traction
        Matrix< DDRMat > tUnperturbedTestTraction =
                this->testTraction( aNormal, aTestDofTypes, aCMFunctionType );
        adtesttractiondu_FD.set_size(
                tUnperturbedTestTraction.n_rows(),
                tDerNumDof,
                0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed test traction contribution to dtesttractiondu
                    adtesttractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tUnperturbedTestTraction / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFIDerivative->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add unperturbed test traction contribution to dtesttractiondu
                    adtesttractiondu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->testTraction( aNormal, aTestDofTypes, aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFIDerivative->set_coeff( tCoeff );

        // set value for storage
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ) = adtesttractiondu_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_ddivfluxdu_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              addivfluxdu_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFIDerivative =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
        uint tDerNumFields = tFIDerivative->get_number_of_fields();

        // compute unperturbed divflux
        Matrix< DDRMat > tDivFlux = this->divflux( aCMFunctionType );

        // set size for derivative
        addivfluxdu_FD.set_size( tDivFlux.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed div flux contribution to ddivfluxdu
                    addivfluxdu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tDivFlux / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFIDerivative->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add contribution of point to ddivfluxdu
                    addivfluxdu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->divflux( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }

        // reset the coefficients values
        tFIDerivative->set_coeff( tCoeff );

        // set value for storage
        mddivfluxdu( tDofIndex ) = addivfluxdu_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_ddivstraindu_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              addivstraindu_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFIDerivative =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
        uint tDerNumFields = tFIDerivative->get_number_of_fields();

        // compute unperturbed div strain
        Matrix< DDRMat > tDivStrain = this->divstrain( aCMFunctionType );

        // set size for derivative
        addivstraindu_FD.set_size( tDivStrain.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed div strain contribution to ddivstraindu
                    addivstraindu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tDivStrain / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFIDerivative->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add contribution of point to ddivstraindu
                    addivstraindu_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->divstrain( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFIDerivative->set_coeff( tCoeff );

        // set value for storage
        mddivstraindu( tDofIndex ) = addivstraindu_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dEnergydDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adEnergydDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute the unperturbed energy
        Matrix< DDRMat > tEnergy = this->Energy( aCMFunctionType );

        // set size for derivative
        adEnergydDOF_FD.set_size( tEnergy.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed energy contribution to denergydu
                    adEnergydDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tEnergy / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficents
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // pertub the coefficent
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add contribution of point to denergydu
                    adEnergydDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->Energy( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mEnergyDof( tDofIndex ) = adEnergydDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dEnergyDotdDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adEnergyDotdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute unperturbed energy dot
        Matrix< DDRMat > tEnergyDot = this->EnergyDot( aCMFunctionType );

        // set size for derivative
        adEnergyDotdDOF_FD.set_size( tEnergyDot.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed energydot contribution to denergydotdu
                    adEnergyDotdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tEnergyDot / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add point contribution to dEnergyDotdu
                    adEnergyDotdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->EnergyDot( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mEnergyDotDof( tDofIndex ) = adEnergyDotdDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dGradEnergydDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adGradEnergydDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute unperturbed grad energy
        Matrix< DDRMat > tGradEnergy = this->gradEnergy( aCMFunctionType );

        // set size for derivative
        adGradEnergydDOF_FD.set_size( tGradEnergy.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed grad energy contribution to dgradenergydu
                    adGradEnergydDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tGradEnergy / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add point contribution to dGradEnergydu
                    adGradEnergydDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->gradEnergy( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mGradEnergyDof( tDofIndex ) = adGradEnergydDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dGradEnergyDotdDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adGradEnergyDotdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute unperturbed grad energy dot
        Matrix< DDRMat > tGradEnergyDot = this->gradEnergyDot( aCMFunctionType );

        // set size for derivative
        adGradEnergyDotdDOF_FD.set_size( tGradEnergyDot.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed grad energy dot contribution to dgradenergydotdu
                    adGradEnergyDotdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tGradEnergyDot / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add point contribution to dGradHdu
                    adGradEnergyDotdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->gradEnergyDot( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mGradEnergyDotDof( tDofIndex ) = adGradEnergyDotdDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dGradDivFluxdDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adGradDivFluxdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute unperturbed grad div flux
        Matrix< DDRMat > tGradDivFlux = this->graddivflux( aCMFunctionType );

        // set size for derivative
        adGradDivFluxdDOF_FD.set_size( tGradDivFlux.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed grad div flux dot contribution to dgraddivfluxdu
                    adGradDivFluxdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tGradDivFlux / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add point contribution to dgraddivfluxdu
                    adGradDivFluxdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->graddivflux( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mGradDivFluxDof( tDofIndex ) = adGradDivFluxdDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dStraindDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adStraindDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the derivative dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of leader dofs wrt which derivative is computed
        uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // compute unperturbed strain
        Matrix< DDRMat > tStrain = this->strain( aCMFunctionType );

        // set size for derivative
        adStraindDOF_FD.set_size( tStrain.n_rows(), tDerNumDof, 0.0 );

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dof counter
        uint tDofCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed strain dot contribution to dstraindu
                    adStraindDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tStrain / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // add point contribution to dstraindu
                    adStraindDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->strain( aCMFunctionType ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mdStraindDof( tDofIndex ) = adStraindDOF_FD;
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dFluxdDV_FD(
            const Vector< gen::PDV_Type >& aDvTypes,
            Matrix< DDRMat >&              adFluxdDV_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDvTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDv     = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // set size for derivative
        uint tNumRow = this->flux( aCMFunctionType ).n_rows();
        adFluxdDV_FD.set_size( tNumRow, tDerNumDv, 0.0 );

        // coefficients for dv type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dv counter
        uint tDvCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // perturbation of the coefficient
                Matrix< DDRMat > tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset constitutive model
                this->reset_eval_flags();

                // evaluate the residual
                Matrix< DDRMat > tFlux_Plus = this->flux( aCMFunctionType );

                // perturbation of the coefficient
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset constitutive model
                this->reset_eval_flags();

                // evaluate the residual
                Matrix< DDRMat > tFlux_Minus = this->flux( aCMFunctionType );

                // evaluate Jacobian
                adFluxdDV_FD.get_column( tDvCounter ) =
                        ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                // update dv counter
                tDvCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );
    }

    //------------------------------------------------------------------------------

    void Constitutive_Model::eval_dStraindDV_FD(
            const Vector< gen::PDV_Type >& aDvTypes,
            Matrix< DDRMat >&              adStraindDV_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );

        // get the field interpolator for type
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDvTypes( 0 ) );

        // get number of coefficients, fields and bases for the considered FI
        uint tDerNumDv     = tFI->get_number_of_space_time_coefficients();
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // set size for derivative
        uint tNumRow = this->strain( aCMFunctionType ).n_rows();
        adStraindDV_FD.set_size( tNumRow, tDerNumDv, 0.0 );

        // coefficients for dv type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // initialize dv counter
        uint tDvCounter = 0;

        // loop over coefficients columns
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // perturbation of the coefficient
                Matrix< DDRMat > tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset constitutive model
                this->reset_eval_flags();

                // evaluate the residual
                Matrix< DDRMat > tStrain_Plus = this->strain( aCMFunctionType );

                // perturbation of the coefficient
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset constitutive model
                this->reset_eval_flags();

                // evaluate the residual
                Matrix< DDRMat > tStrain_Minus = this->strain( aCMFunctionType );

                // evaluate Jacobian
                adStraindDV_FD.get_column( tDvCounter ) =
                        ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                // update dv counter
                tDvCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::flux( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::flux - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mFluxEval )
        {
            // evaluate the flux
            this->eval_flux();

            // set bool for evaluation
            mFluxEval = false;
        }
        // return the flux value
        return mFlux;
    }

    const Matrix< DDRMat >&
    Constitutive_Model::flux(
            int                   aFlatType,
            enum CM_Function_Type aCMFunctionType )
    {
        // no parent implementation
        MORIS_ERROR( false, "Constitutive_Model::flux - Not implemented in parent class." );

        // return the flux value
        return mFlux;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Constitutive_Model::dTestStraindDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Nonlinear_Isotropic::dTestStraindDOF - Only DEFAULT CM function type known in base class." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative of the teststrain wrt. DOF was not evaluated
        if ( mdTestStraindDofEval( tDofIndex ) )
        {
            // evaluate the derivative of the teststrain wrt. DOF
            this->eval_dTestStraindDOF( aDofTypes );

            // set bool for evaluation
            mdTestStraindDofEval( tDofIndex ) = false;
        }
        // return the derivative of the teststrain wrt. DOF values
        return mdTestStraindDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::divflux( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::divflux - Only DEFAULT CM function type known in base class." );

        // if the divergence of the flux was not evaluated
        if ( mDivFluxEval )
        {
            // evaluate the divergence of the flux
            this->eval_divflux();

            // set bool for evaluation
            mDivFluxEval = false;
        }
        // return the divergence of the flux value
        return mDivFlux;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::Energy( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::Energy - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mEnergyEval )
        {
            // evaluate the flux
            this->eval_Energy();

            // set bool for evaluation
            mEnergyEval = false;
        }

        // return the flux value
        return mEnergy;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::EnergyDot( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::EnergyDot - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mEnergyDotEval )
        {
            // evaluate the flux
            this->eval_EnergyDot();

            // set bool for evaluation
            mEnergyDotEval = false;
        }

        // return the flux value
        return mEnergyDot;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::gradEnergyDot( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::gradEnergyDot - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mGradEnergyDotEval )
        {
            // evaluate the flux
            this->eval_gradEnergyDot();

            // set bool for evaluation
            mGradEnergyDotEval = false;
        }

        // return the flux value
        return mGradEnergyDot;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::gradEnergy( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::gradEnergy - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mGradEnergyEval )
        {
            // evaluate the flux
            this->eval_gradEnergy();

            // set bool for evaluation
            mGradEnergyEval = false;
        }

        // return the flux value
        return mGradEnergy;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::graddivflux( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::graddivflux - Only DEFAULT CM function type known in base class." );

        // if the flux was not evaluated
        if ( mGradDivFluxEval )
        {
            // evaluate the flux
            this->eval_graddivflux();

            // set bool for evaluation
            mGradDivFluxEval = false;
        }

        // return the flux value
        return mGradDivFlux;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::ddivfluxdu(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::ddivfluxdu - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::ddivfluxdu - no dependency in this dof type." );

        // get the dof index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative of the divergence of the flux was not evaluated
        if ( mddivfluxduEval( tDofIndex ) )
        {
            // evaluate the derivative of the divergence of the flux
            this->eval_ddivfluxdu( aDofType );

            // set bool for evaluation
            mddivfluxduEval( tDofIndex ) = false;
        }

        // return the divergence of the flux value
        return mddivfluxdu( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::traction(
            const Matrix< DDRMat >& aNormal,
            enum CM_Function_Type   aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::traction - Only DEFAULT CM function type known in base class." );

        // if the traction was not evaluated
        if ( mTractionEval )
        {
            // evaluate the traction
            this->eval_traction( aNormal );

            // set bool for evaluation
            mTractionEval = false;
        }
        // return the traction value
        return mTraction;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::testTraction - Only DEFAULT CM function type known in base class." );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if the test traction was not evaluated
        if ( mTestTractionEval( tTestDofIndex ) )
        {
            // evaluate the test traction
            this->eval_testTraction( aNormal, aTestDofTypes );

            // set bool for evaluation
            mTestTractionEval( tTestDofIndex ) = false;
        }

        // return the test traction value
        return mTestTraction( tTestDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Constitutive_Model::testTraction_trans(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )

    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if the transpose of the test traction was not evaluated
        if ( mTestTractionTransEval( tTestDofIndex ) )
        {
            // evaluate the transpose of the test traction
            mTestTractionTrans( tTestDofIndex ) = trans( this->testTraction( aNormal, aTestDofTypes, aCMFunctionType ) );

            // set bool for evaluation
            mTestTractionTransEval( tTestDofIndex ) = false;
        }

        // return the transpose of the test traction value
        return mTestTractionTrans( tTestDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::stress( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::stress - Only DEFAULT CM function type known in base class." );

        // if the strain was not evaluated
        if ( mStressEval )
        {
            // evaluate the strain
            this->eval_stress();

            // set bool for evaluation
            mStressEval = false;
        }

        // return the strain value
        return mStress;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::strain( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::strain - Only DEFAULT CM function type known in base class." );

        // if the strain was not evaluated
        if ( mStrainEval )
        {
            // evaluate the strain
            this->eval_strain();

            // set bool for evaluation
            mStrainEval = false;
        }

        // return the strain value
        return mStrain;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::divstrain( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::divstrain - Only DEFAULT CM function type known in base class." );

        // if the divergence of the strain was not evaluated
        if ( mDivStrainEval )
        {
            // evaluate the divergence of the strain
            this->eval_divstrain();

            // set bool for evaluation
            mDivStrainEval = false;
        }

        // return the divergence of the strain value
        return mDivStrain;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::ddivstraindu(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::ddivstraindu - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::ddivstraindu - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative of the divergence of the strain was not evaluated
        if ( mddivstrainduEval( tDofIndex ) )
        {
            // evaluate the derivative of the divergence of the strain
            this->eval_ddivstraindu( aDofType );

            // set bool for evaluation
            mddivstrainduEval( tDofIndex ) = false;
        }

        // return the divergence of the strain value
        return mddivstraindu( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::testStrain( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::testStrain - Only DEFAULT CM function type known in base class." );

        // if the test strain was not evaluated
        if ( mTestStrainEval )
        {
            // evaluate the test strain
            this->eval_testStrain();

            // set bool for evaluation
            mTestStrainEval = false;
        }

        // return the test strain value
        return mTestStrain;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::constitutive( enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::constitutive - Only DEFAULT CM function type known in base class." );

        // if the constitutive matrix was not evaluated
        if ( mConstEval )
        {
            // evaluate the constitutive matrix
            this->eval_const();

            // set bool for evaluation
            mConstEval = false;
        }

        // return the constitutive matrix value
        return mConst;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::testStrain_trans( enum CM_Function_Type aCMFunctionType )
    {
        // if the test strain was not evaluated
        if ( mTestStrainTransEval )
        {
            // evaluate the test strain
            mTestStrainTrans = trans( this->testStrain( aCMFunctionType ) );

            // set bool for evaluation
            mTestStrainTransEval = false;
        }

        // return the test strain value
        return mTestStrainTrans;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dFluxdDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),    //
                "Constitutive_Model::dFluxdDOF - "         //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFluxdDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dFluxdDOF( aDofType );

            // set bool for evaluation
            mdFluxdDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdFluxdDof( tDofIndex );
    }

    //-----------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dEnergydDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dEnergydDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dEnergydDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mEnergyDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dEnergydDOF( aDofType );

            // set bool for evaluation
            mEnergyDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mEnergyDof( tDofIndex );
    }

    //-----------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dEnergyDotdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dEnergyDotdDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dEnergyDotdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mEnergyDotDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dEnergyDotdDOF( aDofType );

            // set bool for evaluation
            mEnergyDotDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mEnergyDotDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dGradEnergydDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dGradEnergydDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dGradEnergydDOF - no dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mGradEnergyDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dGradEnergydDOF( aDofType );

            // set bool for evaluation
            mGradEnergyDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mGradEnergyDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dGradEnergyDotdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dGradEnergydDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dGradEnergydDOF - no dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mGradEnergyDotDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dGradEnergyDotdDOF( aDofType );

            // set bool for evaluation
            mGradEnergyDotDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mGradEnergyDotDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dGradDivFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dGradDivFluxdDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dGradDivFluxdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mGradDivFluxDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dGradDivFluxdDOF( aDofType );

            // set bool for evaluation
            mGradDivFluxDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mGradDivFluxDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dTractiondDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dTractiondDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTractiondDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dTractiondDOF( aDofType, aNormal );

            // set bool for evaluation
            mdTractiondDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdTractiondDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dTestTractiondDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dTestTractiondDOF - no dependency in this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dTestTractiondDOF( aDofType, aNormal, aTestDofTypes );

            // set bool for evaluation
            mdTestTractiondDofEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdTestTractiondDof( tTestDofIndex )( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dTestTractiondDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dTestTractiondDOF - no dependency in this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

            // set bool for evaluation
            mdTestTractiondDofEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdTestTractiondDof( tTestDofIndex )( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dstraindx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dstraindx - Only DEFAULT CM function type known in base class." );

        MORIS_ERROR(
                aOrder == 1,
                "Constitutive_Model::dstraindx - Works only for 1st order derivative for now." );

        // if the derivative has not been evaluated yet
        if ( mdStraindxEval( aOrder - 1 ) )
        {
            // evaluate the derivative
            this->eval_dstraindx( aOrder );

            // set bool for evaluation
            mdStraindxEval( aOrder - 1 ) = false;
        }

        // return the derivative
        return mdStraindx( aOrder - 1 );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dStressdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dStressdDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dStressdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdStressdDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dStressdDOF( aDofType );

            // set bool for evaluation
            mdStressdDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdStressdDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dStraindDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dStraindDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dStraindDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdStraindDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dStraindDOF( aDofType );

            // set bool for evaluation
            mdStraindDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdStraindDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dConstdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dConstdDOF - Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "Constitutive_Model::dConstdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdConstdDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dConstdDOF( aDofType );

            // set bool for evaluation
            mdConstdDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdConstdDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dFluxdDV(
            const Vector< gen::PDV_Type >& aDvType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dFluxdDV - Only DEFAULT CM function type known in base class." );

        // if aDvType is not an active dv type
        MORIS_ERROR(
                this->check_dv_dependency( aDvType ),
                "Constitutive_Model::dFluxdDV - no dependency in this dv type." );

        // get the dv index
        uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFluxdDvEval( tDvIndex ) )
        {
            // evaluate the derivative
            this->eval_dFluxdDV( aDvType );

            // set bool for evaluation
            mdFluxdDvEval( tDvIndex ) = false;
        }

        // return the derivative
        return mdFluxdDv( tDvIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dStraindDV(
            const Vector< gen::PDV_Type >& aDvType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dStraindDV - Only DEFAULT CM function type known in base class." );

        // if aDvType is not an active dv type for the property
        MORIS_ERROR(
                this->check_dv_dependency( aDvType ),
                "Constitutive_Model::dStraindDV - no dependency in this dv type." );

        // get the dv index
        uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdStraindDvEval( tDvIndex ) )
        {
            // evaluate the derivative
            this->eval_dStraindDV( aDvType );

            // set bool for evaluation
            mdStraindDvEval( tDvIndex ) = false;
        }

        // return the derivative
        return mdStraindDv( tDvIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >& Constitutive_Model::dConstdDV(
            const Vector< gen::PDV_Type >& aDvType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "Constitutive_Model::dConstdDV - Only DEFAULT CM function type known in base class." );

        // if aDvType is not an active dv type for the property
        MORIS_ERROR(
                this->check_dv_dependency( aDvType ),
                "Constitutive_Model::dConstdDV - no dependency in this dv type." );

        // get the dv index
        uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdConstdDvEval( tDvIndex ) )
        {
            // evaluate the derivative
            this->eval_dConstdDV( aDvType );

            // set bool for evaluation
            mdConstdDvEval( tDvIndex ) = false;
        }

        // return the derivative
        return mdConstdDv( tDvIndex );
    }

}    // namespace moris::fem
