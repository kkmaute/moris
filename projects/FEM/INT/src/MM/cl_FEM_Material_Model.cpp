/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Material_Model.cpp
 *
 */

#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Material_Model::Material_Model(){}

        //------------------------------------------------------------------------------

        void Material_Model::print_names()
        {
            std::cout<<"----------"<<std::endl;
            std::cout<<"MM: "<<mName<<std::endl;

            // properties
            for( uint iProp = 0; iProp < mProperties.size(); iProp++ )
            {
                if( mProperties( iProp ) != nullptr )
                {
                    std::cout<<"Property: "<<mProperties( iProp )->get_name()<<std::endl;
                }
            }
            std::cout<<"----------"<<std::endl;
        }

        //------------------------------------------------------------------------------

        void Material_Model::reset_eval_flags()
        {

            // reset eval flags for internal energy
            mEintEval = true;
            mEintDotEval = true;
            mdEintdxEval = true;
            md2Eintdx2Eval = true;

            mEintDofEval.fill( true );
            mEintDotDofEval.fill( true );
            mdEintdxDofEval.fill( true );
            md2Eintdx2DofEval.fill( true );

            // reset eval flags for density, if it is a dependent variable
            mDensityEval = true;
            mDensityDotEval = true;
            mdDensitydxEval = true;
            md2Densitydx2Eval = true;

            mDensityDofEval.fill( true );
            mDensityDotDofEval.fill( true );
            mdDensitydxDofEval.fill( true );
            md2Densitydx2DofEval.fill( true );

            // reset eval flags for pressure, if it is a dependent variable
            mPressureEval = true;
            mPressureDotEval = true;
            mdPressuredxEval = true;
            md2Pressuredx2Eval = true;

            mPressureDofEval.fill( true );
            mPressureDotDofEval.fill( true );
            mdPressuredxDofEval.fill( true );
            md2Pressuredx2DofEval.fill( true );

            // reset eval flags for temperature, if it is a dependent variable
            mTemperatureEval = true;
            mTemperatureDotEval = true;
            mdTemperaturedxEval = true;
            md2Temperaturedx2Eval = true;

            mTemperatureDofEval.fill( true );
            mTemperatureDotDofEval.fill( true );
            mdTemperaturedxDofEval.fill( true );
            md2Temperaturedx2DofEval.fill( true );

            // reset eval flags for thermodynamic quantities
            mAlphaPEval = true;
            mBetaTEval = true;
            mCvEval = true;
            mCpEval = true;
            mGammaEval = true;

            mAlphaPDofEval.fill( true );
            mBetaTDofEval.fill( true );
            mCvDofEval.fill( true );
            mCpDofEval.fill( true );
            mGammaDofEval.fill( true );

            // reset underlying properties
            for( const std::shared_ptr< Property > & tProp : mProperties )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset evaluation flags for specific material model
            this->reset_specific_eval_flags();
        }

        //------------------------------------------------------------------------------

        void Material_Model::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > aDofTypes )
        {
            // check list of DoF types and set flags for dependent variables accordingly
            // get number of dof types
            uint tNumDofTypes = aDofTypes.size();

            // check list of DoF types and set flags for dependent variables accordingly
            for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // get DoF type in list
                //MSI::Dof_Type tCurrentDofType = aDofTypes( iDOF )( 0 );
                int tCurrentDofType = static_cast< int >( aDofTypes( iDOF )( 0 ) );

                // check if DoF type corresponds to one of thermodynamic state variables
                switch ( tCurrentDofType )
                {
                    // Density DoF type found
                    case static_cast< int >( MSI::Dof_Type::RHO ):
                    {
                        // set trivial flag on density dof type
                        mDensityIsDependent = false;

                        // set function pointers
                        m_get_Density = &Material_Model::density_triv;
                        m_get_DensityDot = &Material_Model::DensityDot_triv;
                        m_get_dnDensitydxn = &Material_Model::dnDensitydxn_triv;
                        m_get_DensityDof = &Material_Model::DensityDOF_triv;
                        m_get_DensityDotDof = &Material_Model::DensityDotDOF_triv;
                        m_get_dnDensitydxnDof = &Material_Model::dnDensitydxnDOF_triv;

                        break;
                    }

                    // Pressure DoF Type Found
                    case static_cast< int >( MSI::Dof_Type::P ):
                    {
                        // set trivial flag on density dof type
                        mPressureIsDependent = false;

                        // set function pointers
                        m_get_Pressure = &Material_Model::pressure_triv;
                        m_get_PressureDot = &Material_Model::PressureDot_triv;
                        m_get_dnPressuredxn = &Material_Model::dnPressuredxn_triv;
                        m_get_PressureDof = &Material_Model::PressureDOF_triv;
                        m_get_PressureDotDof = &Material_Model::PressureDotDOF_triv;
                        m_get_dnPressuredxnDof = &Material_Model::dnPressuredxnDOF_triv;

                        break;
                    }

                    // Temperature DoF Type Found
                    case static_cast< int >( MSI::Dof_Type::TEMP ):
                    {
                        // set trivial flag on density dof type
                        mTemperatureIsDependent = false;

                        // set function pointers
                        m_get_Temperature = &Material_Model::temperature_triv;
                        m_get_TemperatureDot = &Material_Model::TemperatureDot_triv;
                        m_get_dnTemperaturedxn = &Material_Model::dnTemperaturedxn_triv;
                        m_get_TemperatureDof = &Material_Model::TemperatureDOF_triv;
                        m_get_TemperatureDotDof = &Material_Model::TemperatureDotDOF_triv;
                        m_get_dnTemperaturedxnDof = &Material_Model::dnTemperaturedxnDOF_triv;

                        break;
                    }

                    default:
                    {
                        // do nothing
                        break;
                    }
                }
            }

            // FIXME: check for primitive variables, 2 out of 3 required
            uint tNumTDVarsSet = 0;
            if ( !mDensityIsDependent ) tNumTDVarsSet++;
            if ( !mPressureIsDependent ) tNumTDVarsSet++;
            if ( !mTemperatureIsDependent ) tNumTDVarsSet++;
            MORIS_ERROR( tNumTDVarsSet == 2, "Material_Model::set_dof_type_list - 2 out of 3 primitive variables required." );

            // set the dof types
            mDofTypes = aDofTypes;

            // build a map for the dof types
            this->build_dof_type_map();

        }

        //------------------------------------------------------------------------------

        void Material_Model::set_property(
                std::shared_ptr< fem::Property > aProperty,
                std::string                      aPropertyString )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "Material_Model::set_property - MM %s - Unknown aPropertyString : %s \n",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // set the property in the property cell
            mProperties( mPropertyMap[ aPropertyString ] ) = aProperty;
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< fem::Property > & Material_Model::get_property( std::string aPropertyString )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "Material_Model::get_property - MM %s - Unknown aPropertyString : %s \n",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // get the property in the property cell
            return  mProperties( mPropertyMap[ aPropertyString ] );
        }

        //------------------------------------------------------------------------------

        void Material_Model::build_dof_type_map()
        {
            // get number of dof types
            uint tNumDofTypes = mDofTypes.size();

            // determine the max Dof_Type enum
            sint tMaxEnum = 0;
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDOF )( 0 ) ) );
            }
            tMaxEnum++;

            // set map size
            mDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // loop over the dof types
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // fill the dof type map
                mDofTypeMap( static_cast< int >( mDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }
        }

        //------------------------------------------------------------------------------

        void Material_Model::build_global_dof_type_list()
        {
            // get number of dof types
            uint tNumDofTypes = mDofTypes.size();

            // set the size of the dof type list
            uint tCounterMax = tNumDofTypes;

            for ( const std::shared_ptr< Property > & tProperty : mProperties )
            {
                if( tProperty != nullptr )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }
            }
            mGlobalDofTypes.resize( tCounterMax );
            Vector< sint > tCheckList( tCounterMax, -1 );

            // initialize total dof counter
            uint tCounter = 0;

            // get active dof type for material model
            for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                tCheckList( tCounter ) = static_cast< uint >( mDofTypes( iDOF )( 0 ) );
                mGlobalDofTypes( tCounter ) = mDofTypes( iDOF );
                tCounter++;
            }

            for ( const std::shared_ptr< Property > & tProperty : mProperties )
            {
                if( tProperty != nullptr )
                {
                    // get active dof types
                    const Vector< Vector< MSI::Dof_Type > > & tActiveDofType =
                            tProperty->get_dof_type_list();

                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
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
            //uint tNumDirectDofTypes = mDofTypes.size();

            // set flags for evaluation
            mEintDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mEintDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEintdxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            md2Eintdx2DofEval.set_size( tNumGlobalDofTypes, 1, true );

            // set storage for evaluation
            mEintDof.resize( tNumGlobalDofTypes );
            mEintDotDof.resize( tNumGlobalDofTypes );
            mdEintdxDof.resize( tNumGlobalDofTypes );
            md2Eintdx2Dof.resize( tNumGlobalDofTypes );

            // set flags for evaluation
            mDensityDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mDensityDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mdDensitydxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            md2Densitydx2DofEval.set_size( tNumGlobalDofTypes, 1, true );

            // set storage for evaluation
            mDensityDof.resize( tNumGlobalDofTypes );
            mDensityDotDof.resize( tNumGlobalDofTypes );
            mdDensitydxDof.resize( tNumGlobalDofTypes );
            md2Densitydx2Dof.resize( tNumGlobalDofTypes );

            // set flags for evaluation
            mPressureDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mPressureDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mdPressuredxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            md2Pressuredx2DofEval.set_size( tNumGlobalDofTypes, 1, true );

            // set storage for evaluation
            mPressureDof.resize( tNumGlobalDofTypes );
            mPressureDotDof.resize( tNumGlobalDofTypes );
            mdPressuredxDof.resize( tNumGlobalDofTypes );
            md2Pressuredx2Dof.resize( tNumGlobalDofTypes );

            // set flags for evaluation
            mTemperatureDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mTemperatureDotDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mdTemperaturedxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            md2Temperaturedx2DofEval.set_size( tNumGlobalDofTypes, 1, true );

            // set storage for evaluation
            mTemperatureDof.resize( tNumGlobalDofTypes );
            mTemperatureDotDof.resize( tNumGlobalDofTypes );
            mdTemperaturedxDof.resize( tNumGlobalDofTypes );
            md2Temperaturedx2Dof.resize( tNumGlobalDofTypes );

            // set flags for evaluation
            mAlphaPDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mBetaTDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mCvDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mCpDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mGammaDofEval.set_size( tNumGlobalDofTypes, 1, true );

            // set storage for evaluation
            mAlphaPDof.resize( tNumGlobalDofTypes );
            mBetaTDof.resize( tNumGlobalDofTypes );
            mCvDof.resize( tNumGlobalDofTypes );
            mCpDof.resize( tNumGlobalDofTypes );
            mGammaDof.resize( tNumGlobalDofTypes );

            // initialize storage variables specific to child MMs
            this->initialize_spec_storage_vars_and_eval_flags();
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< MSI::Dof_Type > > & Material_Model::get_global_dof_type_list()
        {
            if( mGlobalDofBuild )
            {
                // build the stabilization parameter global dof type list
                this->build_global_dof_type_list();

                // update build flag
                mGlobalDofBuild = false;
            }

            if( mGlobalDofMapBuild )
            {
                // build the stabilization parameter global dof type map
                this->build_global_dof_type_map();

                // update build flag
                mGlobalDofMapBuild = false;
            }

            return mGlobalDofTypes;
        }

        //------------------------------------------------------------------------------

        void Material_Model::build_global_dof_type_map()
        {
            if( mGlobalDofBuild )
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
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ) );
            }
            tMaxEnum++;

            // set the Dof_Type map size
            mGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // fill the Dof_Type map
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // fill the property map
                mGlobalDofTypeMap( static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Material_Model::get_global_dof_type_map()
        {
            if( mGlobalDofMapBuild )
            {
                // build the stabilization parameter global dof type map
                this->build_global_dof_type_map();

                // update build flag
                mGlobalDofMapBuild = false;
            }

            return mGlobalDofTypeMap;
        }

        //------------------------------------------------------------------------------

        bool Material_Model::check_dof_dependency(
                const Vector< MSI::Dof_Type > & aDofType )
        {
            // set bool for dependency
            bool tDofDependency = false;

            // get dof type index
            uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

            // if aDofType is an active dv type for the material model
            if( tDofIndex < this->get_global_dof_type_map().numel() && this->get_global_dof_type_map()( tDofIndex ) != -1 )
            {
                // bool is set to true
                tDofDependency = true;
            }

            // return bool for dependency
            return tDofDependency;
        }

        //------------------------------------------------------------------------------

        void Material_Model::set_field_interpolator_manager(
                Field_Interpolator_Manager * aFieldInterpolatorManager )
        {
            // set the field interpolator manager for the material model
            mFIManager = aFieldInterpolatorManager;

            // loop over the underlying properties
            for( const std::shared_ptr< Property > & tProp : this->get_properties() )
            {
                if (tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( mFIManager );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Material_Model::get_non_unique_dof_types(
                Vector< MSI::Dof_Type > & aDofTypes )
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
            for ( const std::shared_ptr< Property > & tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof type list
                    const Vector< Vector< MSI::Dof_Type > > & tActiveDofType =
                            tProperty->get_dof_type_list();

                    // loop over property dof types
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
            for ( const std::shared_ptr< Property > & tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof type list
                    const Vector< Vector< MSI::Dof_Type > > & tActiveDofType =
                            tProperty->get_dof_type_list();

                    // loop over property dof types
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // populate te dof type list
                        aDofTypes.append( tActiveDofType( iDOF ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        // RETURN FUNCTIONS FOR FIRST EQUATION OF STATE
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::Eint()
        {
            // if the flux was not evaluated
            if( mEintEval )
            {
                // evaluate the flux
                this->eval_Eint();

                // set bool for evaluation
                mEintEval = false;
            }

            // return the flux value
            return mEint;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::EintDot()
        {
            // if the flux was not evaluated
            if( mEintDotEval )
            {
                // evaluate the flux
                this->eval_EintDot();

                // set bool for evaluation
                mEintDotEval = false;
            }

            // return the flux value
            return mEintDot;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::dnEintdxn( uint aOrder )
        {
            switch ( aOrder )
            {
            case 1: // first derivative

                if ( mdEintdxEval )
                {
                    // evaluate the flux
                    this->eval_dEintdx();

                    // set bool for evaluation
                    mdEintdxEval = false;
                }

                // return the flux value
                return mdEintdx;

            case 2: // second derivative

                if ( md2Eintdx2Eval )
                {
                    // evaluate the flux
                    this->eval_d2Eintdx2();

                    // set bool for evaluation
                    md2Eintdx2Eval = false;
                }

                // return the flux value
                return md2Eintdx2;

            default:
                MORIS_ERROR( false, "Material_Model::dnEintdxn - aOrder unknown, only 1 and 2 supported." );
                return mdEintdx;
            }
        }

        //-----------------------------------------------------------------------------
        //-----------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::EintDOF( const Vector< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ASSERT(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::EintDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mEintDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_EintDOF( aDofType );

                // set bool for evaluation
                mEintDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mEintDof( tDofIndex );
        }

        //-----------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::EintDotDOF( const Vector< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ASSERT(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::EintDotDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mEintDotDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_EintDotDOF( aDofType );

                // set bool for evaluation
                mEintDotDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mEintDotDof( tDofIndex );
        }

        //-----------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::dnEintdxnDOF( const Vector< MSI::Dof_Type > & aDofType, uint aOrder )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ASSERT(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::dnEintdxnDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            switch ( aOrder )
            {
            case 1: // first derivative

                // if the derivative has not been evaluated yet
                if( mdEintdxDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dEintdxDOF( aDofType );

                    // set bool for evaluation
                    mdEintdxDofEval( tDofIndex ) = false;
                }

                // return the derivative
                return mdEintdxDof( tDofIndex );

            case 2: // second derivative

                // if the derivative has not been evaluated yet
                if( md2Eintdx2DofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_d2Eintdx2DOF( aDofType );

                    // set bool for evaluation
                    md2Eintdx2DofEval( tDofIndex ) = false;
                }

                // return the derivative
                return md2Eintdx2Dof( tDofIndex );

            default:
                MORIS_ERROR( false, "Material_Model::dnEintdxnDOF - aOrder unknown, only 1 and 2 supported." );
                return mdEintdxDof( 0 );
            }
        }

    //-----------------------------------------------------------------------------

    }/* end_fem_namespace */
}/* end_moris_namespace */

