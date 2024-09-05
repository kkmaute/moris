/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_MM_Perfect_Gas.cpp
 *
 */

#include "cl_FEM_MM_Perfect_Gas.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    MM_Perfect_Gas::MM_Perfect_Gas()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( MM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "IsochoricHeatCapacity" ] = static_cast< uint >( MM_Property_Type::ISOCHORIC_HEAT_CAPACITY );    // constant property
        mPropertyMap[ "SpecificGasConstant" ]   = static_cast< uint >( MM_Property_Type::SPECIFIC_GAS_CONSTANT );      // constant property
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > > &aDofTypes,
            const Vector< std::string >             &aDofStrings )
    {
        // set dof type list
        Material_Model::set_dof_type_list( aDofTypes );

        // loop over the provided dof type
        for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
        {
            // get dof type string
            const std::string &tDofString = aDofStrings( iDof );

            // get dof type
            MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

            // switch on dof type string
            if ( tDofString == "Density" )
            {
                mDofDensity = tDofType;
            }
            else if ( tDofString == "Pressure" )
            {
                mDofPressure = tDofType;
            }
            else if ( tDofString == "Temperature" )
            {
                mDofTemperature = tDofType;
            }
            else
            {
                MORIS_ERROR( false, "MM_Perfect_Gas::set_dof_type_list - Unknown aDofString : %s", tDofString.c_str() );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::set_local_properties()
    {
        // get the isochoric heat capacity properties
        mPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the specific gas constant properties
        mPropSpecificGasConstant = get_property( "SpecificGasConstant" );
    }

    //------------------------------------------------------------------------------
    // SPECIFIC INTERNAL ENERGY (FIRST EQUATION OF STATE)
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_Eint()
    {
        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute internal energy
        mEint = tIsochoricHeatCapacity * this->temperature();
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_EintDot()
    {
        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute internal energy
        mEintDot = tIsochoricHeatCapacity * this->TemperatureDot();
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dEintdx()
    {
        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute internal energy
        mdEintdx = tIsochoricHeatCapacity * this->dnTemperaturedxn( 1 );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Eintdx2()
    {
        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute internal energy
        md2Eintdx2 = tIsochoricHeatCapacity * this->dnTemperaturedxn( 2 );
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_EintDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the heat capacity
        std::shared_ptr< Property > tPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute internal energy
        mEintDof( tDofIndex ) = tPropIsochoricHeatCapacity->val()( 0 ) * this->TemperatureDOF( aDofTypes );

        // dof dependency of the property on the DoF type
        if ( tPropIsochoricHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            mEintDof( tDofIndex ) += this->temperature()( 0 ) * tPropIsochoricHeatCapacity->dPropdDOF( aDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_EintDotDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the heat capacity
        std::shared_ptr< Property > tPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute internal energy
        mEintDotDof( tDofIndex ) = tPropIsochoricHeatCapacity->val()( 0 ) * this->TemperatureDotDOF( aDofTypes );

        // dof dependency of the property on the DoF type
        if ( tPropIsochoricHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            mEintDotDof( tDofIndex ) += this->TemperatureDot()( 0 ) * tPropIsochoricHeatCapacity->dPropdDOF( aDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dEintdxDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the heat capacity
        std::shared_ptr< Property > tPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute internal energy
        mdEintdxDof( tDofIndex ) = tPropIsochoricHeatCapacity->val()( 0 ) * this->dnTemperaturedxnDOF( aDofTypes, 1 );

        // dof dependency of the property on the DoF type
        if ( tPropIsochoricHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            mdEintdxDof( tDofIndex ) += this->dnTemperaturedxn( 1 ) * tPropIsochoricHeatCapacity->dPropdDOF( aDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Eintdx2DOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the heat capacity
        std::shared_ptr< Property > tPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute internal energy
        md2Eintdx2Dof( tDofIndex ) = tPropIsochoricHeatCapacity->val()( 0 ) * this->dnTemperaturedxnDOF( aDofTypes, 2 );

        // dof dependency of the property on the DoF type
        if ( tPropIsochoricHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            md2Eintdx2Dof( tDofIndex ) += this->dnTemperaturedxn( 2 ) * tPropIsochoricHeatCapacity->dPropdDOF( aDofTypes );
        }
    }

    //------------------------------------------------------------------------------
    // DENSITY (SECOND EQUATION OF STATE)
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_density()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // compute density as function of pressure and temperature
        mDensity = this->pressure() / ( tSpecificGasConstant * this->temperature()( 0 ) );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_DensityDot()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get temperature value
        real tTemperature = this->temperature()( 0 );

        // compute density as function of pressure and temperature
        mDensityDot = ( 1.0 / ( tSpecificGasConstant * tTemperature ) ) * (    //
                              this->PressureDot() -                            //
                              this->pressure()( 0 ) * this->TemperatureDot() / tTemperature );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dDensitydx()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get temperature value
        real tTemperature = this->temperature()( 0 );

        // compute density as function of pressure and temperature
        mdDensitydx = ( 1.0 / ( tSpecificGasConstant * tTemperature ) ) * (    //
                              this->dnPressuredxn( 1 ) -                       //
                              this->pressure()( 0 ) * this->dnTemperaturedxn( 1 ) / tTemperature );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Densitydx2()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Densitydx2 - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_DensityDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get temperature value
        real tTemperature = this->temperature()( 0 );

        // compute density as function of pressure and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            mDensityDof( tDofIndex ) = ( -this->pressure()( 0 ) / ( tSpecificGasConstant * tTemperature * tTemperature ) ) *    //
                                       this->TemperatureDOF( aDofTypes );
        }

        if ( aDofTypes( 0 ) == mDofPressure )
        {
            mDensityDof( tDofIndex ) = ( 1.0 / ( tSpecificGasConstant * tTemperature ) ) *    //
                                       this->PressureDOF( aDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_DensityDotDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get temperature value
        real tTemperature = this->temperature()( 0 );

        // compute density as function of pressure and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            mDensityDotDof( tDofIndex ) = 1.0 / ( tSpecificGasConstant * tTemperature * tTemperature ) * (                                                              //
                                                  ( 2.0 * this->pressure()( 0 ) * this->TemperatureDot()( 0 ) / tTemperature ) * this->TemperatureDOF( aDofTypes ) -    //
                                                  this->PressureDot()( 0 ) * this->TemperatureDOF( aDofTypes ) -                                                        //
                                                  this->pressure()( 0 ) * this->TemperatureDotDOF( aDofTypes ) );
        }

        if ( aDofTypes( 0 ) == mDofPressure )
        {
            mDensityDotDof( tDofIndex ) = ( 1.0 / ( tSpecificGasConstant * tTemperature ) ) * (    //
                                                  this->PressureDotDOF( aDofTypes ) -              //
                                                  ( this->TemperatureDot()( 0 ) / tTemperature ) * this->PressureDOF( aDofTypes ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dDensitydxDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get temperature value
        real tTemperature = this->temperature()( 0 );

        // compute density as function of pressure and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            mdDensitydxDof( tDofIndex ) = 1.0 / ( tSpecificGasConstant * tTemperature * tTemperature ) * (                                                              //
                                                  ( 2.0 * this->pressure()( 0 ) * this->dnTemperaturedxn( 1 ) / tTemperature ) * this->TemperatureDOF( aDofTypes ) -    //
                                                  this->dnPressuredxn( 1 ) * this->TemperatureDOF( aDofTypes ) -                                                        //
                                                  this->pressure()( 0 ) * dnTemperaturedxnDOF( aDofTypes, 1 ) );
        }

        if ( aDofTypes( 0 ) == mDofPressure )
        {
            mdDensitydxDof( tDofIndex ) = ( 1.0 / ( tSpecificGasConstant * tTemperature ) ) * (    //
                                                  this->dnPressuredxnDOF( aDofTypes, 1 ) -         //
                                                  ( this->dnTemperaturedxn( 1 ) / tTemperature ) * this->PressureDOF( aDofTypes ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Densitydx2DOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Densitydx2DOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    // PRESSURE (SECOND EQUATION OF STATE)
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_pressure()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // compute pressure as function of density and temperature
        mPressure = tSpecificGasConstant * this->density() * this->temperature();
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_PressureDot()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // compute pressure as function of density and temperature
        mPressureDot = tSpecificGasConstant * (                              //
                               this->DensityDot() * this->temperature() +    //
                               this->density() * this->TemperatureDot() );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dPressuredx()
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // compute pressure as function of density and temperature
        mdPressuredx = tSpecificGasConstant * (                                        //
                               this->temperature()( 0 ) * this->dnDensitydxn( 1 ) +    //
                               this->density()( 0 ) * this->dnTemperaturedxn( 1 ) );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Pressuredx2()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Pressuredx2 - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_PressureDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute pressure as function of density and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
            mPressureDof( tDofIndex ) = tSpecificGasConstant * this->density() * this->TemperatureDOF( aDofTypes );

        if ( aDofTypes( 0 ) == mDofDensity )
            mPressureDof( tDofIndex ) = tSpecificGasConstant * this->temperature() * this->DensityDOF( aDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_PressureDotDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute pressure as function of density and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
            mPressureDotDof( tDofIndex ) = tSpecificGasConstant * (                                                 //
                                                   this->DensityDot()( 0 ) * this->TemperatureDOF( aDofTypes ) +    //
                                                   this->density()( 0 ) * this->TemperatureDotDOF( aDofTypes ) );

        if ( aDofTypes( 0 ) == mDofDensity )
            mPressureDotDof( tDofIndex ) = tSpecificGasConstant * (                                            //
                                                   this->temperature() * this->DensityDotDOF( aDofTypes ) +    //
                                                   this->TemperatureDot() * this->DensityDOF( aDofTypes ) );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dPressuredxDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute pressure as function of density and temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
            mdPressuredxDof( tDofIndex ) = tSpecificGasConstant * (                                                 //
                                                   this->dnDensitydxn( 1 ) * this->TemperatureDOF( aDofTypes ) +    //
                                                   this->density()( 0 ) * this->dnTemperaturedxnDOF( aDofTypes, 1 ) );

        if ( aDofTypes( 0 ) == mDofDensity )
            mdPressuredxDof( tDofIndex ) = tSpecificGasConstant * (                                                      //
                                                   this->temperature()( 0 ) * this->dnDensitydxnDOF( aDofTypes, 1 ) +    //
                                                   this->dnTemperaturedxn( 1 ) * this->DensityDOF( aDofTypes ) );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Pressuredx2DOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Pressuredx2DOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    // TEMPERATURE (SECOND EQUATION OF STATE)
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_temperature()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_temperature - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_TemperatureDot()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDot - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dTemperaturedx()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_dTemperaturedx - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Temperaturedx2()
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Temperaturedx2 - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_TemperatureDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_TemperatureDotDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDotDOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_dTemperaturedxDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_dTemperaturedxDOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_d2Temperaturedx2DOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // FIXME: skip for now as not needed
        MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Temperaturedx2DOF - Not implemented yet." );
    }

    //------------------------------------------------------------------------------
    // THERMODYNAMIC QUANTITIES
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_VolumeExpansivity()
    {
        // compute value 1/T
        mAlphaP = { 1.0 / this->temperature()( 0 ) };
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_VolumeExpansivityDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute value -N/T^2
        mAlphaPDof( tDofIndex ) = -1.0 * this->TemperatureDOF( aDofTypes ) / std::pow( this->temperature()( 0 ), 2.0 );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_IsothermalCompressibility()
    {
        // compute value 1/p
        mBetaT = { 1.0 / this->pressure()( 0 ) };
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_IsothermalCompressibilityDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // compute value -N/p^2
        mBetaTDof( tDofIndex ) = -1.0 * this->PressureDOF( aDofTypes ) / std::pow( this->pressure()( 0 ), 2.0 );
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_Cv()
    {
        MORIS_ERROR( false, " MM_Perfect_Gas::eval_Cv - This function is not implemented yet. " );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_CvDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        MORIS_ERROR( false, " MM_Perfect_Gas::eval_CvDOF - This function is not implemented yet. " );
    }

    //------------------------------------------------------------------------------

    // Cp as function of Cv
    void
    MM_Perfect_Gas::eval_Cp()
    {
        // evaluate Cp
        mCp = this->Cv() + get_property( "SpecificGasConstant" )->val();
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_CpDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // evaluate Cp
        mCpDof( tDofIndex ) = this->CvDOF( aDofTypes );
    }

    //------------------------------------------------------------------------------

    // Cp as function of Cv
    void
    MM_Perfect_Gas::eval_Gamma()
    {
        // evaluate Gamma
        mGamma = ( this->Cv() + get_property( "SpecificGasConstant" )->val() ) / ( this->Cv()( 0 ) );
    }

    //------------------------------------------------------------------------------

    void
    MM_Perfect_Gas::eval_GammaDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get cv value
        real tCv = this->Cv()( 0 );

        // compute the DoF derivative
        mGammaDof( tDofIndex ) = ( -1.0 * get_property( "SpecificGasConstant" )->val()( 0 ) / ( tCv * tCv ) ) * this->CvDOF( aDofTypes );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
