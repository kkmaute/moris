/*
 * cl_FEM_MM_Perfect_Gas.cpp
 *
 *  Created on: Feb 2, 2020
 *  Author: wunsch
 */

#include "cl_FEM_MM_Perfect_Gas.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------

        MM_Perfect_Gas::MM_Perfect_Gas()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( MM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "IsochoricHeatCapacity" ] = static_cast< uint >( MM_Property_Type::ISOCHORIC_HEAT_CAPACITY ); // constant property
            mPropertyMap[ "SpecificGasConstant" ]   = static_cast< uint >( MM_Property_Type::SPECIFIC_GAS_CONSTANT );   // constant property
            
        }

        //--------------------------------------------------------------------------------------------------------------

        void MM_Perfect_Gas::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Material_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if( tDofString == "Density" )
                {
                    mDofDensity = tDofType;
                }
                else if( tDofString == "Pressure" )
                {
                    mDofPressure = tDofType;
                }
                else if( tDofString == "Temperature" )
                {
                    mDofTemperature = tDofType;
                }
                else
                {
                    std::string tErrMsg =
                            std::string( "MM_Perfect_Gas::set_dof_type_list - Unknown aDofString : ") +
                            tDofString;
                    MORIS_ERROR( false , tErrMsg.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::set_local_properties()
        {
            // get the isochoric heat capacity properties
            mPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

            // get the specific gas constant properties
            mPropSpecificGasConstant = get_property( "SpecificGasConstant" );
        }

        //------------------------------------------------------------------------------
        // SPECIFIC INTERNAL ENERGY (FIRST EQUATION OF STATE)
        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_Eint()
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // compute internal energy
            mEint = tIsochoricHeatCapacity * this->temperature()( 0 );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_EintDot()
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // compute internal energy
            mEintDot = tIsochoricHeatCapacity * this->TemperatureDot()( 0 );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_dEintdx()
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // compute internal energy
            mdEintdx = tIsochoricHeatCapacity * this->dnTemperaturedxn( 1 )( 0 );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Eintdx2()
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // compute internal energy
            md2Eintdx2 = tIsochoricHeatCapacity * this->dnTemperaturedxn( 2 )( 0 );
        }        

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_EintDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // compute internal energy
            mEintDof( tDofIndex ) = tIsochoricHeatCapacity * this->TemperatureDOF( aDofTypes );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_EintDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );            

            // compute internal energy
            mEintDotDof( tDofIndex ) = tIsochoricHeatCapacity * this->TemperatureDotDOF( aDofTypes );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_dEintdxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );            

            // compute internal energy
            mdEintdxDof( tDofIndex ) = tIsochoricHeatCapacity * this->dnTemperaturedxnDOF( aDofTypes, 1 );
        }       

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Eintdx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the heat capacity
            real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );            

            // compute internal energy
            md2Eintdx2Dof( tDofIndex ) = tIsochoricHeatCapacity * this->dnTemperaturedxnDOF( aDofTypes, 2 );
        }        

        //------------------------------------------------------------------------------
        // DENSITY (SECOND EQUATION OF STATE)
        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_density()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // compute density as function of pressure and temperature
            mDensity = tSpecificGasConstant * this->temperature() / this->pressure()( 0 );
        }

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_DensityDot()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get pressure value
            real tPressure = this->pressure()( 0 );

            // compute density as function of pressure and temperature
            mDensityDot = tSpecificGasConstant * ( 
                this->TemperatureDot() / tPressure +
                this->temperature()( 0 ) * this->PressureDot() / tPressure / tPressure );
        }     

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_dDensitydx()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get pressure value
            real tPressure = this->pressure()( 0 );

            // compute density as function of pressure and temperature
            mdDensitydx = tSpecificGasConstant * ( 
                this->dnTemperaturedxn( 1 ) / tPressure +
                this->temperature()( 0 ) * this->dnPressuredxn( 1 ) / tPressure / tPressure );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Densitydx2()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Densitydx2 - Not implemented yet." );
        }  

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_DensityDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get pressure value
            real tPressure = this->pressure()( 0 );

            // compute density as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mDensityDof( tDofIndex ) = tSpecificGasConstant * this->TemperatureDOF( aDofTypes ) / tPressure;

            if ( aDofTypes( 0 ) == mDofPressure )
                mDensityDof( tDofIndex ) = tSpecificGasConstant * this->temperature()( 0 ) * this->PressureDOF( aDofTypes ) / tPressure / tPressure;
        }   

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_DensityDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get pressure value
            real tPressure = this->pressure()( 0 );

            // compute density as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mDensityDof( tDofIndex ) = tSpecificGasConstant / tPressure * ( 
                    this->TemperatureDotDOF( aDofTypes ) - this->PressureDot()( 0 ) * this->TemperatureDOF( aDofTypes ) / tPressure );

            if ( aDofTypes( 0 ) == mDofPressure )
                mDensityDof( tDofIndex ) = tSpecificGasConstant / tPressure / tPressure * ( 
                    ( 2.0 * this->temperature()( 0 ) / tPressure - this->TemperatureDot()( 0 ) ) * this->PressureDOF( aDofTypes ) - 
                    this->temperature()( 0 ) * this->PressureDotDOF( aDofTypes ) );
        }       

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_dDensitydxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get pressure value
            real tPressure = this->pressure()( 0 );

            // compute density as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mDensityDof( tDofIndex ) = tSpecificGasConstant / tPressure * ( 
                    this->dnTemperaturedxnDOF( aDofTypes, 1 ) - this->dnPressuredxn( 1 ) * this->TemperatureDOF( aDofTypes ) / tPressure );

            if ( aDofTypes( 0 ) == mDofPressure )
                mDensityDof( tDofIndex ) = tSpecificGasConstant / tPressure / tPressure * ( 
                    2.0 * this->temperature()( 0 ) / tPressure * this->dnPressuredxn( 1 ) * this->PressureDOF( aDofTypes ) - 
                    this->dnTemperaturedxn( 1 ) * this->PressureDOF( aDofTypes ) - 
                    this->temperature()( 0 ) * this->dnPressuredxnDOF( aDofTypes, 1 ) );
        }

        //------------------------------------------------------------------------------
        // PRESSURE (SECOND EQUATION OF STATE)
        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_pressure()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // compute pressure as function of pressure and temperature
            mPressure = tSpecificGasConstant * this->density() * this->temperature();
        }

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_PressureDot()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // compute pressure as function of pressure and temperature
            mPressureDot = tSpecificGasConstant * ( 
                this->DensityDot() * this->temperature() + 
                this->density() * this->TemperatureDot() );
        }     

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_dPressuredx()
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // compute pressure as function of pressure and temperature
            mdPressuredx = tSpecificGasConstant * ( 
                this->temperature() * this->dnDensitydxn( 1 ) + 
                this->density() * this->dnTemperaturedxn( 1 ) );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Pressuredx2()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Pressuredx2 - Not implemented yet." );
        }  

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_PressureDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // compute pressure as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * this->density() * this->TemperatureDOF( aDofTypes );

            if ( aDofTypes( 0 ) == mDofDensity )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * this->temperature() * this->DensityDOF( aDofTypes );
        }   

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_PressureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // compute pressure as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * ( 
                    this->DensityDot()( 0 ) * this->TemperatureDOF( aDofTypes ) + 
                    this->density()( 0 ) * this->TemperatureDotDOF( aDofTypes ) );

            if ( aDofTypes( 0 ) == mDofDensity )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * ( 
                    this->temperature() * this->DensityDotDOF( aDofTypes ) + 
                    this->TemperatureDot() * this->DensityDOF( aDofTypes ) );
        }       

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_dPressuredxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // compute pressure as function of pressure and temperature
            if ( aDofTypes( 0 ) == mDofTemperature )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * ( 
                    this->dnDensitydxn( 1 ) * this->TemperatureDOF( aDofTypes ) + 
                    this->density() * this->dnTemperaturedxnDOF( aDofTypes, 1 ) );

            if ( aDofTypes( 0 ) == mDofDensity )
                mPressureDof( tDofIndex ) = tSpecificGasConstant * ( 
                    this->temperature() * this->dnDensitydxnDOF( aDofTypes, 1 ) + 
                    this->dnTemperaturedxn( 1 ) * this->DensityDOF( aDofTypes ) );
        }        

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Pressuredx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Pressuredx2DOF - Not implemented yet." );
        }    

        //------------------------------------------------------------------------------
        // TEMPERATURE (SECOND EQUATION OF STATE)
        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_temperature()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_temperature - Not implemented yet." );
        }

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_TemperatureDot()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDot - Not implemented yet." );
        }     

        //------------------------------------------------------------------------------          

        void MM_Perfect_Gas::eval_dTemperaturedx()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_dTemperaturedx - Not implemented yet." );
        }

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Temperaturedx2()
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Temperaturedx2 - Not implemented yet." );
        }  

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_TemperatureDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDOF - Not implemented yet." );
        }   

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_TemperatureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_TemperatureDotDOF - Not implemented yet." );
        }       

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_dTemperaturedxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_dTemperaturedxDOF - Not implemented yet." );
        }        

        //------------------------------------------------------------------------------

        void MM_Perfect_Gas::eval_d2Temperaturedx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME: skip for now as not needed
            MORIS_ERROR( false, "MM_Perfect_Gas::eval_d2Temperaturedx2DOF - Not implemented yet." );
        }          

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
