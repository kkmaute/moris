/*
 * cl_FEM_CM_Fluid_Compressible_Ideal.cpp
 *
 *  Created on: Jul 24, 2020
 *  Author: wunsch
 */

#include "cl_FEM_CM_Fluid_Compressible_Ideal.hpp"
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

        CM_Fluid_Compressible_Ideal::CM_Fluid_Compressible_Ideal()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "IsochoricHeatCapacity" ] = CM_Property_Type::ISOCHORIC_HEAT_CAPACITY; // constant property
            mPropertyMap[ "SpecificGasConstant" ]   = CM_Property_Type::SPECIFIC_GAS_CONSTANT;   // constant property
            mPropertyMap[ "DynamicViscosity" ]      = CM_Property_Type::DYNAMIC_VISCOSITY;       // may be a fnct. of T
            mPropertyMap[ "ThermalConductivity" ]   = CM_Property_Type::THERMAL_CONDUCTIVITY;    // may be a fnct. of T
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::set_function_pointers()
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                        {
                    m_eval_strain            = &CM_Fluid_Compressible_Ideal::eval_strain_2d;
                    m_eval_teststrain        = &CM_Fluid_Compressible_Ideal::eval_teststrain_2d;
                    m_eval_velocitymatrix    = &CM_Fluid_Compressible_Ideal::eval_velocitymatrix_2d;
                    m_unfold_tensor          = &CM_Fluid_Compressible_Ideal::unfold_2d;
                    mFlatIdentity = { { 1.0 }, { 1.0 }, { 0.0 } };
                    break;
                        }
                case ( 3 ):
                        {
                    m_eval_strain            = &CM_Fluid_Compressible_Ideal::eval_strain_3d;
                    m_eval_teststrain        = &CM_Fluid_Compressible_Ideal::eval_teststrain_3d;
                    m_eval_velocitymatrix    = &CM_Fluid_Compressible_Ideal::eval_velocitymatrix_3d;
                    m_unfold_tensor          = &CM_Fluid_Compressible_Ideal::unfold_3d;
                    mFlatIdentity = { { 1.0 }, { 1.0 }, { 1.0 }, { 0.0 }, { 0.0 }, { 0.0 } };
                    break;
                        }
                default :
                {
                    MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::set_function_pointers - this function is currently unused, might be used in the future." );
                    break;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if( tDofString == "Velocity" )
                {
                    mDofVelocity = tDofType;
                }
                else if( tDofString == "Density" )
                {
                    mDofDensity = tDofType;
                }
                else if( tDofString == "Temperature" )
                {
                    mDofTemperature = tDofType;
                }
                else
                {
                    std::string tErrMsg =
                            std::string( "CM_Fluid_Compressible_Ideal::set_dof_type_list - Unknown aDofString : ") +
                            tDofString;
                    MORIS_ERROR( false , tErrMsg.c_str() );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::set_property(
                std::shared_ptr< fem::Property > aProperty,
                std::string                      aPropertyString )
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Fluid_Compressible_Ideal::set_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // set the property in the property cell
            mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr< Property > CM_Fluid_Compressible_Ideal::get_property(
                std::string aPropertyString )
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Fluid_Compressible_Ideal::get_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // get the property in the property cell
            return  mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_flux()
        {
            // get the velocity
            Matrix< DDRMat > tVelocity =  mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            // evaluate the thermal flux
            this->eval_thermal_flux();

            // evaluate the velocity matrix
            this->eval_velocityMatrix();

            // compute contribution
            mFlux = mVelocityMatrix * this->stress() -
                    this->Energy() * tVelocity -
                    mThermalFlux;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the velocity
            Field_Interpolator * tFIVelocity  =  mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // evaluate the velocity matrix
            this->eval_velocityMatrix();

            // unfold the flattened stress tensor
            Matrix< DDRMat > tStressTensor;
            this->unfold( this->stress() , tStressTensor );


            // initialize the matrix
            mdFluxdDof( tDofIndex ).set_size( mSpaceDim,
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                // evaluate thermal flux DoF derivative
                this->eval_thermal_dFluxdDOF( aDofTypes );

                // compute contribution
                mdFluxdDof( tDofIndex ).matrix_data() +=
                        mVelocityMatrix * this->dStressdDOF( aDofTypes ) -
                        this->dEnergydDOF( aDofTypes ) * tFIVelocity->val() -
                        mThermalFluxDof;
            }

            // direct dependency on the velocity dof type
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mdFluxdDof( tDofIndex ).matrix_data() +=
                        mVelocityMatrix * this->dStressdDOF( aDofTypes ) +
                        tStressTensor * tFIVelocity->dnNdxn( 1 ) -
                        this->dEnergydDOF( aDofTypes ) * tFIVelocity->val() -
                        this->Energy() * tFIVelocity->dnNdxn( 1 );
            }

            // direct dependency on the temperature dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // evaluate thermal flux DoF derivative
                this->eval_thermal_dFluxdDOF( aDofTypes );

                // compute contribution
                mdFluxdDof( tDofIndex ).matrix_data() +=
                        mVelocityMatrix * this->dStressdDOF( aDofTypes ) -
                        this->dEnergydDOF( aDofTypes ) * tFIVelocity->val() -
                        mThermalFluxDof;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_thermal_flux()
        {
            // get the temperature FI
            Field_Interpolator * tTempFI =  mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the thermal conductivity property
            std::shared_ptr< Property > tThermalConductivity = get_property( "ThermalConductivity" );

            // compute thermal flux q = - k * grad(T)
            mThermalFlux = -1.0 * tThermalConductivity->val()( 0 ) * tTempFI->gradx( 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_thermal_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // initialize the matrix
            mThermalFluxDof.set_size( mSpaceDim,
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // get temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // compute derivative with direct dependency
                mThermalFluxDof.matrix_data() +=
                        -1.0 * tPropThermalConductivity->val()( 0 ) * tFITemp->dnNdxn( 1 );
            }

            // if indirect dependency of conductivity on the dof type
            if ( tPropThermalConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mThermalFluxDof.matrix_data() +=
                        -1.0 * tFITemp->gradx( 1 ) * tPropThermalConductivity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_Energy()
        {
            // get field interpolator values
            Matrix< DDRMat >  tDensity = mFIManager->get_field_interpolators_for_type( mDofDensity )->val();
            Matrix< DDRMat >  tVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();
            Matrix< DDRMat >  tTemperature = mFIManager->get_field_interpolators_for_type( mDofTemperature )->val();

            // get the heat capacity
            Matrix< DDRMat > tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val();

            // compute thermal flux q = - k * grad(T)
            mEnergy = tIsochoricHeatCapacity * tDensity * tTemperature + 0.5 * trans( tVelocity ) * tVelocity * tDensity;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dEnergydDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the FIs
            Field_Interpolator * tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the heat capacity
            Matrix< DDRMat > tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val();

            // initialize the matrix
            mEnergyDof( tDofIndex ).set_size( 1,
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                // compute contribution
                mEnergyDof( tDofIndex ).matrix_data() +=
                        tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->N() +
                        0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->N();
            }

// FIXME: check derivative of norm of velocity
            // direct dependency on the velocity dof type
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mEnergyDof( tDofIndex ).matrix_data() +=
                        tFIDensity->val() * trans( tFIVelocity->val() ) * tFIVelocity->N();
            }

            // direct dependency on the temperature dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // compute contribution
                mEnergyDof( tDofIndex ).matrix_data() +=
                        tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->N();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::CM_Fluid_Compressible_Ideal::eval_EnergyDot()
        {
            // get the FIs
            Field_Interpolator * tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the heat capacity
            Matrix< DDRMat > tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val();

            // compute total energy density
            mEnergyDot =
                    tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->gradt( 1 ) +
                    tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->gradt( 1 ) +
                    0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->gradt( 1 ) +
                    tFIDensity->val() * trans( tFIVelocity->val() ) * tFIVelocity->gradt( 1 ) ;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the FIs
            Field_Interpolator * tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the heat capacity
            Matrix< DDRMat > tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val();

            // initialize the matrix
            mEnergyDotDof( tDofIndex ).set_size( 1,
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                // compute contribution
                mEnergyDotDof( tDofIndex ).matrix_data() +=
                        tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->dnNdtn( 1 ) +
                        tIsochoricHeatCapacity * tFITemp->gradt( 1 ) * tFIDensity->N() +
                        0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->dnNdtn( 1 ) +
                        trans( tFIVelocity->val() ) * tFIVelocity->gradt( 1 ) * tFIDensity->N() ;
            }

            // direct dependency on the velocity dof type
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mEnergyDotDof( tDofIndex ).matrix_data() +=
                        tFIDensity->gradt( 1 ) * trans( tFIVelocity->val() ) * tFIVelocity->N()   +
                        tFIDensity->val() * trans( tFIVelocity->gradt( 1 ) ) * tFIVelocity->N()  +
                        tFIDensity->val() * trans( tFIVelocity->val() ) * tFIVelocity->dnNdtn( 1 ) ;
            }

            // direct dependency on the temperature dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // compute contribution
                mEnergyDotDof( tDofIndex ).matrix_data() +=
                        tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->dnNdtn( 1 ) +
                        tIsochoricHeatCapacity * tFIDensity->gradt( 1 ) * tFITemp->N() ;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_stress()
        {
            // get velocity FI
            Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // compute Stress
            mStress = 2.0 * tPropDynamicViscosity->val() *
                    ( this->strain() -  ( 1 / 3 ) * tFIVelocity->div() * mFlatIdentity  ) -
                    this->pressure() * mFlatIdentity ;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dStressdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the FIs
            Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // initialize the matrix
            mdStressdDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3,
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                // compute contribution
                mdStressdDof( tDofIndex ).matrix_data() +=
                        mFlatIdentity * this->dPressuredDOF( aDofTypes );
            }

            // direct dependency on the velocity dof type
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mdStressdDof( tDofIndex ).matrix_data() +=
                        2.0 * tPropDynamicViscosity->val() *
                        ( this->dStraindDOF( aDofTypes ) - ( 1 / 3 ) * mFlatIdentity * tFIVelocity->div_operator() );
            }

            // direct dependency on the temperature dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // compute contribution
                mdStressdDof( tDofIndex ).matrix_data() +=
                        mFlatIdentity * this->dPressuredDOF( aDofTypes );
            }

            // if indirect dependency of viscosity
            if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdStressdDof( tDofIndex ).matrix_data() +=
                        2.0 * ( this->strain() - ( 1 / 3 ) * tFIVelocity->div() * mFlatIdentity ) *
                        tPropDynamicViscosity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // Function is not implemented
            MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::eval_dTractiondDOF - not implemented yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // Function is not implemented
            MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::eval_traction - not implemented yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::pressure()
        {
            // get field interpolator values
            Matrix< DDRMat >  tDensity = mFIManager->get_field_interpolators_for_type( mDofDensity )->val();
            Matrix< DDRMat >  tTemperature = mFIManager->get_field_interpolators_for_type( mDofTemperature )->val();

            // get the specific gas constant
            Matrix< DDRMat >  tSpecificGasConstant = get_property( "SpecificGasConstant" )->val();

            // return the pressure
            mPressure = tDensity * tSpecificGasConstant * tTemperature;
            return mPressure;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::dPressuredDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get field interpolators
            Field_Interpolator * tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the specific gas constant
            Matrix< DDRMat >  tSpecificGasConstant = get_property( "SpecificGasConstant" )->val();

            // if Density DoF
            if( aDofTypes( 0 ) == mDofDensity )
            {
                mPressureDof = tSpecificGasConstant * tFITemp->val() * tFIDensity->dnNdxn( 1 );
            }

            // if Temperature DoF
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // compute derivative with direct dependency
                mPressureDof = tSpecificGasConstant * tFIDensity->val() * tFITemp->dnNdxn( 1 );
            }

            // return the pressure DoF deriv
            return mPressureDof;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_strain_2d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 3, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = 0.5 * ( tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 ) );
        }

        void CM_Fluid_Compressible_Ideal::eval_strain_3d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 6, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
            mStrain( 3, 0 ) = 0.5 * ( tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 ) );
            mStrain( 4, 0 ) = 0.5 * ( tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 ) );
            mStrain( 5, 0 ) = 0.5 * ( tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_teststrain_2d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn =  mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        void CM_Fluid_Compressible_Ideal::eval_teststrain_3d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 6, tNumBases * 3, 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

            mTestStrain( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init mdStraindDof
            mdStraindDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute derivative
                mdStraindDof( tDofIndex ).matrix_data() += this->testStrain().matrix_data();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Compressible_Ideal::eval_velocitymatrix_2d()
        {
            Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            mVelocityMatrix = {
                    { tU( 0 ),    0.0 , tU( 1 ) },
                    {    0.0 , tU( 1 ), tU( 0 ) } };
        }

        void CM_Fluid_Compressible_Ideal::eval_velocitymatrix_3d()
        {
            Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            mVelocityMatrix = {
                    { tU( 0 ),    0.0 ,  0.0 ,    0.0 , tU( 2 ), tU( 1 ) },
                    {    0.0 , tU( 1 ),  0.0 , tU( 2 ),    0.0 , tU( 0 ) },
                    {    0.0 ,  0.0 , tU( 2 ), tU( 1 ), tU( 0 ),    0.0  } };
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Compressible_Ideal::unfold_2d(
                const Matrix< DDRMat > & aFlattenedTensor,
                Matrix< DDRMat > & aExpandedTensor)
        {
            aExpandedTensor = {
                    { aFlattenedTensor( 0 ), aFlattenedTensor( 2 ) },
                    { aFlattenedTensor( 2 ), aFlattenedTensor( 1 ) } };
        }

        void CM_Fluid_Compressible_Ideal::unfold_3d(
                const Matrix< DDRMat > & aFlattenedTensor,
                Matrix< DDRMat > & aExpandedTensor)
        {
            aExpandedTensor = {
                    { aFlattenedTensor( 0 ), aFlattenedTensor( 5 ), aFlattenedTensor( 4 ) },
                    { aFlattenedTensor( 5 ), aFlattenedTensor( 1 ), aFlattenedTensor( 3 ) },
                    { aFlattenedTensor( 4 ), aFlattenedTensor( 3 ), aFlattenedTensor( 2 ) } };
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */



































