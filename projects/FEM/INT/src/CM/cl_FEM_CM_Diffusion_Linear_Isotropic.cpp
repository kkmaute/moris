/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic.cpp
 *
 */

#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        CM_Diffusion_Linear_Isotropic::CM_Diffusion_Linear_Isotropic()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ] = static_cast< uint >( CM_Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ] = static_cast< uint >( CM_Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "EigenStrain" ]  = static_cast< uint >( CM_Property_Type::EIGEN_STRAIN );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof types
            for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // if temperature dof type string
                if ( tDofString == "Temperature" )
                {
                    mTempDof = tDofType;
                }
                else if ( tDofString == "Theta" )
                {
                    mThetaDof = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Diffusion_Linear_Isotropic::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::set_local_properties()
        {
            // set the conductivity property
            mPropConductivity = get_property( "Conductivity" );

            // set the heat capacity property
            mPropHeatCapacity = get_property( "HeatCapacity" );

            // set the density property
            mPropDensity = get_property( "Density" );

            // set the eigenstrain property
            mPropEigenStrain = get_property( "EigenStrain" );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_flux()
        {
            // compute flux
            // Note: this is a numerical flux and not a physical flux which is the negative
            //       of the numerical flux
            mFlux = mPropConductivity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );

            if ( mPropEigenStrain != nullptr )
            {
                // Get field interpolator for theta
                Field_Interpolator* tFITheta = mFIManager->get_field_interpolators_for_type( mThetaDof );

                // compute normalized gradient of theta
                const Matrix< DDRMat >& tGradTheta = tFITheta->gradx( 1 );

                // compute norm of spatial gradient of theta
                const real tNorm = norm( tGradTheta );

                // add eigen strain contribution
                if ( tNorm > MORIS_REAL_EPS )
                {
                    mFlux += mPropEigenStrain->val()( 0 ) * mPropConductivity->val()( 0 ) * tGradTheta / tNorm;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_Energy()
        {
            // if density and heat capacity properties are set
            if ( mPropDensity != nullptr && mPropHeatCapacity != nullptr )
            {
                // compute enthalpy
                mEnergy = mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->val();
            }
            else
            {
                // if no capacity or density is given, set Energy to zero
                mEnergy = 0.0;
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_EnergyDot()
        {
            // if density and heat capacity properties are set
            if ( mPropDensity != nullptr && mPropHeatCapacity != nullptr )
            {
                // compute rate of enthalpy
                mEnergyDot = mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradt( 1 );
            }
            else
            {
                // if no capacity or density is given, set EnergyDot to zero
                mEnergyDot = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradt( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_gradEnergy()
        {
            // if density and heat capacity properties are set
            if ( mPropDensity != nullptr && mPropHeatCapacity != nullptr )
            {
                // compute rate of gradient of enthalpy
                mGradEnergy = mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
            }
            else
            {
                // if no capacity or density is given, set gradEnergy to zero
                mGradEnergy = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_gradEnergyDot()
        {
            // if density and heat capacity properties are set
            if ( mPropDensity != nullptr && mPropHeatCapacity != nullptr )
            {
                // compute rate of gradient of enthalpy
                mGradEnergyDot = mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradxt();
            }
            else
            {
                // if no capacity or density is given, set gradEnergyDot to zero
                mGradEnergyDot = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradxt();
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_divflux()
        {
            // compute the divergence of the flux
            mDivFlux = mPropConductivity->val() * this->divstrain();
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_graddivflux()
        {
            // get conductivity property value
            moris::real tK = mPropConductivity->val()( 0 );

            // get spatial interpolation order
            mtk::Interpolation_Order tInterpOrder =
                    mFIManager->get_field_interpolators_for_type( mTempDof )->get_space_interpolation_order();

            // only compute if interpolation is not 1 or 0, in that case simply keep the vector of zeros
            if ( tInterpOrder == mtk::Interpolation_Order::LINEAR || tInterpOrder == mtk::Interpolation_Order::CONSTANT )
            {
                // initialize results with zeros
                mGradDivFlux.set_size( mSpaceDim, 1, 0.0 );
            }
            else
            {
                // matrix for purely isotropic case
                // FIXME: this implementation is slow and needs to be improved
                Matrix< DDRMat > tKijIsotropic;
                switch ( mSpaceDim )
                {
                    case 2:
                    {
                        tKijIsotropic = {
                            { tK, 0, 0, tK },
                            { 0, tK, tK, 0 }
                        };
                        break;
                    }
                    case 3:
                    {
                        tKijIsotropic = {
                            { tK, 0, 0, 0, 0, tK, 0, tK, 0, 0 },
                            { 0, tK, 0, tK, 0, 0, 0, 0, tK, 0 },
                            { 0, 0, tK, 0, tK, 0, tK, 0, 0, 0 }
                        };
                        break;
                    }
                    default:
                        MORIS_ASSERT( false, "CM_Diffusion_Linear_Isotropic::eval_graddivflux: Number of spatial dimensions must be 2 or 3" );
                }

                // compute grad div flux
                mGradDivFlux = tKijIsotropic * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 3 );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_traction( const Matrix< DDRMat >& aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_testTraction(
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // compute test traction
            mTestTraction( tTestDofIndex ) =
                    trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) * mPropConductivity->val()( 0 ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_strain()
        {
            // compute strain
            mStrain = mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_divstrain()
        {
            // get the temperature gradient
            Matrix< DDRMat > tTempGrad =
                    mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 2 );

            // evaluate the divergence of the strain
            mDivStrain = sum( tTempGrad( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_testStrain()
        {
            // compute test strain
            mTestStrain = mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = mPropConductivity->val()( 0 ) * I;
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dFluxdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mdFluxdDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get temperature FI
            Field_Interpolator* tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the mTempDof type
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) +=
                        mPropConductivity->val()( 0 ) * tFITemp->dnNdxn( 1 );
            }

            // if indirect dependency on the dof type
            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ) +=
                        tFITemp->gradx( 1 ) * mPropConductivity->dPropdDOF( aDofTypes );
            }

            // if direct dependency on the mThetaDof type
            if ( aDofTypes( 0 ) == mThetaDof )
            {
                // Get field interpolator for theta
                Field_Interpolator* tFITheta = mFIManager->get_field_interpolators_for_type( mThetaDof );

                // compute normalized gradient of Theta
                const Matrix< DDRMat >& tBTheta    = tFITheta->dnNdxn( 1 );
                const Matrix< DDRMat >& tGradTheta = tFITheta->gradx( 1 );

                // compute norm of spatial gradient of theta
                real tNorm = norm( tGradTheta );

                // add eigen strain contribution
                if ( tNorm > MORIS_REAL_EPS )
                {
                    Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;
                    Matrix< DDRMat > tNormBTheta    = tBTheta / tNorm;

                    // compute derivative with direct dependency
                    mdFluxdDof( tDofIndex ) +=
                            mPropEigenStrain->val()( 0 ) * mPropConductivity->val()( 0 ) * ( tNormBTheta - tNormGradTheta * trans( tNormGradTheta ) * tNormBTheta );
                }
            }
        }
        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dEnergydDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mEnergyDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( mPropDensity == nullptr || mPropHeatCapacity == nullptr )
            {
                return;
            }

            // get the temperature FI
            Field_Interpolator* tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if the dof type is temperature
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * tFITemp->N();
            }

            // if density property depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDof( tDofIndex ) +=
                        mPropHeatCapacity->val()( 0 ) * tFITemp->val() * mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * tFITemp->val() * mPropHeatCapacity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dEnergyDotdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mEnergyDotDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( mPropDensity == nullptr || mPropHeatCapacity == nullptr )
            {
                return;
            }

            // get the temperature FI
            Field_Interpolator* tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if the dof type is temperature
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * tFITemp->dnNdtn( 1 );
            }

            // if density property depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        mPropHeatCapacity->val()( 0 ) * tFITemp->gradt( 1 ) * mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * tFITemp->gradt( 1 ) * mPropHeatCapacity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dGradEnergydDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( mPropDensity == nullptr || mPropHeatCapacity == nullptr )
            {
                return;
            }

            // temperature dof type
            Field_Interpolator* tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if the dof type is temperature
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mGradEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * tFITemp->dnNdxn( 1 );
            }

            // if density depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        mPropHeatCapacity->val()( 0 ) * tFITemp->gradx( 1 ) * mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * tFITemp->gradx( 1 ) * mPropHeatCapacity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dGradEnergyDotdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDotDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( mPropDensity == nullptr || mPropHeatCapacity == nullptr )
            {
                return;
            }

            // temperature dof type
            Field_Interpolator* tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the dof type
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mGradEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * mPropHeatCapacity->val()( 0 ) * tFITemp->d2Ndxt();
            }

            // if density depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        mPropHeatCapacity->val()( 0 ) * tFITemp->gradxt() * mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * tFITemp->gradxt() * mPropHeatCapacity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dGradDivFluxdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get conductivity property value
            moris::real tK = mPropConductivity->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize results with zeros
            mGradDivFluxDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get spatial interpolation order
            mtk::Interpolation_Order tInterpOrder = tFIDer->get_space_interpolation_order();

            // only compute if interpolation is not 1 or 0, in that case simply keep the vector of zeros
            if ( tInterpOrder != mtk::Interpolation_Order::LINEAR && tInterpOrder != mtk::Interpolation_Order::CONSTANT )
            {
                // matrix for purely isotropic case
                // FIXME: this implementation is slow and doesn't support indirect DoF dependency of conductivity
                Matrix< DDRMat > tKijIsotropic;
                switch ( mSpaceDim )
                {
                    case 2:
                    {
                        tKijIsotropic = {
                            { tK, 0, 0, tK },
                            { 0, tK, tK, 0 }
                        };
                        break;
                    }
                    case 3:
                    {
                        tKijIsotropic = {
                            { tK, 0, 0, 0, 0, tK, 0, tK, 0, 0 },
                            { 0, tK, 0, tK, 0, 0, 0, 0, tK, 0 },
                            { 0, 0, tK, 0, tK, 0, tK, 0, 0, 0 }
                        };
                        break;
                    }
                    default:
                        MORIS_ASSERT( false, "CM_Diffusion_Linear_Isotropic::eval_dGradDivFluxdDOF: Number of spatial dimensions must be 2 or 3" );
                }

                // FIXME: indirect dependencies missing, spatial derivatives of properties needed
                // if direct dependency on the dof type
                if ( aDofTypes( 0 ) == mTempDof )
                {
                    // compute derivative with direct dependency
                    mGradDivFluxDof( tDofIndex ) +=
                            tKijIsotropic * mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 3 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_ddivfluxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator* tFIDer = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if temperature dof
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ) += mPropConductivity->val()( 0 ) * this->ddivstraindu( aDofTypes );
            }

            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mddivfluxdu( tDofIndex ) += this->divstrain() * mPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_ddivstraindu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if temperature dof type
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tTempd2Ndx2 = tFIDer->dnNdxn( 2 );

                // fill ddivstrain/du
                mddivstraindu( tDofIndex ) = tTempd2Ndx2.get_row( 0 ) + tTempd2Ndx2.get_row( 1 );

                if ( tTempd2Ndx2.n_rows() == 6 )
                {
                    mddivstraindu( tDofIndex ) += tTempd2Ndx2.get_row( 2 );
                }
            }
        }

        //------------------------------------------------------------------------------
        void
        CM_Diffusion_Linear_Isotropic::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute derivative
            mdTractiondDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(), mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // if conductivity depends on dof type
            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) * aNormal * mPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const Matrix< DDRMat >&             aJump,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the test dof type FI
            Field_Interpolator* tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( tFITest->get_number_of_space_time_coefficients(), tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if conductivity depends on dof type
            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) * aNormal * mPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dStraindDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the matrix for dStraindDof
            mdStraindDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if direct dependency on the dof type
            if ( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = tFIDer->dnNdxn( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dConstdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if conductivity depends on the dof type
            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = mPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dFluxdDV(
                const moris::Cell< PDV_Type >& aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dFluxdDV - This function is not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        CM_Diffusion_Linear_Isotropic::eval_dStraindDV(
                const moris::Cell< PDV_Type >& aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dStraindDV - This function is not implemented." );
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

