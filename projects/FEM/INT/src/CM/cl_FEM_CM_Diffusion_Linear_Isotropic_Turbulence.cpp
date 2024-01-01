/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.cpp
 *
 */

#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

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

        CM_Diffusion_Linear_Isotropic_Turbulence::CM_Diffusion_Linear_Isotropic_Turbulence()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ] = static_cast< uint >( CM_Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ] = static_cast< uint >( CM_Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "EigenStrain" ]  = static_cast< uint >( CM_Property_Type::EIGEN_STRAIN );

            mPropertyMap[ "KinematicViscosity" ] = static_cast< uint >( CM_Property_Type::KIN_VISCOSITY );
            mPropertyMap[ "TurbulentPrandtl" ]   = static_cast< uint >( CM_Property_Type::TURBULENT_PRANDTL );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;

            // init storage for evaluation
            mdChidx.resize( tOrder );
            mdFv1dx.resize( tOrder );
            mdTurbDynViscdx.resize( tOrder );
            mdEffConddx.resize( tOrder );

            // init flag for evaluation
            mdChidxEval.set_size( tOrder, 1, true );
            mdFv1dxEval.set_size( tOrder, 1, true );
            mdTurbDynViscdxEval.set_size( tOrder, 1, true );
            mdEffConddxEval.set_size( tOrder, 1, true );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::reset_eval_flags()
        {
            // call parent implementation
            Constitutive_Model::reset_eval_flags();

            // reset child specific eval flags for chi
            mChiEval = true;
            mdChiduEval.fill( true );
            mdChidxEval.fill( true );
            mdChidxduEval.fill( true );

            // reset child specific eval flags for fv1
            mFv1Eval = true;
            mdFv1duEval.fill( true );
            mdFv1dxEval.fill( true );
            mdFv1dxduEval.fill( true );

            // reset child specific eval flags for turbulence dynamic viscosity
            mTurbDynViscEval = true;
            mdTurbDynViscduEval.fill( true );
            mdTurbDynViscdxEval.fill( true );
            mdTurbDynViscdxduEval.fill( true );

            // reset child specific eval flags for effective conductivity
            mEffCondEval = true;
            mdEffCondduEval.fill( true );
            mdEffConddxEval.fill( true );
            mdEffConddxduEval.fill( true );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::build_global_dof_type_list()
        {
            // call parent implementation
            Constitutive_Model::build_global_dof_type_list();

            // get number of dof types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();

            // init child specific eval flags
            mdChiduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFv1duEval.set_size( tNumGlobalDofTypes, 1, true );
            mdTurbDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEffCondduEval.set_size( tNumGlobalDofTypes, 1, true );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;
            mdChidxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdFv1dxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdEffConddxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdTurbDynViscdxduEval.set_size( tOrder, tNumGlobalDofTypes, true );

            // init child specific storage
            mdChidu.resize( tNumGlobalDofTypes );
            mdFv1du.resize( tNumGlobalDofTypes );
            mdTurbDynViscdu.resize( tNumGlobalDofTypes );
            mdEffConddu.resize( tNumGlobalDofTypes );

            mdChidxdu.resize( tOrder );
            mdFv1dxdu.resize( tOrder );
            mdTurbDynViscdxdu.resize( tOrder );
            mdEffConddxdu.resize( tOrder );
            for( uint iOrder = 0; iOrder < tOrder; iOrder++ )
            {
                mdChidxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdFv1dxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdEffConddxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdTurbDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::set_dof_type_list(
                moris::Vector< moris::Vector< MSI::Dof_Type > > aDofTypes,
                moris::Vector< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof types
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // if temperature dof type string
                if( tDofString == "Temperature" )
                {
                    mTempDof = tDofType;
                }
                else if ( tDofString == "Viscosity" )
                {
                    mDofViscosity  = tDofType;
                }
                else if ( tDofString == "Theta" )
                {
                    mThetaDof  = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Diffusion_Linear_Isotropic_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::set_local_properties()
        {
            // set the conductivity property
            mPropConductivity = get_property( "Conductivity" );

            // set the heat capacity property
            mPropHeatCapacity = get_property( "HeatCapacity" );

            // set the density property
            mPropDensity = get_property( "Density" );

            // set the eigenstrain property
            mPropEigenStrain = get_property( "EigenStrain" );

            // set the kinematic viscosity
            mPropKinViscosity = get_property( "KinematicViscosity" );

            // set the Prandtl turbulence property
            mPropPrandtlT = get_property( "TurbulentPrandtl" );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_flux()
        {
            // compute flux
            // Note: this is a numerical flux and not a physical flux which is the negative
            //       of the numerical flux
            mFlux = this->effective_conductivity()( 0 ) *
                    mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );

            if ( mPropEigenStrain != nullptr )
            {
                MORIS_LOG_INFO( "All theta dof dependencies unchanged" );

                // Get field interpolator for theta
                Field_Interpolator * tFITheta = mFIManager->get_field_interpolators_for_type( mThetaDof );

                // compute normalized gradient of theta
                const Matrix< DDRMat > & tGradTheta = tFITheta->gradx( 1 );

                // compute norm of spatial gradient of theta
                const real tNorm = norm( tGradTheta );

                // add eigen strain contribution
                if ( tNorm > MORIS_REAL_EPS )
                {
                    mFlux += mPropConductivity->val()( 0 ) * tGradTheta / tNorm;
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dFluxdDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mdFluxdDof( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative with indirect dependency through effective conductivity
            mdFluxdDof( tDofIndex ) = tFITemp->gradx( 1 ) * this->deffconddu( aDofTypes );

            // if direct dependency on the mTempDof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) +=
                        this->effective_conductivity()( 0 ) * tFITemp->dnNdxn( 1 );
            }

            // if direct dependency on the mThetaDof type
            if( aDofTypes( 0 ) == mThetaDof )
            {
                MORIS_LOG_INFO( "All theta dof dependencies unchanged" );

                // get field interpolator for theta
                Field_Interpolator * tFITheta =
                        mFIManager->get_field_interpolators_for_type( mThetaDof );

                // compute normalized gradient of Theta
                const Matrix< DDRMat > & tBTheta    = tFITheta->dnNdxn( 1 );
                const Matrix< DDRMat > & tGradTheta = tFITheta->gradx( 1 );

                // compute norm of spatial gradient of theta
                real tNorm = norm( tGradTheta );

                // add eigen strain contribution
                if ( tNorm > MORIS_REAL_EPS )
                {
                    Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;
                    Matrix< DDRMat > tNormBTheta    = tBTheta / tNorm;

                    // compute derivative with direct dependency
                    mdFluxdDof( tDofIndex ) +=
                            mPropConductivity->val()( 0 ) *  ( tNormBTheta - tNormGradTheta * trans( tNormGradTheta ) * tNormBTheta );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_divflux()
        {
            // compute the divergence of the flux
            mDivFlux = this->effective_conductivity() * this->divstrain() +
                    trans( this->deffconddx( 1 ) ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_ddivfluxdu(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex =
                    mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // add contributions of the derivatives of the effective conductivity and grad effective conductivity
            mddivfluxdu( tDofIndex ) =
                    this->effective_conductivity() * this->ddivstraindu( aDofTypes ) +
                    trans( mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 ) ) * this->deffconddxdu( aDofTypes, 1 ) +
                    this->divstrain() * this->deffconddu( aDofTypes );

            // if temperature dof
            if( aDofTypes( 0 ) == mTempDof )
            {
                // add contributions of the derivatives of the strain and div strain
                mddivfluxdu( tDofIndex ) +=
                        trans( this->deffconddx( 1 ) ) * mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTractiondDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute derivative
            mdTractiondDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // add contribution from derivative of effective conductivity wrt test dof
            mTestTraction( tTestDofIndex ) =
                    trans( aNormal * this->deffconddu( aTestDofTypes ) ) * this->strain() +
                    this->effective_conductivity()( 0 ) * trans( this->dStraindDOF( aTestDofTypes ) ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTestTractiondDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const moris::Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the test dof FI
            Field_Interpolator * tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients() );

            // if viscosity is the test dof
            if( aTestDofTypes( 0 ) == mDofViscosity && aDofTypes( 0 ) == mDofViscosity )
            {
                // FIXME: Missing second order derivative of effective dynamic viscosity - FD for now

                Constitutive_Model::eval_dtesttractiondu_FD(
                        aDofTypes,
                        aTestDofTypes,
                        mdTestTractiondDof( tTestDofIndex )( tDofIndex ),
                        1e-6,
                        aNormal,
                        fem::FDScheme_Type::POINT_3_CENTRAL );
            }
            else
            {
                // if effective dynamic viscosity depends on test or derivative dof type
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                        trans( this->dStraindDOF( aTestDofTypes ) ) * aNormal *
                        this->deffconddu( aDofTypes ) +
                        trans( this->deffconddu( aTestDofTypes ) ) *
                        trans( aNormal ) * this->dStraindDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTestTractiondDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the test dof FI
            Field_Interpolator * tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients() );

            // if viscosity is the test dof
            if( aTestDofTypes( 0 ) == mDofViscosity && aDofTypes( 0 ) == mDofViscosity )
            {
                // FIXME: Missing second order derivative of effective dynamic viscosity - FD for now

                Constitutive_Model::eval_dtesttractiondu_FD(
                        aDofTypes,
                        aTestDofTypes,
                        mdTestTractiondDof( tTestDofIndex )( tDofIndex ),
                        1e-6,
                        aNormal,
                        aJump,
                        fem::FDScheme_Type::POINT_1_FORWARD );
            }
            else
            {
                // if effective dynamic viscosity depends on test or derivative dof type
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                        trans( this->dStraindDOF( aTestDofTypes ) ) * aNormal *
                        aJump * this->deffconddu( aDofTypes ) +
                        trans( this->deffconddu( aTestDofTypes ) ) * aJump *
                        trans( aNormal ) * this->dStraindDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = this->effective_conductivity() * I;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dConstdDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // compute derivative with indirect dependency through properties
            mdConstdDof( tDofIndex ) = this->deffconddu( aDofTypes );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_effective_conductivity()
        {
            // compute the effective conductivity
            mEffCond = mPropConductivity->val()( 0 ) +
                    mPropHeatCapacity->val()( 0 ) * this->turbulent_dynamic_viscosity()( 0 ) / mPropPrandtlT->val()( 0 );
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::effective_conductivity(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::effective_conductivity - Only DEFAULT CM function type known in base class." );

            // if the effective conductivity was not evaluated
            if( mEffCondEval )
            {
                // evaluate the effective conductivity
                this->eval_effective_conductivity();

                // set bool for evaluation
                mEffCondEval = false;
            }
            // return the effective conductivity value
            return mEffCond;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddu( const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdEffConddu( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // add contribution from the derivative of the dynamic viscosity wrt aDofTypes
            mdEffConddu( tDofIndex ) =
                    mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdu( aDofTypes ) / mPropPrandtlT->val()( 0 );

            // if conductivity depends on the dof type
            if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                mdEffConddu( tDofIndex ) += mPropConductivity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                mdEffConddu( tDofIndex ) +=
                        mPropHeatCapacity->dPropdDOF( aDofTypes ) * this->turbulent_dynamic_viscosity()( 0 ) / mPropPrandtlT->val()( 0 );
            }

            // if turbulent Prandtl number depends on the dof type
            if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
            {
                mdEffConddu( tDofIndex ) -=
                        mPropHeatCapacity->val()( 0 ) * this->turbulent_dynamic_viscosity()( 0 ) * mPropPrandtlT->dPropdDOF( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdEffCondduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_deffconddu( aDofType );

                // set bool for evaluation
                mdEffCondduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdEffConddu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddx( uint aOrder )
        {
            // FIXME work only for 1st order
            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx - Works only for 1st order derivative for now." );

            // set matrix size
            mdEffConddx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // add contribution from the derivative of the turbulent dynamic viscosity wrt x
            mdEffConddx( aOrder - 1 ) =
                    mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdx( 1 ) / mPropPrandtlT->val()( 0 );

            // if conductivity depends on space
            if( mPropConductivity->check_space_dependency( 1 ) )
            {
                mdEffConddx( aOrder - 1 ) += mPropConductivity->dnPropdxn( 1 );
            }

            // if heat capacity depends on space
            if( mPropHeatCapacity->check_space_dependency( 1 ) )
            {
                mdEffConddx( aOrder - 1 ) +=
                        mPropHeatCapacity->dnPropdxn( 1 ) * this->turbulent_dynamic_viscosity()( 0 ) / mPropPrandtlT->val()( 0 );
            }

            // if turbulent Prandtl depends on space
            if( mPropPrandtlT->check_space_dependency( 1 ) )
            {
                mdEffConddx( aOrder - 1 ) -=
                        mPropHeatCapacity->val()( 0 ) * this->turbulent_dynamic_viscosity()( 0 ) * mPropPrandtlT->dnPropdxn( 1 ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx(
                uint aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if( mdEffConddxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_deffconddx( aOrder );

                // set bool for evaluation
                mdEffConddxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdEffConddx( aOrder - 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu(
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                uint                                 aOrder)
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdEffConddxdu( aOrder - 1 )( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients() );

            // compute the derivative of the effective conductivity wrt aDofTypes and x
            mdEffConddxdu( aOrder - 1 )( tDofIndex ) =
                    mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdxdu( aDofTypes, 1 ) / mPropPrandtlT->val()( 0 );

            // if conductivity depends on space
            if( mPropConductivity->check_space_dependency( 1 ) )
            {
                // if conductivity depends on the dof type
                if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
                {
                    // FIXME dPropdxdu
                    //mdEffConddxdu( aOrder - 1 )( tDofIndex ) += mPropConductivity->dPropdxdDOF( aDofTypes );
                }
            }

            // if heat capacity depends on space
            if( mPropHeatCapacity->check_space_dependency( 1 ) )
            {
                mdEffConddxdu( aOrder - 1 )( tDofIndex ) +=
                        mPropHeatCapacity->dnPropdxn( 1 ) * this->dturbdynviscdu( aDofTypes ) / mPropPrandtlT->val()( 0 );

                // if turbulent Prandtl number depends on the dof type
                if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
                {
                    mdEffConddxdu( aOrder - 1 )( tDofIndex ) -=
                            mPropHeatCapacity->dnPropdxn( 1 ) * this->turbulent_dynamic_viscosity()( 0 ) * mPropPrandtlT->dPropdDOF( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);
                }

                // FIXME dPropdxdu
                // if heat capacity depends on the dof type
                //if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
                //{
                    //mdEffConddxdu( aOrder - 1 )( tDofIndex ) +=
                            //mPropHeatCapacity->dPropdxdDOF( aDofTypes ) * this->turbulent_dynamic_viscosity()( 0 ) / mPropPrandtlT->val()( 0 )
                //}
            }

            // if turbulent Prandtl number depends on space
            if( mPropPrandtlT->check_space_dependency( 1 ) )
            {
                mdEffConddxdu( aOrder - 1 )( tDofIndex ) -=
                        mPropHeatCapacity->val()( 0 ) * mPropPrandtlT->dnPropdxn( 1 ) * this->dturbdynviscdu( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);

                // if heat capacity depends on the dof type
                if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
                {
                    mdEffConddxdu( aOrder - 1 )( tDofIndex ) -=
                            mPropHeatCapacity->dPropdDOF( aDofTypes ) * this->turbulent_dynamic_viscosity()( 0 ) * mPropPrandtlT->dnPropdxn( 1 ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);
                }

                // if turbulent Prandtl number depends on the dof type
                if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
                {
                    // FIXME dPropdxdu
                    //mdEffConddxdu( aOrder - 1 )( tDofIndex ) -=
                    //        mPropHeatCapacity->val()( 0 ) * tTurbDynVisc * mPropPrandtlT->dPropdxdDOF( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);

                    mdEffConddxdu( aOrder - 1 )( tDofIndex ) +=
                            2.0 * mPropHeatCapacity->val()( 0 ) * this->turbulent_dynamic_viscosity()( 0 ) * mPropPrandtlT->dnPropdxn( 1 ) * mPropPrandtlT->dPropdDOF( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),3.0);
                }
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                mdEffConddxdu( aOrder - 1 )( tDofIndex ) +=
                        mPropHeatCapacity->dPropdDOF( aDofTypes ) * this->dturbdynviscdx( 1 ) / mPropPrandtlT->val()( 0 );
            }

            // if turbulent Prandtl number depends on the dof type
            if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
            {
                mdEffConddxdu( aOrder - 1 )( tDofIndex ) -=
                        mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdx( 1 ) * mPropPrandtlT->dPropdDOF( aDofTypes ) / std::pow(mPropPrandtlT->val()( 0 ),2.0);
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                uint                                 aOrder,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdEffConddxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_deffconddxdu( aDofType, aOrder );

                // set bool for evaluation
                mdEffConddxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdEffConddxdu( aOrder - 1 )( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_turbulent_dynamic_viscosity()
        {
            // init mTurbDynVisc
            mTurbDynVisc = 0.0;

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute turbulent viscosity
                mTurbDynVisc = mPropDensity->val()( 0 ) * tModViscosity * this->fv1();
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::turbulent_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::turbulent_dynamic_viscosity - Only DEFAULT CM function type known in base class." );

            // if the turbulent dynamic viscosity was not evaluated
            if( mTurbDynViscEval )
            {
                // evaluate the turbulent dynamic viscosity
                this->eval_turbulent_dynamic_viscosity();

                // set bool for evaluation
                mTurbDynViscEval = false;
            }
            // return the turbulent dynamic viscosity value
            return mTurbDynVisc;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdu(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdTurbDynViscdu( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // add contribution from dfv1du
                mdTurbDynViscdu( tDofIndex ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->val() * this->dfv1du( aDofTypes );

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dSPdu
                    mdTurbDynViscdu( tDofIndex ) +=
                            mPropDensity->val()( 0 ) * this->fv1() * tFIModViscosity->N();
                }

                // if density depends on dof
                if( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution from drhodu
                    mdTurbDynViscdu( tDofIndex ) +=
                            tFIModViscosity->val() * this->fv1() * mPropDensity->dPropdDOF( aDofTypes );
                }
            }
            else
            {
                mdTurbDynViscdu( tDofIndex ).fill( 0.0 );
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dturbdynviscdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Turbulence::dturbdynviscdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdTurbDynViscduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdu( aDofType );

                // set bool for evaluation
                mdTurbDynViscduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdTurbDynViscdu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdx( uint aOrder )
        {
            // set matrix size
            mdTurbDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute dTurbDynViscdx
                mdTurbDynViscdx( aOrder - 1 ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * this->fv1() +
                        mPropDensity->val()( 0 ) * this->dfv1dx( 1 ) * tModViscosity;

                // if density depends on space
                if( mPropDensity->check_space_dependency( 1 ) )
                {
                    // add contribution from density space derivative
                    mdTurbDynViscdx( aOrder - 1 ) += this->fv1() * tModViscosity * mPropDensity->dnPropdxn( 1 );
                }
            }
            else
            {
                mdTurbDynViscdx( aOrder - 1 ).fill( 0.0 );
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dturbdynviscdx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dturbdynviscdx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if( mdTurbDynViscdxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdx( aOrder );

                // set bool for evaluation
                mdTurbDynViscdxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdTurbDynViscdx( aOrder - 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdxdu(
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                uint                                 aOrder )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // add contribution from dfv1du
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * this->dfv1du( aDofTypes ) +
                        mPropDensity->val()( 0 ) * tFIModViscosity->val()( 0 ) * this->dfv1dxdu( aDofTypes, 1 );

                // if density depends on space
                if( mPropDensity->check_space_dependency( 1 ) )
                {
                    // add contribution from density space derivative
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            mPropDensity->dnPropdxn( 1 ) * tFIModViscosity->val()( 0 ) * this->dfv1du( aDofTypes );
                }

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dviscositytdxdu
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            mPropDensity->val()( 0 ) * this->fv1() * tFIModViscosity->dnNdxn( 1 ) +
                            mPropDensity->val()( 0 ) * this->dfv1dx( 1 ) * tFIModViscosity->N();

                    // if density depends on space
                    if( mPropDensity->check_space_dependency( 1 ) )
                    {
                        // add contribution from density space derivative
                        mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                                this->fv1() * mPropDensity->dnPropdxn( 1 ) * tFIModViscosity->N();
                    }
                }

                // if density depends on dof
                if( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            tFIModViscosity->gradx( 1 ) * this->fv1() * mPropDensity->dPropdDOF( aDofTypes ) +
                            this->dfv1dx( 1 ) * tFIModViscosity->val()( 0 ) * mPropDensity->dPropdDOF( aDofTypes );

                    // if density depends on space
                    if( mPropDensity->check_space_dependency( 1 ) )
                    {
                        // FIXME dPropdxdu
                        //mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                        //tFIModViscosity->val()( 0 ) * tFv1 * mPropDensity->dPropdxdDOF( aDofTypes )
                    }
                }
            }
            else
            {
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
            }
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                uint                                 aOrder,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdTurbDynViscdxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdxdu( aDofType, aOrder );

                // set bool for evaluation
                mdTurbDynViscdxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real CM_Diffusion_Linear_Isotropic_Turbulence::chi(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::chi - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if( mChiEval )
            {
                // evaluate chi
                mChi = compute_chi(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mChiEval = false;
            }
            // return the diffusion coefficient
            return mChi;
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dchidu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdChiduEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidu( tDofIndex ) );

                // set bool for evaluation
                mdChiduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdChidu( tDofIndex );
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if( mdChidxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                compute_dchidx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        mdChidx( aOrder - 1 ) );

                // set bool for evaluation
                mdChidxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdChidx( aOrder - 1 );
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dchidxdu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                uint                                 aOrder,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dchidxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdChidxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidxdu( aOrder - 1 )( tDofIndex ) );

                // set bool for evaluation
                mdChidxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdChidxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real CM_Diffusion_Linear_Isotropic_Turbulence::fv1(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::fv1 - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if( mFv1Eval )
            {
                // evaluate
                mFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mFv1Eval = false;
            }
            // return
            return mFv1;
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dfv1du(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1du - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1du - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdFv1duEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dfv1du(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdFv1du( tDofIndex ) );

                // set bool for evaluation
                mdFv1duEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFv1du( tDofIndex );
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if( mdFv1dxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                compute_dfv1dx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        mdFv1dx( aOrder - 1 ) );

                // set bool for evaluation
                mdFv1dxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdFv1dx( aOrder - 1 );
        }

        const Matrix< DDRMat > & CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dxdu(
                const moris::Vector< MSI::Dof_Type > & aDofType,
                uint                                 aOrder,
                enum CM_Function_Type                aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::dfv1dxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdFv1dxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                compute_dfv1dxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdFv1dxdu( aOrder - 1 )( tDofIndex ) );

                // set bool for evaluation
                mdFv1dxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdFv1dxdu( aOrder - 1 )( tDofIndex );
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

