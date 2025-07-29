/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence_Smagorinsky.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Incompressible_Turbulence_Smagorinsky.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Fluid_Incompressible_Turbulence_Smagorinsky::CM_Fluid_Incompressible_Turbulence_Smagorinsky()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "Viscosity" ]    = static_cast< uint >( CM_Property_Type::VISCOSITY );
        mPropertyMap[ "WallDistance" ] = static_cast< uint >( CM_Property_Type::WALL_DISTANCE );
    }

    //-------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::set_local_properties()
    {
        // set the density property
        mPropDensity = get_property( "Density" );

        // set the dynamic viscosity property
        mPropViscosity = get_property( "Viscosity" );
        
        // set the wall distance property
        mPropWallDistance = get_property( "WallDistance" );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::reset_specific_eval_flags()
    {
        // call the parent implementation
        CM_Fluid_Incompressible_Turbulence::reset_specific_eval_flags();

        // reset child specific eval flags
        mFluidStrainRateEval = true;
        mdFluidStrainRateduEval.fill( true );
        mdFluidStrainRatedxEval.fill( true );
        mdFluidStrainRatedxduEval.fill( true );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::initialize_spec_storage_vars_and_eval_flags()
    {
        // call the parent implementation
        CM_Fluid_Incompressible_Turbulence::initialize_spec_storage_vars_and_eval_flags();

        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();

        // init child specific eval flags
        mdFluidStrainRateduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdFluidStrainRatedxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdFluidStrainRatedxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );

        // init child specific storage
        mdFluidStrainRatedu.resize( tNumGlobalDofTypes );
        mdFluidStrainRatedx.resize( mMaxSpaceDerOrder );
        mdFluidStrainRatedxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdFluidStrainRatedxdu( iOrder ).resize( tNumGlobalDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_turbulent_dynamic_viscosity()
    {
        // compute the strain rate form
        real tStrainRateForm = 2.0 * dot( this->strain(), this->strain() );

        // compute turbulent viscosity
        mTurbDynVisc = mPropDensity->val()( 0 )                                   //
                     * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )    //
                     * std::pow( tStrainRateForm, 0.5 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // add contribution to derivative of turbulent dynamic viscosity
        mdTurbDynViscdu( tDofIndex ) =                                       //
                mPropDensity->val()( 0 )                                     //
                * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )    //
                * this->dfluidstrainratedu( aDofTypes ) * 0.5 / std::pow( this->fluid_strain_rate()( 0 ), 0.5 );

        // if wall distance depends on dof
        if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
        {
            // add contribution from dwallDistancedu
            mdTurbDynViscdu( tDofIndex ) +=                              //
                    mPropDensity->val()( 0 )                             //
                    * std::pow( this->fluid_strain_rate()( 0 ), 0.5 )    //
                    * 2.0 * std::pow( mKappa, 2.0 ) * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes );
        }

        // if density depends on dof
        if ( mPropDensity->check_dof_dependency( aDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdu -"    //
                    "Dependence of density on dof not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdx( uint aOrder )
    {
        // set matrix size
        mdTurbDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // add contribution from space derivative of strain rate
        mdTurbDynViscdx( aOrder - 1 ) =                                      //
                mPropDensity->val()( 0 )                                     //
                * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )    //
                * this->dfluidstrainratedx( aOrder ) * 0.5 / std::pow( this->fluid_strain_rate()( 0 ), 0.5 );

        // if wall distance has a space dependency
        if ( mPropWallDistance->check_space_dependency() )
        {
            // add contribution from density space derivative
            mdTurbDynViscdx( aOrder - 1 ) +=                             //
                    mPropDensity->val()( 0 )                             //
                    * std::pow( this->fluid_strain_rate()( 0 ), 0.5 )    //
                    * 2.0 * std::pow( mKappa, 2.0 ) * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdx( aOrder );
        }

        // if density has a space dependency
        if ( mPropDensity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                    "Dependence of density on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // add contribution from strain rate form
        mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) =                                                                             //
                mPropDensity->val()( 0 )                                                                                           //
                * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )                                                          //
                * ( ( 0.5 * this->dfluidstrainratedxdu( aDofTypes, aOrder ) / std::pow( this->fluid_strain_rate()( 0 ), 0.5 ) )    //
                        - ( 0.25 * this->dfluidstrainratedx( aOrder ) * this->dfluidstrainratedu( aDofTypes ) / std::pow( this->fluid_strain_rate()( 0 ), 1.5 ) ) );

        // if wall distance depends on dof
        if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
        {
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=                                                             //
                    mPropDensity->val()( 0 )                                                                            //
                    * ( 0.5 * this->dfluidstrainratedx( aOrder ) / std::pow( this->fluid_strain_rate()( 0 ), 0.5 ) )    //
                    * std::pow( mKappa, 2.0 ) * 2.0 * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes );
        }

        // if density depends on dof
        if ( mPropDensity->check_dof_dependency( aDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,                                                                  //
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                    "Dependence of density on dof not accounted for." );
        }

        // if wall distance has a space dependency
        if ( mPropWallDistance->check_space_dependency() )
        {
            // add contribution for strain rate form
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=                                                                   //
                    mPropDensity->val()( 0 )                                                                                  //
                    * std::pow( mKappa, 2.0 ) * 2.0 * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdx( aOrder )    //
                    * 0.5 * this->dfluidstrainratedu( aDofTypes ) / std::pow( this->fluid_strain_rate()( 0 ), 0.5 );

            // if wall distance depends on dof
            if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
            {
                // add contribution from wall distance
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=                                                 //
                        mPropDensity->val()( 0 )                                                                //
                        * std::pow( this->fluid_strain_rate()( 0 ), 0.5 )                                       //
                        * std::pow( mKappa, 2.0 ) * 2.0                                                         //
                        * ( mPropWallDistance->dPropdx( aOrder ) * mPropWallDistance->dPropdDOF( aDofTypes )    //
                                + mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdxdDOF( aDofTypes, aOrder ) );
            }

            // if density depends on dof
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // assume that density prop does not depend on dof
                MORIS_ERROR( false,                                                                  //
                        "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                        "Dependence of density on dof not accounted for." );
            }
        }

        // if density has a space dependency
        if ( mPropDensity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,                                                                  //
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                    "Dependence of density on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_testturbdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // compute test turbulent dyn. viscosity
        mTestTurbDynVisc( tTestDofIndex ) = this->dturbdynviscdu( aTestDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dtestturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        //        // get the test dof FI
        //        Field_Interpolator* tFITest =
        //                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );
        //
        //        // get the derivative dof FI
        //        Field_Interpolator* tFIDer =
        //                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
        //
        //        // this term is not zero for Smagorinsky
        //        mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ).set_size(    //
        //                tFITest->get_number_of_space_time_coefficients(),      //
        //                tFIDer->get_number_of_space_time_coefficients() );


        // add contribution to derivative of strain rate form
        mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ) =                                                                                                 //
                mPropDensity->val()( 0 )                                                                                                                    //
                * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )                                                                                   //
                * ( 2.0 * trans( this->dStraindDOF( aTestDofTypes ) ) * this->dStraindDOF( aDofTypes ) / std::pow( this->fluid_strain_rate()( 0 ), 0.5 )    //
                        - 0.25 * trans( this->dfluidstrainratedu( aTestDofTypes ) ) * this->dfluidstrainratedu( aDofTypes ) / std::pow( this->fluid_strain_rate()( 0 ), 1.5 ) );

        // if wall distance depends on dof
        if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
        {
            mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ) +=                                                                    //
                    mPropDensity->val()( 0 )                                                                                        //
                    * trans( this->dfluidstrainratedu( aTestDofTypes ) ) * 0.5 / std::pow( this->fluid_strain_rate()( 0 ), 0.5 )    //
                    * 2.0 * std::pow( mKappa, 2.0 ) * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes );
        }

        // if density depends on dof
        if ( mPropDensity->check_dof_dependency( aDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible::eval_dtestturbdynviscdu -"    //
                    "Dependence of density on dof not accounted for." );
        }

        // if wall distance depends on test dof
        if ( mPropWallDistance->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,                                                                    //
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dtestturbdynviscdu -"    //
                    "Dependence of wall distance on test dof not accounted for." );
        }

        // if density depends on test dof
        if ( mPropDensity->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,                                                                    //
                    "CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dtestturbdynviscdu -"    //
                    "Dependence of density on test dof not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::fluid_strain_rate(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                       //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::fluid_strain_rate - "    //
                "Only DEFAULT CM function type known in base class." );

        // if not evaluated
        if ( mFluidStrainRateEval )
        {
            // evaluate
            this->eval_fluid_strain_rate();

            // set bool for evaluation
            mFluidStrainRateEval = false;
        }

        // return value
        return mFluidStrainRate;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_fluid_strain_rate()
    {
        mFluidStrainRate = std::max( 2.0 * dot( this->strain(), this->strain() ), mFluidStrainRateTol );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                        //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                              //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedu - "    //
                "No dependency on this dof type." );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        if ( mdFluidStrainRateduEval( tDofIndex ) )
        {
            // compute derivative
            this->eval_dfluidstrainratedu( aDofTypes );

            // set bool for evaluation
            mdFluidStrainRateduEval( tDofIndex ) = false;
        }
        // return value
        return mdFluidStrainRatedu( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dfluidstrainratedu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init storage
        mdFluidStrainRatedu( tDofIndex ).set_size(    //
                1,                                    //
                tFI->get_number_of_space_time_coefficients() );

        // check for strain rate value
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol && aDofTypes( 0 ) == mDofVelocity )
        {
            // fill the derivative
            mdFluidStrainRatedu( tDofIndex ) =    //
                    4.0 * trans( this->strain() ) * this->dStraindDOF( aDofTypes );
        }
        else
        {
            mdFluidStrainRatedu( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                        //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedx - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for space order
        MORIS_ERROR( aOrder == 1,                                                          //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedx - "    //
                "Supported only for 1st order space derivative." );

        // if the derivative has not been evaluated yet
        if ( mdFluidStrainRatedxEval( aOrder - 1 ) )
        {
            // evaluate the derivative
            this->eval_dfluidstrainratedx( aOrder );

            // set bool for evaluation
            mdFluidStrainRatedxEval( aOrder - 1 ) = false;
        }

        // return the derivative
        return mdFluidStrainRatedx( aOrder - 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dfluidstrainratedx(
            uint aOrder )
    {
        // init storage
        mdFluidStrainRatedx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // check for strain rate value
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol )
        {
            // fill the derivative
            mdFluidStrainRatedx( aOrder - 1 ) = 4.0 * trans( this->dstraindx( aOrder ) ) * this->strain();
        }
        else
        {
            mdFluidStrainRatedx( aOrder - 1 ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                          //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                            //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
                "Works only for 1st order space derivative." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                                //
                "CM_Fluid_Incompressible_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFluidStrainRatedxduEval( aOrder - 1, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dfluidstrainratedxdu( aDofTypes, aOrder );

            // set bool for evaluation
            mdFluidStrainRatedxduEval( aOrder - 1, tDofIndex ) = false;
        }

        // return the derivative
        return mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::eval_dfluidstrainratedxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init storage
        mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex ).set_size(    //
                mSpaceDim,                                            //
                tFI->get_number_of_space_time_coefficients() );

        // check for strain rate value
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol && aDofTypes( 0 ) == mDofVelocity )
        {
            // compute the derivative
            mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex ) =                                     //
                    4.0 * ( trans( this->dstraindx( aOrder ) ) * this->dStraindDOF( aDofTypes )    //
                            + this->dstraindxdu( aDofTypes, aOrder, this->strain() ) );
        }
        else
        {
            mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::select_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMRequestType )
        {
            case CM_Request_Type::TURB_DYN_VISC:
            {
                return this->turbulent_dynamic_viscosity();
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                return this->dturbdynviscdx( 1 );
                break;
            }
            case CM_Request_Type::TEST_TURB_DYN_VISC:
            {
                //  get the test dof index
                uint tTestDofIndex                = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );
                mTestTurbDynVisc( tTestDofIndex ) = trans( this->testturbdynvisc( aTestDofTypes ) );
                return mTestTurbDynVisc( tTestDofIndex );
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_RATE:
            {
                return this->fluid_strain_rate();
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_RATE_SPACE_DER:
            {
                return this->dfluidstrainratedx( 1 );
                break;
            }
            default:
                MORIS_ERROR( false,                                                                  //
                        "CM_Fluid_Incompressible_Turbulence_Smagorinsky::select_derivative_FD - "    //
                        "aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Smagorinsky::set_derivative_FD(
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
            case CM_Request_Type::TURB_DYN_VISC:
            {
                // set value to storage
                mdTurbDynViscdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                // set value to storage
                mdTurbDynViscdxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TEST_TURB_DYN_VISC:
            {
                // get the test dof index
                uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

                // set value to storage
                mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_RATE:
            {
                // set value to storage
                mdFluidStrainRatedu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_RATE_SPACE_DER:
            {
                // set value to storage
                mdFluidStrainRatedxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            default:
                MORIS_ERROR( false,                                                               //
                        "CM_Fluid_Incompressible_Turbulence_Smagorinsky::set_derivative_FD - "    //
                        "aCMRequestType undefined" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
