/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky.cpp
 *
 */

#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Conductivity" ] = static_cast< uint >( CM_Property_Type::CONDUCTIVITY );
        mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "HeatCapacity" ] = static_cast< uint >( CM_Property_Type::HEAT_CAPACITY );
        mPropertyMap[ "EigenStrain" ]  = static_cast< uint >( CM_Property_Type::EIGEN_STRAIN );

        mPropertyMap[ "TurbulentPrandtl" ]   = static_cast< uint >( CM_Property_Type::TURBULENT_PRANDTL );
        mPropertyMap[ "WallDistance" ] = static_cast< uint >( CM_Property_Type::WALL_DISTANCE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_function_pointers()
    {
        switch ( mSpaceDim )
        {
            case 2:
            {
                m_eval_fluid_strain   = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_fluid_strain_2d;
                m_eval_dfluidstraindu = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindu_2d;
                m_eval_dfluidstraindx   = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindx_2d;
                m_eval_dfluidstraindxdu = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindxdu_2d;
                break;
            }
            case 3:
            {
                m_eval_fluid_strain   = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_fluid_strain_3d;
                m_eval_dfluidstraindu = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindu_3d;
                m_eval_dfluidstraindx   = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindx_3d;
                m_eval_dfluidstraindxdu = &CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindxdu_3d;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_function_pointers - only works for 2d and 3d." );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_space_dim( uint aSpaceDim )
    {
        // check that space dimension is 1, 2, 3
        MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4,                                        //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_space_dim - "    //
                "Wrong space dimension." );

        // set space dimension
        mSpaceDim = aSpaceDim;

        // set function pointers
        this->set_function_pointers();
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::reset_specific_eval_flags()
    {
        // call the parent implementation
        CM_Diffusion_Linear_Isotropic_Turbulence::reset_specific_eval_flags();

        // reset child specific eval flags
        mFluidStrainEval = true;
        mdFluidStrainduEval.fill( true );
        mdFluidStraindxEval.fill( true );
        mdFluidStraindxduEval.fill( true );

        mFluidStrainRateEval = true;
        mdFluidStrainRateduEval.fill( true );
        mdFluidStrainRatedxEval.fill( true );
        mdFluidStrainRatedxduEval.fill( true );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::initialize_spec_storage_vars_and_eval_flags()
    {
        // call the parent implementation
        CM_Diffusion_Linear_Isotropic_Turbulence::initialize_spec_storage_vars_and_eval_flags();

        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();

        // init child specific eval flags
        mdFluidStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdFluidStraindxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdFluidStraindxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );

        mdFluidStrainRateduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdFluidStrainRatedxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdFluidStrainRatedxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );

        // init child specific storage
        mdFluidStraindu.resize( tNumGlobalDofTypes );
        mdFluidStraindx.resize( mMaxSpaceDerOrder );
        mdFluidStraindxdu.resize( mMaxSpaceDerOrder );
        mdFluidStrainRatedu.resize( tNumGlobalDofTypes );
        mdFluidStrainRatedx.resize( mMaxSpaceDerOrder );
        mdFluidStrainRatedxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdFluidStraindxdu( iOrder ).resize( tNumGlobalDofTypes );
            mdFluidStrainRatedxdu( iOrder ).resize( tNumGlobalDofTypes );
        }

    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
            const Vector< std::string >&             aDofStrings )
    {
        // set dof type list
        Constitutive_Model::set_dof_type_list( aDofTypes );

        // loop over the provided dof types
        for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
        {
            // get dof type string
            const std::string& tDofString = aDofStrings( iDof );

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
            else if ( tDofString == "Velocity" )
            {
                mDofVelocity = tDofType;
            }
            else
            {
                // error unknown dof string
                MORIS_ERROR( false,                                                                     //
                        "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_dof_type_list - "    //
                        "Unknown aDofString : %s \n",
                        tDofString.c_str() );
            }
        }
    }

    //-------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_local_properties()
    {
        // call parent implementation
        CM_Diffusion_Linear_Isotropic_Turbulence::set_local_properties();

        // set the wall distance property
        mPropWallDistance = get_property( "WallDistance" );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_turbulent_dynamic_viscosity()
    {

        // compute turbulent viscosity
        mTurbDynVisc = mPropDensity->val()( 0 )                                   //
                     * std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 )    //
                     * std::pow( this->fluid_strain_rate()( 0 ), 0.5 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // add contribution to derivative of turbulent dynamic viscosity
        mdTurbDynViscdu( tDofIndex ) =
                mPropDensity->val()( 0 ) *                                   //
                std::pow( mKappa * mPropWallDistance->val()( 0 ), 2.0 ) *    //
                this->dfluidstrainratedu( aDofTypes ) / 2.0 / std::pow( this->fluid_strain_rate()( 0 ), 0.5 );

        // if wall distance depends on dof
        if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
        {
            // add contribution from dwallDistancedu
            mdTurbDynViscdu( tDofIndex ) +=                              //
                    mPropDensity->val()( 0 )                             //
                    * std::pow( this->fluid_strain_rate()( 0 ), 0.5 )    //
                    * std::pow( mKappa, 2.0 ) * 2.0 * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes );
        }

        // if density depends on dof
        if ( mPropDensity->check_dof_dependency( aDofTypes ) )
        {
            // assume that density prop does not depend on dof
            MORIS_ERROR( false,                                                                      //
                    "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdu -"    //
                    "Dependence of density on dof not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdx(
            uint aOrder )
    {
        // add contribution from strain rate form space dependency
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
                    * std::pow( mKappa, 2.0 ) * 2.0 * mPropWallDistance->val()( 0 ) * mPropWallDistance->dPropdx( aOrder );
        }

        // if density has a space dependency
        if ( mPropDensity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,                                                                        //
                    "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                    "Dependence of density on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdxdu(
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
            MORIS_ERROR( false,                                                                        //
                    "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
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
                MORIS_ERROR( false,                                                                        //
                        "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                        "Dependence of density on dof not accounted for." );
            }
        }

        // if density has a space dependency
        if ( mPropDensity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,                                                                        //
                    "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dturbdynviscdxdu -"    //
                    "Dependence of density on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_testturbdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // compute test turbulent dyn. viscosity
        mTestTurbDynVisc( tTestDofIndex ) = this->dturbdynviscdu( aTestDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dtestturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
         uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

         // get the test dof FI
         Field_Interpolator* tFITest =
                 mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

         // get the derivative dof FI
         Field_Interpolator* tFIDer =
                 mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

         // if wall distance depends on test dof
         if ( mPropWallDistance->check_dof_dependency( aTestDofTypes ) )
         {
             // assume that density prop does not depend on dof
             MORIS_ERROR( false,                                                                          //
                     "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dtestturbdynviscdu -"    //
                     "Dependence of wall distance on test dof not accounted for." );
         }

         // if density depends on test dof
         if ( mPropDensity->check_dof_dependency( aTestDofTypes ) )
         {
             // assume that density prop does not depend on dof
             MORIS_ERROR( false,                                                                          //
                     "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dtestturbdynviscdu -"    //
                     "Dependence of density on test dof not accounted for." );
         }

         // this term is zero for Smagorinsky
         mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ).set_size(    //
                 tFITest->get_number_of_space_time_coefficients(),      //
                 tFIDer->get_number_of_space_time_coefficients(),
                 0.0 );
    }

    //------------------------------------------------------------------------------
    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::fluid_strain(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                     //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::fluid_strain"    //
                "- Only DEFAULT CM function type known in base class." );

        // if the fluid strain was not evaluated
        if ( mFluidStrainEval )
        {
            // evaluate the fluid strain
            this->eval_fluid_strain();

            // set bool for evaluation
            mFluidStrainEval = false;
        }

        // return the fluid strain value
        return mFluidStrain;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_fluid_strain_2d()
    {
        // get the velocity spatial gradient from velocity FI
        const Matrix< DDRMat >& tVelocityGradx =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mFluidStrain.set_size( 3, 1 );
        mFluidStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mFluidStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mFluidStrain( 2, 0 ) = 0.5 * ( tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_fluid_strain_3d()
    {
        // get the velocity spatial gradient from velocity FI
        const Matrix< DDRMat >& tVelocityGradx =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mFluidStrain.set_size( 6, 1 );
        mFluidStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mFluidStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mFluidStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
        mFluidStrain( 3, 0 ) = 0.5 * ( tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 ) );
        mFluidStrain( 4, 0 ) = 0.5 * ( tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 ) );
        mFluidStrain( 5, 0 ) = 0.5 * ( tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 ) );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                          //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dFluidStraindu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                                //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindu - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFluidStrainduEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dfluidstraindu( aDofTypes );

            // set bool for evaluation
            mdFluidStrainduEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdFluidStraindu( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init storage
        mdFluidStraindu( tDofIndex ).set_size(    //
                ( mSpaceDim - 1 ) * 3,            //
                tFI->get_number_of_space_time_coefficients() );

        // if velocity dof
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute derivative
            ( this->*m_eval_dfluidstraindu )();
        }
        else
        {
            mdFluidStraindu( tDofIndex ).fill( 0.0 );
        }

        // set bool for evaluation
        mdFluidStrainduEval( tDofIndex ) = false;
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindu_2d()
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( mDofVelocity ) );

        // get velocity field interpolator
        Field_Interpolator* tFIVelocity =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // compute velocity gradient
        const Matrix< DDRMat >& tdnNdxn = tFIVelocity->dnNdxn( 1 );

        // get number of bases for velocity
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // build the test strain
        mdFluidStraindu( tDofIndex ).set_size( 3, tNumBases * 2, 0.0 );

        mdFluidStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mdFluidStraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindu_3d()
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( mDofVelocity ) );

        // get velocity field interpolator
        Field_Interpolator* tFIVelocity =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // compute displacement gradient
        const Matrix< DDRMat >& tdnNdxn = tFIVelocity->dnNdxn( 1 );

        // get number of bases for displacement
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // build the test strain
        mdFluidStraindu( tDofIndex ).set_size( 6, tNumBases * 3, 0.0 );

        mdFluidStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 4, 4 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 5, 5 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mdFluidStraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

        mdFluidStraindu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mdFluidStraindu( tDofIndex )( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                          //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindx - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for space order
        MORIS_ERROR( aOrder == 1,                                                       //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dstraindx - "    //
                "Supported only for 1st order space derivative." );

        // if the derivative has not been evaluated yet
        if ( mdFluidStraindxEval( aOrder - 1 ) )
        {
            // evaluate the derivative
            this->eval_dfluidstraindx( aOrder );

            // set bool for evaluation
            mdFluidStraindxEval( aOrder - 1 ) = false;
        }

        // return the derivative
        return mdFluidStraindx( aOrder - 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindx_2d(
            uint aOrder )
    {
        // set size for dfluidstraindx
        mdFluidStraindx( aOrder - 1 ).set_size( 3, 2 );

        // get the velocity gradient
        const Matrix< DDRMat >& tVelocityGrad =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

        // fill dfluidstraindx
        mdFluidStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
        mdFluidStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 2, 0 );
        mdFluidStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 2, 1 );
        mdFluidStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
        mdFluidStraindx( aOrder - 1 )( 2, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
        mdFluidStraindx( aOrder - 1 )( 2, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindx_3d(
            uint aOrder )
    {
        // set size for dfluidstraindx
        mdFluidStraindx( aOrder - 1 ).set_size( 6, 3 );

        // get the velocity gradient
        const Matrix< DDRMat >& tVelocityGrad =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

        // fill dfluidstraindx
        mdFluidStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
        mdFluidStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 5, 0 );
        mdFluidStraindx( aOrder - 1 )( 0, 2 ) = tVelocityGrad( 4, 0 );
        mdFluidStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 5, 1 );
        mdFluidStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
        mdFluidStraindx( aOrder - 1 )( 1, 2 ) = tVelocityGrad( 3, 1 );
        mdFluidStraindx( aOrder - 1 )( 2, 0 ) = tVelocityGrad( 4, 2 );
        mdFluidStraindx( aOrder - 1 )( 2, 1 ) = tVelocityGrad( 3, 2 );
        mdFluidStraindx( aOrder - 1 )( 2, 2 ) = tVelocityGrad( 2, 2 );
        mdFluidStraindx( aOrder - 1 )( 3, 0 ) = 0.5 * ( tVelocityGrad( 4, 1 ) + tVelocityGrad( 5, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 3, 1 ) = 0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 3, 2 ) = 0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 4, 0 ) = 0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 4, 1 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 5, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 4, 2 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
        mdFluidStraindx( aOrder - 1 )( 5, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
        mdFluidStraindx( aOrder - 1 )( 5, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) );
        mdFluidStraindx( aOrder - 1 )( 5, 2 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 4, 1 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            const Matrix< DDRMat >&        aJump,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                            //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                              //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindxdu - "    //
                "Works only for 1st order derivative for now." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                                  //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstraindxdu - "    //
                "No dependency on this dof type." );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFluidStraindxduEval( aOrder - 1, tDofIndex ) )
        {
            // get the dof FI
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init storage
            mdFluidStraindxdu( aOrder - 1 )( tDofIndex ).set_size(    //
                    mSpaceDim,                                        //
                    tFI->get_number_of_space_time_coefficients() );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute derivative
                ( this->*m_eval_dfluidstraindxdu )( aOrder, aJump );
            }
            else
            {
                mdFluidStraindxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
            }

            // set bool for evaluation
            mdFluidStraindxduEval( aOrder - 1, tDofIndex ) = false;
        }

        // return value
        return mdFluidStraindxdu( aOrder - 1 )( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindxdu_2d(
            uint                    aOrder,
            const Matrix< DDRMat >& aJump )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( mDofVelocity ) );

        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // get number of bases for velocity
        uint tNumBases = tFI->get_number_of_space_time_bases();

        // get the velocity gradient
        const Matrix< DDRMat >& tdnNdxn = tFI->dnNdxn( 2 );

        // fill dstraindxdu
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =    //
                aJump( 0 ) * tdnNdxn.get_row( 0 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 2 );
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =    //
                aJump( 1 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 0 );

        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =    //
                aJump( 0 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 1 );
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =    //
                aJump( 1 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstraindxdu_3d(
            uint                           aOrder,
            const Matrix< DDRMat >&        aJump )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( mDofVelocity ) );

        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // get number of bases for velocity
        uint tNumBases = tFI->get_number_of_space_time_bases();

        // get the velocity gradient
        const Matrix< DDRMat >& tdnNdxn = tFI->dnNdxn( 2 );

        // fill dstraindxdu
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =                                                       //
                aJump( 0 ) * tdnNdxn.get_row( 0 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 2 );         //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =                                      //
                aJump( 1 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 0 );    //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                  //
                aJump( 2 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 0 );    //

        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =                                                       //
                aJump( 0 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 1 );         //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =                                      //
                aJump( 1 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 5 );    //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                  //
                aJump( 2 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 5 );    //

        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) =                                                       //
                aJump( 0 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 3 );         //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =                                      //
                aJump( 1 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 4 );    //
        mdFluidStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                  //
                aJump( 2 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 4 );    //
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::fluid_strain_rate(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                             //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::fluid_strain_rate - "    //
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
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_fluid_strain_rate()
    {
        mFluidStrainRate = std::max( 2.0 * dot( this->fluid_strain(), this->fluid_strain() ), mFluidStrainRateTol );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                              //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                                    //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedu - "    //
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
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstrainratedu(
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
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol )
        {
            // fill the derivative
            mdFluidStrainRatedu( tDofIndex ) =    //
                    4.0 * trans( this->fluid_strain() ) * this->dfluidstraindu( aDofTypes );
        }
        else
        {
            mdFluidStrainRatedu( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                              //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedx - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for space order
        MORIS_ERROR( aOrder == 1,                                                                //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedx - "    //
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
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstrainratedx(
            uint aOrder )
    {
        // init storage
        mdFluidStrainRatedx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // check for strain rate value
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol )
        {
            // fill the derivative
            mdFluidStrainRatedx( aOrder - 1 ) = 4.0 * trans( this->dfluidstraindx( aOrder ) ) * this->fluid_strain();
        }
        else
        {
            mdFluidStrainRatedx( aOrder - 1 ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                                //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                                  //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
                "Works only for 1st order space derivative." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                                      //
                "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::dfluidstrainratedxdu - "    //
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
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::eval_dfluidstrainratedxdu(
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
        if ( this->fluid_strain_rate()( 0 ) > mFluidStrainRateTol )
        {
            // compute the derivative
            mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex ) =                                             //
                    4.0 * ( trans( this->dfluidstraindx( aOrder ) ) * this->dfluidstraindu( aDofTypes )    //
                            + this->dfluidstraindxdu( aDofTypes, aOrder, this->fluid_strain() ) );
        }
        else
        {
            mdFluidStrainRatedxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::select_derivative_FD(
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
            case CM_Request_Type::FLUID_STRAIN:
            {
                return this->fluid_strain();
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_SPACE_DER:
            {
                mFluidStrain = trans( this->dfluidstraindx( 1 ) ) * aJump;
                return mFluidStrain;
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
                MORIS_ERROR( false,                                                                        //
                        "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::select_derivative_FD - "    //
                        "aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_derivative_FD(
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
            case CM_Request_Type::FLUID_STRAIN:
            {
                // set value to storage
                mdFluidStraindu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::FLUID_STRAIN_SPACE_DER:
            {
                // set value to storage
                mdFluidStraindxdu( 0 )( tDofIndex ) = aDerivativeFD;
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
                MORIS_ERROR( false,                                                                     //
                        "CM_Diffusion_Linear_Isotropic_Turbulence_Smagorinsky::set_derivative_FD - "    //
                        "aCMRequestType undefined" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
