/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

// LINALG/src
#include "fn_trans.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "Viscosity" ]    = static_cast< uint >( CM_Property_Type::VISCOSITY );
        mPropertyMap[ "KinViscosity" ] = static_cast< uint >( CM_Property_Type::KIN_VISCOSITY );

        // init storage for evaluation
        mdChidx.resize( mMaxSpaceDerOrder );
        mdFv1dx.resize( mMaxSpaceDerOrder );
        mdTurbDynViscdx.resize( mMaxSpaceDerOrder );
        mdEffDynViscdx.resize( mMaxSpaceDerOrder );

        // init flag for evaluation
        mdChidxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdFv1dxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdTurbDynViscdxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdEffDynViscdxEval.set_size( mMaxSpaceDerOrder, 1, true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::reset_specific_eval_flags()
    {
        // call parent implementation
        CM_Fluid_Incompressible_Turbulence::reset_specific_eval_flags();

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
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::initialize_spec_storage_vars_and_eval_flags()
    {
        // call parent implementation
        CM_Fluid_Incompressible_Turbulence::initialize_spec_storage_vars_and_eval_flags();

        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();

        // init child specific eval flags
        mdChiduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdFv1duEval.set_size( tNumGlobalDofTypes, 1, true );

        mdChidxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );
        mdFv1dxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );

        // init child specific storage
        mdChidu.resize( tNumGlobalDofTypes );
        mdFv1du.resize( tNumGlobalDofTypes );

        mdChidxdu.resize( mMaxSpaceDerOrder );
        mdFv1dxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdChidxdu( iOrder ).resize( tNumGlobalDofTypes );
            mdFv1dxdu( iOrder ).resize( tNumGlobalDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
            const Vector< std::string >&             aDofStrings )
    {
        // set dof type list
        Constitutive_Model::set_dof_type_list( aDofTypes );

        // loop over the provided dof type
        for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
        {
            // get dof type string
            const std::string& tDofString = aDofStrings( iDof );

            // get dof type
            MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

            // switch on dof type string
            if ( tDofString == "Velocity" )
            {
                mDofVelocity = tDofType;
            }
            else if ( tDofString == "Pressure" )
            {
                mDofPressure = tDofType;
            }
            else if ( tDofString == "Viscosity" )
            {
                mDofViscosity = tDofType;
            }
            else
            {
                // error unknown dof string
                MORIS_ERROR( false,
                        "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::set_dof_type_list - Unknown aDofString : %s \n",
                        tDofString.c_str() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::set_local_properties()
    {
        // set the density property
        mPropDensity = get_property( "Density" );

        // set the dynamic viscosity property
        mPropViscosity = get_property( "Viscosity" );

        // set the kinematic viscosity property
        mPropKinViscosity = get_property( "KinViscosity" );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_turbulent_dynamic_viscosity()
    {
        // init mTurbDynVisc
        mTurbDynVisc = 0.0;

        // get the viscosity dof type FI
        Field_Interpolator* tFIModViscosity =
                mFIManager->get_field_interpolators_for_type( mDofViscosity );

        // get the modified viscosity value
        real tModViscosity = tFIModViscosity->val()( 0 );

        // if modified viscosity is positive
        if ( tModViscosity >= 0.0 )
        {
            // compute turbulent viscosity
            mTurbDynVisc = mPropDensity->val()( 0 ) * tModViscosity * this->fv1();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get derivative dof type FI
        Field_Interpolator* tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set matrix size
        mdTurbDynViscdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

        // get the viscosity dof type FI
        Field_Interpolator* tFIModViscosity =
                mFIManager->get_field_interpolators_for_type( mDofViscosity );

        // get the modified viscosity value
        real tModViscosity = tFIModViscosity->val()( 0 );

        // if modified viscosity is positive
        if ( tModViscosity >= 0.0 )
        {
            // add contribution from dfv1du
            mdTurbDynViscdu( tDofIndex ) =
                    mPropDensity->val()( 0 ) * tFIModViscosity->val() * this->dfv1du( aDofTypes );

            // if dof type is viscosity
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to dSPdu
                mdTurbDynViscdu( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * this->fv1() * tFIModViscosity->N();
            }

            // if density depends on dof
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
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

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dturbdynviscdx(
            uint aOrder )
    {
        // set matrix size
        mdTurbDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // get the viscosity dof type FI
        Field_Interpolator* tFIModViscosity =
                mFIManager->get_field_interpolators_for_type( mDofViscosity );

        // get the modified viscosity value
        real tModViscosity = tFIModViscosity->val()( 0 );

        // if modified viscosity is positive
        if ( tModViscosity >= 0.0 )
        {
            // compute dTurbDynViscdx
            mdTurbDynViscdx( aOrder - 1 ) =                              //
                    mPropDensity->val()( 0 ) * (                         //
                            tFIModViscosity->gradx( 1 ) * this->fv1()    //
                            + this->dfv1dx( 1 ) * tModViscosity );

            // if density depends on space
            if ( mPropDensity->check_space_dependency() )
            {
                // assume that density prop does not depend on x
                MORIS_ERROR( false,                                                          //
                        "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdx -"    //
                        "Dependence of density on space not accounted for." );
            }
        }
        else
        {
            mdTurbDynViscdx( aOrder - 1 ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dturbdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get derivative dof type FI
        Field_Interpolator* tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set matrix size
        mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).set_size(    //
                mSpaceDim,                                        //
                tFIDer->get_number_of_space_time_coefficients() );

        // get the viscosity dof type FI
        Field_Interpolator* tFIModViscosity =
                mFIManager->get_field_interpolators_for_type( mDofViscosity );

        // get the modified viscosity value
        real tModViscosity = tFIModViscosity->val()( 0 );

        // if modified viscosity is positive
        if ( tModViscosity >= 0.0 )
        {
            // add contribution from dfv1du
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) =                             //
                    mPropDensity->val()( 0 ) * (                                       //
                            tFIModViscosity->gradx( 1 ) * this->dfv1du( aDofTypes )    //
                            + tFIModViscosity->val()( 0 ) * this->dfv1dxdu( aDofTypes, 1 ) );

            // if dof type is viscosity
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to dviscositytdxdu
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=               //
                        mPropDensity->val()( 0 ) * (                          //
                                this->fv1() * tFIModViscosity->dnNdxn( 1 )    //
                                + this->dfv1dx( 1 ) * tFIModViscosity->N() );
            }

            // if density depends on dof
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // assume that density prop does not depend on dof
                MORIS_ERROR( false,                                                                        //
                        "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dturbdynviscdxdu - "    //
                        "Dependence of density on dof not accounted for." );
            }

            // if density depends on space
            if ( mPropDensity->check_space_dependency() )
            {
                // assume that density prop does not depend on x
                MORIS_ERROR( false,                                                                        //
                        "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dturbdynviscdxdu - "    //
                        "Dependence of density on space not accounted for." );
            }
        }
        else
        {
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_testturbdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        MORIS_ERROR( false,                                                                       //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_testturbdynvisc - "    //
                "Not implemented as not needed." );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dtestturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // check for acceptable test dof
        MORIS_ERROR( aTestDofTypes( 0 ) == mDofVelocity || aTestDofTypes( 0 ) == mDofPressure,    //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::eval_dtestturbdynviscdu"    //
                " - Only admissible test dofs are velocity or pressure" );

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

        // this term is always zero for Spalart-Allmaras
        mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ).set_size(    //
                tFITest->get_number_of_space_time_coefficients(),      //
                tFIDer->get_number_of_space_time_coefficients(),
                0.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    real CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::chi(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,              //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::chi - "    //
                "Only DEFAULT CM function type known in base class." );

        // if the diffusion coefficient was not evaluated
        if ( mChiEval )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidu(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                 //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                        //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidu - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdChiduEval( tDofIndex ) )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                 //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidx - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                   //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidx - "    //
                "Works only for 1st order derivative for now." );

        // if the derivative has not been evaluated yet
        if ( mdChidxEval( aOrder - 1 ) )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidxdu(
            const Vector< MSI::Dof_Type >& aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                   //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                     //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidxdu - "    //
                "Works only for 1st order space derivative." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                          //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dchidxdu - "    //
                "No dependency on this dof type." );


        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdChidxduEval( aOrder - 1, tDofIndex ) )
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

    real CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::fv1(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,              //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::fv1 - "    //
                "Only DEFAULT CM function type known in base class." );

        // if not evaluated
        if ( mFv1Eval )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1du(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                 //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1du - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                        //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1du - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFv1duEval( tDofIndex ) )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                 //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dx - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                   //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dx - "    //
                "Works only for 1st order derivative for now." );

        // if the derivative has not been evaluated yet
        if ( mdFv1dxEval( aOrder - 1 ) )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dxdu(
            const Vector< MSI::Dof_Type >& aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                   //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                                     //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dxdu - "    //
                "Works only for 1st order derivative for now." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                          //
                "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::dfv1dxdu - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdFv1dxduEval( aOrder - 1, tDofIndex ) )
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

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::select_derivative_FD(
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
            case CM_Request_Type::STRAIN_SPACE_DER:
            {
                mStrain = trans( this->dstraindx( 1, aCMFunctionType ) ) * aJump;
                return mStrain;
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
            case CM_Request_Type::DIV_FLUX:
            {
                return this->divflux();
                break;
            }
            case CM_Request_Type::DIV_STRAIN:
            {
                return this->divstrain();
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC:
            {
                return this->turbulent_dynamic_viscosity();
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
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                return this->dturbdynviscdx( 1 );
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC:
            {
                mEffDynVisc = this->effective_dynamic_viscosity();
                return mEffDynVisc;
                break;
            }
            case CM_Request_Type::TEST_EFF_DYN_VISC:
            {
                // get the test dof index
                uint tTestDofIndex               = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );
                mTestEffDynVisc( tTestDofIndex ) = trans( this->testeffdynvisc( aTestDofTypes ) );
                return mTestEffDynVisc( tTestDofIndex );
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC_SPACE_DER:
            {
                return this->deffdynviscdx( 1 );
                break;
            }
            default:
                MORIS_ERROR( false,                                                      //
                        "CM_Fluid_Incompressible_Turbulence::select_derivative_FD - "    //
                        "aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::set_derivative_FD(
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
            case CM_Request_Type::STRAIN_SPACE_DER:
            {
                // set value to storage
                mdStraindxdu( 0 )( tDofIndex ) = aDerivativeFD;
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
            case CM_Request_Type::DIV_FLUX:
            {
                // set value to storage
                mddivfluxdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::DIV_STRAIN:
            {
                // set value to storage
                mddivstraindu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC:
            {
                // set value to storage
                mdTurbDynViscdu( tDofIndex ) = aDerivativeFD;
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
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                // set value to storage
                mdTurbDynViscdxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC:
            {
                // set value to storage
                mdEffDynViscdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TEST_EFF_DYN_VISC:
            {
                // get the test dof index
                uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

                // set value to storage
                mdTestEffDynViscdu( tTestDofIndex )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC_SPACE_DER:
            {
                // set value to storage
                mdEffDynViscdxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            default:
                MORIS_ERROR( false,                                                                    //
                        "CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras::set_derivative_FD - "    //
                        "aCMRequestType undefined" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
