/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Fluid_Incompressible::CM_Fluid_Incompressible()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "Density" ]   = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "Viscosity" ] = static_cast< uint >( CM_Property_Type::VISCOSITY );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::set_function_pointers()
    {
        switch ( mSpaceDim )
        {
            case 2:
            {
                m_eval_strain       = &CM_Fluid_Incompressible::eval_strain_2d;
                m_eval_divstrain    = &CM_Fluid_Incompressible::eval_divstrain_2d;
                m_eval_teststrain   = &CM_Fluid_Incompressible::eval_teststrain_2d;
                m_eval_dstraindx    = &CM_Fluid_Incompressible::eval_dstraindx_2d;
                m_eval_dstraindxdu  = &CM_Fluid_Incompressible::eval_dstraindxdu_2d;
                m_eval_ddivstraindu = &CM_Fluid_Incompressible::eval_ddivstraindu_2d;
                m_flatten_normal    = &CM_Fluid_Incompressible::flatten_normal_2d;
                break;
            }
            case 3:
            {
                m_eval_strain       = &CM_Fluid_Incompressible::eval_strain_3d;
                m_eval_divstrain    = &CM_Fluid_Incompressible::eval_divstrain_3d;
                m_eval_teststrain   = &CM_Fluid_Incompressible::eval_teststrain_3d;
                m_eval_dstraindx    = &CM_Fluid_Incompressible::eval_dstraindx_3d;
                m_eval_dstraindxdu  = &CM_Fluid_Incompressible::eval_dstraindxdu_3d;
                m_eval_ddivstraindu = &CM_Fluid_Incompressible::eval_ddivstraindu_3d;
                m_flatten_normal    = &CM_Fluid_Incompressible::flatten_normal_3d;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "CM_Fluid_Incompressible::set_function_pointers - only works for 2d and 3d." );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::set_local_properties()
    {
        // set the density property
        mPropDensity = get_property( "Density" );

        // set the viscosity property
        mPropViscosity = get_property( "Viscosity" );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::set_dof_type_list(
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
            else
            {
                // error unknown dof string
                MORIS_ERROR( false,
                        "CM_Fluid_Incompressible::set_dof_type_list - Unknown aDofString : %s \n",
                        tDofString.c_str() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::reset_specific_eval_flags()
    {
        // reset child specific eval flags
        mdStraindxduEval.fill( true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::initialize_spec_storage_vars_and_eval_flags()
    {
        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();

        // init child specific flags
        mdStraindxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );

        // init child specific storage
        mdStraindxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdStraindxdu( iOrder ).resize( tNumGlobalDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_flux()
    {
        // get the pressure FI
        Field_Interpolator* tPressureFI =
                mFIManager->get_field_interpolators_for_type( mDofPressure );

        // create identity matrix
        Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );

        // evaluate pressure contribution to flux
        Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );

        tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * tPressureFI->val();

        // compute flux
        mFlux = -1.0 * tP + 2.0 * mPropViscosity->val()( 0 ) * this->strain();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_divflux()
    {
        // get the pressure FI
        Field_Interpolator* tPressureFI =
                mFIManager->get_field_interpolators_for_type( mDofPressure );

        // compute flux
        mDivFlux = -1.0 * tPressureFI->gradx( 1 ) + 2.0 * mPropViscosity->val()( 0 ) * this->divstrain();

        // check for dependence of dyn. viscosity on space
        if ( mPropViscosity->check_space_dependency() )
        {
            // assume that viscosity prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible::eval_divflux -"    //
                    "Dependence of dynamic viscosity on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_ddivfluxdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the corresponding FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set size for ddivflux/du
        mddivfluxdu( tDofIndex ).set_size(    //
                mSpaceDim,                    //
                tFI->get_number_of_space_time_coefficients() );

        // if velocity dof
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // fill ddivstrain/dv
            mddivfluxdu( tDofIndex ) =
                    2.0 * mPropViscosity->val()( 0 ) * this->ddivstraindu( aDofTypes );
        }
        // if pressure dof
        else if ( aDofTypes( 0 ) == mDofPressure )
        {
            // fill ddivstrain/dp
            mddivfluxdu( tDofIndex ) = -1.0 * tFI->dnNdxn( 1 );
        }
        else
        {
            mddivfluxdu( tDofIndex ).fill( 0.0 );
        }

        if ( mPropViscosity->check_dof_dependency( aDofTypes ) )
        {
            // fill ddivstrain/du
            mddivfluxdu( tDofIndex ) +=
                    2.0 * this->divstrain() * mPropViscosity->dPropdDOF( aDofTypes );
        }

        // check for dependence of dyn. viscosity on space
        if ( mPropViscosity->check_space_dependency() )
        {
            // assume that viscosity prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible::eval_ddivfluxdu - "    //
                    "Dependence of dynamic viscosity on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_traction(
            const Matrix< DDRMat >& aNormal )
    {
        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute the traction
        mTraction = tFlatNormal * this->flux();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // get the dof FI
        Field_Interpolator* tFITest =
                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

        // init mdFluxdDof
        mTestTraction( tTestDofIndex ).set_size(    //
                mSpaceDim,                          //
                tFITest->get_number_of_space_time_coefficients() );

        // if test traction wrt velocity
        if ( aTestDofTypes( 0 ) == mDofVelocity )
        {
            // compute test traction wrt velocity
            mTestTraction( tTestDofIndex ) =
                    2.0 * mPropViscosity->val()( 0 ) * tFlatNormal * this->testStrain();
        }
        // if test traction wrt pressure
        else if ( aTestDofTypes( 0 ) == mDofPressure )
        {
            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
            Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
            tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

            // build the dtesttractiondP
            mTestTraction( tTestDofIndex ) = -1.0 * tFlatNormal * tII * tFITest->N();
        }
        else
        {
            mTestTraction( tTestDofIndex ).fill( 0.0 );
        }

        // if viscosity property depends on test dof type
        if ( mPropViscosity->check_dof_dependency( aTestDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mTestTraction( tTestDofIndex ) +=
                    2.0 * tFlatNormal * this->strain() *    //
                    mPropViscosity->dPropdDOF( aTestDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_strain_2d()
    {
        // get the velocity spatial gradient from velocity FI
        const Matrix< DDRMat >& tVelocityGradx =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mStrain.set_size( 3, 1 );
        mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mStrain( 2, 0 ) = 0.5 * ( tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_strain_3d()
    {
        // get the velocity spatial gradient from velocity FI
        const Matrix< DDRMat >& tVelocityGradx =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mStrain.set_size( 6, 1 );
        mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
        mStrain( 3, 0 ) = 0.5 * ( tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 ) );
        mStrain( 4, 0 ) = 0.5 * ( tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 ) );
        mStrain( 5, 0 ) = 0.5 * ( tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_divstrain_2d()
    {
        // set size for div strain
        mDivStrain.set_size( 2, 1 );

        // get the velocity gradient
        const Matrix< DDRMat >& tVelocityGrad =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

        // fill div strain
        mDivStrain( 0 ) = tVelocityGrad( 0, 0 ) + 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
        mDivStrain( 1 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) ) + tVelocityGrad( 1, 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_divstrain_3d()
    {
        // set size for div strain
        mDivStrain.set_size( 3, 1 );

        // get the velocity gradient
        const Matrix< DDRMat >& tVelocityGrad =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

        // fill div strain
        mDivStrain( 0 ) = tVelocityGrad( 0, 0 )
                        + 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) )
                        + 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
        mDivStrain( 1 ) = 0.5 * ( tVelocityGrad( 5, 0 ) + tVelocityGrad( 0, 1 ) )
                        + tVelocityGrad( 1, 1 )
                        + 0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
        mDivStrain( 2 ) = 0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) )
                        + 0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) )
                        + tVelocityGrad( 2, 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_ddivstraindu_2d( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof type FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set size for ddivstrain/du
        mddivstraindu( tDofIndex ).set_size(                     //
                2,                                               //
                tFI->get_number_of_space_time_coefficients(),    //
                0.0 );

        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // get the 2nd order derivative of the shape functions d2Ndx2
            const Matrix< DDRMat >& tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

            // get number of space time bases
            uint tNumBases = tFI->get_number_of_space_time_bases();

            // fill ddivstrain/du
            mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )             = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 );
            mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 2 );
            mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = 0.5 * tVelocityd2Ndx2.get_row( 2 );
            mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_ddivstraindu_3d( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof type FI
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set size for ddivstrain/du
        mddivstraindu( tDofIndex ).set_size( 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // get the 2nd order derivative of the shape functions d2Ndx2
            const Matrix< DDRMat >& tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

            // get number of space time bases
            const uint tNumBases = tFI->get_number_of_space_time_bases();

            // fill ddivstrain/du
            mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
            mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 5 );
            mddivstraindu( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 4 );

            mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 5 );
            mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
            mddivstraindu( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 3 );

            mddivstraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 4 );
            mddivstraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 3 );
            mddivstraindu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_dstraindx_2d( uint aOrder )
    {
        switch ( aOrder )
        {
            case 1:
            {
                // set size for dstraindx
                mdStraindx( aOrder - 1 ).set_size( 3, 2 );

                // get the velocity gradient
                const Matrix< DDRMat >& tVelocityGrad =
                        mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

                // fill dstraindx
                mdStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
                mdStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 2, 0 );
                mdStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 2, 1 );
                mdStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
                mdStraindx( aOrder - 1 )( 2, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
                mdStraindx( aOrder - 1 )( 2, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
                break;
            }
            default:
                MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dstraindx_2d - order not supported." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_dstraindx_3d( uint aOrder )
    {
        switch ( aOrder )
        {
            case 1:
            {
                // set size for dstraindx
                mdStraindx( aOrder - 1 ).set_size( 6, 3 );

                // get the velocity gradient
                const Matrix< DDRMat >& tVelocityGrad =
                        mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

                // fill dstraindx
                mdStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
                mdStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 5, 0 );
                mdStraindx( aOrder - 1 )( 0, 2 ) = tVelocityGrad( 4, 0 );
                mdStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 5, 1 );
                mdStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
                mdStraindx( aOrder - 1 )( 1, 2 ) = tVelocityGrad( 3, 1 );
                mdStraindx( aOrder - 1 )( 2, 0 ) = tVelocityGrad( 4, 2 );
                mdStraindx( aOrder - 1 )( 2, 1 ) = tVelocityGrad( 3, 2 );
                mdStraindx( aOrder - 1 )( 2, 2 ) = tVelocityGrad( 2, 2 );
                mdStraindx( aOrder - 1 )( 3, 0 ) = 0.5 * ( tVelocityGrad( 4, 1 ) + tVelocityGrad( 5, 2 ) );
                mdStraindx( aOrder - 1 )( 3, 1 ) = 0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) );
                mdStraindx( aOrder - 1 )( 3, 2 ) = 0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
                mdStraindx( aOrder - 1 )( 4, 0 ) = 0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) );
                mdStraindx( aOrder - 1 )( 4, 1 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 5, 2 ) );
                mdStraindx( aOrder - 1 )( 4, 2 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
                mdStraindx( aOrder - 1 )( 5, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
                mdStraindx( aOrder - 1 )( 5, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) );
                mdStraindx( aOrder - 1 )( 5, 2 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 4, 1 ) );

                break;
            }
            default:
            {
                MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dstraindx_3d - order not supported." );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible::dstraindxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            const Matrix< DDRMat >&        aJump,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                "CM_Fluid_Incompressible::dstraindxdu - Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,
                "CM_Fluid_Incompressible::dstraindxdu - Works only for 1st order derivative for now." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdStraindxduEval( aOrder - 1, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dstraindxdu( aDofTypes, aOrder, aJump );

            // set bool for evaluation
            mdStraindxduEval( aOrder - 1, tDofIndex ) = false;
        }

        // return the derivative
        return mdStraindxdu( aOrder - 1 )( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void CM_Fluid_Incompressible::eval_dstraindxdu_2d(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            const Matrix< DDRMat >&        aJump )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init storage
        mdStraindxdu( aOrder - 1 )( tDofIndex ).set_size(    //
                mSpaceDim,                                   //
                tFI->get_number_of_space_time_coefficients() );

        // if dof is velocity
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // get number of bases for velocity
            uint tNumBases = tFI->get_number_of_space_time_bases();

            // get the velocity gradient
            const Matrix< DDRMat >& tdnNdxn = tFI->dnNdxn( 2 );

            // fill dstraindxdu
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =    //
                    aJump( 0 ) * tdnNdxn.get_row( 0 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 2 );
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =    //
                    aJump( 1 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 0 );

            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =    //
                    aJump( 0 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 1 );
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =    //
                    aJump( 1 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 2 ) * tdnNdxn.get_row( 2 );
        }
        else
        {
            mdStraindxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void CM_Fluid_Incompressible::eval_dstraindxdu_3d(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder,
            const Matrix< DDRMat >&        aJump )
    {
        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init storage
        mdStraindxdu( aOrder - 1 )( tDofIndex ).set_size(    //
                mSpaceDim,                                   //
                tFI->get_number_of_space_time_coefficients() );

        // if dof is velocity
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // get number of bases for velocity
            uint tNumBases = tFI->get_number_of_space_time_bases();

            // get the velocity gradient
            const Matrix< DDRMat >& tdnNdxn = tFI->dnNdxn( 2 );

            // fill dstraindxdu
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =                                                       //
                    aJump( 0 ) * tdnNdxn.get_row( 0 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 2 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =                                           //
                    aJump( 1 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 0 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                       //
                    aJump( 2 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 0 );

            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =                                                       //
                    aJump( 0 ) * tdnNdxn.get_row( 5 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 1 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =                                           //
                    aJump( 1 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 5 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                       //
                    aJump( 2 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 1 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 5 );

            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) =                                                       //
                    aJump( 0 ) * tdnNdxn.get_row( 4 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 3 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =                                           //
                    aJump( 1 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 5 ) * tdnNdxn.get_row( 4 );    //
            mdStraindxdu( aOrder - 1 )( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =                                       //
                    aJump( 2 ) * tdnNdxn.get_row( 2 ) + 0.5 * aJump( 3 ) * tdnNdxn.get_row( 3 ) + 0.5 * aJump( 4 ) * tdnNdxn.get_row( 4 );
        }
        else
        {
            mdStraindxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_teststrain_2d()
    {
        // get velocity field interpolator
        Field_Interpolator* tFIVelocity =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // compute velocity gradient
        const Matrix< DDRMat >& tdnNdxn = tFIVelocity->dnNdxn( 1 );

        // get number of bases for velocity
        const uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // build the test strain
        mTestStrain.set_size( 3, tNumBases * 2, 0.0 );

        mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_teststrain_3d()
    {
        // get velocity field interpolator
        Field_Interpolator* tFIVelocity =
                mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // compute displacement gradient
        const Matrix< DDRMat >& tdnNdxn = tFIVelocity->dnNdxn( 1 );

        // get number of bases for displacement
        const uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

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

    void
    CM_Fluid_Incompressible::eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init mdFluxdDof
        mdFluxdDof( tDofIndex ).set_size(    //
                ( mSpaceDim - 1 ) * 3,       //
                tFI->get_number_of_space_time_coefficients() );

        // if velocity dof
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // build dfluxdv
            mdFluxdDof( tDofIndex ) =
                    2.0 * mPropViscosity->val()( 0 ) * this->dStraindDOF( aDofTypes );
        }
        // if pressure dof
        else if ( aDofTypes( 0 ) == mDofPressure )
        {
            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
            Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );

            tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

            // build the dfluxdp
            mdFluxdDof( tDofIndex ) = -1.0 * tII * tFI->N();
        }
        else
        {
            mdFluxdDof( tDofIndex ).fill( 0.0 );
        }

        // if indirect dependency on the dof type through viscosity
        if ( mPropViscosity->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mdFluxdDof( tDofIndex ) +=
                    2.0 * this->strain() * mPropViscosity->dPropdDOF( aDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // get the dof type as a uint
        const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // flatten normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute dtractiondu
        mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the test dof FI
        Field_Interpolator* tFITest =
                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

        // get the derivative dof FI
        Field_Interpolator* tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init the dTestTractiondDof
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(    //
                tFITest->get_number_of_space_time_coefficients(),     //
                tFIDer->get_number_of_space_time_coefficients() );

        // flatten normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // if viscosity property depends on test or derivative dof type
        if ( mPropViscosity->check_dof_dependency( aDofTypes ) )
        {
            // compute contribution to dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                    2.0 * trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *    //
                    aJump * mPropViscosity->dPropdDOF( aDofTypes );
        }
        else
        {
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).fill( 0.0 );
        }

        // if viscosity property depends on test or derivative dof type
        if ( mPropViscosity->check_dof_dependency( aTestDofTypes ) )
        {
            // compute contribution to dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                    2.0 * trans( mPropViscosity->dPropdDOF( aTestDofTypes ) ) * trans( aJump ) *    //
                    tFlatNormal * this->dStraindDOF( aDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::eval_dStraindDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator* tFI =    //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init mdStraindDof
        mdStraindDof( tDofIndex ).set_size(    //
                ( mSpaceDim - 1 ) * 3,         //
                tFI->get_number_of_space_time_coefficients() );

        // if velocity dof
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute derivative
            mdStraindDof( tDofIndex ) = this->testStrain();
        }
        else
        {
            mdStraindDof( tDofIndex ).fill( 0.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::flatten_normal_2d(
            const Matrix< DDRMat >& aNormal,
            Matrix< DDRMat >&       aFlatNormal )
    {
        aFlatNormal.set_size( 2, 3, 0.0 );
        aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
        aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::flatten_normal_3d(
            const Matrix< DDRMat >& aNormal,
            Matrix< DDRMat >&       aFlatNormal )
    {
        aFlatNormal.set_size( 3, 6, 0.0 );
        aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
        aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
        aFlatNormal( 2, 2 ) = aNormal( 2, 0 );
        aFlatNormal( 0, 4 ) = aNormal( 2, 0 );
        aFlatNormal( 0, 5 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 3 ) = aNormal( 2, 0 );
        aFlatNormal( 1, 5 ) = aNormal( 0, 0 );
        aFlatNormal( 2, 3 ) = aNormal( 1, 0 );
        aFlatNormal( 2, 4 ) = aNormal( 0, 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible::select_derivative_FD(
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
            default:
                MORIS_ERROR( false, "CM_Fluid_Incompressible::select_derivative_FD: aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible::set_derivative_FD(
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
            default:
                MORIS_ERROR( false, "CM_Fluid_Incompressible::set_derivative_FD: aCMRequestType undefined" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
