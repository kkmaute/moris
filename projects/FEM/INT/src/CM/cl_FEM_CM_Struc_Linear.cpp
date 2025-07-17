/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear.cpp
 *
 */

#include "cl_FEM_CM_Struc_Linear.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_dot.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    CM_Struc_Linear::CM_Struc_Linear()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type_Lin::MAX_ENUM ), nullptr );

        mPropertyMap[ "CTE" ]                  = static_cast< uint >( CM_Property_Type_Lin::CTE );
        mPropertyMap[ "PropertyTemperature" ]  = static_cast< uint >( CM_Property_Type_Lin::TEMP_PROP );
        mPropertyMap[ "ReferenceTemperature" ] = static_cast< uint >( CM_Property_Type_Lin::TEMP_REF );
        mPropertyMap[ "AxisymRotationAxis" ]   = static_cast< uint >( CM_Property_Type_Lin::ROT_AXI );
        mPropertyMap[ "EigenStrain" ]          = static_cast< uint >( CM_Property_Type_Lin::EIGEN_STRAIN );
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > > &aDofTypes,
            const Vector< std::string >             &aDofStrings )
    {
        // set dof type list
        Constitutive_Model::set_dof_type_list( aDofTypes );

        // loop over the provided dof type
        for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
        {
            // get dof type string
            const std::string &tDofString = aDofStrings( iDof );

            // get dof type
            MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

            // if displacement dof type string
            if ( tDofString == "Displacement" )
            {
                mDofDispl = tDofType;
            }
            // if temperature dof type string
            else if ( tDofString == "Temperature" )
            {
                mDofTemp = tDofType;
            }
            // if pressure dof type string
            else if ( tDofString == "Pressure" )
            {
                mDofPressure = tDofType;
            }
            else
            {
                // error unknown dof string
                MORIS_ERROR( false,
                        "CM_Struc_Linear::set_dof_type_list - Unknown aDofString : %s \n",
                        tDofString.c_str() );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear::initialize_spec_storage_vars_and_eval_flags()
    {
        mGeometricStiffnessEval.set_size( mGlobalDofTypes.size(), 1, true );
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear::set_local_properties()
    {
        // set the CTE property
        mPropCTE = get_property( "CTE" );

        // set the given temperature property
        mPropTemp = get_property( "PropertyTemperature" );

        // set the reference temperature property
        mPropTRef = get_property( "ReferenceTemperature" );

        // set the reference temperature property
        mPropRotAxis = get_property( "AxisymRotationAxis" );

        // set the eigen-strain property
        mPropEigenStrain = get_property( "EigenStrain" );

        // check that properties needed for thermal strains exist
        if ( mPropCTE )
        {
            MORIS_ASSERT( mPropTRef,
                    "CM_Struc_Linear_Isotropic::set_local_properties - ReferenceTemperature property does not exist.\n" );
        }

        if ( mPlaneType == Model_Type::AXISYMMETRIC )
        {
            MORIS_ASSERT( mPropRotAxis,
                    "CM_Struc_Linear_Isotropic::set_local_properties - Rotation Axis property not defined.\n" );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear::set_function_pointers()
    {
        switch ( mSpaceDim )
        {
            case 2:
            {
                m_eval_strain     = &CM_Struc_Linear::eval_strain_2d;
                m_eval_teststrain = &CM_Struc_Linear::eval_teststrain_2d;
                m_flatten_normal  = &CM_Struc_Linear::flatten_normal_2d;

                switch ( mPlaneType )
                {
                    case Model_Type::PLANE_STRESS:
                    {
                        mStrain.set_size( 3, 1, 0.0 );
                        mConst.set_size( 3, 3, 0.0 );
                        m_eval_inv_bulk_modulus = &CM_Struc_Linear::eval_inv_bulk_modulus_plane_stress;

                        // list number of normal stresses and strains
                        mNumNormalStress = 2;
                        mNumNormalStrain = 2;

                        switch ( mTensorType )
                        {
                            case Model_Type::FULL:
                            {
                                mConstFunc = &CM_Struc_Linear::full_plane_stress;
                                break;
                            }
                            case Model_Type::DEVIATORIC:
                            {
                                mConstFunc = &CM_Struc_Linear::deviatoric_plane_stress;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false, "Only full and deviatoric tensors implemented for plane stress" );
                            }
                        }
                        break;
                    }
                    case Model_Type::PLANE_STRAIN:
                    {
                        mStrain.set_size( 4, 1, 0.0 );
                        mConst.set_size( 4, 4, 0.0 );

                        // list number of normal stresses and strains
                        mNumNormalStress = 3;
                        mNumNormalStrain = 3;

                        switch ( mTensorType )
                        {
                            case Model_Type::FULL:
                            {
                                mConstFunc = &CM_Struc_Linear::full_plane_strain;
                                break;
                            }
                            case Model_Type::DEVIATORIC:
                            {
                                mConstFunc = &CM_Struc_Linear::deviatoric_plane_strain;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false, "Only full and deviatoric tensors implemented for plane strain" );
                            }
                        }
                        break;
                    }
                    case Model_Type::AXISYMMETRIC:
                    {
                        mStrain.set_size( 4, 1, 0.0 );
                        mConst.set_size( 4, 4, 0.0 );

                        // list number of normal stresses and strains
                        mNumNormalStress = 3;
                        mNumNormalStrain = 3;

                        switch ( mTensorType )
                        {
                            case Model_Type::FULL:
                            {
                                mConstFunc = &CM_Struc_Linear::full_axisymmetric;
                                break;
                            }
                            case Model_Type::DEVIATORIC:
                            {
                                mConstFunc = &CM_Struc_Linear::deviatoric_axisymmetric;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false, "Only full and deviatoric tensors implemented for axisymmetric" );
                            }
                        }
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false,
                                "Linear isotropic elasticity in 2d requires "
                                "plane stress, plane strain, or axisymmetric models" );
                    }
                }
                break;
            }
            case 3:
            {
                m_eval_strain     = &CM_Struc_Linear::eval_strain_3d;
                m_eval_teststrain = &CM_Struc_Linear::eval_teststrain_3d;
                m_flatten_normal  = &CM_Struc_Linear::flatten_normal_3d;

                mStrain.set_size( 6, 1, 0.0 );
                mConst.set_size( 6, 6, 0.0 );

                // list number of normal stresses and strains
                mNumNormalStress = 3;
                mNumNormalStrain = 3;

                switch ( mTensorType )
                {
                    case Model_Type::FULL:
                    {
                        mConstFunc = &CM_Struc_Linear::full_3d;
                        break;
                    }
                    case Model_Type::DEVIATORIC:
                    {
                        mConstFunc = &CM_Struc_Linear::deviatoric_3d;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "Only full and deviatoric tensors implemented for plane strain" );
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Linear isotropic elasticity implemented only for 2d and 3d" );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear::reset_specific_eval_flags()
    {
        mGeometricStiffnessEval.fill( true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_flux()
    {
        // compute flux
        mFlux = this->constitutive() * this->strain();

        // if pressure dof
        if ( mDofPressure != MSI::Dof_Type::UNDEFINED )
        {
            // get the pressure FI
            Field_Interpolator *tPressureFI =
                    mFIManager->get_field_interpolators_for_type( mDofPressure );

            // create identity matrix
            const Matrix< DDRMat > tI( mNumNormalStress, 1, 1.0 );

            // evaluate pressure contribution to flux
            Matrix< DDRMat > tP( mConst.n_rows(), 1, 0.0 );
            tP( { 0, mNumNormalStress - 1 }, { 0, 0 } ) = tI * tPressureFI->val();

            // add contribution to the flux
            mFlux -= tP;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_traction( const Matrix< DDRMat > &aNormal )
    {
        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute the traction
        mTraction = tFlatNormal * this->flux();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_testTraction(
            const Matrix< DDRMat >        &aNormal,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // if test traction wrt displacement
        if ( aTestDofTypes( 0 ) == mDofDispl )
        {
            // compute test traction wrt displacement
            mTestTraction( tTestDofIndex ) = tFlatNormal * this->constitutive() * this->testStrain();
        }

        // if test traction wrt pressure
        if ( aTestDofTypes( 0 ) == mDofPressure )
        {
            // compute test traction wrt pressure
            mTestTraction( tTestDofIndex ) = tFlatNormal * this->dFluxdDOF( aTestDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_strain_2d()
    {
        // get the displacement spatial gradient from displacement FI
        const Matrix< DDRMat > &tDisplGradx =
                mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

        // evaluate the strain
        mStrain.fill( 0.0 );

        // x and y normal strains
        mStrain( 0 ) = tDisplGradx( 0, 0 );
        mStrain( 1 ) = tDisplGradx( 1, 1 );

        // assign normal strain in azimuthal direction
        if ( mPlaneType == Model_Type::AXISYMMETRIC )
        {
            // get the displacements and outward radial vector for azimuthal strain
            const Matrix< DDRMat > &tDispl      = mFIManager->get_field_interpolators_for_type( mDofDispl )->val();
            const Matrix< DDRMat > &tOtbdRadVec = mPropRotAxis->val();

            // adjust radius to be larger than MORIS_REAL_EPS
            const real tRadius = std::max( tOtbdRadVec( 1 ), MORIS_REAL_EPS );

            // normal strain in azimuthal direction u_r / r
            // here {u}.*{n_r} / (r) where {n_r} = unit outward radial vector from line to point
            mStrain( 2 ) = dot( tDispl, tOtbdRadVec( { 2, 3 } ) ) / tRadius;

            // 12 shear stress
            mStrain( 3 ) = tDisplGradx( 1, 0 ) + tDisplGradx( 0, 1 );
        }

        // plane strain
        else if ( mPlaneType == Model_Type::PLANE_STRAIN )
        {
            mStrain( 2 ) = 0.0;
            mStrain( 3 ) = tDisplGradx( 1, 0 ) + tDisplGradx( 0, 1 );
        }

        // 12 shear stress for plane stress or plane strain
        else
        {
            mStrain( 2 ) = tDisplGradx( 1, 0 ) + tDisplGradx( 0, 1 );
        }

        // if thermal expansion
        if ( mPropCTE )
        {
            // build thermal expansion vector
            Matrix< DDRMat > tThermalExpansionVector( mStrain.numel(), 1, 0.0 );
            Matrix< DDRMat > tI( mNumNormalStrain, 1, 1.0 );

            tThermalExpansionVector( { 0, mNumNormalStrain - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

            // get temperature field interpolator
            Field_Interpolator *tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemp );

            // check that a unique definition of temperature is provided
            MORIS_ASSERT( !mPropTemp || !tFITemp,
                    "CM_Struc_Linear::eval_strain_2d - cannot specify both temperature as dof and as a property.\n" );

            // check if temperature as a property is defined
            if ( mPropTemp )
            {
                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - mPropTemp->val() );
            }
            else
            {
                // check that temperature interpolator exists
                MORIS_ASSERT( tFITemp,
                        "CM_Struc_Linear::eval_strain_2d - temperature interpolator does not exist.\n" );

                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - tFITemp->val() );
            }
        }

        // if eigen-strain is defined
        if ( mPropEigenStrain != nullptr )
        {
            // add eigen strain contribution; note the sign
            mStrain += mPropEigenStrain->val();
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mStrain ),
                "CM_Struc_Linear::eval_strain_2d - Strain contains NAN or INF, exiting!" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_strain_3d()
    {
        // get the displacement spatial gradient from displacement FI
        const Matrix< DDRMat > &tDisplGradx =
                mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

        // evaluate the strain
        mStrain( 0, 0 ) = tDisplGradx( 0, 0 );
        mStrain( 1, 0 ) = tDisplGradx( 1, 1 );
        mStrain( 2, 0 ) = tDisplGradx( 2, 2 );
        mStrain( 3, 0 ) = tDisplGradx( 1, 2 ) + tDisplGradx( 2, 1 );
        mStrain( 4, 0 ) = tDisplGradx( 0, 2 ) + tDisplGradx( 2, 0 );
        mStrain( 5, 0 ) = tDisplGradx( 0, 1 ) + tDisplGradx( 1, 0 );

        // if thermal expansion
        if ( mPropCTE != nullptr )
        {
            // build thermal expansion vector
            Matrix< DDRMat > tThermalExpansionVector( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );

            tThermalExpansionVector( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

            // get temperature field interpolator
            Field_Interpolator *tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemp );

            // check that a unique definition of temperature is provided
            MORIS_ASSERT( !mPropTemp || !tFITemp,
                    "CM_Struc_Linear::eval_strain_2d - cannot specify both temperature as dof and as a property.\n" );

            // check if temperature as a property is defined
            if ( mPropTemp )
            {
                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - mPropTemp->val() );
            }
            else
            {
                // check that temperature interpolator exists
                MORIS_ASSERT( tFITemp,
                        "CM_Struc_Linear::eval_strain_2d - temperature interpolator does not exist.\n" );

                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - tFITemp->val() );
            }
        }

        // if eigen-strain is defined
        if ( mPropEigenStrain != nullptr )
        {
            // add eigen strain contribution; note the sign
            mStrain += mPropEigenStrain->val();
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mStrain ),
                "CM_Struc_Linear::eval_strain_3d - Strain contains NAN or INF, exiting!" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_teststrain_2d()
    {
        // get displacement field interpolator
        Field_Interpolator *tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

        // compute displacement gradient
        const Matrix< DDRMat > &tdnNdxn = tFIDispl->dnNdxn( 1 );

        // get number of bases for displacement
        uint tNumBases = tFIDispl->get_number_of_space_time_bases();

        // build the test strain
        mTestStrain.set_size( mStrain.numel(), tNumBases * 2, 0.0 );

        // [dN/dx1]
        mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

        // [dN/dx2]
        mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        if ( mPlaneType == Model_Type::AXISYMMETRIC )
        {
            // get the displacements and outward radial vector for azimuthal strain
            // in the form tOtbdRadVec = {{2*pi*r},{r},{n1},{n2}}
            // where n1 and n2 are components of the unit outboard normal
            const Matrix< DDRMat > &tOtbdRadVec = mPropRotAxis->val();

            // compute interpolation function and location
            const Matrix< DDRMat > &tN = tFIDispl->NBuild();

            /*
             * Axisymmetric strain using u_r and radial location
             * This is essentially [N]*{u}.*{n_r}/r.
             * Since {u} = {{u1},{0}} for the u1 vector,  {u}.*{n_r} = u1*n_r1
             * Same goes for u2 direction
             */
            mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } )             = tN * tOtbdRadVec( 2 ) / tOtbdRadVec( 1 );
            mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tN * tOtbdRadVec( 3 ) / tOtbdRadVec( 1 );

            // [ dN/dX2   dN/dX1 ]
            mTestStrain( { 3, 3 }, { 0, tNumBases - 1 } )             = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        // plane strain (note: row index 2 will be all 0s)
        else if ( mPlaneType == Model_Type::PLANE_STRAIN )
        {
            // [ dN/dX2   dN/dX1 ]
            mTestStrain( { 3, 3 }, { 0, tNumBases - 1 } )             = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        // plane stress
        else
        {
            // [ dN/dX2   dN/dX1 ]
            mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } )             = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_teststrain_3d()
    {
        // get displacement field interpolator
        Field_Interpolator *tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

        // compute displacement gradient
        const Matrix< DDRMat > &tdnNdxn = tFIDispl->dnNdxn( 1 );

        // get number of bases for displacement
        uint tNumBases = tFIDispl->get_number_of_space_time_bases();

        // build the test strain
        mTestStrain.set_size( 6, tNumBases * 3, 0.0 );
        mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mTestStrain( { 4, 4 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 5, 5 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

        mTestStrain( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_const()
    {
        MORIS_ERROR( 0, "eval_const is not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    CM_Struc_Linear::eval_inv_bulk_modulus()
    {
        MORIS_ERROR( 0, "eval_inv_bulk_modulus is not implemented in the base class CM_Struc_Linear" );
        return 0.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_inv_bulk_modulus_generic(
            const real &aNu,
            const real &aEMod,
            real       &aInvBulkModulus )
    {
        MORIS_ERROR( 0, "eval_inv_bulk_modulus_generic is not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_inv_bulk_modulus_plane_stress(
            const real &aNu,
            const real &aEMod,
            real       &aInvBulkModulus )
    {
        MORIS_ERROR( 0, "eval_inv_bulk_modulus_plane_stress is not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    CM_Struc_Linear::eval_dInvBulkModulusdDOF(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        MORIS_ERROR( 0, "eval_dInvBulkModulusdDOF is not implemented in the base class CM_Struc_Linear" );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type as a uint
        const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the dof FI
        Field_Interpolator *tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // initialize mdFluxdDof
        mdFluxdDof( tDofIndex ).set_size( mConst.n_rows(), tFI->get_number_of_space_time_coefficients() );

        // if displacements or temperature
        if ( aDofTypes( 0 ) == mDofDispl || aDofTypes( 0 ) == mDofTemp )
        {
            mdFluxdDof( tDofIndex ) =
                    this->constitutive() * this->dStraindDOF( aDofTypes );
        }
        else
        {
            mdFluxdDof( tDofIndex ).fill( 0.0 );
        }

        // if pressure dof
        if ( aDofTypes( 0 ) == mDofPressure )
        {
            // create identity matrix
            Matrix< DDRMat > tI( mNumNormalStress, 1, 1.0 );
            Matrix< DDRMat > tII( mConst.n_rows(), 1, 0.0 );
            tII( { 0, mNumNormalStress - 1 }, { 0, 0 } ) = tI.matrix_data();

            // get shape function for presure field
            Matrix< DDRMat > tPressureN = tFI->N();

            // build the dfluxdp
            mdFluxdDof( tDofIndex ) -= tII * tPressureN;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_dTractiondDOF(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Matrix< DDRMat >        &aNormal )
    {
        // get the dof type as a uint
        const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // flatten normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute derivative
        mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Matrix< DDRMat >        &aNormal,
            const Matrix< DDRMat >        &aJump,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get test dof type index
        const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // initialize the dTestTractiondDof
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(), mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients() );

        mdTestTractiondDof( tTestDofIndex )( tDofIndex ).fill( 0.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_dStraindDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the dof FI
        Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // initialize mdStraindDof
        mdStraindDof( tDofIndex ).set_size( mStrain.numel(), tFI->get_number_of_space_time_coefficients() );

        // if displacement dof
        if ( aDofTypes( 0 ) == mDofDispl )
        {
            // compute derivative
            mdStraindDof( tDofIndex ) = this->testStrain();
        }
        else
        {
            mdStraindDof( tDofIndex ).fill( 0.0 );
        }

        // if temperature dof
        if ( mPropCTE && aDofTypes( 0 ) == mDofTemp )
        {
            // build thermal expansion vector
            Matrix< DDRMat > tThermalExpansionVector( mStrain.numel(), 1, 0.0 );
            Matrix< DDRMat > tI( mNumNormalStrain, 1, 1.0 );
            tThermalExpansionVector( { 0, mNumNormalStrain - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

            // compute derivatives
            mdStraindDof( tDofIndex ) -= tThermalExpansionVector * tFI->N();
        }

        // if properties depend on dofs
        if ( mPropCTE )
        {
            // create identity matrix
            Matrix< DDRMat > tI( mNumNormalStrain, 1, 1.0 );
            Matrix< DDRMat > tII( mStrain.numel(), 1, 0.0 );

            tII( { 0, mNumNormalStrain - 1 }, { 0, 0 } ) = tI.matrix_data();

            //  dof dependency of CTE
            if ( mPropCTE->check_dof_dependency( aDofTypes ) )
            {
                if ( mPropTemp )
                {
                    // compute derivatives
                    mdStraindDof( tDofIndex ) +=
                            tII * mPropCTE->dPropdDOF( aDofTypes ) * ( mPropTRef->val()( 0 ) - mPropTemp->val()( 0 ) );
                }
                else
                {
                    // get temperature field interpolator
                    Field_Interpolator *tFIT =
                            mFIManager->get_field_interpolators_for_type( mDofTemp );

                    // compute derivatives
                    mdStraindDof( tDofIndex ) +=
                            tII * mPropCTE->dPropdDOF( aDofTypes ) * ( mPropTRef->val()( 0 ) - tFIT->val()( 0 ) );
                }
            }

            //  dof dependency of temperature property
            if ( mPropTemp && mPropTemp->check_dof_dependency( aDofTypes ) )
            {
                // compute derivatives
                mdStraindDof( tDofIndex ) -=
                        tII * mPropCTE->val() * mPropTemp->dPropdDOF( aDofTypes );
            }

            //  dof dependency of reference temperature
            if ( mPropTRef->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR( false,
                        "CM_Struc_Linear::eval_dStraindDOF - %s",
                        "dof dependency of reference temperature not implemented.\n" );
            }
        }

        // if eigen-strain is defined
        if ( mPropEigenStrain != nullptr )
        {
            if ( mPropEigenStrain->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR( false,
                        "CM_Struc_Linear::eval_dStraindDOF - %s",
                        "dof dependency of eigen strain not implemented.\n" );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_dConstdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        MORIS_ERROR( false, "CM_Struc_Linear::eval_dConstdDOF - Not implemented." );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::flatten_normal_2d(
            const Matrix< DDRMat > &aNormal,
            Matrix< DDRMat >       &aFlatNormal )
    {
        // num cols based on number of flux terms
        aFlatNormal.set_size( 2, mConst.n_rows(), 0.0 );
        aFlatNormal( 0, 0 )                   = aNormal( 0, 0 );
        aFlatNormal( 0, mConst.n_rows() - 1 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 1 )                   = aNormal( 1, 0 );
        aFlatNormal( 1, mConst.n_rows() - 1 ) = aNormal( 0, 0 );
    }

    void
    CM_Struc_Linear::flatten_normal_3d(
            const Matrix< DDRMat > &aNormal,
            Matrix< DDRMat >       &aFlatNormal )
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

    void
    CM_Struc_Linear::set_space_dim( uint aSpaceDim )
    {
        // check that space dimension is 1, 2, 3
        MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4,
                "Constitutive_Model::set_space_dim - wrong space dimension." );

        // set space dimension
        mSpaceDim = aSpaceDim;

        // set function pointers
        this->set_function_pointers();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::set_model_type( Model_Type aModelType )
    {
        // fixme: currently cannot set plane type and a tensor type at the same time from an input file
        // store model type based on input
        if ( aModelType == Model_Type::PLANE_STRESS or aModelType == Model_Type::PLANE_STRAIN or aModelType == Model_Type::AXISYMMETRIC )
        {
            mPlaneType = aModelType;
        }
        else if ( aModelType == Model_Type::FULL or aModelType == Model_Type::HYDROSTATIC or aModelType == Model_Type::DEVIATORIC )
        {
            mTensorType = aModelType;
        }
        else
        {
            MORIS_ASSERT( false,
                    "CM_Struc_Linear::set_model_type - Specified linear isotropic elasticity model type doesn't exist." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::full_plane_stress( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "full_plane_stress function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::deviatoric_plane_stress( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "deviatoric_plane_stress function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::full_plane_strain( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "full_plane_strain function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::deviatoric_plane_strain( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "deviatoric_plane_strain function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::full_axisymmetric( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "full_axisymmetric function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::deviatoric_axisymmetric( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "deviatoric_axisymmetric function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::full_3d( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "full_3d function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::deviatoric_3d( std::initializer_list< const real > &&tParams )
    {
        MORIS_ERROR( 0, "deviatoric_3d function not implemented in the base class CM_Struc_Linear" );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Struc_Linear::GeometricStiffness(
            const Vector< MSI::Dof_Type > &aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        switch ( aCMFunctionType )
        {
            // second Piola-Kirchhoff stress
            case CM_Function_Type::DEFAULT:
            {
                // if the derivative of the test traction based on second Piola_Kirchhoff stress was not evaluated
                if ( mGeometricStiffnessEval( tDofIndex ) )
                {
                    // evaluate the traction based on second Piola_Kirchhoff stress
                    this->eval_geometric_stiffness( aDofTypes );

                    // set bool for evaluation
                    mGeometricStiffnessEval( tDofIndex ) = false;
                }
                // return the traction based on second Piola_Kirchhoff stress
                return mGeometricStiffness;
            }
            default:
            {
                MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::GeometricStiffness - Only PK2 implemented." );
                return mGeometricStiffness;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear::eval_geometric_stiffness( const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // call to evaluate stress
        this->eval_flux();

        // convert to stress matrix
        Matrix< DDRMat > tStressMatrix( mSpaceDim, mSpaceDim );

        if ( mSpaceDim == 2 )
        {
            tStressMatrix.set_size( mSpaceDim, mSpaceDim );
            tStressMatrix( 0, 0 ) = mFlux( 0 );
            tStressMatrix( 1, 1 ) = mFlux( 1 );
            tStressMatrix( 0, 1 ) = mFlux( 2 );
            tStressMatrix( 1, 0 ) = mFlux( 2 );
        }
        else
        {
            tStressMatrix( 0, 0 ) = mFlux( 0 );
            tStressMatrix( 1, 1 ) = mFlux( 1 );
            tStressMatrix( 2, 2 ) = mFlux( 2 );
            tStressMatrix( 1, 2 ) = mFlux( 3 );
            tStressMatrix( 2, 1 ) = mFlux( 3 );
            tStressMatrix( 0, 2 ) = mFlux( 4 );
            tStressMatrix( 2, 0 ) = mFlux( 4 );
            tStressMatrix( 0, 1 ) = mFlux( 5 );
            tStressMatrix( 1, 0 ) = mFlux( 5 );
        }

        // get the dof FI
        Field_Interpolator *tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get the derivative of the displacement shape functions
        const Matrix< DDRMat > &tdnNdxn =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->dnNdxn( 1 );

        // init mdFluxdDof
        mGeometricStiffness.set_size( tFI->get_number_of_space_time_coefficients(), tFI->get_number_of_space_time_coefficients(), 0.0 );

        // get number of coefficients
        uint tNumSpaceBases = tFI->get_number_of_space_bases();

        MORIS_ASSERT( tFI->get_number_of_space_time_coefficients() == tNumSpaceBases * mSpaceDim,
                "CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_geometric_stiffness_second_piola_kirchhoff_2d -"
                "number of spatial bases does not match; check for time interpolation order; should be constant" );

        // see MORIS - theory on OneDrive for derivation leading to geometric stiffness below
        uint tOffset = 0;

        for ( uint q = 0; q < mSpaceDim; ++q )
        {
            for ( uint p = 0; p < tNumSpaceBases; ++p )
            {
                for ( uint s = 0; s < tNumSpaceBases; ++s )
                {
                    for ( uint i = 0; i < mSpaceDim; ++i )
                    {
                        for ( uint j = 0; j < mSpaceDim; ++j )
                        {
                            mGeometricStiffness( tOffset + p, tOffset + s ) += tdnNdxn( i, p ) * tdnNdxn( j, s ) * tStressMatrix( i, j );
                        }
                    }
                }
            }
            tOffset += tNumSpaceBases;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
