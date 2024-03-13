/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff.cpp
 *
 */

#include "cl_FEM_CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

#include "fn_clip_value.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff(
                enum Constitutive_Type aConstType )
        {
            // set the property pointer cell size
            mProperties.resize( mProperties.size() + static_cast< uint >( CM_Property_Type_Saint_Venant_Kirchhoff::MAX_ENUM ), nullptr );

            // This index is due to the existence of properties in the parent class
            uint tCurrentIndexOffSet = static_cast< uint >( CM_Property_Type::MAX_ENUM );

            // populate the map
            mPropertyMap[ "YoungsModulus" ] = static_cast< uint >( CM_Property_Type_Saint_Venant_Kirchhoff::EMOD ) + tCurrentIndexOffSet;
            mPropertyMap[ "PoissonRatio" ]  = static_cast< uint >( CM_Property_Type_Saint_Venant_Kirchhoff::NU ) + tCurrentIndexOffSet;
        }

        //------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_local_properties()
        {
            // set the parent properties
            CM_Struc_Nonlinear_Isotropic::set_local_properties();

            // set the Young's modulus property
            mPropEMod = get_property( "YoungsModulus" );

            // set the Poisson ratio property
            mPropPoisson = get_property( "PoissonRatio" );

            // check that essential properties exist
            MORIS_ASSERT( mPropEMod,
                    "CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_local_properties - Young's modulus property does not exist.\n");

            MORIS_ASSERT( mPropPoisson,
                    "CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_local_properties - Poisson ratio property does not exist.\n");
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_function_pointers()
        {
            // set function pointer from the parent class
            CM_Struc_Nonlinear_Isotropic::set_function_pointers();

            // switch on space dimension
            switch( mSpaceDim )
            {
                case 2 :
                {
                    m_eval_d1PKStressdDOF = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_d1PKStressdDOF_2d;
                    m_d1PKNdFdF           = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_symbolic_d1PKNdFdF_2d;

                    switch( mPlaneType )
                    {
                        case Model_Type::PLANE_STRESS :
                        {
                            mConst.set_size( 3, 3, 0.0 );

                            switch( mTensorType )
                            {
                                case Model_Type::FULL :
                                {
                                    mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_plane_stress;
                                    break;
                                }
                                case Model_Type::DEVIATORIC :
                                {
                                    mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_plane_stress;
                                    break;
                                }
                                default:
                                {
                                    MORIS_ERROR(false, "Only full and deviatoric tensors implemented for plane stress");
                                }
                            }
                            break;
                        }
                        case Model_Type::PLANE_STRAIN :
                        {
                            mConst.set_size( 3, 3, 0.0 );

                            switch( mTensorType )
                            {
                                case Model_Type::FULL :
                                {
                                    mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_plane_strain;
                                    break;
                                }
                                case Model_Type::DEVIATORIC :
                                {
                                    mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_plane_strain;
                                    break;
                                }
                                default:
                                {
                                    MORIS_ERROR(false, "Only full and deviatoric tensors implemented for plane strain");
                                }
                            }
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "Nonlinear isotropic elasticity in 2d requires "
                                    "plane stress, plane strain models");
                        }
                    }
                    break;
                }
                case 3 :
                {
                    m_eval_d1PKStressdDOF = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_d1PKStressdDOF_3d;
                    m_d1PKNdFdF           = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_symbolic_d1PKNdFdF_3d;

                    mConst.set_size( 6, 6, 0.0 );

                    switch(mTensorType)
                    {
                        case Model_Type::FULL :
                        {
                            mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_3d;
                            break;
                        }
                        case Model_Type::DEVIATORIC :
                        {
                            mConstFunc = &CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_3d;
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "Only full and deviatoric tensors implemented for plane strain");
                        }
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Nonlinear isotropic elasticity implemented only for 2d and 3d" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_const()
        {
            // get the Poisson's ratio value
            const real tNu = mPropPoisson->val()( 0 );

            // get the Youngs' modulus value
            const real tEmod = mPropEMod->val()( 0 );

            // evaluate the constitutive matrix
            ( this->*mConstFunc )( tEmod, tNu );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_flux_first_piola_kirchhoff()
        {
            // get the second Piola-Kirchhoff stress in full
            Matrix< DDRMat > t2PKStressFull;
            this->voigt_to_full_sym_stress( this->flux( CM_Function_Type::PK2 ), t2PKStressFull );

            // evaluate 1PK stress for Saint Venant-Kirchhoff from 2PK stress
            this->full_to_voigt_nonsym( this->deformation_gradient() * t2PKStressFull, m1PKStress );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_flux_second_piola_kirchhoff()
        {
            // evaluate 2PK stress for Saint Venant-Kirchhoff based on constitutive model
            m2PKStress = this->constitutive() * this->strain( CM_Function_Type::LAGRANGIAN );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_flux_cauchy()
        {
            // get the 2nd Piola-Kirchhoff stress in full form
            Matrix< DDRMat > t2PKStressFull;
            this->voigt_to_full_sym_stress( this->flux( CM_Function_Type::PK2 ), t2PKStressFull );

            // evaluate the Cauchy stress
            Matrix< DDRMat > tCauchyStressFull =
                    this->deformation_gradient() * t2PKStressFull * trans( this->deformation_gradient() ) / this->volume_change_jacobian();

            // cast from full to voigt notation
            this->full_to_voigt_sym_stress( tCauchyStressFull, mCauchyStress );
        }

        //------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_d1PKStressdDOF_2d(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            md1PKStressdu( tDofIndex ).set_size(
                    4,
                    tFI->get_number_of_space_time_coefficients() );

            // if displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // get the PK2 stress
                const Matrix< DDRMat > & t2PKStress = this->flux( CM_Function_Type::PK2 );

                // get the derivative of the PK2 stress
                const Matrix< DDRMat > & td2PKdu = this->dFluxdDOF( aDofTypes, CM_Function_Type::PK2 );

                // get the deformation gradient
                const Matrix< DDRMat > & tF = this->deformation_gradient();

                // get the derivative of the deformation gradient
                const Matrix< DDRMat > & tdFdu = this->test_deformation_gradient( aDofTypes );

                // compute derivative
                md1PKStressdu( tDofIndex ).get_row( 0 ) =
                          tF( 0, 0 ) * td2PKdu.get_row( 0 ) + tF( 0, 1 ) * td2PKdu.get_row( 2 )
                        + t2PKStress( 0 ) * tdFdu.get_row( 0 ) + t2PKStress( 2 ) * tdFdu.get_row( 1 );

                md1PKStressdu( tDofIndex ).get_row( 1 ) =
                        tF( 0, 0 ) * td2PKdu.get_row( 2 ) + tF( 0, 1 ) * td2PKdu.get_row( 1 )
                        + t2PKStress( 2 ) * tdFdu.get_row( 0 ) + t2PKStress( 1 ) * tdFdu.get_row( 1 );

                md1PKStressdu( tDofIndex ).get_row( 2 ) =
                        tF( 1, 0 ) * td2PKdu.get_row( 0 ) + tF( 1, 1 ) * td2PKdu.get_row( 2 )
                        + t2PKStress( 0 ) * tdFdu.get_row( 2 ) + t2PKStress( 2 ) * tdFdu.get_row( 3 );

                md1PKStressdu( tDofIndex ).get_row( 3 ) =
                        tF( 1, 0 ) * td2PKdu.get_row( 2 ) + tF( 1, 1 ) * td2PKdu.get_row( 1 )
                        + t2PKStress( 2 ) * tdFdu.get_row( 2 ) + t2PKStress( 1 ) * tdFdu.get_row( 3 );
            }

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - elastic modulus dependency on DOF not implemented yet" );
            }

            // if Poisson ratio depends on dof type
            if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - Poisson's ratio dependency on DOF not implemented yet" );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_d1PKStressdDOF_3d(
        const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            md1PKStressdu( tDofIndex ).set_size(
                    9,
                    tFI->get_number_of_space_time_coefficients() );

            // if displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // get the PK2 stress
                const Matrix< DDRMat > & t2PKStress = this->flux( CM_Function_Type::PK2 );

                // get the derivative of the PK2 stress
                const Matrix< DDRMat > & td2PKdu = this->dFluxdDOF( aDofTypes, CM_Function_Type::PK2 );

                // get the deformation gradient
                const Matrix< DDRMat > & tF = this->deformation_gradient();

                // get the derivative of the deformation gradient
                const Matrix< DDRMat > & tdFdu = this->test_deformation_gradient( aDofTypes );

                // compute derivative
                md1PKStressdu( tDofIndex ).get_row( 0 ) =
                          tF( 0, 0 ) * td2PKdu.get_row( 0 ) + tF( 0, 1 ) * td2PKdu.get_row( 5 ) + tF( 0, 2 ) * td2PKdu.get_row( 4 )
                        + t2PKStress( 0 ) * tdFdu.get_row( 0 ) + t2PKStress( 5 ) * tdFdu.get_row( 1 ) + t2PKStress( 4 ) * tdFdu.get_row( 2 );

                md1PKStressdu( tDofIndex ).get_row( 1 ) =
                        tF( 0, 0 ) * td2PKdu.get_row( 5 ) + tF( 0, 1 ) * td2PKdu.get_row( 1 )  + tF( 0, 2 ) * td2PKdu.get_row( 3 )
                        + t2PKStress( 5 ) * tdFdu.get_row( 0 ) + t2PKStress( 1 ) * tdFdu.get_row( 1 ) + t2PKStress( 3 ) * tdFdu.get_row( 2 );

                md1PKStressdu( tDofIndex ).get_row( 2 ) =
                        tF( 0, 0 ) * td2PKdu.get_row( 4 ) + tF( 0, 1 ) * td2PKdu.get_row( 3 )  + tF( 0, 2 ) * td2PKdu.get_row( 2 )
                        + t2PKStress( 4 ) * tdFdu.get_row( 0 ) + t2PKStress( 3 ) * tdFdu.get_row( 1 ) + t2PKStress( 2 ) * tdFdu.get_row( 2 );

                md1PKStressdu( tDofIndex ).get_row( 3 ) =
                          tF( 1, 0 ) * td2PKdu.get_row( 0 ) + tF( 1, 1 ) * td2PKdu.get_row( 5 ) + tF( 1, 2 ) * td2PKdu.get_row( 4 )
                        + t2PKStress( 0 ) * tdFdu.get_row( 3 ) + t2PKStress( 5 ) * tdFdu.get_row( 4 ) + t2PKStress( 4 ) * tdFdu.get_row( 5 );

                md1PKStressdu( tDofIndex ).get_row( 4 ) =
                        tF( 1, 0 ) * td2PKdu.get_row( 5 ) + tF( 1, 1 ) * td2PKdu.get_row( 1 )  + tF( 1, 2 ) * td2PKdu.get_row( 3 )
                        + t2PKStress( 5 ) * tdFdu.get_row( 3 ) + t2PKStress( 1 ) * tdFdu.get_row( 4 ) + t2PKStress( 3 ) * tdFdu.get_row( 5 );

                md1PKStressdu( tDofIndex ).get_row( 5 ) =
                        tF( 1, 0 ) * td2PKdu.get_row( 4 ) + tF( 1, 1 ) * td2PKdu.get_row( 3 )  + tF( 1, 2 ) * td2PKdu.get_row( 2 )
                        + t2PKStress( 4 ) * tdFdu.get_row( 3 ) + t2PKStress( 3 ) * tdFdu.get_row( 4 ) + t2PKStress( 2 ) * tdFdu.get_row( 5 );

                md1PKStressdu( tDofIndex ).get_row( 6 ) =
                        tF( 2, 0 ) * td2PKdu.get_row( 0 ) + tF( 2, 1 ) * td2PKdu.get_row( 5 ) + tF( 2, 2 ) * td2PKdu.get_row( 4 )
                        + t2PKStress( 0 ) * tdFdu.get_row( 6 ) + t2PKStress( 5 ) * tdFdu.get_row( 7 ) + t2PKStress( 4 ) * tdFdu.get_row( 8 );

                md1PKStressdu( tDofIndex ).get_row( 7 ) =
                        tF( 2, 0 ) * td2PKdu.get_row( 5 ) + tF( 2, 1 ) * td2PKdu.get_row( 1 )  + tF( 2, 2 ) * td2PKdu.get_row( 3 )
                        + t2PKStress( 5 ) * tdFdu.get_row( 6 ) + t2PKStress( 1 ) * tdFdu.get_row( 7 ) + t2PKStress( 3 ) * tdFdu.get_row( 8 );

                md1PKStressdu( tDofIndex ).get_row( 8 ) =
                        tF( 2, 0 ) * td2PKdu.get_row( 4 ) + tF( 2, 1 ) * td2PKdu.get_row( 3 )  + tF( 2, 2 ) * td2PKdu.get_row( 2 )
                        + t2PKStress( 4 ) * tdFdu.get_row( 6 ) + t2PKStress( 3 ) * tdFdu.get_row( 7 ) + t2PKStress( 2 ) * tdFdu.get_row( 8 );

            }

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - elastic modulus dependency on DOF not implemented yet" );
            }

            // if Poisson ratio depends on dof type
            if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - Poisson's ratio dependency on DOF not implemented yet" );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_d2PKStressdDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            md2PKStressdu( tDofIndex ).set_size(
                    mConst.n_rows(),
                    tFI->get_number_of_space_time_coefficients() );

            // if displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // compute derivative
                md2PKStressdu( tDofIndex ) =
                        this->constitutive() * this->dStraindDOF( aDofTypes, CM_Function_Type::LAGRANGIAN );
            }

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - elastic modulus dependency on DOF not implemented yet" );
            }

            // if Poisson ratio depends on dof type
            if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - Poisson's ratio dependency on DOF not implemented yet" );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dCauchyStressdDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            mdCauchyStressdu( tDofIndex ).set_size(
                    mConst.n_rows(),
                    tFI->get_number_of_space_time_coefficients() );

            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dCauchyStressdDOF - Not implemented yet.");
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_traction_first_piola_kirchhoff(
                const Matrix< DDRMat > & aNormal )
        {
            // get first Piola-Kirchhoff in full notation
            Matrix< DDRMat > t1PKStressFull;
            this->voigt_to_full_nonsym( this->flux( CM_Function_Type::PK1 ), t1PKStressFull );

            // evaluate traction based on 1PK stress
            m1PKTraction = t1PKStressFull * aNormal;
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_traction_cauchy(
                const Matrix< DDRMat > & aNormal )
        {
            // FIXME need to implement traction based on cauchy stress with normal in current config
            MORIS_ASSERT( false, "CM_Struc_Nonlinear_Isotropic::eval_traction_cauchy - Not implemented yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTractiondDOF_first_piola_kirchhoff(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

           // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal_nonsym( aNormal, tFlatNormal );

            // compute test traction wrt dof
            md1PKTractiondu( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes, CM_Function_Type::PK1 );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTractiondDOF_cauchy(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // FIXME need to implement traction based on cauchy stress with normal in current configuration
            MORIS_ASSERT( false, "CM_Struc_Nonlinear_Isotropic::eval_dTractiondDOF_cauchy - Not implemented yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_testTraction_first_piola_kirchhoff(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal_nonsym( aNormal, tFlatNormal );

            // compute test traction wrt dof
            m1PKTestTraction( tTestDofIndex ) =  this->dTractiondDOF( aTestDofTypes, aNormal, CM_Function_Type::PK1 );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_testTraction_cauchy(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // FIXME need to implement traction based on cauchy stress with normal in current configuration
            MORIS_ASSERT( false, "CM_Struc_Nonlinear_Isotropic::eval_testTraction_cauchy - Not implemented yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTestTractiondDOF_first_piola_kirchhoff(
                const Vector< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >      & aNormal,
                const Matrix< DDRMat >      & aJump,
                const Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // init the dTestTractiondDof
            md1PKTestTractiondu( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients());

            // get the Poisson's ratio value
            const real tNu = mPropPoisson->val()( 0 );

            // get the Youngs' modulus value
            const real tEmod = mPropEMod->val()( 0 );

            // compute the Lame's constants
            real tLame1 = tEmod * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
            real tLame2 = tEmod / ( 2.0 * ( 1.0 + tNu ) );

            // compute the second order derivative of the traction based on the first Piola-Kirchhoff stress wrt the deformation gradient
            // d2P_ij n_j/dF_kl dF_mn
            Matrix< DDRMat > td1PKNdFdF;
            this->d1PKNdFdF( td1PKNdFdF, tLame1, tLame2, aNormal );

            // flatten the jump
            Matrix< DDRMat > tJumpFlat;
            proj_jump( aJump, tJumpFlat );

            // compute the second order derivative of the traction based on the first Piola-Kirchhoff stress wrt the required dof
            // d2P_ij n_j/du_k du_l
            md1PKTestTractiondu( tTestDofIndex )( tDofIndex ) =
                    trans( this->dStraindDOF( aDofTypes, CM_Function_Type::DEFORMATION_GRADIENT ) ) *    //
                    trans( td1PKNdFdF ) * tJumpFlat *                                                    //
                    this->dStraindDOF( aDofTypes, CM_Function_Type::DEFORMATION_GRADIENT );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTestTractiondDOF_cauchy(
                const Vector< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >      & aNormal,
                const Matrix< DDRMat >      & aJump,
                const Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dTestTractiondDOF_cauchy - Not implemented." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dConstdDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dConstdDOF - Not implemented." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_space_dim(
                uint aSpaceDim )
        {
            // check that space dimension is 1, 2, 3
            MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4,
                    "Constitutive_Model::set_space_dim - wrong space dimension.");

            // set space dimension
            mSpaceDim = aSpaceDim;

            // set function pointers
            this->set_function_pointers();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_plane_stress(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / ( 1 - std::pow( aNu, 2 ) );

            mConst( 0, 0 ) = tPre;
            mConst( 1, 1 ) = tPre;
            mConst( 0, 1 ) = tPre * aNu;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 2, 2 ) = tPre * 0.5 * (1.0 - aNu );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_plane_stress(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / ((1 + aNu) * 2.0);

            mConst( 0, 0 ) =  tPre;
            mConst( 1, 1 ) =  tPre;
            mConst( 0, 1 ) = -tPre;
            mConst( 1, 0 ) = -tPre;
            mConst( 2, 2 ) =  tPre;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_plane_strain(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / (1.0 + aNu ) / (1.0 - 2.0 * aNu ) ;

            mConst( 0, 0 ) = tPre * ( 1.0 - aNu );
            mConst( 0, 1 ) = tPre * aNu;
            mConst( 0, 2 ) = 0.0;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 1, 1 ) = tPre * ( 1.0 - aNu );
            mConst( 1, 2 ) = 0.0;
            mConst( 2, 0 ) = 0.0;
            mConst( 2, 1 ) = 0.0;
            mConst( 2, 2 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_plane_strain(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / (3.0 * (1.0 + aNu ) );

            mConst( 0, 0 ) = tPre * 4.0;
            mConst( 0, 1 ) = tPre;
            mConst( 0, 2 ) = tPre;
            mConst( 1, 0 ) = tPre;
            mConst( 1, 1 ) = tPre * 4.0;
            mConst( 1, 2 ) = tPre;
            mConst( 2, 0 ) = tPre;
            mConst( 2, 1 ) = tPre;
            mConst( 2, 2 ) = tPre * 4.0;
            mConst( 3, 3 ) = tPre * 3.0 / 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::full_3d(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / (1.0 + aNu ) / (1.0 - 2.0 * aNu );

            mConst( 0, 0 ) = tPre * ( 1.0 - aNu );
            mConst( 0, 1 ) = tPre * aNu;
            mConst( 0, 2 ) = tPre * aNu;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 1, 1 ) = tPre * ( 1.0 - aNu );
            mConst( 1, 2 ) = tPre * aNu;
            mConst( 2, 0 ) = tPre * aNu;
            mConst( 2, 1 ) = tPre * aNu;
            mConst( 2, 2 ) = tPre * ( 1.0 - aNu );
            mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
            mConst( 4, 4 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
            mConst( 5, 5 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::deviatoric_3d(
                const real & aEmod,
                const real & aNu )
        {
            const real tPre = aEmod / ( 3.0 * ( 1.0 + aNu ) );

            mConst( 0, 0 ) = tPre * 4.0;
            mConst( 0, 1 ) = tPre;
            mConst( 0, 2 ) = tPre;
            mConst( 1, 0 ) = tPre;
            mConst( 1, 1 ) = tPre * 4.0;
            mConst( 1, 2 ) = tPre;
            mConst( 2, 0 ) = tPre;
            mConst( 2, 1 ) = tPre;
            mConst( 2, 2 ) = tPre * 4.0;
            mConst( 3, 3 ) = tPre * 3.0 / 2.0;
            mConst( 4, 4 ) = tPre * 3.0 / 2.0;
            mConst( 5, 5 ) = tPre * 3.0 / 2.0;
        }

        // Function to compute the second derivative of the traction PK1*N wrt F
        // based on a routine symbolically generated in MATLAB
        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_symbolic_d1PKNdFdF_2d(
                Matrix< DDRMat >&       ad1PKNdFdF,
                const real&             aLame1,    // lambda
                const real&             aLame2,    // mu
                const Matrix< DDRMat >& aNorm )
        {
            //  unpack the components of the deformation gradient
            real tF11 = this->strain(CM_Function_Type::DEFORMATION_GRADIENT)( 0, 0 );
            real tF12 = this->strain(CM_Function_Type::DEFORMATION_GRADIENT)( 1, 0 );
            real tF21 = this->strain(CM_Function_Type::DEFORMATION_GRADIENT)( 2, 0 );
            real tF22 = this->strain(CM_Function_Type::DEFORMATION_GRADIENT)( 3, 0 );

            // compute combination of material parameters
            real tPre           = aLame1 + 2.0 * aLame2;

            // compute all necessary components of dPdFdF
            real td2P11dF11dF11 = 3.0 * tPre * tF11;
            real td2P11dF11dF12 = tPre * tF12;
            real td2P11dF11dF21 = tPre * tF21;
            real td2P11dF11dF22 = aLame1 * tF22;
            real td2P11dF12dF12 = tPre * tF11;
            real td2P11dF12dF21 = aLame2 * tF22;
            real td2P11dF12dF22 = aLame2 * tF21;
            real td2P11dF21dF21 = tPre * tF11;
            real td2P11dF21dF22 = aLame2 * tF12;
            real td2P11dF22dF22 = aLame1 * tF11;

            real td2P12dF11dF11 = tPre * tF12;
            real td2P12dF11dF12 = tPre * tF11;
            real td2P12dF11dF21 = aLame2 * tF22;
            real td2P12dF11dF22 = aLame2 * tF21;
            real td2P12dF12dF12 = 3.0 * tPre * tF12;
            real td2P12dF12dF21 = aLame1 * tF21;
            real td2P12dF12dF22 = tPre * tF22;
            real td2P12dF21dF21 = aLame1 * tF12;
            real td2P12dF21dF22 = aLame2 * tF11;
            real td2P12dF22dF22 = tPre * tF12;

            real td2P21dF11dF11 = tPre * tF21;
            real td2P21dF11dF12 = aLame2 * tF22;
            real td2P21dF11dF21 = tPre * tF11;
            real td2P21dF11dF22 = aLame2 * tF12;
            real td2P21dF12dF12 = aLame1 * tF21;
            real td2P21dF12dF21 = aLame1 * tF12;
            real td2P21dF12dF22 = aLame2 * tF11;
            real td2P21dF21dF21 = 3.0 * tPre * tF21;
            real td2P21dF21dF22 = tPre * tF22;
            real td2P21dF22dF22 = tPre * tF21;

            real td2P22dF11dF11 = aLame1 * tF22;
            real td2P22dF11dF12 = aLame2 * tF21;
            real td2P22dF11dF21 = aLame2 * tF12;
            real td2P22dF11dF22 = aLame1 * tF11;
            real td2P22dF12dF12 = tPre * tF22;
            real td2P22dF12dF21 = aLame2 * tF11;
            real td2P22dF12dF22 = tPre * tF12;
            real td2P22dF21dF21 = tPre * tF22;
            real td2P22dF21dF22 = tPre * tF21;
            real td2P22dF22dF22 = 3.0 * tPre * tF22;

            // unpack the components of the normal
            real tN1 = aNorm( 0, 0 );
            real tN2 = aNorm( 1, 0 );

            // compute all necessary components of dPNdFdF
            ad1PKNdFdF.set_size( 8, 4 );
            ad1PKNdFdF( 0, 0 ) = td2P11dF11dF11 * tN1 + td2P12dF11dF11 * tN2;
            ad1PKNdFdF( 0, 1 ) = td2P11dF11dF12 * tN1 + td2P12dF11dF12 * tN2;
            ad1PKNdFdF( 0, 2 ) = td2P11dF11dF21 * tN1 + td2P12dF11dF21 * tN2;
            ad1PKNdFdF( 0, 3 ) = td2P11dF11dF22 * tN1 + td2P12dF11dF22 * tN2;

            ad1PKNdFdF( 1, 0 ) = ad1PKNdFdF( 0, 1 );
            ad1PKNdFdF( 1, 1 ) = td2P11dF12dF12 * tN1 + td2P12dF12dF12 * tN2;
            ad1PKNdFdF( 1, 2 ) = td2P11dF12dF21 * tN1 + td2P12dF12dF21 * tN2;
            ad1PKNdFdF( 1, 3 ) = td2P11dF12dF22 * tN1 + td2P12dF12dF22 * tN2;

            ad1PKNdFdF( 2, 0 ) = ad1PKNdFdF( 0, 2 );
            ad1PKNdFdF( 2, 1 ) = ad1PKNdFdF( 1, 2 );
            ad1PKNdFdF( 2, 2 ) = td2P11dF21dF21 * tN1 + td2P12dF21dF21 * tN2;
            ad1PKNdFdF( 2, 3 ) = td2P11dF21dF22 * tN1 + td2P12dF21dF22 * tN2;

            ad1PKNdFdF( 3, 0 ) = ad1PKNdFdF( 0, 3 );
            ad1PKNdFdF( 3, 1 ) = ad1PKNdFdF( 1, 3 );
            ad1PKNdFdF( 3, 2 ) = ad1PKNdFdF( 2, 3 );
            ad1PKNdFdF( 3, 3 ) = td2P11dF22dF22 * tN1 + td2P12dF22dF22 * tN2;

            ad1PKNdFdF( 4, 0 ) = td2P21dF11dF11 * tN1 + td2P22dF11dF11 * tN2;
            ad1PKNdFdF( 4, 1 ) = td2P21dF11dF12 * tN1 + td2P22dF11dF12 * tN2;
            ad1PKNdFdF( 4, 2 ) = td2P21dF11dF21 * tN1 + td2P22dF11dF21 * tN2;
            ad1PKNdFdF( 4, 3 ) = td2P21dF11dF22 * tN1 + td2P22dF11dF22 * tN2;

            ad1PKNdFdF( 5, 0 ) = ad1PKNdFdF( 4, 1 );
            ad1PKNdFdF( 5, 1 ) = td2P21dF12dF12 * tN1 + td2P22dF12dF12 * tN2;
            ad1PKNdFdF( 5, 2 ) = td2P21dF12dF21 * tN1 + td2P22dF12dF21 * tN2;
            ad1PKNdFdF( 5, 3 ) = td2P21dF12dF22 * tN1 + td2P22dF12dF22 * tN2;

            ad1PKNdFdF( 6, 0 ) = ad1PKNdFdF( 4, 2 );
            ad1PKNdFdF( 6, 1 ) = ad1PKNdFdF( 5, 2 );
            ad1PKNdFdF( 6, 2 ) = td2P21dF21dF21 * tN1 + td2P22dF21dF21 * tN2;
            ad1PKNdFdF( 6, 3 ) = td2P21dF21dF22 * tN1 + td2P22dF21dF22 * tN2;

            ad1PKNdFdF( 7, 0 ) = ad1PKNdFdF( 4, 3 );
            ad1PKNdFdF( 7, 1 ) = ad1PKNdFdF( 5, 3 );
            ad1PKNdFdF( 7, 2 ) = ad1PKNdFdF( 6, 3 );
            ad1PKNdFdF( 7, 3 ) = td2P21dF22dF22 * tN1 + td2P22dF22dF22 * tN2;
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_symbolic_d1PKNdFdF_3d(
                Matrix< DDRMat >&       ad1PKNdFdF,
                const real&             aLame1,    // lambda
                const real&             aLame2,    // mu
                const Matrix< DDRMat >& aNorm )
        {
            //  unpack the components of the deformation gradient
            real tF11 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 0, 0 );
            real tF12 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 1, 0 );
            real tF13 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 2, 0 );
            real tF21 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 3, 0 );
            real tF22 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 4, 0 );
            real tF23 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 5, 0 );
            real tF31 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 6, 0 );
            real tF32 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 7, 0 );
            real tF33 = this->strain( CM_Function_Type::DEFORMATION_GRADIENT )( 8, 0 );

            // compute combination of material parameters
            real tPre = aLame1 + 2.0 * aLame2;

            // compute all necessary components of dPdFdF
            real td2P11dF11dF11 = 3.0 * tPre * tF11;
            real td2P11dF11dF12 = tPre * tF12;
            real td2P11dF11dF13 = tPre * tF13;
            real td2P11dF11dF21 = tPre * tF21;
            real td2P11dF11dF22 = aLame1 * tF22;
            real td2P11dF11dF23 = aLame1 * tF23;
            real td2P11dF11dF31 = tPre * tF31;
            real td2P11dF11dF32 = aLame1 * tF32;
            real td2P11dF11dF33 = aLame1 * tF33;
            real td2P11dF12dF12 = tPre * tF11;
            real td2P11dF12dF13 = 0.0;
            real td2P11dF12dF21 = aLame2 * tF22;
            real td2P11dF12dF22 = aLame2 * tF21;
            real td2P11dF12dF23 = 0.0;
            real td2P11dF12dF31 = aLame2 * tF32;
            real td2P11dF12dF32 = aLame2 * tF31;
            real td2P11dF12dF33 = 0.0;
            real td2P11dF13dF13 = tPre * tF11;
            real td2P11dF13dF21 = aLame2 * tF23;
            real td2P11dF13dF22 = 0.0;
            real td2P11dF13dF23 = aLame2 * tF21;
            real td2P11dF13dF31 = aLame2 * tF33;
            real td2P11dF13dF32 = 0.0;
            real td2P11dF13dF33 = aLame2 * tF31;
            real td2P11dF21dF21 = tPre * tF11;
            real td2P11dF21dF22 = aLame2 * tF12;
            real td2P11dF21dF23 = aLame2 * tF13;
            real td2P11dF21dF31 = 0.0;
            real td2P11dF21dF32 = 0.0;
            real td2P11dF21dF33 = 0.0;
            real td2P11dF22dF22 = aLame1 * tF11;
            real td2P11dF22dF23 = 0.0;
            real td2P11dF22dF31 = 0.0;
            real td2P11dF22dF32 = 0.0;
            real td2P11dF22dF33 = 0.0;
            real td2P11dF23dF23 = aLame1 * tF11;
            real td2P11dF23dF31 = 0.0;
            real td2P11dF23dF32 = 0.0;
            real td2P11dF23dF33 = 0.0;
            real td2P11dF31dF31 = tPre * tF11;
            real td2P11dF31dF32 = aLame2 * tF12;
            real td2P11dF31dF33 = aLame2 * tF13;
            real td2P11dF32dF32 = aLame1 * tF11;
            real td2P11dF32dF33 = 0.0;
            real td2P11dF33dF33 = aLame1 * tF11;

            real td2P12dF11dF11 = tPre * tF12;
            real td2P12dF11dF12 = tPre * tF11;
            real td2P12dF11dF13 = 0.0;
            real td2P12dF11dF21 = aLame2 * tF22;
            real td2P12dF11dF22 = aLame2 * tF21;
            real td2P12dF11dF23 = 0.0;
            real td2P12dF11dF31 = aLame2 * tF32;
            real td2P12dF11dF32 = aLame2 * tF31;
            real td2P12dF11dF33 = 0.0;
            real td2P12dF12dF12 = 3.0 * tPre * tF12;
            real td2P12dF12dF13 = tPre * tF13;
            real td2P12dF12dF21 = aLame1 * tF21;
            real td2P12dF12dF22 = tPre * tF22;
            real td2P12dF12dF23 = aLame1 * tF23;
            real td2P12dF12dF31 = aLame1 * tF31;
            real td2P12dF12dF32 = tPre * tF32;
            real td2P12dF12dF33 = aLame1 * tF33;
            real td2P12dF13dF13 = tPre * tF12;
            real td2P12dF13dF21 = 0.0;
            real td2P12dF13dF22 = aLame2 * tF23;
            real td2P12dF13dF23 = aLame2 * tF22;
            real td2P12dF13dF31 = 0.0;
            real td2P12dF13dF32 = aLame2 * tF33;
            real td2P12dF13dF33 = aLame2 * tF32;
            real td2P12dF21dF21 = aLame1 * tF12;
            real td2P12dF21dF22 = aLame2 * tF11;
            real td2P12dF21dF23 = 0.0;
            real td2P12dF21dF31 = 0.0;
            real td2P12dF21dF32 = 0.0;
            real td2P12dF21dF33 = 0.0;
            real td2P12dF22dF22 = tPre * tF12;
            real td2P12dF22dF23 = aLame2 * tF13;
            real td2P12dF22dF31 = 0.0;
            real td2P12dF22dF32 = 0.0;
            real td2P12dF22dF33 = 0.0;
            real td2P12dF23dF23 = aLame1 * tF12;
            real td2P12dF23dF31 = 0.0;
            real td2P12dF23dF32 = 0.0;
            real td2P12dF23dF33 = 0.0;
            real td2P12dF31dF31 = aLame1 * tF12;
            real td2P12dF31dF32 = aLame2 * tF11;
            real td2P12dF31dF33 = 0.0;
            real td2P12dF32dF32 = tPre * tF12;
            real td2P12dF32dF33 = aLame2 * tF13;
            real td2P12dF33dF33 = aLame1 * tF12;

            real td2P13dF11dF11 = tPre * tF13;
            real td2P13dF11dF12 = 0.0;
            real td2P13dF11dF13 = tPre * tF11;
            real td2P13dF11dF21 = aLame2 * tF23;
            real td2P13dF11dF22 = 0.0;
            real td2P13dF11dF23 = aLame2 * tF21;
            real td2P13dF11dF31 = aLame2 * tF33;
            real td2P13dF11dF32 = 0.0;
            real td2P13dF11dF33 = aLame2 * tF31;
            real td2P13dF12dF12 = tPre * tF13;
            real td2P13dF12dF13 = tPre * tF12;
            real td2P13dF12dF21 = 0.0;
            real td2P13dF12dF22 = aLame2 * tF23;
            real td2P13dF12dF23 = aLame2 * tF22;
            real td2P13dF12dF31 = 0.0;
            real td2P13dF12dF32 = aLame2 * tF33;
            real td2P13dF12dF33 = aLame2 * tF32;
            real td2P13dF13dF13 = 3.0 * tPre * tF13;
            real td2P13dF13dF21 = aLame1 * tF21;
            real td2P13dF13dF22 = aLame1 * tF22;
            real td2P13dF13dF23 = tPre * tF23;
            real td2P13dF13dF31 = aLame1 * tF31;
            real td2P13dF13dF32 = aLame1 * tF32;
            real td2P13dF13dF33 = tPre * tF33;
            real td2P13dF21dF21 = aLame1 * tF13;
            real td2P13dF21dF22 = 0.0;
            real td2P13dF21dF23 = aLame2 * tF11;
            real td2P13dF21dF31 = 0.0;
            real td2P13dF21dF32 = 0.0;
            real td2P13dF21dF33 = 0.0;
            real td2P13dF22dF22 = aLame1 * tF13;
            real td2P13dF22dF23 = aLame2 * tF12;
            real td2P13dF22dF31 = 0.0;
            real td2P13dF22dF32 = 0.0;
            real td2P13dF22dF33 = 0.0;
            real td2P13dF23dF23 = tPre * tF13;
            real td2P13dF23dF31 = 0.0;
            real td2P13dF23dF32 = 0.0;
            real td2P13dF23dF33 = 0.0;
            real td2P13dF31dF31 = aLame1 * tF13;
            real td2P13dF31dF32 = 0.0;
            real td2P13dF31dF33 = aLame2 * tF11;
            real td2P13dF32dF32 = aLame1 * tF13;
            real td2P13dF32dF33 = aLame2 * tF12;
            real td2P13dF33dF33 = tPre * tF13;

            real td2P21dF11dF11 = tPre * tF21;
            real td2P21dF11dF12 = aLame2 * tF22;
            real td2P21dF11dF13 = aLame2 * tF23;
            real td2P21dF11dF21 = tPre * tF11;
            real td2P21dF11dF22 = aLame2 * tF12;
            real td2P21dF11dF23 = aLame2 * tF13;
            real td2P21dF11dF31 = 0.0;
            real td2P21dF11dF32 = 0.0;
            real td2P21dF11dF33 = 0.0;
            real td2P21dF12dF12 = aLame1 * tF21;
            real td2P21dF12dF13 = 0.0;
            real td2P21dF12dF21 = aLame1 * tF12;
            real td2P21dF12dF22 = aLame2 * tF11;
            real td2P21dF12dF23 = 0.0;
            real td2P21dF12dF31 = 0.0;
            real td2P21dF12dF32 = 0.0;
            real td2P21dF12dF33 = 0.0;
            real td2P21dF13dF13 = aLame1 * tF21;
            real td2P21dF13dF21 = aLame1 * tF13;
            real td2P21dF13dF22 = 0.0;
            real td2P21dF13dF23 = aLame2 * tF11;
            real td2P21dF13dF31 = 0.0;
            real td2P21dF13dF32 = 0.0;
            real td2P21dF13dF33 = 0.0;
            real td2P21dF21dF21 = 3.0 * tPre * tF21;
            real td2P21dF21dF22 = tPre * tF22;
            real td2P21dF21dF23 = tPre * tF23;
            real td2P21dF21dF31 = tPre * tF31;
            real td2P21dF21dF32 = aLame1 * tF32;
            real td2P21dF21dF33 = aLame1 * tF33;
            real td2P21dF22dF22 = tPre * tF21;
            real td2P21dF22dF23 = 0.0;
            real td2P21dF22dF31 = aLame2 * tF32;
            real td2P21dF22dF32 = aLame2 * tF31;
            real td2P21dF22dF33 = 0.0;
            real td2P21dF23dF23 = tPre * tF21;
            real td2P21dF23dF31 = aLame2 * tF33;
            real td2P21dF23dF32 = 0.0;
            real td2P21dF23dF33 = aLame2 * tF31;
            real td2P21dF31dF31 = tPre * tF21;
            real td2P21dF31dF32 = aLame2 * tF22;
            real td2P21dF31dF33 = aLame2 * tF23;
            real td2P21dF32dF32 = aLame1 * tF21;
            real td2P21dF32dF33 = 0.0;
            real td2P21dF33dF33 = aLame1 * tF21;

            real td2P22dF11dF11 = aLame1 * tF22;
            real td2P22dF11dF12 = aLame2 * tF21;
            real td2P22dF11dF13 = 0.0;
            real td2P22dF11dF21 = aLame2 * tF12;
            real td2P22dF11dF22 = aLame1 * tF11;
            real td2P22dF11dF23 = 0.0;
            real td2P22dF11dF31 = 0.0;
            real td2P22dF11dF32 = 0.0;
            real td2P22dF11dF33 = 0.0;
            real td2P22dF12dF12 = tPre * tF22;
            real td2P22dF12dF13 = aLame2 * tF23;
            real td2P22dF12dF21 = aLame2 * tF11;
            real td2P22dF12dF22 = tPre * tF12;
            real td2P22dF12dF23 = aLame2 * tF13;
            real td2P22dF12dF31 = 0.0;
            real td2P22dF12dF32 = 0.0;
            real td2P22dF12dF33 = 0.0;
            real td2P22dF13dF13 = aLame1 * tF22;
            real td2P22dF13dF21 = 0.0;
            real td2P22dF13dF22 = aLame1 * tF13;
            real td2P22dF13dF23 = aLame2 * tF12;
            real td2P22dF13dF31 = 0.0;
            real td2P22dF13dF32 = 0.0;
            real td2P22dF13dF33 = 0.0;
            real td2P22dF21dF21 = tPre * tF22;
            real td2P22dF21dF22 = tPre * tF21;
            real td2P22dF21dF23 = 0.0;
            real td2P22dF21dF31 = aLame2 * tF32;
            real td2P22dF21dF32 = aLame2 * tF31;
            real td2P22dF21dF33 = 0.0;
            real td2P22dF22dF22 = 3.0 * tPre * tF22;
            real td2P22dF22dF23 = tPre * tF23;
            real td2P22dF22dF31 = aLame1 * tF31;
            real td2P22dF22dF32 = tPre * tF32;
            real td2P22dF22dF33 = aLame1 * tF33;
            real td2P22dF23dF23 = tPre * tF22;
            real td2P22dF23dF31 = 0.0;
            real td2P22dF23dF32 = aLame2 * tF33;
            real td2P22dF23dF33 = aLame2 * tF32;
            real td2P22dF31dF31 = aLame1 * tF22;
            real td2P22dF31dF32 = aLame2 * tF21;
            real td2P22dF31dF33 = 0.0;
            real td2P22dF32dF32 = tPre * tF22;
            real td2P22dF32dF33 = aLame2 * tF23;
            real td2P22dF33dF33 = aLame1 * tF22;

            real td2P23dF11dF11 = aLame1 * tF23;
            real td2P23dF11dF12 = 0.0;
            real td2P23dF11dF13 = aLame2 * tF21;
            real td2P23dF11dF21 = aLame2 * tF13;
            real td2P23dF11dF22 = 0.0;
            real td2P23dF11dF23 = aLame1 * tF11;
            real td2P23dF11dF31 = 0.0;
            real td2P23dF11dF32 = 0.0;
            real td2P23dF11dF33 = 0.0;
            real td2P23dF12dF12 = aLame1 * tF23;
            real td2P23dF12dF13 = aLame2 * tF22;
            real td2P23dF12dF21 = 0.0;
            real td2P23dF12dF22 = aLame2 * tF13;
            real td2P23dF12dF23 = aLame1 * tF12;
            real td2P23dF12dF31 = 0.0;
            real td2P23dF12dF32 = 0.0;
            real td2P23dF12dF33 = 0.0;
            real td2P23dF13dF13 = tPre * tF23;
            real td2P23dF13dF21 = aLame2 * tF11;
            real td2P23dF13dF22 = aLame2 * tF12;
            real td2P23dF13dF23 = tPre * tF13;
            real td2P23dF13dF31 = 0.0;
            real td2P23dF13dF32 = 0.0;
            real td2P23dF13dF33 = 0.0;
            real td2P23dF21dF21 = tPre * tF23;
            real td2P23dF21dF22 = 0.0;
            real td2P23dF21dF23 = tPre * tF21;
            real td2P23dF21dF31 = aLame2 * tF33;
            real td2P23dF21dF32 = 0.0;
            real td2P23dF21dF33 = aLame2 * tF31;
            real td2P23dF22dF22 = tPre * tF23;
            real td2P23dF22dF23 = tPre * tF22;
            real td2P23dF22dF31 = 0.0;
            real td2P23dF22dF32 = aLame2 * tF33;
            real td2P23dF22dF33 = aLame2 * tF32;
            real td2P23dF23dF23 = 3.0 * tPre * tF23;
            real td2P23dF23dF31 = aLame1 * tF31;
            real td2P23dF23dF32 = aLame1 * tF32;
            real td2P23dF23dF33 = tPre * tF33;
            real td2P23dF31dF31 = aLame1 * tF23;
            real td2P23dF31dF32 = 0.0;
            real td2P23dF31dF33 = aLame2 * tF21;
            real td2P23dF32dF32 = aLame1 * tF23;
            real td2P23dF32dF33 = aLame2 * tF22;
            real td2P23dF33dF33 = tPre * tF23;

            real td2P31dF11dF11 = tPre * tF31;
            real td2P31dF11dF12 = aLame2 * tF32;
            real td2P31dF11dF13 = aLame2 * tF33;
            real td2P31dF11dF21 = 0.0;
            real td2P31dF11dF22 = 0.0;
            real td2P31dF11dF23 = 0.0;
            real td2P31dF11dF31 = tPre * tF11;
            real td2P31dF11dF32 = aLame2 * tF12;
            real td2P31dF11dF33 = aLame2 * tF13;
            real td2P31dF12dF12 = aLame1 * tF31;
            real td2P31dF12dF13 = 0.0;
            real td2P31dF12dF21 = 0.0;
            real td2P31dF12dF22 = 0.0;
            real td2P31dF12dF23 = 0.0;
            real td2P31dF12dF31 = aLame1 * tF12;
            real td2P31dF12dF32 = aLame2 * tF11;
            real td2P31dF12dF33 = 0.0;
            real td2P31dF13dF13 = aLame1 * tF31;
            real td2P31dF13dF21 = 0.0;
            real td2P31dF13dF22 = 0.0;
            real td2P31dF13dF23 = 0.0;
            real td2P31dF13dF31 = aLame1 * tF13;
            real td2P31dF13dF32 = 0.0;
            real td2P31dF13dF33 = aLame2 * tF11;
            real td2P31dF21dF21 = tPre * tF31;
            real td2P31dF21dF22 = aLame2 * tF32;
            real td2P31dF21dF23 = aLame2 * tF33;
            real td2P31dF21dF31 = tPre * tF21;
            real td2P31dF21dF32 = aLame2 * tF22;
            real td2P31dF21dF33 = aLame2 * tF23;
            real td2P31dF22dF22 = aLame1 * tF31;
            real td2P31dF22dF23 = 0.0;
            real td2P31dF22dF31 = aLame1 * tF22;
            real td2P31dF22dF32 = aLame2 * tF21;
            real td2P31dF22dF33 = 0.0;
            real td2P31dF23dF23 = aLame1 * tF31;
            real td2P31dF23dF31 = aLame1 * tF23;
            real td2P31dF23dF32 = 0.0;
            real td2P31dF23dF33 = aLame2 * tF21;
            real td2P31dF31dF31 = 3.0 * tPre * tF31;
            real td2P31dF31dF32 = tPre * tF32;
            real td2P31dF31dF33 = tPre * tF33;
            real td2P31dF32dF32 = tPre * tF31;
            real td2P31dF32dF33 = 0.0;
            real td2P31dF33dF33 = tPre * tF31;

            real td2P32dF11dF11 = aLame1 * tF32;
            real td2P32dF11dF12 = aLame2 * tF31;
            real td2P32dF11dF13 = 0.0;
            real td2P32dF11dF21 = 0.0;
            real td2P32dF11dF22 = 0.0;
            real td2P32dF11dF23 = 0.0;
            real td2P32dF11dF31 = aLame2 * tF12;
            real td2P32dF11dF32 = aLame1 * tF11;
            real td2P32dF11dF33 = 0.0;
            real td2P32dF12dF12 = tPre * tF32;
            real td2P32dF12dF13 = aLame2 * tF33;
            real td2P32dF12dF21 = 0.0;
            real td2P32dF12dF22 = 0.0;
            real td2P32dF12dF23 = 0.0;
            real td2P32dF12dF31 = aLame2 * tF11;
            real td2P32dF12dF32 = tPre * tF12;
            real td2P32dF12dF33 = aLame2 * tF13;
            real td2P32dF13dF13 = aLame1 * tF32;
            real td2P32dF13dF21 = 0.0;
            real td2P32dF13dF22 = 0.0;
            real td2P32dF13dF23 = 0.0;
            real td2P32dF13dF31 = 0.0;
            real td2P32dF13dF32 = aLame1 * tF13;
            real td2P32dF13dF33 = aLame2 * tF12;
            real td2P32dF21dF21 = aLame1 * tF32;
            real td2P32dF21dF22 = aLame2 * tF31;
            real td2P32dF21dF23 = 0.0;
            real td2P32dF21dF31 = aLame2 * tF22;
            real td2P32dF21dF32 = aLame1 * tF21;
            real td2P32dF21dF33 = 0.0;
            real td2P32dF22dF22 = tPre * tF32;
            real td2P32dF22dF23 = aLame2 * tF33;
            real td2P32dF22dF31 = aLame2 * tF21;
            real td2P32dF22dF32 = tPre * tF22;
            real td2P32dF22dF33 = aLame2 * tF23;
            real td2P32dF23dF23 = aLame1 * tF32;
            real td2P32dF23dF31 = 0.0;
            real td2P32dF23dF32 = aLame1 * tF23;
            real td2P32dF23dF33 = aLame2 * tF22;
            real td2P32dF31dF31 = tPre * tF32;
            real td2P32dF31dF32 = tPre * tF31;
            real td2P32dF31dF33 = 0.0;
            real td2P32dF32dF32 = 3.0 * tPre * tF32;
            real td2P32dF32dF33 = tPre * tF33;
            real td2P32dF33dF33 = tPre * tF32;

            real td2P33dF11dF11 = aLame1 * tF33;
            real td2P33dF11dF12 = 0.0;
            real td2P33dF11dF13 = aLame2 * tF31;
            real td2P33dF11dF21 = 0.0;
            real td2P33dF11dF22 = 0.0;
            real td2P33dF11dF23 = 0.0;
            real td2P33dF11dF31 = aLame2 * tF13;
            real td2P33dF11dF32 = 0.0;
            real td2P33dF11dF33 = aLame1 * tF11;
            real td2P33dF12dF12 = aLame1 * tF33;
            real td2P33dF12dF13 = aLame2 * tF32;
            real td2P33dF12dF21 = 0.0;
            real td2P33dF12dF22 = 0.0;
            real td2P33dF12dF23 = 0.0;
            real td2P33dF12dF31 = 0.0;
            real td2P33dF12dF32 = aLame2 * tF13;
            real td2P33dF12dF33 = aLame1 * tF12;
            real td2P33dF13dF13 = tPre * tF33;
            real td2P33dF13dF21 = 0.0;
            real td2P33dF13dF22 = 0.0;
            real td2P33dF13dF23 = 0.0;
            real td2P33dF13dF31 = aLame2 * tF11;
            real td2P33dF13dF32 = aLame2 * tF12;
            real td2P33dF13dF33 = tPre * tF13;
            real td2P33dF21dF21 = aLame1 * tF33;
            real td2P33dF21dF22 = 0.0;
            real td2P33dF21dF23 = aLame2 * tF31;
            real td2P33dF21dF31 = aLame2 * tF23;
            real td2P33dF21dF32 = 0.0;
            real td2P33dF21dF33 = aLame1 * tF21;
            real td2P33dF22dF22 = aLame1 * tF33;
            real td2P33dF22dF23 = aLame2 * tF32;
            real td2P33dF22dF31 = 0.0;
            real td2P33dF22dF32 = aLame2 * tF23;
            real td2P33dF22dF33 = aLame1 * tF22;
            real td2P33dF23dF23 = tPre * tF33;
            real td2P33dF23dF31 = aLame2 * tF21;
            real td2P33dF23dF32 = aLame2 * tF22;
            real td2P33dF23dF33 = tPre * tF23;
            real td2P33dF31dF31 = tPre * tF33;
            real td2P33dF31dF32 = 0.0;
            real td2P33dF31dF33 = tPre * tF31;
            real td2P33dF32dF32 = tPre * tF33;
            real td2P33dF32dF33 = tPre * tF32;
            real td2P33dF33dF33 = 3.0 * tPre * tF33;

            // unpack the components of the normal
            real tN1 = aNorm( 0, 0 );
            real tN2 = aNorm( 1, 0 );
            real tN3 = aNorm( 2, 0 );

            // compute all necessary components of dPNdFdF
            ad1PKNdFdF.set_size( 27, 9 );
            ad1PKNdFdF( 0, 0 ) = td2P11dF11dF11 * tN1 + td2P12dF11dF11 * tN2 + td2P13dF11dF11 * tN3;
            ad1PKNdFdF( 0, 1 ) = td2P11dF11dF12 * tN1 + td2P12dF11dF12 * tN2 + td2P13dF11dF12 * tN3;
            ad1PKNdFdF( 0, 2 ) = td2P11dF11dF13 * tN1 + td2P12dF11dF13 * tN2 + td2P13dF11dF13 * tN3;
            ad1PKNdFdF( 0, 3 ) = td2P11dF11dF21 * tN1 + td2P12dF11dF21 * tN2 + td2P13dF11dF21 * tN3;
            ad1PKNdFdF( 0, 4 ) = td2P11dF11dF22 * tN1 + td2P12dF11dF22 * tN2 + td2P13dF11dF22 * tN3;
            ad1PKNdFdF( 0, 5 ) = td2P11dF11dF23 * tN1 + td2P12dF11dF23 * tN2 + td2P13dF11dF23 * tN3;
            ad1PKNdFdF( 0, 6 ) = td2P11dF11dF31 * tN1 + td2P12dF11dF31 * tN2 + td2P13dF11dF31 * tN3;
            ad1PKNdFdF( 0, 7 ) = td2P11dF11dF32 * tN1 + td2P12dF11dF32 * tN2 + td2P13dF11dF32 * tN3;
            ad1PKNdFdF( 0, 8 ) = td2P11dF11dF33 * tN1 + td2P12dF11dF33 * tN2 + td2P13dF11dF33 * tN3;
            ad1PKNdFdF( 1, 0 ) = ad1PKNdFdF( 0, 1 );
            ad1PKNdFdF( 1, 1 ) = td2P11dF12dF12 * tN1 + td2P12dF12dF12 * tN2 + td2P13dF12dF12 * tN3;
            ad1PKNdFdF( 1, 2 ) = td2P11dF12dF13 * tN1 + td2P12dF12dF13 * tN2 + td2P13dF12dF13 * tN3;
            ad1PKNdFdF( 1, 3 ) = td2P11dF12dF21 * tN1 + td2P12dF12dF21 * tN2 + td2P13dF12dF21 * tN3;
            ad1PKNdFdF( 1, 4 ) = td2P11dF12dF22 * tN1 + td2P12dF12dF22 * tN2 + td2P13dF12dF22 * tN3;
            ad1PKNdFdF( 1, 5 ) = td2P11dF12dF23 * tN1 + td2P12dF12dF23 * tN2 + td2P13dF12dF23 * tN3;
            ad1PKNdFdF( 1, 6 ) = td2P11dF12dF31 * tN1 + td2P12dF12dF31 * tN2 + td2P13dF12dF31 * tN3;
            ad1PKNdFdF( 1, 7 ) = td2P11dF12dF32 * tN1 + td2P12dF12dF32 * tN2 + td2P13dF12dF32 * tN3;
            ad1PKNdFdF( 1, 8 ) = td2P11dF12dF33 * tN1 + td2P12dF12dF33 * tN2 + td2P13dF12dF33 * tN3;
            ad1PKNdFdF( 2, 0 ) = ad1PKNdFdF( 0, 2 );
            ad1PKNdFdF( 2, 1 ) = ad1PKNdFdF( 1, 2 );
            ad1PKNdFdF( 2, 2 ) = td2P11dF13dF13 * tN1 + td2P12dF13dF13 * tN2 + td2P13dF13dF13 * tN3;
            ad1PKNdFdF( 2, 3 ) = td2P11dF13dF21 * tN1 + td2P12dF13dF21 * tN2 + td2P13dF13dF21 * tN3;
            ad1PKNdFdF( 2, 4 ) = td2P11dF13dF22 * tN1 + td2P12dF13dF22 * tN2 + td2P13dF13dF22 * tN3;
            ad1PKNdFdF( 2, 5 ) = td2P11dF13dF23 * tN1 + td2P12dF13dF23 * tN2 + td2P13dF13dF23 * tN3;
            ad1PKNdFdF( 2, 6 ) = td2P11dF13dF31 * tN1 + td2P12dF13dF31 * tN2 + td2P13dF13dF31 * tN3;
            ad1PKNdFdF( 2, 7 ) = td2P11dF13dF32 * tN1 + td2P12dF13dF32 * tN2 + td2P13dF13dF32 * tN3;
            ad1PKNdFdF( 2, 8 ) = td2P11dF13dF33 * tN1 + td2P12dF13dF33 * tN2 + td2P13dF13dF33 * tN3;
            ad1PKNdFdF( 3, 0 ) = ad1PKNdFdF( 0, 3 );
            ad1PKNdFdF( 3, 1 ) = ad1PKNdFdF( 1, 3 );
            ad1PKNdFdF( 3, 2 ) = ad1PKNdFdF( 2, 3 );
            ad1PKNdFdF( 3, 3 ) = td2P11dF21dF21 * tN1 + td2P12dF21dF21 * tN2 + td2P13dF21dF21 * tN3;
            ad1PKNdFdF( 3, 4 ) = td2P11dF21dF22 * tN1 + td2P12dF21dF22 * tN2 + td2P13dF21dF22 * tN3;
            ad1PKNdFdF( 3, 5 ) = td2P11dF21dF23 * tN1 + td2P12dF21dF23 * tN2 + td2P13dF21dF23 * tN3;
            ad1PKNdFdF( 3, 6 ) = td2P11dF21dF31 * tN1 + td2P12dF21dF31 * tN2 + td2P13dF21dF31 * tN3;
            ad1PKNdFdF( 3, 7 ) = td2P11dF21dF32 * tN1 + td2P12dF21dF32 * tN2 + td2P13dF21dF32 * tN3;
            ad1PKNdFdF( 3, 8 ) = td2P11dF21dF33 * tN1 + td2P12dF21dF33 * tN2 + td2P13dF21dF33 * tN3;
            ad1PKNdFdF( 4, 0 ) = ad1PKNdFdF( 0, 4 );
            ad1PKNdFdF( 4, 1 ) = ad1PKNdFdF( 1, 4 );
            ad1PKNdFdF( 4, 2 ) = ad1PKNdFdF( 2, 4 );
            ad1PKNdFdF( 4, 3 ) = ad1PKNdFdF( 3, 4 );
            ad1PKNdFdF( 4, 4 ) = td2P11dF22dF22 * tN1 + td2P12dF22dF22 * tN2 + td2P13dF22dF22 * tN3;
            ad1PKNdFdF( 4, 5 ) = td2P11dF22dF23 * tN1 + td2P12dF22dF23 * tN2 + td2P13dF22dF23 * tN3;
            ad1PKNdFdF( 4, 6 ) = td2P11dF22dF31 * tN1 + td2P12dF22dF31 * tN2 + td2P13dF22dF31 * tN3;
            ad1PKNdFdF( 4, 7 ) = td2P11dF22dF32 * tN1 + td2P12dF22dF32 * tN2 + td2P13dF22dF32 * tN3;
            ad1PKNdFdF( 4, 8 ) = td2P11dF22dF33 * tN1 + td2P12dF22dF33 * tN2 + td2P13dF22dF33 * tN3;
            ad1PKNdFdF( 5, 0 ) = ad1PKNdFdF( 0, 5 );
            ad1PKNdFdF( 5, 1 ) = ad1PKNdFdF( 1, 5 );
            ad1PKNdFdF( 5, 2 ) = ad1PKNdFdF( 2, 5 );
            ad1PKNdFdF( 5, 3 ) = ad1PKNdFdF( 3, 5 );
            ad1PKNdFdF( 5, 4 ) = ad1PKNdFdF( 4, 5 );
            ad1PKNdFdF( 5, 5 ) = td2P11dF23dF23 * tN1 + td2P12dF23dF23 * tN2 + td2P13dF23dF23 * tN3;
            ad1PKNdFdF( 5, 6 ) = td2P11dF23dF31 * tN1 + td2P12dF23dF31 * tN2 + td2P13dF23dF31 * tN3;
            ad1PKNdFdF( 5, 7 ) = td2P11dF23dF32 * tN1 + td2P12dF23dF32 * tN2 + td2P13dF23dF32 * tN3;
            ad1PKNdFdF( 5, 8 ) = td2P11dF23dF33 * tN1 + td2P12dF23dF33 * tN2 + td2P13dF23dF33 * tN3;
            ad1PKNdFdF( 6, 0 ) = ad1PKNdFdF( 0, 6 );
            ad1PKNdFdF( 6, 1 ) = ad1PKNdFdF( 1, 6 );
            ad1PKNdFdF( 6, 2 ) = ad1PKNdFdF( 2, 6 );
            ad1PKNdFdF( 6, 3 ) = ad1PKNdFdF( 3, 6 );
            ad1PKNdFdF( 6, 4 ) = ad1PKNdFdF( 4, 6 );
            ad1PKNdFdF( 6, 5 ) = ad1PKNdFdF( 5, 6 );
            ad1PKNdFdF( 6, 6 ) = td2P11dF31dF31 * tN1 + td2P12dF31dF31 * tN2 + td2P13dF31dF31 * tN3;
            ad1PKNdFdF( 6, 7 ) = td2P11dF31dF32 * tN1 + td2P12dF31dF32 * tN2 + td2P13dF31dF32 * tN3;
            ad1PKNdFdF( 6, 8 ) = td2P11dF31dF33 * tN1 + td2P12dF31dF33 * tN2 + td2P13dF31dF33 * tN3;
            ad1PKNdFdF( 7, 0 ) = ad1PKNdFdF( 0, 7 );
            ad1PKNdFdF( 7, 1 ) = ad1PKNdFdF( 1, 7 );
            ad1PKNdFdF( 7, 2 ) = ad1PKNdFdF( 2, 7 );
            ad1PKNdFdF( 7, 3 ) = ad1PKNdFdF( 3, 7 );
            ad1PKNdFdF( 7, 4 ) = ad1PKNdFdF( 4, 7 );
            ad1PKNdFdF( 7, 5 ) = ad1PKNdFdF( 5, 7 );
            ad1PKNdFdF( 7, 6 ) = ad1PKNdFdF( 6, 7 );
            ad1PKNdFdF( 7, 7 ) = td2P11dF32dF32 * tN1 + td2P12dF32dF32 * tN2 + td2P13dF32dF32 * tN3;
            ad1PKNdFdF( 7, 8 ) = td2P11dF32dF33 * tN1 + td2P12dF32dF33 * tN2 + td2P13dF32dF33 * tN3;
            ad1PKNdFdF( 8, 0 ) = ad1PKNdFdF( 0, 8 );
            ad1PKNdFdF( 8, 1 ) = ad1PKNdFdF( 1, 8 );
            ad1PKNdFdF( 8, 2 ) = ad1PKNdFdF( 2, 8 );
            ad1PKNdFdF( 8, 3 ) = ad1PKNdFdF( 3, 8 );
            ad1PKNdFdF( 8, 4 ) = ad1PKNdFdF( 4, 8 );
            ad1PKNdFdF( 8, 5 ) = ad1PKNdFdF( 5, 8 );
            ad1PKNdFdF( 8, 6 ) = ad1PKNdFdF( 6, 8 );
            ad1PKNdFdF( 8, 7 ) = ad1PKNdFdF( 7, 8 );
            ad1PKNdFdF( 8, 8 ) = td2P11dF33dF33 * tN1 + td2P12dF33dF33 * tN2 + td2P13dF33dF33 * tN3;

            ad1PKNdFdF( 9, 0 )  = td2P21dF11dF11 * tN1 + td2P22dF11dF11 * tN2 + td2P23dF11dF11 * tN3;
            ad1PKNdFdF( 9, 1 )  = td2P21dF11dF12 * tN1 + td2P22dF11dF12 * tN2 + td2P23dF11dF12 * tN3;
            ad1PKNdFdF( 9, 2 )  = td2P21dF11dF13 * tN1 + td2P22dF11dF13 * tN2 + td2P23dF11dF13 * tN3;
            ad1PKNdFdF( 9, 3 )  = td2P21dF11dF21 * tN1 + td2P22dF11dF21 * tN2 + td2P23dF11dF21 * tN3;
            ad1PKNdFdF( 9, 4 )  = td2P21dF11dF22 * tN1 + td2P22dF11dF22 * tN2 + td2P23dF11dF22 * tN3;
            ad1PKNdFdF( 9, 5 )  = td2P21dF11dF23 * tN1 + td2P22dF11dF23 * tN2 + td2P23dF11dF23 * tN3;
            ad1PKNdFdF( 9, 6 )  = td2P21dF11dF31 * tN1 + td2P22dF11dF31 * tN2 + td2P23dF11dF31 * tN3;
            ad1PKNdFdF( 9, 7 )  = td2P21dF11dF32 * tN1 + td2P22dF11dF32 * tN2 + td2P23dF11dF32 * tN3;
            ad1PKNdFdF( 9, 8 )  = td2P21dF11dF33 * tN1 + td2P22dF11dF33 * tN2 + td2P23dF11dF33 * tN3;
            ad1PKNdFdF( 10, 0 ) = ad1PKNdFdF( 9, 1 );
            ad1PKNdFdF( 10, 1 ) = td2P21dF12dF12 * tN1 + td2P22dF12dF12 * tN2 + td2P23dF12dF12 * tN3;
            ad1PKNdFdF( 10, 2 ) = td2P21dF12dF13 * tN1 + td2P22dF12dF13 * tN2 + td2P23dF12dF13 * tN3;
            ad1PKNdFdF( 10, 3 ) = td2P21dF12dF21 * tN1 + td2P22dF12dF21 * tN2 + td2P23dF12dF21 * tN3;
            ad1PKNdFdF( 10, 4 ) = td2P21dF12dF22 * tN1 + td2P22dF12dF22 * tN2 + td2P23dF12dF22 * tN3;
            ad1PKNdFdF( 10, 5 ) = td2P21dF12dF23 * tN1 + td2P22dF12dF23 * tN2 + td2P23dF12dF23 * tN3;
            ad1PKNdFdF( 10, 6 ) = td2P21dF12dF31 * tN1 + td2P22dF12dF31 * tN2 + td2P23dF12dF31 * tN3;
            ad1PKNdFdF( 10, 7 ) = td2P21dF12dF32 * tN1 + td2P22dF12dF32 * tN2 + td2P23dF12dF32 * tN3;
            ad1PKNdFdF( 10, 8 ) = td2P21dF12dF33 * tN1 + td2P22dF12dF33 * tN2 + td2P23dF12dF33 * tN3;
            ad1PKNdFdF( 11, 0 ) = ad1PKNdFdF( 9, 2 );
            ad1PKNdFdF( 11, 1 ) = ad1PKNdFdF( 10, 2 );
            ad1PKNdFdF( 11, 2 ) = td2P21dF13dF13 * tN1 + td2P22dF13dF13 * tN2 + td2P23dF13dF13 * tN3;
            ad1PKNdFdF( 11, 3 ) = td2P21dF13dF21 * tN1 + td2P22dF13dF21 * tN2 + td2P23dF13dF21 * tN3;
            ad1PKNdFdF( 11, 4 ) = td2P21dF13dF22 * tN1 + td2P22dF13dF22 * tN2 + td2P23dF13dF22 * tN3;
            ad1PKNdFdF( 11, 5 ) = td2P21dF13dF23 * tN1 + td2P22dF13dF23 * tN2 + td2P23dF13dF23 * tN3;
            ad1PKNdFdF( 11, 6 ) = td2P21dF13dF31 * tN1 + td2P22dF13dF31 * tN2 + td2P23dF13dF31 * tN3;
            ad1PKNdFdF( 11, 7 ) = td2P21dF13dF32 * tN1 + td2P22dF13dF32 * tN2 + td2P23dF13dF32 * tN3;
            ad1PKNdFdF( 11, 8 ) = td2P21dF13dF33 * tN1 + td2P22dF13dF33 * tN2 + td2P23dF13dF33 * tN3;
            ad1PKNdFdF( 12, 0 ) = ad1PKNdFdF( 9, 3 );
            ad1PKNdFdF( 12, 1 ) = ad1PKNdFdF( 10, 3 );
            ad1PKNdFdF( 12, 2 ) = ad1PKNdFdF( 11, 3 );
            ad1PKNdFdF( 12, 3 ) = td2P21dF21dF21 * tN1 + td2P22dF21dF21 * tN2 + td2P23dF21dF21 * tN3;
            ad1PKNdFdF( 12, 4 ) = td2P21dF21dF22 * tN1 + td2P22dF21dF22 * tN2 + td2P23dF21dF22 * tN3;
            ad1PKNdFdF( 12, 5 ) = td2P21dF21dF23 * tN1 + td2P22dF21dF23 * tN2 + td2P23dF21dF23 * tN3;
            ad1PKNdFdF( 12, 6 ) = td2P21dF21dF31 * tN1 + td2P22dF21dF31 * tN2 + td2P23dF21dF31 * tN3;
            ad1PKNdFdF( 12, 7 ) = td2P21dF21dF32 * tN1 + td2P22dF21dF32 * tN2 + td2P23dF21dF32 * tN3;
            ad1PKNdFdF( 12, 8 ) = td2P21dF21dF33 * tN1 + td2P22dF21dF33 * tN2 + td2P23dF21dF33 * tN3;
            ad1PKNdFdF( 13, 0 ) = ad1PKNdFdF( 9, 4 );
            ad1PKNdFdF( 13, 1 ) = ad1PKNdFdF( 10, 4 );
            ad1PKNdFdF( 13, 2 ) = ad1PKNdFdF( 11, 4 );
            ad1PKNdFdF( 13, 3 ) = ad1PKNdFdF( 12, 4 );
            ad1PKNdFdF( 13, 4 ) = td2P21dF22dF22 * tN1 + td2P22dF22dF22 * tN2 + td2P23dF22dF22 * tN3;
            ad1PKNdFdF( 13, 5 ) = td2P21dF22dF23 * tN1 + td2P22dF22dF23 * tN2 + td2P23dF22dF23 * tN3;
            ad1PKNdFdF( 13, 6 ) = td2P21dF22dF31 * tN1 + td2P22dF22dF31 * tN2 + td2P23dF22dF31 * tN3;
            ad1PKNdFdF( 13, 7 ) = td2P21dF22dF32 * tN1 + td2P22dF22dF32 * tN2 + td2P23dF22dF32 * tN3;
            ad1PKNdFdF( 13, 8 ) = td2P21dF22dF33 * tN1 + td2P22dF22dF33 * tN2 + td2P23dF22dF33 * tN3;
            ad1PKNdFdF( 14, 0 ) = ad1PKNdFdF( 9, 5 );
            ad1PKNdFdF( 14, 1 ) = ad1PKNdFdF( 10, 5 );
            ad1PKNdFdF( 14, 2 ) = ad1PKNdFdF( 11, 5 );
            ad1PKNdFdF( 14, 3 ) = ad1PKNdFdF( 12, 5 );
            ad1PKNdFdF( 14, 4 ) = ad1PKNdFdF( 13, 5 );
            ad1PKNdFdF( 14, 5 ) = td2P21dF23dF23 * tN1 + td2P22dF23dF23 * tN2 + td2P23dF23dF23 * tN3;
            ad1PKNdFdF( 14, 6 ) = td2P21dF23dF31 * tN1 + td2P22dF23dF31 * tN2 + td2P23dF23dF31 * tN3;
            ad1PKNdFdF( 14, 7 ) = td2P21dF23dF32 * tN1 + td2P22dF23dF32 * tN2 + td2P23dF23dF32 * tN3;
            ad1PKNdFdF( 14, 8 ) = td2P21dF23dF33 * tN1 + td2P22dF23dF33 * tN2 + td2P23dF23dF33 * tN3;
            ad1PKNdFdF( 15, 0 ) = ad1PKNdFdF( 9, 6 );
            ad1PKNdFdF( 15, 1 ) = ad1PKNdFdF( 10, 6 );
            ad1PKNdFdF( 15, 2 ) = ad1PKNdFdF( 11, 6 );
            ad1PKNdFdF( 15, 3 ) = ad1PKNdFdF( 12, 6 );
            ad1PKNdFdF( 15, 4 ) = ad1PKNdFdF( 13, 6 );
            ad1PKNdFdF( 15, 5 ) = ad1PKNdFdF( 14, 6 );
            ad1PKNdFdF( 15, 6 ) = td2P21dF31dF31 * tN1 + td2P22dF31dF31 * tN2 + td2P23dF31dF31 * tN3;
            ad1PKNdFdF( 15, 7 ) = td2P21dF31dF32 * tN1 + td2P22dF31dF32 * tN2 + td2P23dF31dF32 * tN3;
            ad1PKNdFdF( 15, 8 ) = td2P21dF31dF33 * tN1 + td2P22dF31dF33 * tN2 + td2P23dF31dF33 * tN3;
            ad1PKNdFdF( 16, 0 ) = ad1PKNdFdF( 9, 7 );
            ad1PKNdFdF( 16, 1 ) = ad1PKNdFdF( 10, 7 );
            ad1PKNdFdF( 16, 2 ) = ad1PKNdFdF( 11, 7 );
            ad1PKNdFdF( 16, 3 ) = ad1PKNdFdF( 12, 7 );
            ad1PKNdFdF( 16, 4 ) = ad1PKNdFdF( 13, 7 );
            ad1PKNdFdF( 16, 5 ) = ad1PKNdFdF( 14, 7 );
            ad1PKNdFdF( 16, 6 ) = ad1PKNdFdF( 15, 7 );
            ad1PKNdFdF( 16, 7 ) = td2P21dF32dF32 * tN1 + td2P22dF32dF32 * tN2 + td2P23dF32dF32 * tN3;
            ad1PKNdFdF( 16, 8 ) = td2P21dF32dF33 * tN1 + td2P22dF32dF33 * tN2 + td2P23dF32dF33 * tN3;
            ad1PKNdFdF( 17, 0 ) = ad1PKNdFdF( 9, 8 );
            ad1PKNdFdF( 17, 1 ) = ad1PKNdFdF( 10, 8 );
            ad1PKNdFdF( 17, 2 ) = ad1PKNdFdF( 11, 8 );
            ad1PKNdFdF( 17, 3 ) = ad1PKNdFdF( 12, 8 );
            ad1PKNdFdF( 17, 4 ) = ad1PKNdFdF( 13, 8 );
            ad1PKNdFdF( 17, 5 ) = ad1PKNdFdF( 14, 8 );
            ad1PKNdFdF( 17, 6 ) = ad1PKNdFdF( 15, 8 );
            ad1PKNdFdF( 17, 7 ) = ad1PKNdFdF( 16, 8 );
            ad1PKNdFdF( 17, 8 ) = td2P21dF33dF33 * tN1 + td2P22dF33dF33 * tN2 + td2P23dF33dF33 * tN3;

            ad1PKNdFdF( 18, 0 ) = td2P31dF11dF11 * tN1 + td2P32dF11dF11 * tN2 + td2P33dF11dF11 * tN3;
            ad1PKNdFdF( 18, 1 ) = td2P31dF11dF12 * tN1 + td2P32dF11dF12 * tN2 + td2P33dF11dF12 * tN3;
            ad1PKNdFdF( 18, 2 ) = td2P31dF11dF13 * tN1 + td2P32dF11dF13 * tN2 + td2P33dF11dF13 * tN3;
            ad1PKNdFdF( 18, 3 ) = td2P31dF11dF21 * tN1 + td2P32dF11dF21 * tN2 + td2P33dF11dF21 * tN3;
            ad1PKNdFdF( 18, 4 ) = td2P31dF11dF22 * tN1 + td2P32dF11dF22 * tN2 + td2P33dF11dF22 * tN3;
            ad1PKNdFdF( 18, 5 ) = td2P31dF11dF23 * tN1 + td2P32dF11dF23 * tN2 + td2P33dF11dF23 * tN3;
            ad1PKNdFdF( 18, 6 ) = td2P31dF11dF31 * tN1 + td2P32dF11dF31 * tN2 + td2P33dF11dF31 * tN3;
            ad1PKNdFdF( 18, 7 ) = td2P31dF11dF32 * tN1 + td2P32dF11dF32 * tN2 + td2P33dF11dF32 * tN3;
            ad1PKNdFdF( 18, 8 ) = td2P31dF11dF33 * tN1 + td2P32dF11dF33 * tN2 + td2P33dF11dF33 * tN3;
            ad1PKNdFdF( 19, 0 ) = ad1PKNdFdF( 18, 1 );
            ad1PKNdFdF( 19, 1 ) = td2P31dF12dF12 * tN1 + td2P32dF12dF12 * tN2 + td2P33dF12dF12 * tN3;
            ad1PKNdFdF( 19, 2 ) = td2P31dF12dF13 * tN1 + td2P32dF12dF13 * tN2 + td2P33dF12dF13 * tN3;
            ad1PKNdFdF( 19, 3 ) = td2P31dF12dF21 * tN1 + td2P32dF12dF21 * tN2 + td2P33dF12dF21 * tN3;
            ad1PKNdFdF( 19, 4 ) = td2P31dF12dF22 * tN1 + td2P32dF12dF22 * tN2 + td2P33dF12dF22 * tN3;
            ad1PKNdFdF( 19, 5 ) = td2P31dF12dF23 * tN1 + td2P32dF12dF23 * tN2 + td2P33dF12dF23 * tN3;
            ad1PKNdFdF( 19, 6 ) = td2P31dF12dF31 * tN1 + td2P32dF12dF31 * tN2 + td2P33dF12dF31 * tN3;
            ad1PKNdFdF( 19, 7 ) = td2P31dF12dF32 * tN1 + td2P32dF12dF32 * tN2 + td2P33dF12dF32 * tN3;
            ad1PKNdFdF( 19, 8 ) = td2P31dF12dF33 * tN1 + td2P32dF12dF33 * tN2 + td2P33dF12dF33 * tN3;
            ad1PKNdFdF( 20, 0 ) = ad1PKNdFdF( 18, 2 );
            ad1PKNdFdF( 20, 1 ) = ad1PKNdFdF( 19, 2 );
            ad1PKNdFdF( 20, 2 ) = td2P31dF13dF13 * tN1 + td2P32dF13dF13 * tN2 + td2P33dF13dF13 * tN3;
            ad1PKNdFdF( 20, 3 ) = td2P31dF13dF21 * tN1 + td2P32dF13dF21 * tN2 + td2P33dF13dF21 * tN3;
            ad1PKNdFdF( 20, 4 ) = td2P31dF13dF22 * tN1 + td2P32dF13dF22 * tN2 + td2P33dF13dF22 * tN3;
            ad1PKNdFdF( 20, 5 ) = td2P31dF13dF23 * tN1 + td2P32dF13dF23 * tN2 + td2P33dF13dF23 * tN3;
            ad1PKNdFdF( 20, 6 ) = td2P31dF13dF31 * tN1 + td2P32dF13dF31 * tN2 + td2P33dF13dF31 * tN3;
            ad1PKNdFdF( 20, 7 ) = td2P31dF13dF32 * tN1 + td2P32dF13dF32 * tN2 + td2P33dF13dF32 * tN3;
            ad1PKNdFdF( 20, 8 ) = td2P31dF13dF33 * tN1 + td2P32dF13dF33 * tN2 + td2P33dF13dF33 * tN3;
            ad1PKNdFdF( 21, 0 ) = ad1PKNdFdF( 18, 3 );
            ad1PKNdFdF( 21, 1 ) = ad1PKNdFdF( 19, 3 );
            ad1PKNdFdF( 21, 2 ) = ad1PKNdFdF( 20, 3 );
            ad1PKNdFdF( 21, 3 ) = td2P31dF21dF21 * tN1 + td2P32dF21dF21 * tN2 + td2P33dF21dF21 * tN3;
            ad1PKNdFdF( 21, 4 ) = td2P31dF21dF22 * tN1 + td2P32dF21dF22 * tN2 + td2P33dF21dF22 * tN3;
            ad1PKNdFdF( 21, 5 ) = td2P31dF21dF23 * tN1 + td2P32dF21dF23 * tN2 + td2P33dF21dF23 * tN3;
            ad1PKNdFdF( 21, 6 ) = td2P31dF21dF31 * tN1 + td2P32dF21dF31 * tN2 + td2P33dF21dF31 * tN3;
            ad1PKNdFdF( 21, 7 ) = td2P31dF21dF32 * tN1 + td2P32dF21dF32 * tN2 + td2P33dF21dF32 * tN3;
            ad1PKNdFdF( 21, 8 ) = td2P31dF21dF33 * tN1 + td2P32dF21dF33 * tN2 + td2P33dF21dF33 * tN3;
            ad1PKNdFdF( 22, 0 ) = ad1PKNdFdF( 18, 4 );
            ad1PKNdFdF( 22, 1 ) = ad1PKNdFdF( 19, 4 );
            ad1PKNdFdF( 22, 2 ) = ad1PKNdFdF( 20, 4 );
            ad1PKNdFdF( 22, 3 ) = ad1PKNdFdF( 21, 4 );
            ad1PKNdFdF( 22, 4 ) = td2P31dF22dF22 * tN1 + td2P32dF22dF22 * tN2 + td2P33dF22dF22 * tN3;
            ad1PKNdFdF( 22, 5 ) = td2P31dF22dF23 * tN1 + td2P32dF22dF23 * tN2 + td2P33dF22dF23 * tN3;
            ad1PKNdFdF( 22, 6 ) = td2P31dF22dF31 * tN1 + td2P32dF22dF31 * tN2 + td2P33dF22dF31 * tN3;
            ad1PKNdFdF( 22, 7 ) = td2P31dF22dF32 * tN1 + td2P32dF22dF32 * tN2 + td2P33dF22dF32 * tN3;
            ad1PKNdFdF( 22, 8 ) = td2P31dF22dF33 * tN1 + td2P32dF22dF33 * tN2 + td2P33dF22dF33 * tN3;
            ad1PKNdFdF( 23, 0 ) = ad1PKNdFdF( 18, 5 );
            ad1PKNdFdF( 23, 1 ) = ad1PKNdFdF( 19, 5 );
            ad1PKNdFdF( 23, 2 ) = ad1PKNdFdF( 20, 5 );
            ad1PKNdFdF( 23, 3 ) = ad1PKNdFdF( 21, 5 );
            ad1PKNdFdF( 23, 4 ) = ad1PKNdFdF( 22, 5 );
            ad1PKNdFdF( 23, 5 ) = td2P31dF23dF23 * tN1 + td2P32dF23dF23 * tN2 + td2P33dF23dF23 * tN3;
            ad1PKNdFdF( 23, 6 ) = td2P31dF23dF31 * tN1 + td2P32dF23dF31 * tN2 + td2P33dF23dF31 * tN3;
            ad1PKNdFdF( 23, 7 ) = td2P31dF23dF32 * tN1 + td2P32dF23dF32 * tN2 + td2P33dF23dF32 * tN3;
            ad1PKNdFdF( 23, 8 ) = td2P31dF23dF33 * tN1 + td2P32dF23dF33 * tN2 + td2P33dF23dF33 * tN3;
            ad1PKNdFdF( 24, 0 ) = ad1PKNdFdF( 18, 6 );
            ad1PKNdFdF( 24, 1 ) = ad1PKNdFdF( 19, 6 );
            ad1PKNdFdF( 24, 2 ) = ad1PKNdFdF( 20, 6 );
            ad1PKNdFdF( 24, 3 ) = ad1PKNdFdF( 21, 6 );
            ad1PKNdFdF( 24, 4 ) = ad1PKNdFdF( 22, 6 );
            ad1PKNdFdF( 24, 5 ) = ad1PKNdFdF( 23, 6 );
            ad1PKNdFdF( 24, 6 ) = td2P31dF31dF31 * tN1 + td2P32dF31dF31 * tN2 + td2P33dF31dF31 * tN3;
            ad1PKNdFdF( 24, 7 ) = td2P31dF31dF32 * tN1 + td2P32dF31dF32 * tN2 + td2P33dF31dF32 * tN3;
            ad1PKNdFdF( 24, 8 ) = td2P31dF31dF33 * tN1 + td2P32dF31dF33 * tN2 + td2P33dF31dF33 * tN3;
            ad1PKNdFdF( 25, 0 ) = ad1PKNdFdF( 18, 7 );
            ad1PKNdFdF( 25, 1 ) = ad1PKNdFdF( 19, 7 );
            ad1PKNdFdF( 25, 2 ) = ad1PKNdFdF( 20, 7 );
            ad1PKNdFdF( 25, 3 ) = ad1PKNdFdF( 21, 7 );
            ad1PKNdFdF( 25, 4 ) = ad1PKNdFdF( 22, 7 );
            ad1PKNdFdF( 25, 5 ) = ad1PKNdFdF( 23, 7 );
            ad1PKNdFdF( 25, 6 ) = ad1PKNdFdF( 24, 7 );
            ad1PKNdFdF( 25, 7 ) = td2P31dF32dF32 * tN1 + td2P32dF32dF32 * tN2 + td2P33dF32dF32 * tN3;
            ad1PKNdFdF( 25, 8 ) = td2P31dF32dF33 * tN1 + td2P32dF32dF33 * tN2 + td2P33dF32dF33 * tN3;
            ad1PKNdFdF( 26, 0 ) = ad1PKNdFdF( 18, 8 );
            ad1PKNdFdF( 26, 1 ) = ad1PKNdFdF( 19, 8 );
            ad1PKNdFdF( 26, 2 ) = ad1PKNdFdF( 20, 8 );
            ad1PKNdFdF( 26, 3 ) = ad1PKNdFdF( 21, 8 );
            ad1PKNdFdF( 26, 4 ) = ad1PKNdFdF( 22, 8 );
            ad1PKNdFdF( 26, 5 ) = ad1PKNdFdF( 23, 8 );
            ad1PKNdFdF( 26, 6 ) = ad1PKNdFdF( 24, 8 );
            ad1PKNdFdF( 26, 7 ) = ad1PKNdFdF( 25, 8 );
            ad1PKNdFdF( 26, 8 ) = td2P31dF33dF33 * tN1 + td2P32dF33dF33 * tN2 + td2P33dF33dF33 * tN3;
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

