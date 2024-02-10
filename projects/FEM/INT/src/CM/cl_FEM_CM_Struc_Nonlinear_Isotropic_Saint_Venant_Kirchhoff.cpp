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
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "YoungsModulus" ] = static_cast< uint >( CM_Property_Type::EMOD );
            mPropertyMap[ "PoissonRatio" ]  = static_cast< uint >( CM_Property_Type::NU );
        }

        //------------------------------------------------------------------------------

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::set_local_properties()
        {
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

            // FIXME cannot be negative
            // check for small volume change Jacobian
            real tVolumeChangeJacobian = clip_value( this->volume_change_jacobian() );

            // evaluate the Cauchy stress
            Matrix< DDRMat > tCauchyStressFull =
                    this->deformation_gradient() * t2PKStressFull * trans( this->deformation_gradient() ) /
                    tVolumeChangeJacobian;

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
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_traction_second_piola_kirchhoff(
                const Matrix< DDRMat > & aNormal )
        {
            // get second Piola-Kirchhoff in full notation
            Matrix< DDRMat > t2PKStressFull;
            this->voigt_to_full_sym_stress( this->flux( CM_Function_Type::PK2 ), t2PKStressFull );

            // evaluate traction based on 2PK stress
            m2PKTraction = t2PKStressFull * aNormal;
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
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTractiondDOF_second_piola_kirchhoff(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

           // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute test traction wrt dof
            md2PKTractiondu( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes, CM_Function_Type::PK2 );
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
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_testTraction_second_piola_kirchhoff(
                const Matrix< DDRMat >      & aNormal,
                const Vector< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute test traction wrt dof
            m2PKTestTraction( tTestDofIndex ) = this->dTractiondDOF( aTestDofTypes, aNormal, CM_Function_Type::PK2 );
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

            // FIXME Nils please provide implementation!!!

//            // get the Poisson's ratio value
//            const real tNu = mPropPoisson->val()( 0 );
//
//            // get the Youngs' modulus value
//            const real tEmod = mPropEMod->val()( 0 );
//
//        	real tLame1 = tEmod*tNu/((1.0+tNu)*(1.0-2.0*tNu)); 	// lam
//        	real tLame2 = tEmod/(2.0*(1.0+tNu)); 				// mu
//
//            const Matrix< DDRMat > & tDisplGradx =
//                    mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );
//
//            Matrix< DDRMat > projTestTraction;
//            projTestTraction.set_size( 8, mSpaceDim*mSpaceDim, 0.0 );
//
//            Matrix< DDRMat > tJump;
//            tJump.set_size( 8, mSpaceDim*mSpaceDim, 0.0 );
//            tJump = {
//            		{ aJump(0,0), 0, 0, 0 },
//					{ 0, aJump(0,0), 0, 0 },
//					{ 0, 0, aJump(0,0), 0 },
//					{ 0, 0, 0, aJump(0,0) },
//					{ aJump(1,0), 0, 0, 0 },
//					{ 0, aJump(1,0), 0, 0 },
//					{ 0, 0, aJump(1,0), 0 },
//					{ 0, 0, 0, aJump(1,0) },
//            };
//
//            projTestTraction( 0, 0 ) = (tLame1 + 2*tLame2)*(3*aNormal( 0, 0 ) + 3*tDisplGradx( 0, 0 )*aNormal( 0, 0 ) + tDisplGradx( 0, 1 )*aNormal( 1, 0 ));
//            projTestTraction( 0, 1 ) = (tLame1 + 2*tLame2)*(aNormal( 1, 0 ) + tDisplGradx( 0, 0 )*aNormal( 1, 0 ) + tDisplGradx( 0, 1 )*aNormal( 0, 0 ));
//            projTestTraction( 0, 2 ) = aNormal( 0, 0 )*(tDisplGradx( 1, 0 )*tLame1 + 2*tDisplGradx( 1, 0 )*tLame2) + 2*tLame2*aNormal( 1, 0 )*(tDisplGradx( 1, 1 )/ 2.0 + 1.0 / 2.0);
//            projTestTraction( 0, 3 ) = tDisplGradx( 1, 0 )*tLame2*aNormal( 1, 0 ) + tLame1*aNormal( 0, 0 )*(tDisplGradx( 1, 1 ) + 1.0);
//
//            projTestTraction( 1, 0 ) = (tLame1 + 2*tLame2)*(aNormal( 1, 0 ) + tDisplGradx( 0, 0 )*aNormal( 1, 0 ) + tDisplGradx( 0, 1 )*aNormal( 0, 0 ));
//            projTestTraction( 1, 1 ) = (tLame1 + 2*tLame2)*(aNormal( 0, 0 ) + tDisplGradx( 0, 0 )*aNormal( 0, 0 ) + 3*tDisplGradx( 0, 1 )*aNormal( 1, 0 ));
//            projTestTraction( 1, 2 ) = 2*tLame2*aNormal( 0, 0 )*(tDisplGradx( 1, 1 )/ 2.0 + 1.0 / 2.0) + tDisplGradx( 1, 0 )*tLame1*aNormal( 1, 0 );
//            projTestTraction( 1, 3 ) = aNormal( 1, 0 )*(tLame1*(tDisplGradx( 1, 1 ) + 1) + 2*tLame2*(tDisplGradx( 1, 1 ) + 1)) + tDisplGradx( 1, 0 )*tLame2*aNormal( 0, 0 );
//
//            projTestTraction( 2, 0 ) = aNormal( 0, 0 )*(tDisplGradx( 1, 0 )*tLame1 + 2*tDisplGradx( 1, 0 )*tLame2) + 2*tLame2*aNormal( 1, 0 )*(tDisplGradx( 1, 1 )/ 2.0 + 1.0 / 2.0);
//            projTestTraction( 2, 1 ) = 2*tLame2*aNormal( 0, 0 )*(tDisplGradx( 1, 1 )/ 2.0 + 1.0 / 2.0) + tDisplGradx( 1, 0 )*tLame1*aNormal( 1, 0 );
//            projTestTraction( 2, 2 ) = aNormal( 0, 0 )*(tLame1 + 2*tLame2)*(tDisplGradx( 0, 0 ) + 1) + tDisplGradx( 0, 1 )*tLame1*aNormal( 1, 0 );
//            projTestTraction( 2, 3 ) = tDisplGradx( 0, 1 )*tLame2*aNormal( 0, 0 ) + tLame2*aNormal( 1, 0 )*(tDisplGradx( 0, 0 ) + 1);
//
//            projTestTraction( 3, 0 ) = tDisplGradx( 1, 0 )*tLame2*aNormal( 1, 0 ) + tLame1*aNormal( 0, 0 )*(tDisplGradx( 1, 1 ) + 1);
//            projTestTraction( 3, 1 ) = aNormal( 1, 0 )*(tLame1*(tDisplGradx( 1, 1 ) + 1) + 2*tLame2*(tDisplGradx( 1, 1 ) + 1)) + tDisplGradx( 1, 0 )*tLame2*aNormal( 0, 0 );
//            projTestTraction( 3, 2 ) = tDisplGradx( 0, 1 )*tLame2*aNormal( 0, 0 ) + tLame2*aNormal( 1, 0 )*(tDisplGradx( 0, 0 ) + 1);
//            projTestTraction( 3, 3 ) = tDisplGradx( 0, 1 )*aNormal( 1, 0 )*(tLame1 + 2*tLame2) + tLame1*aNormal( 0, 0 )*(tDisplGradx( 0, 0 ) + 1);
//
//            projTestTraction( 4, 0 ) = tDisplGradx( 1, 0 )*aNormal( 0, 0 )*(tLame1 + 2*tLame2) + tLame1*aNormal( 1, 0 )*(tDisplGradx( 1, 1 ) + 1);
//            projTestTraction( 4, 1 ) = tDisplGradx( 1, 0 )*tLame2*aNormal( 1, 0 ) + tLame2*aNormal( 0, 0 )*(tDisplGradx( 1, 1 ) + 1);
//            projTestTraction( 4, 2 ) = aNormal( 0, 0 )*(tLame1*(tDisplGradx( 0, 0 ) + 1) + 2*tLame2*(tDisplGradx( 0, 0 ) + 1)) + tDisplGradx( 0, 1 )*tLame2*aNormal( 1, 0 );
//            projTestTraction( 4, 3 ) = tDisplGradx( 0, 1 )*tLame2*aNormal( 0, 0 ) + tLame1*aNormal( 1, 0 )*(tDisplGradx( 0, 0 ) + 1);
//
//            projTestTraction( 5, 0 ) = tDisplGradx( 1, 0 )*tLame2*aNormal( 1, 0 ) + tLame2*aNormal( 0, 0 )*(tDisplGradx( 1, 1 ) + 1);
//            projTestTraction( 5, 1 ) = aNormal( 1, 0 )*(tLame1 + 2*tLame2)*(tDisplGradx( 1, 1 ) + 1) + tDisplGradx( 1, 0 )*tLame1*aNormal( 0, 0 );
//            projTestTraction( 5, 2 ) = 2*tLame2*aNormal( 1, 0 )*(tDisplGradx( 0, 0 )/ 2.0 + 1.0 / 2.0) + tDisplGradx( 0, 1 )*tLame1*aNormal( 0, 0 );
//            projTestTraction( 5, 3 ) = aNormal( 1, 0 )*(tDisplGradx( 0, 1 )*tLame1 + 2*tDisplGradx( 0, 1 )*tLame2) + 2*tLame2*aNormal( 0, 0 )*(tDisplGradx( 0, 0 )/ 2.0 + 1.0 / 2.0);
//
//            projTestTraction( 6, 0 ) = aNormal( 0, 0 )*(tLame1*(tDisplGradx( 0, 0 ) + 1) + 2*tLame2*(tDisplGradx( 0, 0 ) + 1)) + tDisplGradx( 0, 1 )*tLame2*aNormal( 1, 0 );
//            projTestTraction( 6, 1 ) = 2*tLame2*aNormal( 1, 0 )*(tDisplGradx( 0, 0 )/ 2.0 + 1.0 / 2.0) + tDisplGradx( 0, 1 )*tLame1*aNormal( 0, 0 );
//            projTestTraction( 6, 2 ) = (tLame1 + 2*tLame2)*(aNormal( 1, 0 ) + 3*tDisplGradx( 1, 0 )*aNormal( 0, 0 ) + tDisplGradx( 1, 1 )*aNormal( 1, 0 ));
//            projTestTraction( 6, 3 ) = (tLame1 + 2*tLame2)*(aNormal( 0, 0 ) + tDisplGradx( 1, 0 )*aNormal( 1, 0 ) + tDisplGradx( 1, 1 )*aNormal( 0, 0 ));
//
//            projTestTraction( 7, 0 ) = tDisplGradx( 0, 1 )*tLame2*aNormal( 0, 0 ) + tLame1*aNormal( 1, 0 )*(tDisplGradx( 0, 0 ) + 1);
//            projTestTraction( 7, 1 ) = aNormal( 1, 0 )*(tDisplGradx( 0, 1 )*tLame1 + 2*tDisplGradx( 0, 1 )*tLame2) + 2*tLame2*aNormal( 0, 0 )*(tDisplGradx( 0, 0 )/ 2.0 + 1.0 / 2.0);
//            projTestTraction( 7, 2 ) = (tLame1 + 2*tLame2)*(aNormal( 0, 0 ) + tDisplGradx( 1, 0 )*aNormal( 1, 0 ) + tDisplGradx( 1, 1 )*aNormal( 0, 0 ));
//            projTestTraction( 7, 3 ) = (tLame1 + 2*tLame2)*(3*aNormal( 1, 0 ) + tDisplGradx( 1, 0 )*aNormal( 0, 0 ) + 3*tDisplGradx( 1, 1 )*aNormal( 1, 0 ));
//
//
//
//            md1PKTestTractiondu( tTestDofIndex )( tDofIndex ) = trans( this->dTestStraindDOF( aDofTypes, CM_Function_Type::LAGRANGIAN ) ) * trans(projTestTraction) * tJump
//            		* this->dTestStraindDOF( aDofTypes, CM_Function_Type::LAGRANGIAN );

            MORIS_ASSERT( false, "CM_Struc_Nonlinear_Isotropic::eval_dTestTractiondDOF_first_piola_kirchhoff - Not implemented yet." );
        }

        void
        CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff::eval_dTestTractiondDOF_second_piola_kirchhoff(
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
            md2PKTestTractiondu( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients());

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // FIXME flatten the 2nd derivative of (traction*ajump)
            Matrix< DDRMat > projTestTraction;
            projTestTraction.set_size( mSpaceDim * mSpaceDim, mSpaceDim * mSpaceDim, 0.0 );
            proj_sym( this->constitutive() * trans( tFlatNormal ) * aJump, projTestTraction );

            // FIXME compute the derivative
            md2PKTestTractiondu( tTestDofIndex )( tDofIndex ) =
                    trans( this->dTestStraindDOF( aDofTypes, CM_Function_Type::LAGRANGIAN ) ) *
                    projTestTraction *
                    this->dTestStraindDOF( aDofTypes, CM_Function_Type::LAGRANGIAN );

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                MORIS_LOG_INFO( "CM_Struc_Nonlinear_Isotropic::eval_dTestTractiondDOF_second_piola_kirchhoff - dof dependency of elastic modulus not implemented" );
            }
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
            mConst( 0, 2 ) = tPre * aNu;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 1, 1 ) = tPre * ( 1.0 - aNu );
            mConst( 1, 2 ) = tPre * aNu;
            mConst( 2, 0 ) = tPre * aNu;
            mConst( 2, 1 ) = tPre * aNu;
            mConst( 2, 2 ) = tPre * ( 1.0 - aNu );
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

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

