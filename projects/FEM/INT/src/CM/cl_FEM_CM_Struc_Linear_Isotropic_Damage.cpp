/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_Isotropic_Damage.cpp
 *
 */

#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > aDofTypes,
                Vector< std::string >           aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // if displacement dof type string
                if ( tDofString == "Displacement" )
                {
                    mDofDispl = tDofType;
                }
                // if nonlocal equivalent strain dof type string
                else if ( tDofString == "NonlocalEqStrain" )
                {
                    mDofNonlocalEqStrain = tDofType;
                }
                // if nonlocal equivalent strain history dof type string
                else if ( tDofString == "History" )
                {
                    mDofHistory = tDofType;
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
                            "CM_Struc_Linear_Isotropic_Damage::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::set_parameters(
                Vector< Matrix< DDRMat > > aParameters )
        {
            // set parameters
            mParameters = aParameters;

            // check provided parameters for damage
            uint tParamSize = aParameters.size();
            MORIS_ERROR( tParamSize == 3,
                    "CM_Struc_Linear_Isotropic_Damage::set_parameters - 3 set of parameters need to be provided." );

            // unpack the parameters for evaluation of local equivalent strain
            uint tLEqStrainType = mParameters( 0 )( 0, 0 );
            switch ( tLEqStrainType )
            {
                case 0:
                {
                    // check correct number of parameters provided for the local equivalent strain
                    MORIS_ERROR( mParameters( 0 ).numel() - 1 == 0,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 0 parameter need to be provided." );

                    // set function pointer for local equivalent strain
                    m_eval_equivalent_strain = &CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_LemaitreChaboche;
                    m_eval_dEqStraindu       = &CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_LemaitreChaboche;
                    break;
                }
                case 1:
                {
                    // check correct number of parameters provided for the local equivalent strain
                    MORIS_ERROR( mParameters( 0 ).numel() - 1 == 1,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 1 parameter need to be provided." );

                    // unpack the parameters
                    mLEqStrainParam = mParameters( 0 )( { 0, 0 }, { 1, mParameters( 0 ).numel() - 1 } );
                    mK              = mParameters( 0 )( 1 );

                    // set function pointer for evaluation of local equivalent strain
                    if ( mSpaceDim == 2 && mPlaneType == Model_Type::PLANE_STRESS )
                    {
                        m_eval_equivalent_strain = &CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_2d_plane_stress;
                        m_eval_dEqStraindu       = &CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_stress;
                    }
                    else if ( mSpaceDim == 2 && mPlaneType == Model_Type::PLANE_STRAIN )
                    {
                        m_eval_equivalent_strain = &CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_2d_plane_strain;
                        m_eval_dEqStraindu       = &CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_strain;
                    }
                    else if ( mSpaceDim == 3 )
                    {
                        m_eval_equivalent_strain = &CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_3d;
                        m_eval_dEqStraindu       = &CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_3d;
                    }
                    else
                    {
                        MORIS_ERROR( false, "Unknown equivalent strain implementation based on de Vree et al." );
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Unknown equivalent strain implementation." );
                    break;
                }
            }

            // unpack the parameters for evaluation of damage
            uint tDamageType = mParameters( 1 )( 0, 0 );
            switch ( tDamageType )
            {
                case 0:
                {
                    // check correct number of parameters provided for damage law
                    MORIS_ERROR( mParameters( 1 ).numel() - 1 == 2,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 2 parameters need to be provided." );

                    // unpack the parameters
                    mDamageParam = mParameters( 1 )( { 0, 0 }, { 1, mParameters( 1 ).numel() - 1 } );
                    mKappa0      = mParameters( 1 )( 1 );
                    mAlpha       = mParameters( 1 )( 2 );

                    // set function pointer for evaluation of damage
                    m_eval_damage = &CM_Struc_Linear_Isotropic_Damage::eval_damage_linear;
                    m_eval_dDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_linear;
                    break;
                }
                case 1:
                {
                    // check correct number of parameters provided for the local damage law
                    MORIS_ERROR( mParameters( 1 ).numel() - 1 == 3,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 3 parameters need to be provided." );

                    // unpack the parameters
                    mDamageParam = mParameters( 1 )( { 0, 0 }, { 1, mParameters( 1 ).numel() - 1 } );
                    mKappa0      = mParameters( 1 )( 1 );
                    mAlpha       = mParameters( 1 )( 2 );
                    mBeta        = mParameters( 1 )( 3 );

                    // set function pointer for evaluation of damage
                    m_eval_damage = &CM_Struc_Linear_Isotropic_Damage::eval_damage_exponential;
                    m_eval_dDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_exponential;
                    break;
                }
                case 2:
                {
                    // check correct number of parameters provided for the damage law
                    MORIS_ERROR( mParameters( 1 ).numel() - 1 == 2,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 2 parameters need to be provided." );

                    // unpack the parameters
                    mDamageParam = mParameters( 1 )( { 0, 0 }, { 1, mParameters( 1 ).numel() - 1 } );
                    mKappa0      = mParameters( 1 )( 1 );
                    mAlpha       = mParameters( 1 )( 2 );

                    // set function pointer for evaluation of damage
                    m_eval_damage = &CM_Struc_Linear_Isotropic_Damage::eval_damage_smooth_exponential;
                    m_eval_dDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_smooth_exponential;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Unknown damage implementation." );
                }
            }

            // unpack the parameters for evaluation of smooth damage
            uint tSmoothingType = mParameters( 2 )( 0, 0 );

            switch ( tSmoothingType )
            {
                case 0:
                {
                    // check correct number of parameters provided for smoothing
                    MORIS_ERROR( mParameters( 2 ).numel() - 1 == 0,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 0 parameters need to be provided." );

                    // set function pointer for evaluation of damage
                    m_eval_smooth_damage   = &CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_noSmoothing;
                    m_eval_dSmoothDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_noSmoothing;
                    break;
                }
                case 1:
                {
                    // check correct number of parameters provided for smoothing
                    MORIS_ERROR( mParameters( 2 ).numel() - 1 == 1,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 1 parameter needs to be provided." );

                    // unpack the parameters
                    mSmoothParam = mParameters( 2 )( { 0, 0 }, { 1, mParameters( 2 ).numel() - 1 } );
                    mSmoothC     = mParameters( 2 )( 1 );

                    // set function pointer for smoothing
                    m_eval_smooth_damage = &CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_ks;
                    m_eval_dSmoothDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_ks;
                    break;
                }
                case 2:
                {
                    // check correct number of parameters provided for smoothing
                    MORIS_ERROR( mParameters( 2 ).numel() - 1 == 1,
                            "CM_Struc_Linear_Isotropic_Damage::set_parameters - 1 parameter needs to be provided." );

                    // unpack the parameters
                    mSmoothParam = mParameters( 2 )( { 0, 0 }, { 1, mParameters( 2 ).numel() - 1 } );
                    mSmoothC     = mParameters( 2 )( 1 );

                    // set function pointer for smoothing
                    m_eval_smooth_damage = &CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_corrected_ks;
                    m_eval_dSmoothDamagedu = &CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_corrected_ks;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Unknown smoothing implementation." );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::reset_eval_flags()
        {
            // call parent implementation
            Constitutive_Model::reset_eval_flags();

            // reset child specific eval flags for equivalent strain
            mEqStrainEval = true;
            mdEqStrainduEval.fill( true );

            // reset child specific eval flags for nonlocal equivalent strain history
            mHistoryEval = true;
            mdHistoryduEval.fill( true );

            // reset child specific eval flags for nonlocal equivalent strain history
            mHistoryRefEval = true;
            mdHistoryRefduEval.fill( true );

            // reset child specific eval flags for damage
            mDamageEval = true;
            mdDamageduEval.fill( true );

            // reset child specific eval flags for smoothed damage
            mSmoothDamageEval = true;
            mdSmoothDamageduEval.fill( true );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::build_global_dof_type_list()
        {
            // call parent implementation
            Constitutive_Model::build_global_dof_type_list();

            // get number of dof types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();

            // init child specific eval flags
            mdDamageduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdSmoothDamageduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEqStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdHistoryduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdHistoryRefduEval.set_size( tNumGlobalDofTypes, 1, true );

            // init child specific storage
            mdDamagedu.resize( tNumGlobalDofTypes );
            mdSmoothDamagedu.resize( tNumGlobalDofTypes );
            mdEqStraindu.resize( tNumGlobalDofTypes );
            mdHistorydu.resize( tNumGlobalDofTypes );
            mdHistoryRefdu.resize( tNumGlobalDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::damage(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::damage - Only DEFAULT CM function type known in base class." );

            // if the value was not evaluated
            if ( mDamageEval )
            {
                // compute the value
                this->eval_damage();

                // set bool for evaluation
                mDamageEval = false;
            }
            // return the value
            return mDamage;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::dDamagedu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::dDamagedu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Struc_Linear_Isotropic_Damage::dDamagedu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdDamageduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dDamagedu( aDofType );

                // set bool for evaluation
                mdDamageduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdDamagedu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::smooth_damage(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::smooth_damage - Only DEFAULT CM function type known in base class." );

            // if the value was not evaluated
            if ( mSmoothDamageEval )
            {
                // compute the value
                this->eval_smooth_damage();

                // set bool for evaluation
                mSmoothDamageEval = false;
            }
            // return the value
            return mSmoothDamage;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::dSmoothDamagedu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::dSmoothDamagedu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Struc_Linear_Isotropic_Damage::dSmoothDamagedu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdSmoothDamageduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dSmoothDamagedu( aDofType );

                // set bool for evaluation
                mdSmoothDamageduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdSmoothDamagedu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::equivalent_strain(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::equivalent_strain - Only DEFAULT CM function type known in base class." );

            // if the value was not evaluated
            if ( mEqStrainEval )
            {
                // evaluate the value
                this->eval_equivalent_strain();

                // set bool for evaluation
                mEqStrainEval = false;
            }
            // return the value
            return mEqStrain;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::dEqStraindu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::dEqStraindu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Struc_Linear_Isotropic_Damage::dEqStraindu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdEqStrainduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dEqStraindu( aDofType );

                // set bool for evaluation
                mdEqStrainduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdEqStraindu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::history(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::history - Only DEFAULT CM function type known in base class." );

            // if the value was not evaluated
            if ( mHistoryEval )
            {
                // compute the value
                this->eval_history();

                // set bool for evaluation
                mHistoryEval = false;
            }
            // return the value
            return mHistory;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::dHistorydu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::dHistorydu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Struc_Linear_Isotropic_Damage::dHistorydu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdHistoryduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dHistorydu( aDofType );

                // set bool for evaluation
                mdHistoryduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdHistorydu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::history_ref(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::history - Only DEFAULT CM function type known in base class." );

            // if the value was not evaluated
            if ( mHistoryRefEval )
            {
                // compute the value
                this->eval_history_ref();

                // set bool for evaluation
                mHistoryRefEval = false;
            }
            // return the value
            return mHistoryRef;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::dHistoryRefdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Linear_Isotropic_Damage::dHistoryRefdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Struc_Linear_Isotropic_Damage::dHistoryRefdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdHistoryRefduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dHistoryRefdu( aDofType );

                // set bool for evaluation
                mdHistoryRefduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdHistoryRefdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_flux()
        {
            // call the parent contribution
            CM_Struc_Linear_Isotropic::eval_flux();

            // compute the damage value
            real tDamage = this->smooth_damage()( 0 );

            // compute the flux with damage contribution
            mFlux = ( 1.0 - tDamage ) * mFlux;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // call the parent contribution
            CM_Struc_Linear_Isotropic::eval_dFluxdDOF( aDofTypes );

            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute the damage value
            real tDamage = this->smooth_damage()( 0 );

            // compute dFLuxdDOF with damage contribution
            mdFluxdDof( tDofIndex ) = ( 1.0 - tDamage ) * mdFluxdDof( tDofIndex );

            // compute the damage value
            if ( tDamage < 1.0 )
            {
                // compute contribution of damage derivative
                mdFluxdDof( tDofIndex ) -= this->flux() * this->dSmoothDamagedu( aDofTypes ) / ( 1.0 - tDamage );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&      aNormal )
        {
            // get derivative dof type index
            const uint tDerDofIndex = mDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute derivative of traction wrt dof
            mdTractiondDof( tDerDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_testTraction(
                const Matrix< DDRMat >&      aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute test traction wrt test dof
            mTestTraction( tTestDofIndex ) = tFlatNormal * this->dFluxdDOF( aTestDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&      aNormal,
                const Matrix< DDRMat >&      aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get number of test and derivative dof coefficients
            real tNumTestDofCoeff = mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients();
            real tNumDerivativeDofCoeff = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients();

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( tNumTestDofCoeff, tNumDerivativeDofCoeff );

            // if test dof is displacement and derivative dof is nonlocal equivalent strain
            if ( aTestDofTypes( 0 ) == mDofDispl )
            {
                // compute the damage
                real tDamage = this->smooth_damage()( 0 );

                // compute derivative of test traction
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                        -1.0 * this->testTraction_trans( aNormal, aTestDofTypes ) * aJump * this->dSmoothDamagedu( aDofTypes ) / ( 1.0 - tDamage );
            }
            else
            {
                // set contribution set to zero
                // FIXME check that correct
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain()
        {
            // call function pointer for evaluation
            ( this->*m_eval_equivalent_strain )( mLEqStrainParam );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_LemaitreChaboche(
                const Matrix< DDRMat >& aLEqStrainParam )
        {
            // compute the energy release rate
            Matrix< DDRMat > tY =
                    trans( this->strain() ) * this->constitutive() * this->strain() / mPropEMod->val()( 0 );

            // check that the energy release rate is nonzero
            if ( tY( 0 ) > 0.0 )
            {
                // compute the equivalent strain
                mEqStrain = { { std::pow( tY( 0 ), 0.5 ) } };
            }
            else
            {
                // set equivalent strain to zero
                mEqStrain = { { 0.0 } };
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_2d_plane_stress(
                const Matrix< DDRMat >& aLEqStrainParam )
        {
            // grab Poisson's ratio value
            real tNu = mPropPoisson->val()( 0 );

            // unpack the strain components
            real tEpsilon11 = this->strain()( 0 );
            real tEpsilon22 = this->strain()( 1 );
            real tEpsilon33 = -tNu * ( tEpsilon11 + tEpsilon22 ) / ( 1.0 - tNu );
            real tEpsilon12 = this->strain()( 2 );

            // unpack the square of the strain components
            real tEpsilon1122Square = std::pow( tEpsilon11 - tEpsilon22, 2.0 );
            real tEpsilon2233Square = std::pow( tEpsilon22 - tEpsilon33, 2.0 );
            real tEpsilon3311Square = std::pow( tEpsilon33 - tEpsilon11, 2.0 );
            real tEpsilon12Square   = std::pow( tEpsilon12, 2.0 );

            // compute first invariant
            real tI1 = tEpsilon11 + tEpsilon22 + tEpsilon33;

            // compute second deviatoric invariant
            real tJ2 = ( ( tEpsilon1122Square + tEpsilon2233Square + tEpsilon3311Square ) / 6.0 ) + tEpsilon12Square;

            // compute local equivalent strain
            mEqStrain = ( mK - 1.0 ) * tI1 / ( 2.0 * mK * ( 1 - 2.0 * tNu ) )    //
                      + ( 0.5 / mK ) * std::pow( std::pow( ( mK - 1.0 ) * tI1 / ( 1 - 2.0 * tNu ), 2.0 ) + 12.0 * mK * tJ2 / std::pow( 1 + tNu, 2.0 ), 0.5 );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_2d_plane_strain(
                const Matrix< DDRMat >& aLEqStrainParam )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_2d_plane_strain - case not implemented." );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_equivalent_strain_deVree_3d(
                const Matrix< DDRMat >& aLEqStrainParam )
        {
            // unpack the strain components
            real tEpsilon11 = this->strain()( 0 );
            real tEpsilon22 = this->strain()( 1 );
            real tEpsilon33 = this->strain()( 2 );
            real tEpsilon12 = this->strain()( 5 );
            real tEpsilon23 = this->strain()( 3 );
            real tEpsilon31 = this->strain()( 4 );

            // unpack the square of the strain components
            real tEpsilon1122Square = std::pow( tEpsilon11 - tEpsilon22, 2.0 );
            real tEpsilon2233Square = std::pow( tEpsilon22 - tEpsilon33, 2.0 );
            real tEpsilon3311Square = std::pow( tEpsilon33 - tEpsilon11, 2.0 );
            real tEpsilon12Square = std::pow( tEpsilon12, 2.0 );
            real tEpsilon23Square = std::pow( tEpsilon23, 2.0 );
            real tEpsilon31Square = std::pow( tEpsilon31, 2.0 );

            // compute first invariant
            real tI1 = tEpsilon11 + tEpsilon22 + tEpsilon33;

            // compute second deviatoric invariant
            real tJ2 = ( ( tEpsilon1122Square + tEpsilon2233Square + tEpsilon3311Square ) / 6.0 )    //
                     + tEpsilon12Square + tEpsilon23Square + tEpsilon31Square;

            // grab Poisson's ratio value
            real tNu = mPropPoisson->val()( 0 );

            // compute local equivalent strain
            mEqStrain = ( mK - 1.0 ) * tI1 / ( 2.0 * mK * ( 1 - 2.0 * tNu ) )    //
                      + ( 0.5 / mK ) * std::pow( std::pow( ( mK - 1.0 ) * tI1 / ( 1 - 2.0 * tNu ), 2.0 ) + 12.0 * mK * tJ2 / std::pow( 1 + tNu, 2.0 ), 0.5 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // call function pointer for evaluation
            ( this->*m_eval_dEqStraindu )( aDofTypes );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_LemaitreChaboche(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdEqStraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the equivalent strain value
            real tEqStrain = this->equivalent_strain()( 0 );

            // check that the equivalent strain is nonzero
            if ( tEqStrain > 0.0 )
            {
                // if derivative dof type is nonlocal equivalent strain
                if ( aDofTypes( 0 ) == mDofDispl )
                {
                    mdEqStraindu( tDofIndex ) =
                            2.0 * trans( this->strain() ) * this->constitutive() * this->dStraindDOF( aDofTypes ) / mPropEMod->val()( 0 );
                }
                else
                {
                    mdEqStraindu( tDofIndex ).fill( 0.0 );
                }

                // if Young's modulus property depends on dof type
                if ( mPropEMod->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_LemaitreChaboche - Youngs modulus depends on dof, not handled." );
                }

                // if Poisson ratio property depends on dof type
                if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_LemaitreChaboche - Poisson's ratio depends on dof, not handled." );
                }

                mdEqStraindu( tDofIndex ) = mdEqStraindu( tDofIndex ) / ( 2.0 * tEqStrain );
            }
            else
            {
                mdEqStraindu( tDofIndex ).fill( 0.0 );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_stress(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdEqStraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the equivalent strain value
            real tEqStrain = this->equivalent_strain()( 0 );

            // check that the equivalent strain is nonzero
            if ( tEqStrain > 0.0 )
            {
                // if derivative dof type is nonlocal equivalent strain
                if ( aDofTypes( 0 ) == mDofDispl )
                {
                    // grab Poisson's ratio value
                    real tNu = mPropPoisson->val()( 0 );

                    // unpack the strain components
                    real tEpsilon11 = this->strain()( 0 );
                    real tEpsilon22 = this->strain()( 1 );
                    real tEpsilon33 = -tNu * ( tEpsilon11 + tEpsilon22 ) / ( 1.0 - tNu );
                    real tEpsilon12 = this->strain()( 2 );

                    // unpack the strain components derivatives
                    const Matrix< DDRMat >& tdEpsilon11du = this->dStraindDOF( aDofTypes ).get_row( 0 );
                    const Matrix< DDRMat >& tdEpsilon22du = this->dStraindDOF( aDofTypes ).get_row( 1 );
                    const Matrix< DDRMat >& tdEpsilon33du = -tNu * ( tdEpsilon11du + tdEpsilon22du ) / ( 1.0 - tNu );
                    const Matrix< DDRMat >& tdEpsilon12du = this->dStraindDOF( aDofTypes ).get_row( 2 );

                    // unpack the difference of the strain components
                    real tEpsilon1122 = tEpsilon11 - tEpsilon22;
                    real tEpsilon2233 = tEpsilon22 - tEpsilon33;
                    real tEpsilon3311 = tEpsilon33 - tEpsilon11;

                    // unpack the square of the strain components
                    real tEpsilon1122Square = std::pow( tEpsilon1122, 2.0 );
                    real tEpsilon2233Square = std::pow( tEpsilon2233, 2.0 );
                    real tEpsilon3311Square = std::pow( tEpsilon3311, 2.0 );
                    real tEpsilon12Square   = std::pow( tEpsilon12, 2.0 );

                    // compute first invariant
                    real tI1 = tEpsilon11 + tEpsilon22 + tEpsilon33;

                    // compute second deviatoric invariant
                    real tJ2 = ( ( tEpsilon1122Square + tEpsilon2233Square + tEpsilon3311Square ) / 6.0 ) + tEpsilon12Square;

                    // compute the square ro term
                    real tSquareRootTerm = std::pow(                                   //
                            std::pow( ( mK - 1.0 ) * tI1 / ( 1 - 2.0 * tNu ), 2.0 )    //
                                    + 12.0 * mK * tJ2 / std::pow( 1 + tNu, 2.0 ),
                            0.5 );

                    // compute dLEqStraindI1
                    real dLEqStraindI1 = ( mK - 1.0 ) / ( 2.0 * mK * ( 1 - 2.0 * tNu ) )    //
                                       + tI1 * std::pow( mK - 1.0, 2.0 ) / ( 2.0 * mK * std::pow( 1 - 2.0 * tNu, 2.0 ) * tSquareRootTerm );

                    // compute dLEqStraindJ2
                    real dLEqStraindJ2 = 3.0 / ( std::pow( 1.0 + tNu, 2.0 ) * tSquareRootTerm );

                    // compute dI1du
                    Matrix< DDRMat > dI1du = tdEpsilon11du + tdEpsilon22du + tdEpsilon33du;

                    // compute dJ2dEpsilonij
                    real tdJ2dEpsilon11 = ( tEpsilon1122 - tEpsilon3311 ) / 3.0;
                    real tdJ2dEpsilon22 = ( tEpsilon2233 - tEpsilon1122 ) / 3.0;
                    real tdJ2dEpsilon33 = ( tEpsilon3311 - tEpsilon2233 ) / 3.0;
                    real tdJ2dEpsilon12 = 2.0 * tEpsilon12;

                    // compute dJ2du
                    Matrix< DDRMat > dJ2du =                    //
                            tdJ2dEpsilon11 * tdEpsilon11du      //
                            + tdJ2dEpsilon22 * tdEpsilon22du    //
                            + tdJ2dEpsilon33 * tdEpsilon33du    //
                            + tdJ2dEpsilon12 * tdEpsilon12du;

                    // compute the derivative of the local equivalent strain
                    mdEqStraindu( tDofIndex ) = dLEqStraindI1 * dI1du + dLEqStraindJ2 * dJ2du;
                }
                else
                {
                    mdEqStraindu( tDofIndex ).fill( 0.0 );
                }

                // if Poisson ratio property depends on dof type
                if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_stress - Poisson's ratio depends on dof, not handled." );
                }
            }
            else
            {
                mdEqStraindu( tDofIndex ).fill( 0.0 );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_strain(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_2d_plane_strain - case not implemented." );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_3d(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdEqStraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the equivalent strain value
            real tEqStrain = this->equivalent_strain()( 0 );

            // check that the equivalent strain is nonzero
            if ( tEqStrain > 0.0 )
            {
                // if derivative dof type is nonlocal equivalent strain
                if ( aDofTypes( 0 ) == mDofDispl )
                {
                    // unpack the strain components
                    real tEpsilon11 = this->strain()( 0 );
                    real tEpsilon22 = this->strain()( 1 );
                    real tEpsilon33 = this->strain()( 2 );
                    real tEpsilon23 = this->strain()( 3 );
                    real tEpsilon31 = this->strain()( 4 );
                    real tEpsilon12 = this->strain()( 5 );

                    // unpack the strain components derivatives
                    const Matrix< DDRMat >& tdEpsilon11du = this->dStraindDOF( aDofTypes ).get_row( 0 );
                    const Matrix< DDRMat >& tdEpsilon22du = this->dStraindDOF( aDofTypes ).get_row( 1 );
                    const Matrix< DDRMat >& tdEpsilon33du = this->dStraindDOF( aDofTypes ).get_row( 2 );
                    const Matrix< DDRMat >& tdEpsilon23du = this->dStraindDOF( aDofTypes ).get_row( 3 );
                    const Matrix< DDRMat >& tdEpsilon31du = this->dStraindDOF( aDofTypes ).get_row( 4 );
                    const Matrix< DDRMat >& tdEpsilon12du = this->dStraindDOF( aDofTypes ).get_row( 5 );

                    // unpack the difference of the strain components
                    real tEpsilon1122 = tEpsilon11 - tEpsilon22;
                    real tEpsilon2233 = tEpsilon22 - tEpsilon33;
                    real tEpsilon3311 = tEpsilon33 - tEpsilon11;

                    // unpack the square of the strain components
                    real tEpsilon1122Square = std::pow( tEpsilon1122, 2.0 );
                    real tEpsilon2233Square = std::pow( tEpsilon2233, 2.0 );
                    real tEpsilon3311Square = std::pow( tEpsilon3311, 2.0 );
                    real tEpsilon12Square   = std::pow( tEpsilon12, 2.0 );
                    real tEpsilon23Square   = std::pow( tEpsilon23, 2.0 );
                    real tEpsilon31Square   = std::pow( tEpsilon31, 2.0 );

                    // compute first invariant
                    real tI1 = tEpsilon11 + tEpsilon22 + tEpsilon33;

                    // compute second deviatoric invariant
                    real tJ2 = ( ( tEpsilon1122Square + tEpsilon2233Square + tEpsilon3311Square ) / 6.0 )    //
                             + tEpsilon12Square + tEpsilon23Square + tEpsilon31Square;

                    // grab Poisson's ratio value
                    real tNu = mPropPoisson->val()( 0 );

                    // compute the square ro term
                    real tSquareRootTerm = std::pow(                                   //
                            std::pow( ( mK - 1.0 ) * tI1 / ( 1 - 2.0 * tNu ), 2.0 )    //
                                    + 12.0 * mK * tJ2 / std::pow( 1 + tNu, 2.0 ),
                            0.5 );

                    // compute dLEqStraindI1
                    real dLEqStraindI1 = ( mK - 1.0 ) / ( 2.0 * mK * ( 1 - 2.0 * tNu ) )    //
                                       + tI1 * std::pow( mK - 1.0, 2.0 ) / ( 2.0 * mK * std::pow( 1 - 2.0 * tNu, 2.0 ) * tSquareRootTerm );

                    // compute dLEqStraindJ2
                    real dLEqStraindJ2 = 3.0 / ( std::pow( 1.0 + tNu, 2.0 ) * tSquareRootTerm );

                    // compute dI1du
                    Matrix< DDRMat > dI1du = tdEpsilon11du + tdEpsilon22du + tdEpsilon33du;

                    // compute dJ2dEpsilonij
                    real tdJ2dEpsilon11 = ( tEpsilon1122 - tEpsilon3311 ) / 3.0;
                    real tdJ2dEpsilon22 = ( tEpsilon2233 - tEpsilon1122 ) / 3.0;
                    real tdJ2dEpsilon33 = ( tEpsilon3311 - tEpsilon2233 ) / 3.0;
                    real tdJ2dEpsilon12 = 2.0 * tEpsilon12;
                    real tdJ2dEpsilon23 = 2.0 * tEpsilon23;
                    real tdJ2dEpsilon31 = 2.0 * tEpsilon31;

                    // compute dJ2du
                    Matrix< DDRMat > dJ2du =                    //
                            tdJ2dEpsilon11 * tdEpsilon11du      //
                            + tdJ2dEpsilon22 * tdEpsilon22du    //
                            + tdJ2dEpsilon33 * tdEpsilon33du    //
                            + tdJ2dEpsilon12 * tdEpsilon12du    //
                            + tdJ2dEpsilon23 * tdEpsilon23du    //
                            + tdJ2dEpsilon31 * tdEpsilon31du;

                    // compute the derivative of the local equivalent strain
                    mdEqStraindu( tDofIndex ) = dLEqStraindI1 * dI1du + dLEqStraindJ2 * dJ2du;
                }
                else
                {
                    mdEqStraindu( tDofIndex ).fill( 0.0 );
                }

                // if Poisson ratio property depends on dof type
                if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::eval_dEqStraindu_deVree_3d - Poisson's ratio depends on dof, not handled." );
                }
            }
            else
            {
                mdEqStraindu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_damage()
        {
            // call function pointer for evaluation
            ( this->*m_eval_damage )( mDamageParam );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_damage_linear(
                const Matrix< DDRMat >& aDamageParam )
        {
            // get the field interpolator for nonlocal equivalent strain history
            Field_Interpolator* tFIHistory = mFIManager->get_field_interpolators_for_type( mDofHistory );

            // get the nonlocal equivalent strain history value
            real tHistory = tFIHistory->val()( 0 );

            // compute damage value based on history
            mDamage = ( mAlpha / tHistory ) * ( tHistory - mKappa0 ) / ( mAlpha - mKappa0 );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_damage_exponential(
                const Matrix< DDRMat >& aDamageParam )
        {
            // get the field interpolator for nonlocal equivalent strain history
            Field_Interpolator* tFIHistory = mFIManager->get_field_interpolators_for_type( mDofHistory );

            // get the nonlocal equivalent strain history value
            real tHistory = tFIHistory->val()( 0 );

            // FIXME
            // if nonlocal equivalent strain history larger than threshold
            if ( tHistory > 1.0e-12 )
            {
                // compute damage value based on history
                mDamage = 1.0 - ( mKappa0 / tHistory ) * ( 1.0 - mAlpha + mAlpha * std::exp( -mBeta * ( tHistory - mKappa0 ) ) );
            }
            else
            {
                // set damage to zero
                mDamage = 1.0 - ( mKappa0 / 1.0e-12 ) * ( 1.0 - mAlpha + mAlpha * std::exp( -mBeta * ( 1.0e-12 - mKappa0 ) ) );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_damage_smooth_exponential(
                const Matrix< DDRMat >& aDamageParam )
        {
            // get the field interpolator for nonlocal equivalent strain history
            Field_Interpolator* tFIHistory = mFIManager->get_field_interpolators_for_type( mDofHistory );

            // get the nonlocal equivalent strain history value
            real tHistory = tFIHistory->val()( 0 );

            // compute damage value based on history
            mDamage = std::exp( -std::exp( mAlpha * ( mKappa0 - tHistory ) ) );
        }


        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // call function pointer for evaluation
            ( this->*m_eval_dDamagedu )( aDofTypes );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_linear(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdDamagedu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // if derivative dof is history dof
            if ( aDofTypes( 0 ) == mDofHistory )
            {
                // get the nonlocal equivalent strain history value
                real tHistory = tFI->val()( 0 );

                // get the damage value
                real tDamage = this->damage()( 0 );

                // compute the derivative of the damage value based on history
                mdDamagedu( tDofIndex ) = tDamage * mKappa0 * tFI->N() / ( tHistory * ( tHistory - mKappa0 ) );
            }
            else
            {
                mdDamagedu( tDofIndex ).fill( 0.0 );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_exponential(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdDamagedu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // if derivative dof is history dof
            if ( aDofTypes( 0 ) == mDofHistory )
            {
                // get the nonlocal equivalent strain history value
                real tHistory = tFI->val()( 0 );

                // get the damage value
                real tDamage = this->damage()( 0 );

                // if nonlocal equivalent strain history larger than threshold
                if ( tHistory > 1.0e-12 )
                {

                    // compute the derivative of the damage value based on history
                    mdDamagedu( tDofIndex ) =                                                                    //
                            ( 1.0 - tDamage                                                                      //
                                    + mKappa0 * mAlpha * mBeta * std::exp( mBeta * ( mKappa0 - tHistory ) ) )    //
                            * tFI->N() / tHistory;
                }
                else
                {
                    mdDamagedu( tDofIndex ) =                                                                   //
                            ( 1.0 - tDamage                                                                     //
                                    + mKappa0 * mAlpha * mBeta * std::exp( mBeta * ( mKappa0 - 1.0e-12 ) ) )    //
                            * tFI->N() / 1.0e-12;
                }
            }
            else
            {
                mdDamagedu( tDofIndex ).fill( 0.0 );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dDamagedu_smooth_exponential(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set the derivative size
            mdDamagedu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // if derivative dof is history dof
            if ( aDofTypes( 0 ) == mDofHistory )
            {
                // get the nonlocal equivalent strain history value
                real tHistory = tFI->val()( 0 );

                // get the damage value
                real tDamage = this->damage()( 0 );

                // compute ddamagedu based on history
                mdDamagedu( tDofIndex ) =
                        mAlpha * tDamage * std::exp( mAlpha * ( mKappa0 - tHistory ) ) * tFI->N();
            }
            else
            {
                mdDamagedu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage()
        {
            // call function pointer for evaluation
            ( this->*m_eval_smooth_damage )( mDamageParam );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_noSmoothing(
                const Matrix< DDRMat >& aSmoothParam )
        {
            // FIXME check correct implementation
            // no smoothing and smooth damage is equal to damage
            mSmoothDamage = this->damage();

            // include the zero part
            if ( mSmoothDamage( 0, 0 ) < 0.0 )
            {
                mSmoothDamage = { { 0.0 } };
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_ks(
                const Matrix< DDRMat >& aSmoothParam )
        {
            // smoothing with ks function
            mSmoothDamage = std::log( 1.0 + std::exp( mSmoothC * this->damage()( 0 ) ) ) / mSmoothC;
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_smooth_damage_corrected_ks(
                const Matrix< DDRMat >& aSmoothParam )
        {
            // smoothing with corrected ks function
            mSmoothDamage = std::log( 1.0 + std::exp( mSmoothC * this->damage()( 0 ) ) ) / std::log( 1.0 + std::exp( mSmoothC ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // call function pointer for evaluation
            ( this->*m_eval_dSmoothDamagedu )( aDofTypes );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_noSmoothing(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // no smoothing and smooth damage is equal to damage
            real tDamage = this->damage()( 0 );

            // include the zero part
            if ( tDamage < 0.0 )
            {
                // fill with zero as no damage part
                mdSmoothDamagedu( tDofIndex ).fill( 0.0 );
            }
            else
            {
                // no smoothing and smooth damage derivative is equal to damage derivative
                mdSmoothDamagedu( tDofIndex ) = this->dDamagedu( aDofTypes );
            }
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_ks(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the damage value
            real tDamage = this->damage()( 0 );

            // smoothing with ks function
            mdSmoothDamagedu( tDofIndex ) =    //
                    std::exp( mSmoothC * tDamage ) * this->dDamagedu( aDofTypes ) / ( 1.0 + std::exp( mSmoothC * tDamage ) );
        }

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dSmoothDamagedu_corrected_ks(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the damage value
            real tDamage = this->damage()( 0 );

            // smoothing with ks function
            mdSmoothDamagedu( tDofIndex ) =    //
                    mSmoothC * std::exp( mSmoothC * tDamage ) * this->dDamagedu( aDofTypes ) / ( std::log( 1.0 + std::exp( mSmoothC ) ) * ( 1.0 + std::exp( mSmoothC * tDamage ) ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_history()
        {
            // get the nonlocal equivalent strain field interpolator
            Field_Interpolator* tFINlEqStrain = mFIManager->get_field_interpolators_for_type( mDofNonlocalEqStrain );

            // get the nonlocal equivalent strain value
            real tNlEqStrain = tFINlEqStrain->val()( 0 );

            // get the reference equivalent strain value
            mHistory = this->history_ref();

            // check if nonlocal equivalent strain larger than history
            if ( tNlEqStrain > mHistory( 0 ) )
            {
                // update history value
                mHistory( 0 ) = tNlEqStrain;
            }
        }
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dHistorydu(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdHistorydu
            mdHistorydu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the nonlocal equivalent strain field interpolator
            Field_Interpolator* tFINlEqStrain = mFIManager->get_field_interpolators_for_type( mDofNonlocalEqStrain );

            // get the nonlocal equivalent strain value
            real tNlEqStrain = tFINlEqStrain->val()( 0 );

            // get the reference equivalent strain value
            real tHistoryRef = this->history_ref()( 0 );

            // check if nonlocal equivalent strain larger than history
            if ( tNlEqStrain > tHistoryRef && aDofTypes( 0 ) == mDofNonlocalEqStrain )
            {
                // MORIS_LOG_INFO( "CM_Struc_Linear_Isotropic_Damage::eval_dHistorydu - Growth" );
                // compute derivative of updated history value
                mdHistorydu( tDofIndex ) = tFI->N();
            }
            else if ( tNlEqStrain <= tHistoryRef && aDofTypes( 0 ) == mDofHistory )
            {
                // MORIS_LOG_INFO( "CM_Struc_Linear_Isotropic_Damage::eval_dHistorydu - No growth" );
                // compute derivative of kept history value
                mdHistorydu( tDofIndex ) = this->dHistoryRefdu( aDofTypes );
            }
            else
            {
                mdHistorydu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_history_ref()
        {
            // grab and save current evaluation point
            Field_Interpolator* tFIHistory = mFIManager->get_field_interpolators_for_type( mDofHistory );
            const Matrix<DDRMat>& tSpaceLocalParamPoint = tFIHistory->get_space();
            const Matrix<DDRMat> tTimeLocalParamPoint = tFIHistory->get_time();

            // create evaluation point at start of time slab
            Matrix< DDRMat > tStartLocalParamPoint( mSpaceDim + 1, 1 );
            tStartLocalParamPoint( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tSpaceLocalParamPoint.matrix_data();
            tStartLocalParamPoint( mSpaceDim )    = -1.0;

            // pass in evaluation point at start of time slab
            tFIHistory->set_space_time( tStartLocalParamPoint );

            // evaluate the nonlocal equivalent strain history
            mHistoryRef = tFIHistory->val();

            // avoid using zero value for history variable
            mHistoryRef = {{ std::max( mHistoryRef( 0 ), mEpsilon ) }};

            // restore current evaluation point
            tStartLocalParamPoint( mSpaceDim ) = tTimeLocalParamPoint( 0 );
            tFIHistory->set_space_time( tStartLocalParamPoint );
        }
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::eval_dHistoryRefdu(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdHistorydu
            mdHistoryRefdu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // check if nonlocal equivalent strain larger than history
            if ( aDofTypes( 0 ) == mDofHistory && this->history_ref()( 0 ) > mEpsilon )
            {
                // grab and save current evaluation point
                const Matrix<DDRMat>& tSpaceLocalParamPoint = tFI->get_space();
                const Matrix<DDRMat> tTimeLocalParamPoint  = tFI->get_time();

                // create evaluation point at start of time slab
                Matrix< DDRMat > tStartLocalParamPoint( mSpaceDim + 1, 1 );
                tStartLocalParamPoint( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tSpaceLocalParamPoint.matrix_data();
                tStartLocalParamPoint( mSpaceDim )    = -1.0;

                // pass in evaluation point at start of time slab
                tFI->set_space_time( tStartLocalParamPoint );

                // compute derivative of kept history value
                mdHistoryRefdu( tDofIndex ) = tFI->N();

                // restore current evaluation point
                tStartLocalParamPoint( mSpaceDim ) = tTimeLocalParamPoint( 0 );
                tFI->set_space_time( tStartLocalParamPoint );
            }
            else
            {
                mdHistoryRefdu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Linear_Isotropic_Damage::select_derivative_FD(
                enum CM_Request_Type                aCMRequestType,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const Matrix< DDRMat >&             aJump,
                enum CM_Function_Type               aCMFunctionType )
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
                case CM_Request_Type::DAMAGE:
                {
                    return this->damage( aCMFunctionType );
                    break;
                }
                case CM_Request_Type::SMOOTH_DAMAGE:
                {
                    return this->smooth_damage( aCMFunctionType );
                    break;
                }
                case CM_Request_Type::EQSTRAIN:
                {
                    return this->equivalent_strain( aCMFunctionType );
                    break;
                }
                case CM_Request_Type::HISTORY:
                {
                    return this->history( aCMFunctionType );
                    break;
                }
                default:
                    MORIS_ERROR( false, "Constitutive_Model::select_derivative_FD: aCMRequestType undefined" );
                    return this->strain( aCMFunctionType );
            }
        }
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_Isotropic_Damage::set_derivative_FD(
                enum CM_Request_Type                aCMRequestType,
                Matrix< DDRMat >&                   aDerivativeFD,
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type               aCMFunctionType )
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
                case CM_Request_Type::DAMAGE:
                {
                    // set value to storage
                    mdDamagedu( tDofIndex ) = aDerivativeFD;
                    break;
                }
                case CM_Request_Type::SMOOTH_DAMAGE:
                {
                    // set value to storage
                    mdSmoothDamagedu( tDofIndex ) = aDerivativeFD;
                    break;
                }
                case CM_Request_Type::EQSTRAIN:
                {
                    // set value to storage
                    mdEqStraindu( tDofIndex ) = aDerivativeFD;
                    break;
                }
                case CM_Request_Type::HISTORY:
                {
                    // set value to storage
                    mdHistorydu( tDofIndex ) = aDerivativeFD;
                    break;
                }
                default:
                    MORIS_ERROR( false, "CM_Struc_Linear_Isotropic_Damage::set_derivative_FD: aCMRequestType undefined" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

