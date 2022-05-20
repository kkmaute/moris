// FIXME: set_model_type not yet working? -> default set to PLANE_STRAIN
// FIXME:: Cleaning up the projected jump mess in Neo_Hookean eval_dTestTractiondDOF_second_piola_kirchhoff
// FIXME:: Testtraction transpose check - consistent with linear case
// FIXME:: no breaks in different cases in flux()

#include "cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp"
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

        CM_Struc_Nonlinear_Isotropic::CM_Struc_Nonlinear_Isotropic(enum Constitutive_Type aConstType)
                {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "YoungsModulus" ]        = static_cast< uint >( CM_Property_Type::EMOD );
            mPropertyMap[ "PoissonRatio" ]         = static_cast< uint >( CM_Property_Type::NU );
                }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::reset_eval_flags()
        {
            // reset flag from parent class
            Constitutive_Model::reset_eval_flags();

            // reset the deformation related flags
            mDefGradEval  = true;
            mRCGDefEval   = true;
            mInvRCGDefEval   = true;
            mLGStrainEval = true;
            mEAStrainEval = true;
            mDGStrainEval = true;

            mTestDefGradEval     = true;
            mdDefGradduEval.fill(  true );
            mdLGStrainduEval.fill( true );
            mdEAStrainduEval.fill( true );
            mdDGStrainduEval.fill( true );

            mLGTestStrainEval = true;
            mEATestStrainEval = true;
            mDGTestStrainEval = true;
            mdLGTestStrainduEval.fill( true );
            mdEATestStrainduEval.fill( true );
            mdDGTestStrainduEval.fill( true );

            // reset the Jacobian related flags
            mVolumeChangeJEval = true;

            // reset the stress related flags
            m1PKStressEval    = true;
            m2PKStressEval    = true;
            mCauchyStressEval = true;

            mFluxProjEval = true;

            md1PKStressduEval.fill(    true );
            md2PKStressduEval.fill(    true );
            mdCauchyStressduEval.fill( true );

            m1PKTractionEval    = true;
            m2PKTractionEval    = true;
            mCauchyTractionEval = true;

            md1PKTractionduEval.fill(    true );
            md2PKTractionduEval.fill(    true );
            mdCauchyTractionduEval.fill( true );

            m1PKTestTractionEval.fill(    true );
            m2PKTestTractionEval.fill(    true );
            mCauchyTestTractionEval.fill( true );

            md1PKTestTractionduEval.fill(    true );
            md2PKTestTractionduEval.fill(    true );
            mdCauchyTestTractionduEval.fill( true );
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::build_global_dof_type_list()
        {
            // build list from parent class
            Constitutive_Model::build_global_dof_type_list();

            // number of dof types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();
            uint tNumDirectDofTypes = mDofTypes.size();

            // init deformation related flags
            mdDefGradduEval.set_size( tNumGlobalDofTypes, 1,  true );
            mdLGStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEAStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdDGStrainduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdLGTestStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEATestStrainduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdDGTestStrainduEval.set_size( tNumGlobalDofTypes, 1, true );

            // init stress related flags
            md1PKStressduEval.set_size( tNumGlobalDofTypes, 1,    true );
            md2PKStressduEval.set_size( tNumGlobalDofTypes, 1,    true );
            mdCauchyStressduEval.set_size( tNumGlobalDofTypes, 1, true );

            md1PKTractionduEval.set_size( tNumGlobalDofTypes, 1,    true );
            md2PKTractionduEval.set_size( tNumGlobalDofTypes, 1,    true );
            mdCauchyTractionduEval.set_size( tNumGlobalDofTypes, 1, true );

            m1PKTestTractionEval.set_size( tNumGlobalDofTypes, 1,    true );
            m2PKTestTractionEval.set_size( tNumGlobalDofTypes, 1,    true );
            mCauchyTestTractionEval.set_size( tNumGlobalDofTypes, 1, true );

            md1PKTestTractionduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes,    true );
            md2PKTestTractionduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes,    true );
            mdCauchyTestTractionduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );

            // init deformation related storage
            mdDefGraddu.resize( tNumGlobalDofTypes );
            mdLGStraindu.resize( tNumGlobalDofTypes );
            mdEAStraindu.resize( tNumGlobalDofTypes );
            mdDGStraindu.resize( tNumGlobalDofTypes );

            mdLGTestStraindu.resize( tNumGlobalDofTypes );
            mdEATestStraindu.resize( tNumGlobalDofTypes );
            mdDGTestStraindu.resize( tNumGlobalDofTypes );

            // init stress related storage
            md1PKStressdu.resize( tNumGlobalDofTypes );
            md2PKStressdu.resize( tNumGlobalDofTypes );
            mdCauchyStressdu.resize( tNumGlobalDofTypes );

            md1PKTractiondu.resize( tNumDirectDofTypes );
            md2PKTractiondu.resize( tNumDirectDofTypes );
            mdCauchyTractiondu.resize( tNumDirectDofTypes );

            m1PKTestTraction.resize( tNumDirectDofTypes );
            m2PKTestTraction.resize( tNumDirectDofTypes );
            mCauchyTestTraction.resize( tNumDirectDofTypes );

            md1PKTestTractiondu.resize( tNumDirectDofTypes );
            md2PKTestTractiondu.resize( tNumDirectDofTypes );
            mdCauchyTestTractiondu.resize( tNumDirectDofTypes );

            for( uint iDirectDof = 0; iDirectDof < tNumDirectDofTypes; iDirectDof++ )
            {
                md1PKTestTractiondu( iDirectDof ).resize( tNumGlobalDofTypes );
                md2PKTestTractiondu( iDirectDof ).resize( tNumGlobalDofTypes );
                mdCauchyTestTractiondu( iDirectDof ).resize( tNumGlobalDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_dof_type_list(
                Cell< Cell< MSI::Dof_Type > > aDofTypes,
                Cell< std::string >           aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // if displacement dof type string
                if( tDofString == "Displacement" )
                {
                    mDofDispl = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Struc_Nonlinear_Isotropic::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_local_properties()
        {
            // set the Young's modulus property
            mPropEMod = get_property( "YoungsModulus" );

            // set the Poisson ratio property
            mPropPoisson = get_property( "PoissonRatio" );

            // check that essential properties exist
            MORIS_ASSERT( mPropEMod,
                    "CM_Struc_Linear_Isotropic::set_local_properties - Young's modulus property does not exist.\n");

            MORIS_ASSERT( mPropPoisson,
                    "CM_Struc_Linear_Isotropic::set_local_properties - Poisson ratio property does not exist.\n");
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_function_pointers()
        {
            switch( mSpaceDim )
            {
                case 2 :
                {
                    // set Voigt notation map
                    mVoigtNonSymMap = {{ 0, 1 }, { 2, 3 }};
                    mVoigtSymMap    = {{ 0, 2 }, { 2, 1 }};

                    m_eval_test_deformation_gradient = &CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient_2d;

                    m_eval_dLGStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dLGStraindDOF_2d;
                    m_eval_dEAStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_2d;

                    m_eval_dLGTestStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dLGTestStraindDOF_2d;
                    m_eval_dEATestStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_2d;

                    m_flatten_normal        = &CM_Struc_Nonlinear_Isotropic::flatten_normal_2d;
                    m_flatten_normal_nonsym = &CM_Struc_Nonlinear_Isotropic::flatten_normal_nonsym_2d;

                    m_proj_sym  = &CM_Struc_Nonlinear_Isotropic::proj_sym_2d;
                    m_proj_nsym  = &CM_Struc_Nonlinear_Isotropic::proj_nsym_2d;

                    //                    mFluxHead.set_size( 4, 4, 0.0 ); // FIXME
                    //                    mCauchy.set_size( 3, 1, 0.0 );    // FIXME
                    mFluxProj.set_size( 4, 4, 0.0);

                    break;
                }
                case 3 :
                {
                    // set Voigt notation map
                    mVoigtNonSymMap = {{ 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 }};
                    mVoigtSymMap    = {{ 0, 5, 4 }, { 5, 1, 3 }, { 4, 3, 2 }};

                    m_eval_test_deformation_gradient = &CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient_3d;

                    m_eval_dLGStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dLGStraindDOF_3d;
                    m_eval_dEAStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_3d;

                    m_eval_dLGTestStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dLGTestStraindDOF_3d;
                    m_eval_dEATestStraindDOF = &CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_3d;

                    m_flatten_normal        = &CM_Struc_Nonlinear_Isotropic::flatten_normal_3d;
                    m_flatten_normal_nonsym = &CM_Struc_Nonlinear_Isotropic::flatten_normal_nonsym_3d;

                    m_proj_sym  = &CM_Struc_Nonlinear_Isotropic::proj_sym_3d;
                    m_proj_nsym  = &CM_Struc_Nonlinear_Isotropic::proj_nsym_3d;

                    // FIXME
                    //                    mFluxHead.set_size( 9, 9, 0.0 );
                    //                    mStrain.set_size( 6, 1, 0.0 );

                    mFluxProj.set_size( 9, 9, 0.0);

                    // list number of normal stresses and strains
                    mNumNormalStress = 3;
                    mNumNormalStrain = 3;

                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Nonlinear isotropic elasticity implemented only for 2d and 3d" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::deformation_gradient()
        {
            // if the deformation gradient was not evaluated
            if( mDefGradEval )
            {
                // evaluate the deformation gradient
                this->eval_deformation_gradient();

                // set bool for evaluation
                mDefGradEval = false;
            }

            // return the deformation gradient value
            return mDefGrad;
        }

        void CM_Struc_Nonlinear_Isotropic::eval_deformation_gradient()
        {
            // get the displacement gradient
            const Matrix< DDRMat > & tDisplGradx =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

            // evaluate the deformation gradient as F = I + du/dX
            mDefGrad = trans( tDisplGradx ) + eye( mSpaceDim, mSpaceDim );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::test_deformation_gradient(
                const Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if the deformation gradient was not evaluated
            if( mTestDefGradEval )
            {
                // evaluate the deformation gradient
                this->eval_test_deformation_gradient( aTestDofTypes );

                // set bool for evaluation
                mTestDefGradEval = false;
            }

            // return the deformation gradient value
            return mTestDefGrad;
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient_2d(
                const Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if test dof type is displacement
            if( aTestDofTypes( 0 ) == mDofDispl )
            {
                // get the derivative of the displacement shape functions
                const Matrix< DDRMat > & tdnNdxn =
                        mFIManager->get_field_interpolators_for_type( mDofDispl )->dnNdxn( 1 );

                // get number of bases for displacement
                uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofDispl )->get_number_of_space_time_bases();

                // set size
                mTestDefGrad.set_size( 4, 2 * tNumBases, 0.0 );

                // evaluate the test deformation gradient
                mTestDefGrad( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn.get_row( 0 );
                mTestDefGrad( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn.get_row( 1 );
                mTestDefGrad( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn.get_row( 0 );
                mTestDefGrad( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn.get_row( 1 );
            }
            else
            {
                // not implemented for other test dof type than displacement
                MORIS_ERROR( false,
                        "CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient - unsupported test dof type" );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient_3d(
                const Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if test dof type is displacement
            if( aTestDofTypes( 0 ) == mDofDispl )
            {
                // get the derivative of the displacement shape functions
                const Matrix< DDRMat > & tdnNdxn =
                        mFIManager->get_field_interpolators_for_type( mDofDispl )->dnNdxn( 1 );

                // get number of bases for displacement
                uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofDispl )->get_number_of_space_time_bases();

                // set size
                mTestDefGrad.set_size( 9, 3 * tNumBases, 0.0 );

                // evaluate the test deformation gradient
                mTestDefGrad( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn.get_row( 0 );
                mTestDefGrad( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn.get_row( 1 );
                mTestDefGrad( { 2, 2 }, { 0, tNumBases - 1 } ) = tdnNdxn.get_row( 2 );
                mTestDefGrad( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn.get_row( 0 );
                mTestDefGrad( { 4, 4 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn.get_row( 1 );
                mTestDefGrad( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn.get_row( 2 );
                mTestDefGrad( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn.get_row( 0 );
                mTestDefGrad( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn.get_row( 1 );
                mTestDefGrad( { 8, 8 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn.get_row( 2 );
            }
            else
            {
                // not implemented for other test dof type than displacement
                MORIS_ERROR( false,
                        "CM_Struc_Nonlinear_Isotropic::eval_test_deformation_gradient - unsupported test dof type" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const real &
        CM_Struc_Nonlinear_Isotropic::volume_change_jacobian(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Nonlinear_Isotropic::volume_change_jacobian - Only DEFAULT CM function type known in base class." );

            // if the Jacobian was not evaluated
            if( mVolumeChangeJEval )
            {
                // evaluate the Jacobian
                this->eval_volume_change_jacobian();

                // set bool for evaluation
                mVolumeChangeJEval = false;
            }
            // return the Jacobian value
            return mVolumeChangeJ;
        }

        void CM_Struc_Nonlinear_Isotropic::eval_volume_change_jacobian()
        {
            // evaluate the volume change Jacobian
            // FIXME better to use a self implemented version of det
            mVolumeChangeJ = det( this->deformation_gradient() );

            // check for negative Jacobian
            MORIS_ASSERT( mVolumeChangeJ > 0,
                    "CM_Struc_Nonlinear_Isotropic::volume_change_jacobian - Negative volume change jacobian." );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::right_cauchy_green_deformation_tensor(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Nonlinear_Isotropic::right_cauchy_green_deformation_tensor - Only DEFAULT CM function type known in base class." );

            // if the right cauchy green deformation tensor was not evaluated
            if( mRCGDefEval )
            {
                // evaluate the right cauchy green deformation tensor
                this->eval_right_cauchy_green_deformation_tensor();

                // set bool for evaluation
                mRCGDefEval = false;
            }
            // return the right cauchy green deformation tensor
            return mRCGDef;
        }

        void CM_Struc_Nonlinear_Isotropic::eval_right_cauchy_green_deformation_tensor()
        {
            // evaluate the right cauchy green deformation tensor (full)
            Matrix< DDRMat > tRCGDefFull = trans( this->deformation_gradient() ) * this->deformation_gradient();

            // store in Voigt notation
            this->full_to_voigt_sym_strain( tRCGDefFull, mRCGDef );
        }

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::inv_right_cauchy_green_deformation_tensor(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Struc_Nonlinear_Isotropic::inv_right_cauchy_green_deformation_tensor - Only DEFAULT CM function type known in base class." );

            // if the inverse of the right cauchy green deformation tensor was not evaluated
            if( mInvRCGDefEval )
            {
                // evaluate the inverse of the right cauchy green deformation tensor
                this->eval_inv_right_cauchy_green_deformation_tensor();

                // set bool for evaluation
                mInvRCGDefEval = false;
            }
            // return the inverse of the right cauchy green deformation tensor (full)
            return mInvRCGDef;
        }

        void CM_Struc_Nonlinear_Isotropic::eval_inv_right_cauchy_green_deformation_tensor()
        {
            // evaluate the right cauchy green deformation tensor (full)
            Matrix< DDRMat > tRCGDefFull;

            // store in Matrix notation
            this->voigt_to_full_sym_strain( this->right_cauchy_green_deformation_tensor(), tRCGDefFull );

            // Calculate the inverse of the right cauchy green deformation tensor (full)
            mInvRCGDef = inv(tRCGDefFull);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::strain(
                enum CM_Function_Type aCMFunctionType  )
        {
            // switch on CM function type
            switch( aCMFunctionType )
            {
                // lagrangian/Green strain
                case CM_Function_Type::LAGRANGIAN :
                {
                    // if the  lagrangian/Green strain was not evaluated
                    if( mLGStrainEval )
                    {
                        // evaluate the  lagrangian/Green strain
                        this->eval_lagrangian_green_strain_tensor();

                        // set bool for evaluation
                        mLGStrainEval = false;
                    }
                    // return the lagrangian/Green strain
                    return mLGStrain;
                }
                // eulerian/Almansi strain
                case CM_Function_Type::EULERIAN :
                {
                    // if the eulerian/Almansi strain was not evaluated
                    if( mEAStrainEval )
                    {
                        // evaluate the eulerian/Almansi strain
                        this->eval_eulerian_almansi_strain_tensor();

                        // set bool for evaluation
                        mEAStrainEval = false;
                    }
                    // return the eulerian/Almansi strain
                    return mEAStrain;
                }
                case CM_Function_Type::DEFORMATION_GRADIENT:
                {
                    // if the eulerian/Almansi strain was not evaluated
                    if( mDGStrainEval )
                    {
                        // evaluate the eulerian/Almansi strain
                        this->eval_deformation_gradient_strain_tensor();

                        // set bool for evaluation
                        mDGStrainEval = false;
                    }
                    // return the eulerian/Almansi strain
                    return mDGStrain;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::strain - Unknown strain type.");
                    return mStrain;
                }
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_lagrangian_green_strain_tensor()
        {
            // evaluate the lagrangian/green strain tensor
            mLGStrain = this->right_cauchy_green_deformation_tensor();

            // substract identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
            mLGStrain({0, mSpaceDim - 1}, { 0, 0 }) -= tI.matrix_data();

            // apply factor 0.5
            mLGStrain = 0.5 * mLGStrain;
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_eulerian_almansi_strain_tensor()
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_eulerian_almansi_strain_tensor - Not implemented yet.");
        }

        // FIXME: storage for deformation gradient is set twice:
        // Strain - Vector notation -> to be consistent with other work conjugates
        // DefGrad: Matrix notation

        void
        CM_Struc_Nonlinear_Isotropic::eval_deformation_gradient_strain_tensor()
        {
            this->full_to_voigt_sym_strain( this->deformation_gradient(), mDGStrain );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::testStrain(
                enum CM_Function_Type aCMFunctionType  )
        {
            // switch on CM function type
            switch( aCMFunctionType )
            {
                // lagrangian/Green strain
                case CM_Function_Type::LAGRANGIAN :
                {
                    // if the  lagrangian/Green strain was not evaluated
                    if( mLGTestStrainEval )
                    {
                        // evaluate the  lagrangian/Green strain
                        this->eval_testStrain_lagrangian_green_strain();

                        // set bool for evaluation
                        mLGTestStrainEval = false;
                    }
                    // return the lagrangian/Green strain
                    return mLGTestStrain;
                }
                // eulerian/Almansi strain
                case CM_Function_Type::EULERIAN :
                {
                    // if the eulerian/Almansi strain was not evaluated
                    if( mEATestStrainEval )
                    {
                        // evaluate the eulerian/Almansi strain
                        this->eval_testStrain_eulerian_almansi_strain();

                        // set bool for evaluation
                        mEATestStrainEval = false;
                    }
                    // return the eulerian/Almansi strain
                    return mEATestStrain;
                }
                // eulerian/Almansi strain
                case CM_Function_Type::DEFORMATION_GRADIENT:
                {
                    // if the eulerian/Almansi strain was not evaluated
                    if( mDGTestStrainEval )
                    {
                        // evaluate the eulerian/Almansi strain
                        this->eval_testStrain_deformation_gradient_strain();

                        // set bool for evaluation
                        mDGTestStrainEval = false;
                    }
                    // return the eulerian/Almansi strain
                    return mDGTestStrain;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::testStrain - Unknown strain type.");
                    return mTestStrain;
                }
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_testStrain_lagrangian_green_strain()
        {
            Cell< MSI::Dof_Type > tDofTypes;
            tDofTypes.resize(1);
            tDofTypes( 0 ) = mDofDispl;

            mLGTestStrain = this->dStraindDOF( tDofTypes , CM_Function_Type::LAGRANGIAN );
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_testStrain_eulerian_almansi_strain()
        {
            MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::eval_testStrain_eulerian_almansi_strain - Unknown strain type.");
        }

        // FIXME: storage of teststrain based on displacement gradient is set twice to be consistent with general formulation
        void
        CM_Struc_Nonlinear_Isotropic::eval_testStrain_deformation_gradient_strain()
        {
            Cell< MSI::Dof_Type > tDofTypes;
            tDofTypes.resize(1);
            tDofTypes( 0 ) = mDofDispl;

            mDGTestStrain = this->dStraindDOF( tDofTypes, CM_Function_Type::DEFORMATION_GRADIENT );
        }
        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::dStraindDOF(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "Constitutive_Model::dStraindDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // lagrangian/Green strain
                case CM_Function_Type::LAGRANGIAN :
                {
                    // if the derivative has not been evaluated yet
                    if( mdLGStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dLGStraindDOF( aDofType );

                        // set bool for evaluation
                        mdLGStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdLGStraindu( tDofIndex );
                }
                // eulerian/Almansi strain
                case CM_Function_Type::EULERIAN :
                {
                    // if the derivative has not been evaluated yet
                    if( mdEAStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dEAStraindDOF( aDofType );

                        // set bool for evaluation
                        mdEAStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdEAStraindu( tDofIndex );
                }
                // deformation gradient strain
                case CM_Function_Type::DEFORMATION_GRADIENT:
                {
                    // if the derivative has not been evaluated yet
                    if( mdDGStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dDGStraindDOF( aDofType );

                        // set bool for evaluation
                        mdDGStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdDGStraindu( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::dStraindDOF - Unknown strain type.");
                    return mdStraindDof( tDofIndex );
                }
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_dLGStraindDOF_2d(
                const Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdLGstraindDof
            mdLGStraindu( tDofIndex ).set_size(
                    3,
                    tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // compute displacement gradient
                const Matrix< DDRMat > & tdnNdxn = tFI->dnNdxn( 1 );

                // get number of bases for displacement
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // get the deformation gradient
                const Matrix< DDRMat > & tF = this->deformation_gradient();

                // populate the derivative of the lagrangian/green strain
                mdLGStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tF( 0, 0 ) * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = tF( 0, 1 ) * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = tF( 0, 0 ) * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } ) + tF( 0, 1 ) * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

                mdLGStraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 0 ) * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 1 ) * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 0 ) * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } ) + tF( 1, 1 ) * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_dLGStraindDOF_3d(
                const Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdLGstraindDof
            mdLGStraindu( tDofIndex ).set_size(
                    6,
                    tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // compute displacement gradient
                const Matrix< DDRMat > & tdnNdxn = tFI->dnNdxn( 1 );

                // get number of bases for displacement
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // get the deformation gradient
                const Matrix< DDRMat > & tF = this->deformation_gradient();

                // FIXME not checked
                // populate the derivative of the lagrangian/green strain
                mdLGStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tF( 0, 0 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = tF( 0, 1 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = tF( 0, 2 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) = tF( 0, 1 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } ) + tF( 0, 2 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 4, 4 }, { 0, tNumBases - 1 } ) = tF( 0, 0 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } ) + tF( 0, 2 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 5, 5 }, { 0, tNumBases - 1 } ) = tF( 0, 0 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } ) + tF( 0, 1 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

                mdLGStraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 0 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 1 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 2 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 1 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } ) + tF( 1, 2 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 4, 4 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 0 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } ) + tF( 1, 2 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tF( 1, 1 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } ) + tF( 1, 0 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );;

                mdLGStraindu( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 0 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 1 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 2 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 2 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } ) + tF( 2, 1 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 2 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } ) + tF( 2, 0 )*tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
                mdLGStraindu( tDofIndex )( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tF( 2, 1 )*tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } ) + tF( 2, 0 )*tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_2d(
                const Cell< MSI::Dof_Type > & aDofTypes )
        {

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_2d - Not implemented yet.");
            }
        }

        void
        CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_3d(
                const Cell< MSI::Dof_Type > & aDofTypes )
        {

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dEAStraindDOF_3d - Not implemented yet.");
            }

        }

        // FIXME:: Storage for derivative of deformation gradient wrt displacement is set twice to include the possibility of
        //    taking the derivative wrt. to dof
        void
        CM_Struc_Nonlinear_Isotropic::eval_dDGStraindDOF(
                const Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdDGStraindu
            mdDGStraindu( tDofIndex ).set_size(
                    6,
                    tFI->get_number_of_space_time_coefficients(), 0.0 );
            //
            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                mdDGStraindu( tDofIndex ) = this->test_deformation_gradient(aDofTypes) ;
            }

        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Struc_Nonlinear_Isotropic::dTestStraindDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type aCMFunctionType )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "Constitutive_Model::dTestStraindDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // lagrangian/Green strain
                case CM_Function_Type::LAGRANGIAN :
                {
                    // if the derivative has not been evaluated yet
                    if( mdLGTestStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dLGTestStraindDOF( aDofType );

                        // set bool for evaluation
                        mdLGTestStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdLGTestStraindu( tDofIndex );
                }
                // eulerian/Almansi strain
                case CM_Function_Type::EULERIAN :
                {
                    // if the derivative has not been evaluated yet
                    if( mdEATestStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dEATestStraindDOF( aDofType );

                        // set bool for evaluation
                        mdEATestStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdEATestStraindu( tDofIndex );
                }
                // deformation gradient strain
                case CM_Function_Type::DEFORMATION_GRADIENT:
                {
                    // if the derivative has not been evaluated yet
                    if( mdDGTestStrainduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dDGTestStraindDOF( aDofType );

                        // set bool for evaluation
                        mdDGTestStrainduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdDGTestStraindu( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::dTestStraindDOF - Unknown strain type.");
                    return mdTestStraindDof( tDofIndex );
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::eval_dLGTestStraindDOF_2d(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get displacement field interpolator
            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // get the derivative of the displacement shape functions
            const Matrix< DDRMat > & tdnNdxn =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases = tFIDispl->get_number_of_space_time_bases();

            // build the derivative of the test strain
            mdLGTestStraindu( tDofIndex ).set_size( 4, tNumBases * 2, 0.0 );

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // populate the derivative of the teststrain
                mdLGTestStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

                mdLGTestStraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            }
        }

        void CM_Struc_Nonlinear_Isotropic::eval_dLGTestStraindDOF_3d(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get displacement field interpolator
            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // get the derivative of the displacement shape functions
            const Matrix< DDRMat > & tdnNdxn =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases = tFIDispl->get_number_of_space_time_bases();

            // build the derivative of the test strain
            mdLGTestStraindu( tDofIndex ).set_size( 9, tNumBases * 3, 0.0 );

            // if derivative dof is displacements dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // populate the derivative of the teststrain
                mdLGTestStraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );

                mdLGTestStraindu( tDofIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 4, 4 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );

                mdLGTestStraindu( tDofIndex )( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
                mdLGTestStraindu( tDofIndex )( { 8, 8 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            }
        }

        void CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_2d(
                const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_2d - Not implemented yet." );
        }

        void CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_3d(
                const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dEATestStraindDOF_3d - Not implemented yet." );
        }

        // FIXME: set mdDGTestStraindu as an empty matrix of the right size to be consistent with general formulation
        void CM_Struc_Nonlinear_Isotropic::eval_dDGTestStraindDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get displacement field interpolator
            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // get number of bases for displacement
            uint tNumBases = tFIDispl->get_number_of_space_time_bases();

            // build the derivative of the test strain
            mdDGTestStraindu( tDofIndex ).set_size( mSpaceDim * mSpaceDim, tNumBases * mSpaceDim, 0.0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::flux(
                enum CM_Function_Type aCMFunctionType)
        {
            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {

                    // if the first Piola_Kirchhoff stress tensor was not evaluated
                    if( m1PKStressEval )
                    {
                        // evaluate the first Piola_Kirchhoff stress tensor
                        this->eval_flux_first_piola_kirchhoff();

                        // set bool for evaluation
                        m1PKStressEval = false;
                    }
                    // return the first Piola_Kirchhoff stress tensor
                    return m1PKStress;

                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {

                    // if the second Piola_Kirchhoff stress tensor was not evaluated
                    if( m2PKStressEval )
                    {
                        // evaluate the second Piola_Kirchhoff stress tensor
                        this->eval_flux_second_piola_kirchhoff();

                        // set bool for evaluation
                        m2PKStressEval = false;
                    }
                    // return the second Piola_Kirchhoff stress tensor
                    return m2PKStress;

                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {

                    // if the Cauchy stress tensor was not evaluated
                    if( mCauchyStressEval )
                    {
                        // evaluate the Cauchy stress tensor
                        this->eval_flux_cauchy();

                        // set bool for evaluation
                        mCauchyStressEval = false;
                    }
                    // return the Cauchy stress tensor
                    return mCauchyStress;

                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::flux - Unknown stress type.");
                    return mFlux;
                }
            }
        }

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::flux(
                int aFlatType,
                enum CM_Function_Type aCMFunctionType)
        {
            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    switch( aFlatType )
                    {
                        // first Piola-Kirchhoff stress
                        case 1 :
                        {
                            if( mFluxProjEval )
                            {
                                // second Piola_Kirchhoff stress tensor
                                this->eval_flux_proj_nsym(CM_Function_Type::PK1);

                                // set bool for evaluation
                                mFluxProjEval = false;
                            }
                            // return the projected second Piola_Kirchhoff stress tensor
                            return mFluxProj;
                        }
                        case 0:
                        {
                            // if the first Piola_Kirchhoff stress tensor was not evaluated
                            if( m1PKStressEval )
                            {
                                // evaluate the first Piola_Kirchhoff stress tensor
                                this->eval_flux_first_piola_kirchhoff();

                                // set bool for evaluation
                                m1PKStressEval = false;
                            }
                            // return the first Piola_Kirchhoff stress tensor
                            return m1PKStress;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::flux - aFlatType unequal to 0 (no projection) or 1 (projection) .");
                        }
                    }
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    switch( aFlatType )
                    {
                        case 1:
                        {
                            if( mFluxProjEval )
                            {
                                // projected second Piola_Kirchhoff stress tensor
                                this->eval_flux_proj_sym(CM_Function_Type::PK2);

                                // set bool for evaluation
                                mFluxProjEval = false;
                            }
                            // return the projected second Piola_Kirchhoff stress tensor
                            return mFluxProj;
                        }
                        case 0:
                        {
                            // if the second Piola_Kirchhoff stress tensor was not evaluated
                            if( m2PKStressEval )
                            {
                                // evaluate the second Piola_Kirchhoff stress tensor
                                this->eval_flux_second_piola_kirchhoff();

                                // set bool for evaluation
                                m2PKStressEval = false;
                            }
                            // return the second Piola_Kirchhoff stress tensor
                            return m2PKStress;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::flux - aFlatType unequal to 0 (no projection) or 1 (projection) .");
                        }
                    }
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    switch( aFlatType )
                    {
                        case 1 :
                        {
                            if( mFluxProjEval )
                            {
                                // projected cauchy stress tensor
                                this->eval_flux_proj_sym(CM_Function_Type::CAUCHY);

                                // set bool for evaluation
                                mFluxProjEval = false;
                            }
                            // return the projected cauchy stress tensor
                            return mFluxProj;
                        }
                        case 0:
                        {
                            // if the Cauchy stress tensor was not evaluated
                            if( mCauchyStressEval )
                            {
                                // evaluate the Cauchy stress tensor
                                this->eval_flux_cauchy();

                                // set bool for evaluation
                                mCauchyStressEval = false;
                            }
                            // return the Cauchy stress tensor
                            return mCauchyStress;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::flux - aFlatType unequal to 0 (no projection) or 1 (projection) .");
                        }
                    }
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::flux - Unknown stress type.");
                    return mFlux;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::dFluxdDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                enum CM_Function_Type aCMFunctionType )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    // if the derivative of first Piola_Kirchhoff stress has not been evaluated yet
                    if( md1PKStressduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_d1PKStressdDOF( aDofType );

                        // set bool for evaluation
                        md1PKStressduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return md1PKStressdu( tDofIndex );
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    // if the derivative of second Piola_Kirchhoff stress has not been evaluated yet
                    if( md2PKStressduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_d2PKStressdDOF( aDofType );

                        // set bool for evaluation
                        md2PKStressduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return md2PKStressdu( tDofIndex );
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    // if the derivative of cauchy stress has not been evaluated yet
                    if( mdCauchyStressduEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dCauchyStressdDOF( aDofType );

                        // set bool for evaluation
                        mdCauchyStressduEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdCauchyStressdu( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::dfluxdDOF - Unknown stress type.");
                    return mdFluxdDof( tDofIndex );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Struc_Nonlinear_Isotropic::traction(
                const Matrix< DDRMat > & aNormal,
                enum CM_Function_Type aCMFunctionType )
        {
            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    // if the traction based on first Piola_Kirchhoff stress was not evaluated
                    if( m1PKTractionEval )
                    {
                        // evaluate the traction based on first Piola_Kirchhoff stress
                        this->eval_traction_first_piola_kirchhoff( aNormal );

                        // set bool for evaluation
                        m1PKTractionEval = false;
                    }
                    // return the traction based on first Piola_Kirchhoff stress
                    return m1PKTraction;
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    // if the traction based on second Piola_Kirchhoff stress was not evaluated
                    if( m2PKTractionEval )
                    {
                        // evaluate the traction based on second Piola_Kirchhoff stress
                        this->eval_traction_second_piola_kirchhoff( aNormal );

                        // set bool for evaluation
                        m2PKTractionEval = false;
                    }
                    // return the traction based on second Piola_Kirchhoff stress
                    return m2PKTraction;
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    // if the traction based on Cauchy stress was not evaluated
                    if( mCauchyTractionEval )
                    {
                        // evaluate the traction based on Cauchy stress
                        this->eval_traction_cauchy( aNormal );

                        // set bool for evaluation
                        mCauchyTractionEval = false;
                    }
                    // return the Cauchy stress tensor
                    return mCauchyTraction;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::traction - Unknown stress type.");
                    return mTraction;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                enum CM_Function_Type               aCMFunctionType )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    // if the derivative of traction based on first Piola_Kirchhoff stress was not evaluated
                    if( md1PKTractionduEval( tDofIndex ) )
                    {
                        // evaluate the derivative of traction based on first Piola_Kirchhoff stress
                        this->eval_dTractiondDOF_first_piola_kirchhoff( aNormal, aDofType );

                        // set bool for evaluation
                        md1PKTractionduEval( tDofIndex ) = false;
                    }
                    // return the derivative of traction based on first Piola_Kirchhoff stress
                    return md1PKTractiondu( tDofIndex );
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    // if the derivative of traction based on second Piola_Kirchhoff stress was not evaluated
                    if( md2PKTractionduEval( tDofIndex ) )
                    {
                        // evaluate the derivative of traction based on second Piola_Kirchhoff stress
                        this->eval_dTractiondDOF_second_piola_kirchhoff( aNormal, aDofType );

                        // set bool for evaluation
                        md2PKTractionduEval( tDofIndex ) = false;
                    }
                    // return the derivative of traction based on second Piola_Kirchhoff stress
                    return md2PKTractiondu( tDofIndex );
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    // if the derivative of traction based on Cauchy stress was not evaluated
                    if( mdCauchyTractionduEval( tDofIndex ) )
                    {
                        // evaluate the derivative of traction based on Cauchy stress
                        this->eval_dTractiondDOF_cauchy( aNormal, aDofType );

                        // set bool for evaluation
                        mdCauchyTractionduEval( tDofIndex ) = false;
                    }
                    // return the derivative of traction based on Cauchy stress tensor
                    return mdCauchyTractiondu( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::dTractiondDOF - Unknown stress type.");
                    return mdTractiondDof( tDofIndex );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::testTraction(
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type               aCMFunctionType )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    // if the test traction based on first Piola_Kirchhoff stress was not evaluated
                    if( m1PKTestTractionEval( tTestDofIndex ) )
                    {
                        // evaluate the test traction based on first Piola_Kirchhoff stress
                        this->eval_testTraction_first_piola_kirchhoff( aNormal, aTestDofTypes );

                        // set bool for evaluation
                        m1PKTestTractionEval( tTestDofIndex ) = false;
                    }
                    // return the traction based on first Piola_Kirchhoff stress
                    return m1PKTestTraction( tTestDofIndex );
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    // if the traction based on second Piola_Kirchhoff stress was not evaluated
                    if( m2PKTestTractionEval( tTestDofIndex ) )
                    {
                        // evaluate the traction based on second Piola_Kirchhoff stress
                        this->eval_testTraction_second_piola_kirchhoff( aNormal, aTestDofTypes );

                        // set bool for evaluation
                        m2PKTestTractionEval( tTestDofIndex ) = false;
                    }
                    // return the traction based on second Piola_Kirchhoff stress
                    return m2PKTestTraction( tTestDofIndex );
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    // if the test traction based on Cauchy stress was not evaluated
                    if( mCauchyTestTractionEval( tTestDofIndex ) )
                    {
                        // evaluate the test traction based on Cauchy stress
                        this->eval_testTraction_cauchy( aNormal, aTestDofTypes );

                        // set bool for evaluation
                        mCauchyTestTractionEval( tTestDofIndex ) = false;
                    }
                    // return the test Cauchy stress tensor
                    return mCauchyTestTraction( tTestDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::testTraction - Unknown stress type.");
                    return mTestTraction( tTestDofIndex );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Struc_Nonlinear_Isotropic::dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const Matrix< DDRMat >&             aJump,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type               aCMFunctionType )
        {
            // get the test dof index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // switch on CM function type
            switch( aCMFunctionType )
            {
                // first Piola-Kirchhoff stress
                case CM_Function_Type::PK1 :
                {
                    // if the derivative of the test traction based on first Piola_Kirchhoff stress was not evaluated
                    if( md1PKTestTractionduEval( tTestDofIndex, tDofIndex ) )
                    {
                        // evaluate the test traction based on first Piola_Kirchhoff stress
                        this->eval_dTestTractiondDOF_first_piola_kirchhoff( aDofTypes, aNormal, aJump, aTestDofTypes );

                        // set bool for evaluation
                        md1PKTestTractionduEval( tTestDofIndex, tDofIndex ) = false;
                    }
                    // return the traction based on first Piola_Kirchhoff stress
                    return md1PKTestTractiondu( tTestDofIndex )( tDofIndex );
                }
                // second Piola-Kirchhoff stress
                case CM_Function_Type::PK2 :
                {
                    // if the derivative of the test traction based on second Piola_Kirchhoff stress was not evaluated
                    if( md2PKTestTractionduEval( tTestDofIndex, tDofIndex ) )
                    {
                        // evaluate the traction based on second Piola_Kirchhoff stress
                        this->eval_dTestTractiondDOF_second_piola_kirchhoff( aDofTypes, aNormal, aJump, aTestDofTypes );

                        // set bool for evaluation
                        md2PKTestTractionduEval( tTestDofIndex, tDofIndex ) = false;
                    }
                    // return the traction based on second Piola_Kirchhoff stress
                    return md2PKTestTractiondu( tTestDofIndex )( tDofIndex );
                }
                // cauchy stress
                case CM_Function_Type::CAUCHY :
                {
                    // if the test traction based on Cauchy stress was not evaluated
                    if( mdCauchyTestTractionduEval( tTestDofIndex,tDofIndex ) )
                    {
                        // evaluate the test traction based on Cauchy stress
                        this->eval_dTestTractiondDOF_cauchy( aDofTypes, aNormal, aJump, aTestDofTypes );

                        // set bool for evaluation
                        mdCauchyTestTractionduEval( tTestDofIndex,tDofIndex ) = false;
                    }
                    // return the test Cauchy stress tensor
                    return mdCauchyTestTractiondu( tTestDofIndex )( tDofIndex );
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Nonlinear_Isotropic::dTestTractiondDOF - Unknown stress type.");
                    return mdTestTractiondDof( tTestDofIndex )( tDofIndex );
                }
            }
        }

        //        void CM_Struc_Nonlinear_Isotropic::eval_dTeststraindDOF_2D()
        //        {
        //            // get displacement field interpolator
        //            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );
        //
        //            // compute displacement gradient
        //            const Matrix< DDRMat > & tdnNdxn = tFIDispl->dnNdxn( 1 );
        //
        //            // get number of bases for displacement
        //            uint tNumBases = tFIDispl->get_number_of_space_time_bases();
        //
        //            // build the derivative of the test strain
        //            mdTeststraindDof.set_size( 4, tNumBases * 2, 0.0 );
        //            mdTeststraindDof( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        //
        //            mdTeststraindDof( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        //        }
        //
        //        void CM_Struc_Nonlinear_Isotropic::eval_dTeststraindDOF_3D()
        //        {
        //            // get displacement field interpolator
        //            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );
        //
        //            // compute displacement gradient
        //            const Matrix< DDRMat > & tdnNdxn = tFIDispl->dnNdxn( 1 );
        //
        //            // get number of bases for displacement
        //            uint tNumBases = tFIDispl->get_number_of_space_time_bases();
        //
        //            // build the derivative of the test strain
        //            mdTeststraindDof.set_size( 9, tNumBases * 3, 0.0 );
        //            mdTeststraindDof( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 1, 1 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 2, 2 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        //
        //            mdTeststraindDof( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 4, 4 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        //
        //            mdTeststraindDof( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        //            mdTeststraindDof( { 8, 8 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        //        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_flux_proj_sym(
                enum CM_Function_Type aCMFunctionType)
        {
            // get the 2nd PK stress in voigt notation
            Matrix< DDRMat > tSymStressVoigt;
            tSymStressVoigt = this->flux( aCMFunctionType );

            // Project the flux
            this->proj_sym(tSymStressVoigt,mFluxProj);
        }

        void CM_Struc_Nonlinear_Isotropic::eval_flux_proj_nsym(
                enum CM_Function_Type aCMFunctionType)
        {
            // get the 2nd PK stress in voigt notation
            Matrix< DDRMat > tNSymStressVoigt;
            tNSymStressVoigt = this->flux( aCMFunctionType );

            // Project the flux
            this->proj_nsym(tNSymStressVoigt,mFluxProj);
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::flatten_normal_2d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            // num cols based on number of flux terms
            aFlatNormal.set_size( 2, mConst.n_rows(), 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, mConst.n_rows() - 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, mConst.n_rows() - 1 ) = aNormal( 0, 0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::flatten_normal_3d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
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
        CM_Struc_Nonlinear_Isotropic ::flatten_normal_nonsym_2d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            // set size based on dimensions
            aFlatNormal.set_size( mSpaceDim, 4 );

            // populate flattened normal
            aFlatNormal = {
                    { aNormal(0,0), aNormal(1,0), 0, 0 },
                    { 0, 0, aNormal(0,0), aNormal(1,0) }
            };
        }

        void
        CM_Struc_Nonlinear_Isotropic ::flatten_normal_nonsym_3d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            // set size based on dimensions
            aFlatNormal.set_size( mSpaceDim, 9, 0.0 );

            aFlatNormal = {
                    { aNormal(0,0), aNormal(1,0), aNormal(2,0), 0, 0, 0, 0, 0, 0 },
                    { 0, 0, 0, aNormal(0,0), aNormal(1,0), aNormal(2,0), 0, 0, 0 },
                    { 0, 0, 0, 0, 0, 0, aNormal(0,0), aNormal(1,0), aNormal(2,0) }
            };
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::proj_sym_2d(
                const Matrix< DDRMat > & aVector,
                Matrix< DDRMat >       & aProjMatrix )
        {
            // representing in voigt notation vector->matrix
            Matrix< DDRMat > tMatrix;
            this->voigt_to_full_sym_stress( aVector, tMatrix );

            // projecting the matrix as ProjMatrix = [aMatrix 0 ; 0 aMatrix]
            aProjMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( { mSpaceDim, 2 * mSpaceDim - 1 }, { mSpaceDim, 2 * mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
        }

        void CM_Struc_Nonlinear_Isotropic::proj_sym_3d(
                const Matrix< DDRMat > & aVector,
                Matrix< DDRMat >       & aProjMatrix )
        {
            // representing in voigt notation vector->matrix
            Matrix< DDRMat > tMatrix;
            this->voigt_to_full_sym_stress( aVector, tMatrix );

            // projecting the matrix as ProjMatrix = [aMatrix 0 0 ; 0 aMatrix 0 ; 0 0 aMatrix]
            aProjMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( { mSpaceDim, 2 * mSpaceDim - 1 }, { mSpaceDim, 2 * mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( {2*mSpaceDim, 3*mSpaceDim-1}, {2*mSpaceDim, 3*mSpaceDim-1}) = tMatrix({0, mSpaceDim-1},{0, mSpaceDim-1});
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::proj_nsym_2d(
                const Matrix< DDRMat > & aVector,
                Matrix< DDRMat >       & aProjMatrix )
        {
            // representing in voigt notation vector->matrix
            Matrix< DDRMat > tMatrix;
            this->voigt_to_full_nonsym( aVector, tMatrix );

            // projecting the matrix as ProjMatrix = [aMatrix 0 ; 0 aMatrix]
            aProjMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( { mSpaceDim, 2 * mSpaceDim - 1 }, { mSpaceDim, 2 * mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
        }

        void CM_Struc_Nonlinear_Isotropic::proj_nsym_3d(
                const Matrix< DDRMat > & aVector,
                Matrix< DDRMat >       & aProjMatrix )
        {
            // representing in voigt notation vector->matrix
            Matrix< DDRMat > tMatrix;
            this->voigt_to_full_nonsym( aVector, tMatrix );

            // projecting the matrix as ProjMatrix = [aMatrix 0 0 ; 0 aMatrix 0 ; 0 0 aMatrix]
            aProjMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( { mSpaceDim, 2 * mSpaceDim - 1 }, { mSpaceDim, 2 * mSpaceDim - 1 } ) = tMatrix( { 0, mSpaceDim - 1 }, { 0, mSpaceDim - 1 } );
            aProjMatrix( {2*mSpaceDim, 3*mSpaceDim-1}, {2*mSpaceDim, 3*mSpaceDim-1}) = tMatrix({0, mSpaceDim-1},{0, mSpaceDim-1});
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_space_dim( uint aSpaceDim )
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

        void CM_Struc_Nonlinear_Isotropic::set_model_type( Model_Type aModelType )
        {
            // fixme: currently cannot set plane type and a tensor type at the same time from an input file
            // store model type based on input
            if ( aModelType == Model_Type::PLANE_STRESS or aModelType == Model_Type::PLANE_STRAIN )
            {
                mPlaneType = aModelType;
            }
            else if ( aModelType == Model_Type::FULL or aModelType == Model_Type::HYDROSTATIC or aModelType == Model_Type::DEVIATORIC)
            {
                mTensorType = aModelType;
            }
            else
            {
                MORIS_ASSERT( false,
                        "CM_Struc_Nonlinear_Isotropic::set_model_type - Specified Nonlinear isotropic elasticity model type doesn't exist." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::voigt_to_full_nonsym(
                const Matrix< DDRMat > & aVoigtVector,
                Matrix< DDRMat >       & aFullMatrix )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtNonSymMap.n_rows();
            uint tNumCol = mVoigtNonSymMap.n_cols();

            // set size for full matrix
            aFullMatrix.set_size( tNumRow, tNumCol );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < tNumCol; jCol++ )
                {
                    // fill the full matrix
                    aFullMatrix( iRow, jCol ) = aVoigtVector( mVoigtNonSymMap( iRow, jCol ) );
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::voigt_to_full_sym_strain(
                const Matrix< DDRMat > & aVoigtVector,
                Matrix< DDRMat >       & aFullMatrix )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtSymMap.n_rows();
            uint tNumCol = mVoigtSymMap.n_cols();

            // set size for full matrix
            aFullMatrix.set_size( tNumRow, tNumCol );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < tNumCol; jCol++ )
                {
                    if( iRow == jCol )
                    {
                        // fill the full matrix
                        aFullMatrix( iRow, jCol ) = aVoigtVector( mVoigtSymMap( iRow, jCol ) );
                    }
                    else
                    {
                        // fill the full matrix
                        aFullMatrix( iRow, jCol ) = 0.5 * aVoigtVector( mVoigtSymMap( iRow, jCol ) );
                    }
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::voigt_to_full_sym_stress(
                const Matrix< DDRMat > & aVoigtVector,
                Matrix< DDRMat >       & aFullMatrix )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtSymMap.n_rows();
            uint tNumCol = mVoigtSymMap.n_cols();

            // set size for full matrix
            aFullMatrix.set_size( tNumRow, tNumCol );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < tNumCol; jCol++ )
                {
                    // fill the full matrix
                    aFullMatrix( iRow, jCol ) = aVoigtVector( mVoigtSymMap( iRow, jCol ) );
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::full_to_voigt_nonsym(
                const Matrix< DDRMat > & aFullMatrix,
                Matrix< DDRMat > & aVoigtVector )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtNonSymMap.n_rows();
            uint tNumCol = mVoigtNonSymMap.n_cols();

            // set size for Voigt vector
            aVoigtVector.set_size( tNumRow * tNumCol, 1 );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < tNumCol; jCol++ )
                {
                    // fill the Voigt vector
                    aVoigtVector( mVoigtNonSymMap( iRow, jCol ) ) = aFullMatrix( iRow, jCol );
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::full_to_voigt_sym_strain(
                const Matrix< DDRMat > & aFullMatrix,
                Matrix< DDRMat > & aVoigtVector )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtSymMap.n_rows();
            uint tNumCol = mVoigtSymMap.n_cols();

            // set size for Voigt vector
            aVoigtVector.set_size( mConst.n_cols(), 1 );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = iRow; jCol < tNumCol; jCol++ )
                {
                    if( iRow == jCol )
                    {
                        // fill the Voigt vector
                        aVoigtVector( mVoigtSymMap( iRow, jCol ) ) = aFullMatrix( iRow, jCol );
                    }
                    else
                    {
                        // fill the Voigt vector
                        aVoigtVector( mVoigtSymMap( iRow, jCol ) ) = 2.0 * aFullMatrix( iRow, jCol );
                    }
                }
            }
        }

        void CM_Struc_Nonlinear_Isotropic::full_to_voigt_sym_stress(
                const Matrix< DDRMat > & aFullMatrix,
                Matrix< DDRMat > & aVoigtVector )
        {
            // get the matrix dimensions
            uint tNumRow = mVoigtSymMap.n_rows();
            uint tNumCol = mVoigtSymMap.n_cols();

            // set size for Voigt vector
            aVoigtVector.set_size( mConst.n_cols(), 1 );

            // loop over the rows
            for ( uint iRow = 0; iRow < tNumRow; iRow++ )
            {
                // loop over the columns
                for ( uint jCol = iRow; jCol < tNumCol; jCol++ )
                {
                    // fill the Voigt vector
                    aVoigtVector( mVoigtSymMap( iRow, jCol ) ) = aFullMatrix( iRow, jCol );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
