
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        CM_Struc_Nonlinear_Isotropic::CM_Struc_Nonlinear_Isotropic()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "YoungsModulus" ]        = static_cast< uint >( CM_Property_Type::EMOD );
            mPropertyMap[ "PoissonRatio" ]         = static_cast< uint >( CM_Property_Type::NU );
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
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

        }

        //------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::set_function_pointers()
        {
            switch( mSpaceDim )
            {
                case 2 :
                {
                    m_eval_strain       = &CM_Struc_Nonlinear_Isotropic::eval_strain_2d;
                    m_eval_teststrain   = &CM_Struc_Nonlinear_Isotropic::eval_teststrain_2d;

                    mStrain.set_size( 3, 1, 0.0 );

                    switch( mPlaneType )
                    {
                        case Model_Type::PLANE_STRAIN:
                        {
                            switch( mTensorType )
                            {
                                case Model_Type::FULL :
                                {
                                    mConstFunc = &CM_Struc_Nonlinear_Isotropic::full_plane_strain;
                                    mConst.set_size( 3, 3, 0.0 );
                                    break;
                                }
                                default:
                                {
                                    MORIS_ERROR(false, "Only full tensors implemented for plane stress");
                                }
                            }
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "Nonlinear isotropic elasticity in 2d requires plane strain models");
                        }
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Nonlinear isotropic elasticity implemented only for 2d" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_flux()
        {
            // compute flux
            mFlux = this->constitutive() * this->strain();

        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mTraction = tFlatNormal * this->flux();
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // if test traction wrt displacement
            if( aTestDofTypes( 0 ) == mDofDispl )
            {
                // compute test traction wrt displacement
                mTestTraction( tTestDofIndex ) = trans( this->testStrain() ) * this->constitutive() * trans( tFlatNormal );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_strain_2d()
        {
            // get the displacement spatial gradient from displacement FI
            Matrix< DDRMat > tDisplGradx =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

            // evaluate the strain
            mStrain.fill( 0.0 );
            mStrain( 0, 0 ) = tDisplGradx( 0, 0 );
            mStrain( 1, 0 ) = tDisplGradx( 1, 1 );
            mStrain( 2, 0 ) = tDisplGradx( 1, 0 ) + tDisplGradx( 0, 1 );
        }

        void CM_Struc_Nonlinear_Isotropic::eval_strain_3d()
        {
            // get the displacement spatial gradient from displacement FI
            Matrix< DDRMat > tDisplGradx =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

            // evaluate the strain
            mStrain.fill( 0.0 );
            mStrain( 0, 0 ) = tDisplGradx( 0, 0 );
            mStrain( 1, 0 ) = tDisplGradx( 1, 1 );
            mStrain( 2, 0 ) = tDisplGradx( 2, 2 );
            mStrain( 3, 0 ) = tDisplGradx( 1, 2 ) + tDisplGradx( 2, 1 );
            mStrain( 4, 0 ) = tDisplGradx( 0, 2 ) + tDisplGradx( 2, 0 );
            mStrain( 5, 0 ) = tDisplGradx( 0, 1 ) + tDisplGradx( 1, 0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_teststrain_2d()
        {
            // // get displacement field interpolator
            // Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // Geometry_Interpolator* tIgGeomInterp = mFIManager->get_IG_geometry_interpolator();

            // // compute displacement gradient
            // const Matrix< DDRMat > & tdnNdxn = tFIDispl->dnNdxn( 1 );

            // // get number of bases for displacement
            // uint tNumBases = tFIDispl->get_number_of_space_time_bases();

            // // Jacobian 
            // const Matrix< DDRMat > &  tJ = tIgGeomInterp->space_jacobian();

            // //dNdx
            // const Matrix< DDRMat > & tdNdxn = tFIDispl->dNdxn(1);

            // // displacement gradient 
            // auto tH = 


            // // displacement gradient field
            // // \frac{\partial x}{\partial u}
            // auto tH = 

            // moris::print(tJ,"tJ");


            // // build the linear test strain
            // mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
            // mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            // mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            // mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            // mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );




        }

        void CM_Struc_Nonlinear_Isotropic::eval_teststrain_3d()
        {
            // get displacement field interpolator
            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // compute displacement gradient
            const Matrix< DDRMat > & tdnNdxn = tFIDispl->dnNdxn( 1 );

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

        void CM_Struc_Nonlinear_Isotropic::eval_const()
        {
            // get the Poisson's ratio value
            moris::real tNu = mPropPoisson->val()( 0 );

            // get the Youngs' modulus value
            moris::real tEmod = mPropEMod->val()( 0 );

            // evaluate the constitutive matrix
            ( this->*mConstFunc )( tEmod, tNu );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            mdFluxdDof( tDofIndex ).set_size(
                    ( mSpaceDim - 1 ) * 3,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            // if displacements or temperature
            if( aDofTypes( 0 ) == mDofDispl )
            {
                mdFluxdDof( tDofIndex ) +=
                        this->constitutive() * this->dStraindDOF( aDofTypes );
            }

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ) +=
                        this->constitutive() *
                        this->strain() *
                        mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
            }

            // if Poisson ratio depends on dof type
            if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute derivative
            mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    0.0 );

            // if elastic modulus depends on dof type
            if ( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        this->testTraction( aNormal, aTestDofTypes ) * trans( aJump ) *
                        mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
            }

            // if Poisson's ratio depends on dof type
            if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdStraindDof
            mdStraindDof( tDofIndex ).set_size(
                    ( mSpaceDim - 1 ) * 3,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            // if displacement dof
            if( aDofTypes( 0 ) == mDofDispl )
            {
                // compute derivative
                mdStraindDof( tDofIndex ) += this->testStrain();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Nonlinear_Isotropic::eval_dConstdDOF - Not implemented." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Nonlinear_Isotropic::flatten_normal_2d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            aFlatNormal.set_size( 2, 3, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
        }

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
            throw;
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
                        "CM_Struc_Nonlinear_Isotropic::set_model_type - Specified linear isotropic elasticity model type doesn't exist." );
            }

            // set function pointers
            this->set_function_pointers();
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Struc_Nonlinear_Isotropic::full_plane_stress(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / ( 1 - std::pow( aNu, 2 ) );

            mConst( 0, 0 ) = tPre;
            mConst( 1, 1 ) = tPre;
            mConst( 0, 1 ) = tPre * aNu;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 2, 2 ) = tPre * 0.5 * (1.0 - aNu );
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Struc_Nonlinear_Isotropic::full_plane_strain(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / (1.0 + aNu ) / (1.0 - 2.0 * aNu ) ;

            mConst( 0, 0 ) = tPre * ( 1.0 - aNu );
            mConst( 0, 1 ) = tPre * aNu;
            mConst( 1, 0 ) = tPre * aNu;
            mConst( 1, 1 ) = tPre * ( 1.0 - aNu );
            mConst( 2, 0 ) = tPre * aNu;
            mConst( 2, 1 ) = tPre * aNu;
            mConst( 3, 2 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
        }
        //--------------------------------------------------------------------------------------------------------------


    } /* namespace fem */
} /* namespace moris */
