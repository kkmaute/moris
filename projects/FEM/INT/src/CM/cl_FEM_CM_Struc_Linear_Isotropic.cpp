
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        CM_Struc_Linear_Isotropic::CM_Struc_Linear_Isotropic()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "YoungsModulus" ]        = static_cast< uint >( CM_Property_Type::EMOD );
            mPropertyMap[ "PoissonRatio" ]         = static_cast< uint >( CM_Property_Type::NU );
            mPropertyMap[ "CTE" ]                  = static_cast< uint >( CM_Property_Type::CTE );
            mPropertyMap[ "ReferenceTemperature" ] = static_cast< uint >( CM_Property_Type::TEMP_REF );
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::set_dof_type_list(
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
                // if temperature dof type string
                else if( tDofString == "Temperature" )
                {
                    mDofTemp = tDofType;
                }
                // if pressure dof type string
                else if( tDofString == "Pressure" )
                {
                    mDofPressure = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Struc_Linear_Isotropic::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::set_local_properties()
        {
            // set the Young's modulus property
            mPropEMod = get_property( "YoungsModulus" );

            // set the Poisson ratio property
            mPropPoisson = get_property( "PoissonRatio" );

            // set the CTE property
            mPropCTE = get_property( "CTE" );

            // set the reference temperature property
            mPropTRef = get_property( "ReferenceTemperature" );
        }

        //------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::set_function_pointers()
        {
            switch( mSpaceDim )
            {
                case 2 :
                {
                    m_eval_strain       = &CM_Struc_Linear_Isotropic::eval_strain_2d;
                    m_eval_teststrain   = &CM_Struc_Linear_Isotropic::eval_teststrain_2d;
                    m_flatten_normal    = &CM_Struc_Linear_Isotropic::flatten_normal_2d;

                    mStrain.set_size( 3, 1, 0.0 );

                    switch( mPlaneType )
                    {
                        case Model_Type::PLANE_STRESS :
                        {
                            m_eval_inv_bulk_modulus = &CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus_plane_stress;

                            switch( mTensorType )
                            {
                                case Model_Type::FULL :
                                {
                                    mConstFunc = &CM_Struc_Linear_Isotropic::full_plane_stress;
                                    mConst.set_size( 3, 3, 0.0 );
                                    break;
                                }
                                case Model_Type::DEVIATORIC :
                                {
                                    mConstFunc = &CM_Struc_Linear_Isotropic::deviatoric_plane_stress;
                                    mConst.set_size( 3, 3, 0.0 );
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
                            switch( mTensorType )
                            {
                                case Model_Type::FULL :
                                {
                                    mConstFunc = &CM_Struc_Linear_Isotropic::full_plane_strain;
                                    mConst.set_size( 4, 3, 0.0 );
                                    break;
                                }
                                case Model_Type::DEVIATORIC :
                                {
                                    mConstFunc = &CM_Struc_Linear_Isotropic::deviatoric_plane_strain;
                                    mConst.set_size( 4, 3, 0.0 );
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
                            MORIS_ERROR(false, "Linear isotropic elasticity in 2d requires plane stress or plane strain models");
                        }
                    }
                    break;
                }
                case 3 :
                {
                    m_eval_strain       = &CM_Struc_Linear_Isotropic::eval_strain_3d;
                    m_eval_teststrain   = &CM_Struc_Linear_Isotropic::eval_teststrain_3d;
                    m_flatten_normal    = &CM_Struc_Linear_Isotropic::flatten_normal_3d;

                    mStrain.set_size( 6, 1, 0.0 );


                    switch(mTensorType)
                    {
                        case Model_Type::FULL :
                        {
                            mConstFunc = &CM_Struc_Linear_Isotropic::full_3d;
                            mConst.set_size( 6, 6, 0.0 );
                            break;
                        }
                        case Model_Type::DEVIATORIC :
                        {
                            mConstFunc = &CM_Struc_Linear_Isotropic::deviatoric_3d;
                            mConst.set_size( 6, 6, 0.0 );
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
                    MORIS_ERROR( false, "Linear isotropic elasticity implemented only for 2d and 3d" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_flux()
        {
            // compute flux
            mFlux = this->constitutive() * this->strain();

            // if pressure dof
            if ( mDofPressure != MSI::Dof_Type::UNDEFINED )
            {
                // get the pressure FI
                Field_Interpolator * tPressureFI =
                        mFIManager->get_field_interpolators_for_type( mDofPressure );

                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );

                // evaluate pressure contribution to flux
                Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * tPressureFI->val();

                // add contribution to the flux
                mFlux -= tP;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mTraction = tFlatNormal * this->flux();
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_testTraction(
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

            // if test traction wrt pressure
            if( aTestDofTypes( 0 ) == mDofPressure )
            {
                // compute test traction wrt pressure
                mTestTraction( tTestDofIndex ) = tFlatNormal * this->dFluxdDOF( aTestDofTypes ) ;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_strain_2d()
        {
            // get the displacement spatial gradient from displacement FI
            Matrix< DDRMat > tDisplGradx =
                    mFIManager->get_field_interpolators_for_type( mDofDispl )->gradx( 1 );

            // evaluate the strain
            mStrain.fill( 0.0 );
            mStrain( 0, 0 ) = tDisplGradx( 0, 0 );
            mStrain( 1, 0 ) = tDisplGradx( 1, 1 );
            mStrain( 2, 0 ) = tDisplGradx( 1, 0 ) + tDisplGradx( 0, 1 );

            // if thermal expansion
            if ( mPropCTE != nullptr )
            {
                // build thermal expansion vector
                Matrix< DDRMat > tThermalExpansionVector( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                tThermalExpansionVector( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

                // get temperature field interpolator
                Field_Interpolator* tFITemp =
                        mFIManager->get_field_interpolators_for_type( mDofTemp );

                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - tFITemp->val() );
            }
        }

        void CM_Struc_Linear_Isotropic::eval_strain_3d()
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

            // if thermal expansion
            if ( mPropCTE != nullptr )
            {
                // build thermal expansion vector
                Matrix< DDRMat > tThermalExpansionVector( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                tThermalExpansionVector( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

                // get temperature field interpolator
                Field_Interpolator* tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemp );

                // add thermal contribution to the strain
                mStrain += tThermalExpansionVector * ( mPropTRef->val() - tFITemp->val() );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_teststrain_2d()
        {
            // get displacement field interpolator
            Field_Interpolator * tFIDispl = mFIManager->get_field_interpolators_for_type( mDofDispl );

            // compute displacement gradient
            const Matrix< DDRMat > & tdnNdxn = tFIDispl->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases = tFIDispl->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        void CM_Struc_Linear_Isotropic::eval_teststrain_3d()
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

        void CM_Struc_Linear_Isotropic::eval_const()
        {
            // get the Poisson's ratio value
            moris::real tNu = mPropPoisson->val()( 0 );

            // get the Youngs' modulus value
            moris::real tEmod = mPropEMod->val()( 0 );

            // evaluate the constitutive matrix
            ( this->*mConstFunc )( tEmod, tNu );
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::real CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus()
        {
            // get Poisson ratio value
            moris::real tNu = mPropPoisson->val()( 0 );

            // get elasticity modulus value
            moris::real tEMod = mPropEMod->val()( 0 );

            // init inverse of the bulk modulus
            moris::real tInvBulkModulus;

            // evaluate inverse of the bulk modulus
            ( this->*m_eval_inv_bulk_modulus )( tNu, tEMod, tInvBulkModulus );

            // return
            return tInvBulkModulus;
        }

        void CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus_generic(
                moris::real aNu,
                moris::real aEMod,
                moris::real & aInvBulkModulus )
        {
            aInvBulkModulus = 3.0 * ( 1.0 - 2.0 * aNu ) / aEMod;
        }

        void CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus_plane_stress(
                moris::real aNu,
                moris::real aEMod,
                moris::real & aInvBulkModulus )
        {
            aInvBulkModulus = 2.0 * ( 1.0 - aNu ) / aEMod;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::Matrix< DDRMat > CM_Struc_Linear_Isotropic::eval_dInvBulkModulusdDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init inverse of the bulk modulus
            Matrix< DDRMat > tdInvBulkModulusdDOF( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if Young's modulus property depends on dof type
            if( mPropEMod->check_dof_dependency( aDofTypes ) )
            {
                tdInvBulkModulusdDOF -=
                        eval_inv_bulk_modulus() * mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
            }

            // if Poisson ratio property depends on dof type
            if( mPropPoisson->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dInvBulkModulusdDOF - Poisson's ratio depends on dof, not handled." );
            }

            // return
            return tdInvBulkModulusdDOF;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
            if( aDofTypes( 0 ) == mDofDispl || aDofTypes( 0 ) == mDofTemp )
            {
                mdFluxdDof( tDofIndex ) +=
                        this->constitutive() * this->dStraindDOF( aDofTypes );
            }

            // if pressure dof
            if ( aDofTypes( 0 ) == mDofPressure )
            {
                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                // get shape function for presure field
                Matrix< DDRMat > tPressureN = tFI->N();

                // build the dfluxdp
                mdFluxdDof( tDofIndex ) -= tII * tPressureN;
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
                MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_dTractiondDOF(
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

        void CM_Struc_Linear_Isotropic::eval_dTestTractiondDOF(
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
                MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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

            // if temperature dof
            if ( mPropCTE!=nullptr && aDofTypes( 0 ) == mDofTemp )
            {
                // build thermal expansion vector
                Matrix< DDRMat > tThermalExpansionVector( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                tThermalExpansionVector( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * mPropCTE->val();

                // compute derivatives
                mdStraindDof( tDofIndex ) -= tThermalExpansionVector * tFI->N();
            }

            // if thermal expansion
            if ( mPropCTE != nullptr && mPropCTE->check_dof_dependency( aDofTypes ) )
            {
                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                // get temperature field interpolator
                Field_Interpolator* tFIT =
                        mFIManager->get_field_interpolators_for_type( mDofTemp );

                // compute derivatives
                mdStraindDof( tDofIndex ) +=
                        tII * mPropCTE->dPropdDOF( aDofTypes ) *
                        ( mPropTRef->val()( 0 ) - tFIT->val()( 0 ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dConstdDOF - Not implemented." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::flatten_normal_2d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            aFlatNormal.set_size( 2, 3, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
        }

        void CM_Struc_Linear_Isotropic::flatten_normal_3d(
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

        void CM_Struc_Linear_Isotropic::set_space_dim( uint aSpaceDim )
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

        void CM_Struc_Linear_Isotropic::set_model_type( Model_Type aModelType )
        {
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
                        "CM_Struc_Linear_Isotropic::set_model_type - Specified linear isotropic elasticity model type doesn't exist." );
            }

            // set function pointers
            this->set_function_pointers();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::real CM_Struc_Linear_Isotropic::get_e_prime()
        {
            moris::real tEPrime;

            // get elasticity modulus value
            real tE  = mPropEMod->val()( 0 );

            // get Poisson ratio value
            real tNu = mPropPoisson->val()( 0 );

            switch ( mPlaneType )
            {
                case Model_Type::PLANE_STRESS :
                {
                    // Eprime = E
                    tEPrime = tE;
                    break;
                }
                case Model_Type::PLANE_STRAIN :
                {
                    // Eprime = E / ( 1 - nu^2 )
                    tEPrime = tE / ( 1 - ( tNu * tNu ) );
                    break;
                }
                default:
                {
                    tEPrime = tE;
                    MORIS_ASSERT( false, "CM_Struc_Linear_Isotropic::get_e_prime() - unknown model type " );
                    break;
                }
            }
            return tEPrime;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::full_plane_stress(
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

        void CM_Struc_Linear_Isotropic::deviatoric_plane_stress(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / ((1 + aNu) * 2.0);

            mConst( 0, 0 ) =  tPre;
            mConst( 1, 1 ) =  tPre;
            mConst( 0, 1 ) = -tPre;
            mConst( 1, 0 ) = -tPre;
            mConst( 2, 2 ) =  tPre;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::full_plane_strain(
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

        void CM_Struc_Linear_Isotropic::deviatoric_plane_strain(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / (3.0 * (1.0 + aNu ) );

            mConst( 0, 0 ) = tPre * 4.0;
            mConst( 0, 1 ) = tPre;
            mConst( 1, 0 ) = tPre;
            mConst( 1, 1 ) = tPre * 4.0;
            mConst( 2, 0 ) = tPre;
            mConst( 2, 1 ) = tPre;
            mConst( 3, 2 ) = tPre * 3.0 / 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::full_3d(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / (1.0 + aNu ) / (1.0 - 2.0 * aNu );

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

        void CM_Struc_Linear_Isotropic::deviatoric_3d(
                moris::real aEmod,
                moris::real aNu )
        {
            moris::real tPre = aEmod / ( 3.0 * ( 1.0 + aNu ) );

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