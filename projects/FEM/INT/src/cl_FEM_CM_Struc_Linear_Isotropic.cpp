
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_flux()
        {
            // compute flux
            mFlux = this->constitutive() * this->strain();
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            Matrix< DDRMat > tNormal;

            // flatten normal
            this->flatten_normal( aNormal, tNormal );

            // compute traction
            mTraction = tNormal * this->flux();
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_testTraction( const Matrix< DDRMat > & aNormal )
        {
            Matrix< DDRMat > tNormal;

            // flatten normal
            this->flatten_normal( aNormal, tNormal );

            // compute test traction
            mTestTraction = trans( this->testStrain() ) * this->constitutive() * trans( tNormal );
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_strain()
        {
            // compute displacement gradient
            Matrix< DDRMat > tGradx;
            tGradx = mDofFI( 0 )->gradx( 1 );

            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    mStrain.set_size( 3, 1 , 0.0 );
                    mStrain( 0, 0 ) = tGradx( 0, 0 );
                    mStrain( 1, 0 ) = tGradx( 1, 1 );
                    mStrain( 2, 0 ) = tGradx( 1, 0 ) + tGradx( 0, 1 );
                     break;
                }
                case( 3 ):
                {
                    mStrain.set_size( 6, 1 , 0.0 );
                    mStrain( 0, 0 ) = tGradx( 0, 0 );
                    mStrain( 1, 0 ) = tGradx( 1, 1 );
                    mStrain( 2, 0 ) = tGradx( 1, 1 );
                    mStrain( 2, 0 ) = tGradx( 1, 2 ) + tGradx( 2, 1 );
                    mStrain( 2, 0 ) = tGradx( 0, 2 ) + tGradx( 2, 0 );
                    mStrain( 2, 0 ) = tGradx( 0, 1 ) + tGradx( 1, 0 );
                    break;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Linear_Isotropic::eval_strain - Flattening of strain tensor only implemented in 2 and 3 D");
                }
            }
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_testStrain()
        {
            // compute temp gradient
            Matrix< DDRMat > tdnNdxn;
            tdnNdxn = mDofFI( 0 )->dnNdxn( 1 );

            uint tNumSpaceTimeBasis = mDofFI( 0 )->get_number_of_space_time_bases();

            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                     mTestStrain.set_size( 3, tNumSpaceTimeBasis * mSpaceDim, 0.0 );
                     mTestStrain( {0,0},{0,tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,tNumSpaceTimeBasis-1});
                     mTestStrain( {2,2},{0,tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,tNumSpaceTimeBasis-1});

                     mTestStrain( {1,1},{tNumSpaceTimeBasis,2*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,tNumSpaceTimeBasis-1});
                     mTestStrain( {2,2},{tNumSpaceTimeBasis,2*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,tNumSpaceTimeBasis-1});
                     break;
                }
                case( 3 ):
                {
                    mTestStrain.set_size( 6, tNumSpaceTimeBasis * mSpaceDim , 0.0 );
                    mTestStrain( {0,0},{0,tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {3,3},{0,tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {5,5},{0,tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({2,2},{0,tNumSpaceTimeBasis-1});

                    mTestStrain( {1,1},{tNumSpaceTimeBasis,2*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {3,3},{tNumSpaceTimeBasis,2*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {4,4},{tNumSpaceTimeBasis,2*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({2,2},{0,tNumSpaceTimeBasis-1});

                    mTestStrain( {2,2},{2*tNumSpaceTimeBasis,3*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({2,2},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {4,4},{2*tNumSpaceTimeBasis,3*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,tNumSpaceTimeBasis-1});
                    mTestStrain( {5,5},{2*tNumSpaceTimeBasis,3*tNumSpaceTimeBasis-1} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,tNumSpaceTimeBasis-1});
                    break;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Linear_Isotropic::eval_strain - Flattening of strain tensor only implemented in 2 and 3 D");
                }
            }
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_const()
        {
            moris::real tNu = mProperties( static_cast< uint >( Property_Type::NU ) )->val()( 0 );

            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    uint tPlainStrainOrStress = 1;

                    switch ( tPlainStrainOrStress )
                    {
                        case ( 1 ):
                        {
                            moris::real tPre = mProperties( static_cast< uint >( Property_Type::EMOD ) )->val()( 0 ) / (1 - std::pow( tNu, 2));

                            // compute conductivity matrix
                            mConst.set_size( 3, 3, 0.0 );
                            mConst( 0, 0 ) = tPre;
                            mConst( 1, 1 ) = tPre;
                            mConst( 0, 1 ) = tPre * tNu;
                            mConst( 1, 0 ) = tPre * tNu;
                            mConst( 2, 2 ) = tPre * 0.5 * (1.0 - tNu );
                            break;
                        }
                        case ( 2 ):
                        {
                            moris::real tPre = mProperties( static_cast< uint >( Property_Type::EMOD ) )->val()( 0 ) / (1.0 + tNu ) / (1.0 - 2.0 * tNu ) ;

                            mConst.set_size( 4, 4, 0.0 );
                            mConst( 0, 0 ) = tPre * ( 1.0 - tNu );
                            mConst( 0, 1 ) = tPre * tNu;
                            mConst( 0, 2 ) = tPre * tNu;

                            mConst( 1, 0 ) = tPre * tNu;
                            mConst( 1, 1 ) = tPre * ( 1.0 - tNu );
                            mConst( 1, 2 ) = tPre * tNu;

                            mConst( 2, 0 ) = tPre * tNu;
                            mConst( 2, 1 ) = tPre * tNu;
                            mConst( 2, 2 ) = tPre * ( 1.0 - tNu );

                            mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR(false, "CM_Struc_Linear_Isotropic::eval_const - In 2D only plain strain or plain stress implemented");
                        }
                }
                break;
            }
            case( 3 ):
            {
                moris::real tPre = mProperties( static_cast< uint >( Property_Type::EMOD ) )->val()( 0 ) / (1.0 + tNu ) / (1.0 - 2.0 * tNu ) ;

                mConst.set_size( 6, 6, 0.0 );
                mConst( 0, 0 ) = tPre * ( 1.0 - tNu );
                mConst( 0, 1 ) = tPre * tNu;
                mConst( 0, 2 ) = tPre * tNu;

                mConst( 1, 0 ) = tPre * tNu;
                mConst( 1, 1 ) = tPre * ( 1.0 - tNu );
                mConst( 1, 2 ) = tPre * tNu;

                mConst( 2, 0 ) = tPre * tNu;
                mConst( 2, 1 ) = tPre * tNu;
                mConst( 2, 2 ) = tPre * ( 1.0 - tNu );

                mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
                mConst( 4, 4 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
                mConst( 5, 5 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
                break;
            }
            default:
            {
                MORIS_ERROR(false, "CM_Struc_Linear_Isotropic::eval_strain - Flattening of strain tensor only implemented in 2 and 3 D");
            }
        }
    }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) = this->constitutive() * this->testStrain();
            }
            else
            {
                // reset the matrix
                mdFluxdDof( tDofIndex ).set_size( 3, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::EMOD ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ).matrix_data() += mDofFI( 0 )->gradx( 1 ) * mProperties( static_cast< uint >( Property_Type::EMOD ) )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                            const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            Matrix< DDRMat > tNormal;

            // flatten normal
            this->flatten_normal( aNormal, tNormal );

            // compute derivative
            mdTractiondDof( tDofIndex ) = tNormal * this->dFluxdDOF( aDofTypes );
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                const Matrix< DDRMat >             & aNormal )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dTestTractiondDOF - Not implemented.");
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = this->testStrain();
            }
            else
            {
                // reset the matrix
                mdStraindDof( tDofIndex ).set_size( 3, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dConstdDOF - Not implemented." );
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::flatten_normal( const Matrix< DDRMat > & aNormal,
                                                              Matrix< DDRMat > & aFlattenedNormal)
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    aFlattenedNormal.set_size( 2, 3, 0.0 );
                    aFlattenedNormal( 0, 0 ) = aNormal( 0,0 );
                    aFlattenedNormal( 0, 2 ) = aNormal( 1,0 );
                    aFlattenedNormal( 1, 1 ) = aNormal( 1,0 );
                    aFlattenedNormal( 1, 2 ) = aNormal( 0,0 );
                     break;
                }
                case( 3 ):
                {
                    aFlattenedNormal.set_size( 3, 6, 0.0 );
                    aFlattenedNormal( 0, 0 ) = aNormal( 0,0 );
                    aFlattenedNormal( 1, 1 ) = aNormal( 1,0 );
                    aFlattenedNormal( 2, 2 ) = aNormal( 2,0 );
                    aFlattenedNormal( 0, 4 ) = aNormal( 2,0 );
                    aFlattenedNormal( 0, 5 ) = aNormal( 1,0 );
                    aFlattenedNormal( 1, 3 ) = aNormal( 2,0 );
                    aFlattenedNormal( 1, 5 ) = aNormal( 0,0 );
                    aFlattenedNormal( 2, 3 ) = aNormal( 1,0 );
                    aFlattenedNormal( 2, 4 ) = aNormal( 0,0 );
                    break;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Linear_Isotropic::flatten_normal - Flattening for normal only implemented in 2 and 3 D");
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
