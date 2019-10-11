
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute flux
            Matrix< DDRMat > tFlux;
            mMasterCM( 0 )->eval_flux( tFlux );

            Matrix< DDRMat > tTestStrain;
            mMasterCM( 0 )->eval_test_strain( tTestStrain );

            Matrix< DDRMat >tShapeFunc2D( 2, 8, 0.0 );
            tShapeFunc2D( 0, 0 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFunc2D( 0, 1 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFunc2D( 0, 2 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFunc2D( 0, 3 ) = mMasterFI( 0 )->N()( 0, 3 );
            tShapeFunc2D( 1, 4 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFunc2D( 1, 5 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFunc2D( 1, 6 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFunc2D( 1, 7 ) = mMasterFI( 0 )->N()( 0, 3 );

//            print(tShapeFunc2D, "tShapeFunc2D");


            // compute the residual
            aResidual( 0 ) = trans( tTestStrain ) * tFlux;
                           //- trans( tShapeFunc2D ) * mMasterProp( 0 )->val()( 0 );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // if property has dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute flux derivative
                    Matrix< DDRMat > tdStressDOF;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdStressDOF );

                    Matrix< DDRMat > tTestStrain;
                    mMasterCM( 0 )->eval_test_strain( tTestStrain );

//                    print(tdStressDOF,"tStress");

                    // compute the jacobian

                    aJacobians( 0 )( iDOF ).matrix_data() += trans( tTestStrain ) * tdStressDOF;

                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
