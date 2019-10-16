
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
            // flatten normal
            Matrix< DDRMat > tNormal( 2, 3, 0.0 );
            tNormal( 0, 0 ) = aNormal( 0,0 );
            tNormal( 0, 2 ) = aNormal( 1,0 );
            tNormal( 1, 1 ) = aNormal( 1,0 );
            tNormal( 1, 2 ) = aNormal( 0,0 );

            // compute traction
            mTraction = tNormal * this->flux();
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_testTraction( const Matrix< DDRMat > & aNormal )
        {
            // flatten normal
            Matrix< DDRMat > tNormal( 2, 3, 0.0 );
            tNormal( 0, 0 ) = aNormal( 0,0 );
            tNormal( 0, 2 ) = aNormal( 1,0 );
            tNormal( 1, 1 ) = aNormal( 1,0 );
            tNormal( 1, 2 ) = aNormal( 0,0 );

            // compute test traction
            mTestTraction = trans( this->testStrain() ) * this->constitutive() * trans( tNormal );
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_strain()
        {
            // compute displacement gradient
            Matrix< DDRMat > tGradx;
            tGradx = mDofFI( 0 )->gradx( 1 );

            // set strain matrix size
            mStrain.set_size( 3, 1 , 0.0 );

            // fill with strain
            mStrain( 0, 0 ) = tGradx( 0, 0 );
            mStrain( 1, 0 ) = tGradx( 1, 1 );
            mStrain( 2, 0 ) = tGradx( 1, 0 ) + tGradx( 0, 1 );
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_testStrain()
        {
            // compute temp gradient
            Matrix< DDRMat > tdnNdxn;
            tdnNdxn = mDofFI( 0 )->dnNdxn( 1 );

            mTestStrain.set_size( 3, 8 , 0.0 );
            mTestStrain( {0,0},{0,3} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,3});
            mTestStrain( {2,2},{0,3} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,3});

            mTestStrain( {1,1},{4,7} ) = mDofFI( 0 )->dnNdxn( 1 )({1,1},{0,3});
            mTestStrain( {2,2},{4,7} ) = mDofFI( 0 )->dnNdxn( 1 )({0,0},{0,3});
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_const()
        {
            // compute conductivity matrix
            mConst.set_size(mSpaceDim+1, mSpaceDim+1, 0.0 );

            mConst( 0, 0 ) = 1;
            mConst( 1, 1 ) = 1;
            mConst( 0, 1 ) = mProperties( 1 )->val()( 0 );
            mConst( 1, 0 ) = mProperties( 1 )->val()( 0 );
            mConst( 2, 2 ) =  0.5 * (1-mProperties( 1 )->val()( 0 ) );

            mConst = mProperties( 0 )->val()( 0 ) / (1 - std::pow(mProperties( 1 )->val()( 0 ), 2))  * mConst;

//            pre = Em/(1-nu*nu);
//
//            mConstit[0] = pre;    // 1st row
//            mConstit[1] = pre*nu;
//            mConstit[2] = 0.0;
//            mConstit[3] = pre*nu; // 2nd row
//            mConstit[4] = pre;
//            mConstit[5] = 0.0;
//            mConstit[6] = 0.0;    // 3rd row
//            mConstit[7] = 0.0;
//            mConstit[8] = pre*(1-nu)/2;


//            pre      = Em/(1-nu*nu);
//            dpre_dnu = 2*Em*nu/(1-nu*nu)/(1-nu*nu);
//            dpre_dx  = dEm/(1-nu*nu);
//
//            mConstit[0] = dpre_dnu*dnu + dpre_dx;    // 1st row
//            mConstit[1] = (dpre_dnu*nu+pre)*dnu + dpre_dx*nu;
//            mConstit[2] = 0.0;
//            mConstit[3] = mConstit[1]; // 2nd row
//            mConstit[4] = mConstit[0];
//            mConstit[5] = 0.0;
//            mConstit[6] = 0.0;    // 3rd row
//            mConstit[7] = 0.0;
//            mConstit[8] = ((dpre_dnu*(1-nu)/2)-(pre/2))*dnu + dpre_dx*(1-nu)/2;
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
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ).matrix_data() += mDofFI( 0 )->gradx( 1 ) * mProperties( 0 )->dPropdDOF( aDofTypes );
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

            // flatten normal
            Matrix< DDRMat > tNormal( 2, 3, 0.0 );
            tNormal( 0, 0 ) = aNormal( 0,0 );
            tNormal( 0, 2 ) = aNormal( 1,0 );
            tNormal( 1, 1 ) = aNormal( 1,0 );
            tNormal( 1, 2 ) = aNormal( 0,0 );

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
    } /* namespace fem */
} /* namespace moris */
