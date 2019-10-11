
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_flux( Matrix< DDRMat > & aFlux )
        {
            // compute conductivity matrix
            Matrix< DDRMat > K;
            this->eval_const( K );

            Matrix< DDRMat > tStrain;
            this->eval_strain( tStrain );

            // compute flux
            aFlux = K * tStrain;
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_strain( Matrix< DDRMat > & aStrain )
        {
            // compute temp gradient
            Matrix< DDRMat > tGradx;
            tGradx = mFieldInterpolators( 0 )->gradx( 1 );

            aStrain.set_size( 3, 1 , 0.0 );

            aStrain( 0, 0 ) = tGradx( 0, 0 );
            aStrain( 1, 0 ) = tGradx( 1, 1 );
            aStrain( 2, 0 ) = tGradx( 1, 0 ) + tGradx( 0, 1 );
        }

//------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic::eval_test_strain( Matrix< DDRMat > & aTestStrain )
        {
            // compute temp gradient
            Matrix< DDRMat > tdnNdxn;
            tdnNdxn = mFieldInterpolators( 0 )->dnNdxn( 1 );

            aTestStrain.set_size( 3, 8 , 0.0 );
            aTestStrain( {0,0},{0,3} ) = mFieldInterpolators( 0 )->dnNdxn( 1 )({0,0},{0,3});
            aTestStrain( {2,2},{0,3} ) = mFieldInterpolators( 0 )->dnNdxn( 1 )({1,1},{0,3});

            aTestStrain( {1,1},{4,7} ) = mFieldInterpolators( 0 )->dnNdxn( 1 )({1,1},{0,3});
            aTestStrain( {2,2},{4,7} ) = mFieldInterpolators( 0 )->dnNdxn( 1 )({0,0},{0,3});



//            aTestStrain( 0, 0 ) = tdnNdxn( 0, 0 );
//            aTestStrain( 2, 0 ) = tdnNdxn( 1, 0 );
//            aTestStrain( 0, 1 ) = tdnNdxn( 0, 1 );
//            aTestStrain( 2, 1 ) = tdnNdxn( 1, 1 );
//            aTestStrain( 0, 2 ) = tdnNdxn( 0, 2 );
//            aTestStrain( 2, 2 ) = tdnNdxn( 1, 2 );
//            aTestStrain( 0, 3 ) = tdnNdxn( 0, 3 );
//            aTestStrain( 2, 3 ) = tdnNdxn( 1, 3 );
//
//            aTestStrain( 1, 4 ) = tdnNdxn( 1, 0 );
//            aTestStrain( 2, 4 ) = tdnNdxn( 0, 0 );
//            aTestStrain( 1, 5 ) = tdnNdxn( 1, 1 );
//            aTestStrain( 2, 5 ) = tdnNdxn( 0, 1 );
//            aTestStrain( 1, 6 ) = tdnNdxn( 1, 2 );
//            aTestStrain( 2, 6 ) = tdnNdxn( 0, 2 );
//            aTestStrain( 1, 7 ) = tdnNdxn( 1, 3 );
//            aTestStrain( 2, 7 ) = tdnNdxn( 0, 3 );


//            aTestStrain( 0, 0 ) = tdnNdxn( 0, 0 );
//            aTestStrain( 2, 0 ) = tdnNdxn( 1, 0 );
//            aTestStrain( 1, 1 ) = tdnNdxn( 1, 0 );
//            aTestStrain( 2, 1 ) = tdnNdxn( 0, 0 );
//
//            aTestStrain( 0, 2 ) = tdnNdxn( 0, 1 );
//            aTestStrain( 2, 2 ) = tdnNdxn( 1, 1 );
//            aTestStrain( 1, 3 ) = tdnNdxn( 1, 1 );
//            aTestStrain( 2, 3 ) = tdnNdxn( 0, 1 );
//
//            aTestStrain( 0, 4 ) = tdnNdxn( 0, 2 );
//            aTestStrain( 2, 4 ) = tdnNdxn( 1, 2 );
//            aTestStrain( 1, 5 ) = tdnNdxn( 1, 2 );
//            aTestStrain( 2, 5 ) = tdnNdxn( 0, 2 );
//
//            aTestStrain( 0, 6 ) = tdnNdxn( 0, 3 );
//            aTestStrain( 2, 6 ) = tdnNdxn( 1, 3 );
//            aTestStrain( 1, 7 ) = tdnNdxn( 1, 3 );
//            aTestStrain( 2, 7 ) = tdnNdxn( 0, 3 );

//            print (aTestStrain,"aTestStrain");
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_const( Matrix< DDRMat > & aConst )
        {
            // compute conductivity matrix
            aConst.set_size(mSpaceDim+1, mSpaceDim+1, 0.0 );

            aConst( 0, 0 ) = 1;
            aConst( 1, 1 ) = 1;
            aConst( 0, 1 ) = mProperties( 1 )->val()( 0 );
            aConst( 1, 0 ) = mProperties( 1 )->val()( 0 );
            aConst( 2, 2 ) =  0.5 * (1-mProperties( 1 )->val()( 0 ) );

            aConst = mProperties( 0 )->val()( 0 ) / (1 - std::pow(mProperties( 1 )->val()( 0 ), 2))  * aConst;

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
        void CM_Struc_Linear_Isotropic::eval_dFluxdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                        Matrix< DDRMat >             & adFluxdDOF )
        {
            // if direct dependency on the dof type
            if( static_cast< uint >( aDofTypes( 0 ) ) < mDofTypeMap.numel() && mDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) ) != -1 )
            {
//                Matrix< DDRMat > KK;
//                this-> eval_flux(KK);
                // compute conductivity matrix
                Matrix< DDRMat > K;
                this->eval_const( K );

                Matrix< DDRMat > tTestStrain;
                this->eval_test_strain( tTestStrain );

//                print(tTestStrain, "tTestStrain");

                // compute derivative with direct dependency
                adFluxdDOF = K * tTestStrain;
            }

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // init matrix size
                if( adFluxdDOF.numel() < 1 )
                {
                    uint tFIIndex = mProperties( 0 )->get_dof_type_map()( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                    Field_Interpolator* tFI = mProperties( 0 )->get_field_interpolators()( tFIIndex ) ;

                    adFluxdDOF.set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );
                }

                // compute derivative with indirect dependency through properties
                adFluxdDOF.matrix_data() += mFieldInterpolators( 0 )->gradx( 1 ) * mProperties( 0 )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dStraindDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                              Matrix< DDRMat >         & adStraindDOF )
        {
//            // get dof index
//            uint tDOFIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );
//
//            // init dConstdDOF
//            adStraindDOF.set_size( mSpaceDim, mFieldInterpolators( tDOFIndex )->get_number_of_space_time_coefficients(), 0.0 );
//
//            // if direct dependency on the dof type
//            if( static_cast< uint >( aDofTypes( 0 ) ) < mDofTypeMap.numel() && mDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) ) != -1 )
//            {
//                // compute derivative with direct dependency
//                adStraindDOF = mFieldInterpolators( 0 )->dnNdxn( 1 );
//            }
        }

//------------------------------------------------------------------------------
        void CM_Struc_Linear_Isotropic::eval_dConstdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                             Matrix< DDRMat >             & adConstdDOF )
        {
            // get dof index
            uint tDOFIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // init dConstdDOF
            adConstdDOF.set_size( 1, mFieldInterpolators( tDOFIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                adConstdDOF = mProperties( 0 )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
