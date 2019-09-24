
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_stress( Matrix< DDRMat > & aStress )
        {
            // compute conductivity matrix
            Matrix< DDRMat > K;
            this->eval_const( K );

            // compute flux
            aStress = K * mFieldInterpolators( 0 )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_strain( Matrix< DDRMat > & aStrain )
        {
            // compute temp gradient
            aStrain = mFieldInterpolators( 0 )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_const( Matrix< DDRMat > & aConst )
        {
            // compute conductivity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );
            aConst = mProperties( 0 )->val()( 0 ) * I;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStressdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                              Matrix< DDRMat >             & adStressdDOF )
        {

            // if direct dependency on the dof type
            if( static_cast< uint >( aDofTypes( 0 ) ) < mDofTypeMap.numel() && mDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) ) != -1 )
            {
                // compute conductivity matrix
                Matrix< DDRMat > I;
                eye( mSpaceDim, mSpaceDim, I );
                Matrix< DDRMat > K = mProperties( 0 )->val()( 0 ) * I;

                // compute derivative with direct dependency
                adStressdDOF = K * mFieldInterpolators( 0 )->dnNdxn( 1 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // init matrix size
                if( adStressdDOF.numel() < 1 )
                {
                    uint tFIIndex = mProperties( 0 )->get_dof_type_map()( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                    Field_Interpolator* tFI = mProperties( 0 )->get_field_interpolators()( tFIIndex ) ;

                    adStressdDOF.set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );
                }

                // compute derivative with indirect dependency through properties
                adStressdDOF.matrix_data() += mFieldInterpolators( 0 )->gradx( 1 ) * mProperties( 0 )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                              Matrix< DDRMat >             & adStraindDOF )
        {

            // if direct dependency on the dof type
            if( static_cast< uint >( aDofTypes( 0 ) ) < mDofTypeMap.numel() && mDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) ) != -1 )
            {
                // compute derivative with direct dependency
                adStraindDOF = mFieldInterpolators( 0 )->dnNdxn( 1 );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dConstdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                             Matrix< DDRMat >             & adConstdDOF )
        {
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
