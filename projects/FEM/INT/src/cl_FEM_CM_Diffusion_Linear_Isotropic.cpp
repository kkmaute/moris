
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_flux()
        {
            // compute flux
            mFlux = this->constitutive() * mDofFI( 0 )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_strain()
        {
            // compute strain
            mStrain = mDofFI( 0 )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = mProperties( 0 )->val()( 0 ) * I;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mFluxDofDer( tDofIndex ) = this->constitutive() * mDofFI( 0 )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mFluxDofDer( tDofIndex ).set_size( mSpaceDim, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mFluxDofDer( tDofIndex ).matrix_data()
                += mDofFI( 0 )->gradx( 1 ) * mProperties( 0 )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mStrainDofDer( tDofIndex ) = mDofFI( 0 )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mStrainDofDer( tDofIndex ).set_size( mSpaceDim, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dConstdDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // reset matrix
            mConstDofDer( tDofIndex ).set_size( 1, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mConstDofDer( tDofIndex ) = mProperties( 0 )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDV( moris::Cell< MSI::Dv_Type > aDvTypes )
        {
            // get the dv type
            uint tDvType = static_cast< uint >( aDvTypes( 0 ) );

            // get the dv type index
            uint tDvIndex = mGlobalDvTypeMap( tDvType );

            // if direct dependency on the dv type
            if( tDvType < mDofTypeMap.numel() && mDofTypeMap( tDvType ) != -1 )
            {
                // compute derivative with direct dependency
                // fixme
            }
            else
            {
                // reset matrix
                mFluxDvDer( tDvIndex ).set_size( mSpaceDim, mDvFI( tDvIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dv_dependency( aDvTypes ) )
            {
                // compute derivative with indirect dependency through properties
                // fixme
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDV( moris::Cell< MSI::Dv_Type > aDvTypes )
        {
            // get the dv type
            uint tDvType = static_cast< uint >( aDvTypes( 0 ) );

            // get the dv type index
            uint tDvIndex = mGlobalDvTypeMap( tDvType );

            // if direct dependency on the dv type
            if( tDvIndex < mDofTypeMap.numel() && mDofTypeMap( tDvIndex ) != -1 )
            {
                // compute derivative with direct dependency
                // fixme
            }
            else
            {
                // reset matrix
                mStrainDvDer( tDvIndex ).set_size( mSpaceDim, mDvFI( tDvIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dConstdDV( moris::Cell< MSI::Dv_Type > aDvTypes )
        {
            // get the dv type
            uint tDvType = static_cast< uint >( aDvTypes( 0 ) );

            // get the dv type index
            uint tDvIndex = mGlobalDvTypeMap( tDvType );

            // reset matrix
            mConstDvDer( tDvIndex ).set_size( 1, mDvFI( tDvIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( 0 )->check_dv_dependency( aDvTypes ) )
            {
                // compute derivative with indirect dependency through properties
                // fixme
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
