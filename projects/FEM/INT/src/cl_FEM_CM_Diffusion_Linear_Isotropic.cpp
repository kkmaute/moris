
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

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
        void CM_Diffusion_Linear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_testTraction( const Matrix< DDRMat > & aNormal )
        {
            // compute test traction
            mTestTraction = trans( mDofFI( 0 )->dnNdxn( 1 ) ) * this->constitutive() * aNormal;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_strain()
        {
            // compute strain
            mStrain = mDofFI( 0 )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_testStrain()
        {
            // compute test strain
            mTestStrain = mDofFI( 0 )->dnNdxn( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 ) * I;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) = this->constitutive() * mDofFI( 0 )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mdFluxdDof( tDofIndex ).set_size( mSpaceDim, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ).matrix_data() += mDofFI( 0 )->gradx( 1 ) * mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute derivative
            mdTractiondDof( tDofIndex ) =  trans( aNormal ) * this->dFluxdDOF( aDofTypes );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                    const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute derivative
            mdTestTractiondDof( tDofIndex ) = trans( mDofFI( 0 )->dnNdxn( 1 ) ) * aNormal * this->dConstdDOF( aDofTypes );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = mDofFI( 0 )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mdStraindDof( tDofIndex ).set_size( mSpaceDim, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size( 1, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dFluxdDV - This function is not implemented.");
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dStraindDV - This function is not implemented.");
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dConstdDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dConstdDV - This function is not implemented.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
