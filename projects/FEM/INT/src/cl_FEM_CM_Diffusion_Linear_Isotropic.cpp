
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        CM_Diffusion_Linear_Isotropic::CM_Diffusion_Linear_Isotropic()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ] = Property_Type::CONDUCTIVITY;

            // populate dof map
            mDofMap[ "Temp" ] = MSI::Dof_Type::TEMP;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_flux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute flux
            mFlux = tPropConductivity->val()( 0 )
                  * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_divflux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute the divergence of the flux
            mDivFlux = tPropConductivity->val() * this->divstrain();
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_testTraction( const Matrix< DDRMat > & aNormal,
                                                               const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute test traction
            mTestTraction( tTestDofIndex ) = trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) )
                                           * tPropConductivity->val()( 0 ) * aNormal;
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_strain()
        {
            // compute strain
            mStrain = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_divstrain()
        {
            // get the temperature gradient
            Matrix< DDRMat > tTempGrad
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 2 );

            // evaluate the divergence of the strain
            mDivStrain = sum( tTempGrad( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_testStrain()
        {
            // compute test strain
            mTestStrain = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
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

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // init the matrix
            mdFluxdDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ).matrix_data()
                += tPropConductivity->val()( 0 )
                 * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ).matrix_data()
                += mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 )
                 * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if temperature dof
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ).matrix_data() += tPropConductivity->val()( 0 ) * this->ddivstraindu( aDofTypes );
            }

            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mddivfluxdu( tDofIndex ).matrix_data() += this->divstrain() * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tTempd2Ndx2 = tFI->dnNdxn( 2 );

                // fill ddivstrain/du
                mddivstraindu( tDofIndex ) = tTempd2Ndx2.get_row( 0 ) + tTempd2Ndx2.get_row( 1 );

                if( tTempd2Ndx2.n_rows() == 6 )
                {
                    mddivstraindu( tDofIndex ).matrix_data() += tTempd2Ndx2.get_row( 2 );
                }
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
            mdTractiondDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes );
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                    const Matrix< DDRMat >             & aNormal,
                                                                    const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                                                                       mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).matrix_data()
                += trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) )
                 * aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        void CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                    const Matrix< DDRMat >             & aNormal,
                                                                    const Matrix< DDRMat >             & aJump,
                                                                    const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                                                                       mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).matrix_data()
                += trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) )
                 * aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mdStraindDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );
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
            mdConstdDof( tDofIndex ).set_size( 1, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dFluxdDV - This function is not implemented.");
        }

//------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dStraindDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dStraindDV - This function is not implemented.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
