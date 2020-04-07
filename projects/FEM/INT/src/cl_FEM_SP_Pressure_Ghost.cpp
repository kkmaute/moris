
//FEM/INT/src
#include "cl_FEM_SP_Pressure_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void SP_Pressure_Ghost::eval_SP()
        {
            // get the viscosity and density property
            std::shared_ptr< Property > tViscosityProp = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tDensityProp   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the velocity FI
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Velocity" ] );

            // compute infinity norm of u
            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );
            for( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute deltaT
            Matrix< DDRMat > tTimeCoeff = mMasterFIManager->get_IP_geometry_interpolator()->get_time_coeff();
            real tDeltaT = tTimeCoeff.max() - tTimeCoeff.min();

            // compute deltaP
            real tDeltaP = tViscosityProp->val()( 0 ) / mElementSize
                         + tDensityProp->val()( 0 ) * tInfinityNorm / 6.0
                         + tDensityProp->val()( 0 ) * mElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * std::pow( mElementSize, 2 * mOrder ) / tDeltaP;
        }

//------------------------------------------------------------------------------
        void SP_Pressure_Ghost::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity and density property
            std::shared_ptr< Property > tViscosityProp = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tDensityProp   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute infinity norm
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Velocity" ] );
            uint tInfinityNormIndex = 0;
            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );
            for( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNormIndex = iDim;
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute deltaT
            Matrix< DDRMat > tTimeCoeff = mMasterFIManager->get_IP_geometry_interpolator()->get_time_coeff();
            real tDeltaT = tTimeCoeff.max() - tTimeCoeff.min();

            // compute deltaP
            real tDeltaP = tViscosityProp->val()( 0 ) / mElementSize
                         + tDensityProp->val()( 0 ) * tInfinityNorm / 6.0
                         + tDensityProp->val()( 0 ) * mElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT );

            // if dof type == velocity
            if( aDofTypes( 0 ) == mMasterDofMap[ "Velocity" ] )
            {
                // compute derivative of infinity norm
                Matrix< DDRMat > tdInfinityNormdu
                = tFI->N().get_row( tInfinityNormIndex ) * tVelocityFI->val()( tInfinityNormIndex ) / tInfinityNorm;

                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                -= mParameters( 0 )( 0 ) * std::pow( mElementSize, 2 * mOrder )
                 * tDensityProp->val()( 0 ) * tdInfinityNormdu
                 / ( 6.0 * std::pow( tDeltaP, 2.0 ) );
            }

            // if viscosity depends on dof type
            if( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                -= mParameters( 0 )( 0 ) * std::pow( mElementSize, 2 * mOrder )
                 * tViscosityProp->dPropdDOF( aDofTypes )
                 / ( mElementSize * std::pow( tDeltaP, 2.0 ) );
            }

            // if density depends on dof type
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data()
                -= mParameters( 0 )( 0 ) * std::pow( mElementSize, 2 * mOrder )
                 * tDensityProp->dPropdDOF( aDofTypes ) * ( tInfinityNorm / 6.0 + mElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT ) )
                 / std::pow( tDeltaP, 2.0 );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


