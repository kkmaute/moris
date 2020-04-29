
//FEM/INT/src
#include "cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        SP_Velocity_Dirichlet_Nitsche::SP_Velocity_Dirichlet_Nitsche()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ] = Property_Type::VISCOSITY;
            mPropertyMap[ "Density" ]   = Property_Type::DENSITY;

            // populate the dof map (default)
            mMasterDofMap[ "Velocity" ] = MSI::Dof_Type::VX;
        }

//------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::eval_SP()
        {
            // get the viscosity and density property
            std::shared_ptr< Property > tPropViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tPropDensity   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

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

            // compute stabilization parameter value
            mPPVal = mParameters( 0 )
                   * (   tPropViscosity->val()( 0 ) / mElementSize
                       + tPropDensity->val()( 0 ) * tInfinityNorm / 6.0 );
        }

//------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity and density property
            std::shared_ptr< Property > tPropViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tPropDensity   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

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

            // FIXME infinity norm
            // if dof type == velocity
            if( aDofTypes( 0 ) == mMasterDofMap[ "Velocity" ] )
            {
            	Matrix< DDRMat > tdInfinityNormdu;
                // compute derivative of infinity norm
            	if( tInfinityNorm == 0.0 )
            	{
                    tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex );
            	}
            	else
            	{
                    tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex ) * tVelocityFI->val()( tInfinityNormIndex ) / tInfinityNorm;
            	}

                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 )( 0 ) * tPropDensity->val()( 0 ) * tdInfinityNormdu / 6.0;
            }

            // if viscosity depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 )( 0 ) * tPropViscosity->dPropdDOF( aDofTypes ) / mElementSize;
            }

            // if density depends on dof type
            if( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 )( 0 ) * tPropDensity->dPropdDOF( aDofTypes ) * tInfinityNorm / 6.0;
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


