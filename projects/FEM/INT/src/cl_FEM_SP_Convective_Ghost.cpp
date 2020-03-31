
//FEM/INT/src
#include "cl_FEM_SP_Convective_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void SP_Convective_Ghost::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator* tVelocityFI
            = mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Velocity" ] );

            // get the density property
            std::shared_ptr< Property > tDensityProp = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get absolute value of u.n
            Matrix< DDRMat > tAbs = trans( tVelocityFI->val() ) * mNormal;
            real tAbsReal = std::abs( tAbs( 0 ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tDensityProp->val()( 0 ) * tAbsReal * std::pow( mElementSize, 2.0 );
        }

//------------------------------------------------------------------------------
        void SP_Convective_Ghost::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI
            = mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Velocity" ] );

            // get absolute value of u.n
            Matrix< DDRMat > tAbs = trans( tVelocityFI->val() ) * mNormal;
            real tAbsReal = std::abs( tAbs( 0 ) );

            // get the density property
            std::shared_ptr< Property > tDensityProp = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mMasterDofMap[ "Velocity" ] )
            {
                // get the sign of u.n
                Matrix< DDRMat > tSign = tAbs / tAbsReal;

                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 ) * std::pow( mElementSize, 2.0 ) * tDensityProp->val()( 0 ) * tSign * trans( mNormal ) * tVelocityFI->N();
            }

            // if density depends on dof
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 ) * std::pow( mElementSize, 2.0 ) * tAbsReal * tDensityProp->dPropdDOF( aDofTypes );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


