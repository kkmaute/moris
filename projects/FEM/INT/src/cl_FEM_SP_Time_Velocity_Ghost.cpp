
//FEM/INT/src
#include "cl_FEM_SP_Time_Velocity_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void SP_Time_Velocity_Ghost::eval_SP()
        {
            // get the density property
            std::shared_ptr< Property > tDensityProp = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute time increment deltat
            Matrix< DDRMat > tTimeCoeff = mMasterFIManager->get_IP_geometry_interpolator()->get_time_coeff();
            real tDeltaT = tTimeCoeff.max() - tTimeCoeff.min();

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * std::pow( mElementSize, 2.0 * ( mOrder - 1.0 ) + 3.0 )
                   * tDensityProp->val()( 0 )  / ( mParameters( 1 )( 0 ) * tDeltaT );
        }

//------------------------------------------------------------------------------
        void SP_Time_Velocity_Ghost::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the density property
            std::shared_ptr< Property > tDensityProp = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute time increment deltat
            Matrix< DDRMat > tTimeCoeff = mMasterFIManager->get_IP_geometry_interpolator()->get_time_coeff();
            real tDeltaT = tTimeCoeff.max() - tTimeCoeff.min();

            // if density depends on dof type
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += mParameters( 0 ) * std::pow( mElementSize, 2.0 * ( mOrder - 1.0 ) + 3.0 )
                * tDensityProp->dPropdDOF( aDofTypes ) / ( mParameters( 1 )( 0 ) * tDeltaT );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


