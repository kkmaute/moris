/*
 * cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Drag_Lift_Coefficient::IQI_Drag_Lift_Coefficient( sint aBeta )
        {
            // fill master dof map (default)
            mMasterDofMap[ "Velocity" ] = MSI::Dof_Type::VX;
            mMasterDofMap[ "Pressure" ] = MSI::Dof_Type::P;

            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]     = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "Viscosity" ]   = static_cast< uint >( Property_Type::VISCOSITY );
            mPropertyMap[ "VelocityMax" ] = static_cast< uint >( Property_Type::VELOCITY_MAX );
            mPropertyMap[ "Diameter" ]    = static_cast< uint >( Property_Type::DIAMETER );

            // set mBeta
            mBeta = aBeta;
        }

        //------------------------------------------------------------------------------

        void IQI_Drag_Lift_Coefficient::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap["Velocity"] );

            // get the pressure FI
            Field_Interpolator * tFIPressure =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap["Pressure"] );

            // get the density property value
            real tDensity = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );

            // get the viscosity property value
            real tViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) )->val()( 0 );

            // get the maximum velocity property value
            real tUMax = mMasterProp( static_cast< uint >( Property_Type::VELOCITY_MAX ) )->val()( 0 );

            // get the diameter property value
            real tDiameter = mMasterProp( static_cast< uint >( Property_Type::DIAMETER ) )->val()( 0 );

            // compute dvxdy - dvydx
            real tVorticity = tFIVelocity->gradx( 1 )( 1, 0 ) - tFIVelocity->gradx( 1 )( 0, 1 );

            // compute deno = 0.5 rho U² D
            real tDeno = 0.5 * tDensity * std::pow( tUMax, 2.0 ) * tDiameter;

            // compute QI
            aQI = {{ ( mBeta * tViscosity * tVorticity * mNormal( 1 ) - tFIPressure->val()( 0 ) * mNormal( 0 ) ) / tDeno }};
        }

        //------------------------------------------------------------------------------

        void IQI_Drag_Lift_Coefficient::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap["Velocity"] );

            // get the pressure FI
            Field_Interpolator * tFIPressure =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap["Pressure"] );

            // get the density property value
            real tDensity = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );

            // get the viscosity property value
            real tViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) )->val()( 0 );

            // get the maximum velocity property value
            real tUMax = mMasterProp( static_cast< uint >( Property_Type::VELOCITY_MAX ) )->val()( 0 );

            // get the diameter property value
            real tDiameter = mMasterProp( static_cast< uint >( Property_Type::DIAMETER ) )->val()( 0 );

            // compute dvxdy - dvydx
            real tVorticity = tFIVelocity->gradx( 1 )( 1, 0 ) - tFIVelocity->gradx( 1 )( 0, 1 );

            // compute deno = 0.5 rho U² D
            real tDeno = 0.5 * tDensity * std::pow( tUMax, 2.0 ) * tDiameter;

            // compute QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    ( mBeta * tViscosity * tVorticity * mNormal( 1 ) - tFIPressure->val()( 0 ) * mNormal( 0 ) ) / tDeno );
        }

        //------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


