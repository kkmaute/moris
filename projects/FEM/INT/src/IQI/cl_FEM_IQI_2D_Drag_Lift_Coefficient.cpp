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

            // get the pressure
            real tPressure = tFIPressure->val()( 0 );

            // get the density property value
            real tDensity = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );

            // get the viscosity property value
            real tViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) )->val()( 0 );

            // get the maximum velocity property value
            real tUMax = mMasterProp( static_cast< uint >( Property_Type::VELOCITY_MAX ) )->val()( 0 );

            // get the diameter property value
            real tDiameter = mMasterProp( static_cast< uint >( Property_Type::DIAMETER ) )->val()( 0 );

            // compute strain-rate
            Matrix< DDRMat > tStrainRate = 0.5 * ( tFIVelocity->gradx( 1 ) + trans( tFIVelocity->gradx( 1 ) ) );

            // compute traction
            Matrix< DDRMat > tTraction = 2.0 * tViscosity * tStrainRate * mNormal - tPressure * mNormal;

            // compute deno = 0.5 rho U² D
            real tDeno = 0.5 * tDensity * std::pow( tUMax, 2.0 ) * tDiameter;

            // compute QI
            if( mBeta == 1 ) // drag coefficient
            {
                aQI = {{ tTraction( 0 ) / tDeno }};
            }
            else if( mBeta == -1 ) // lift coefficient
            {
                aQI = {{ tTraction( 1 ) / tDeno }};
            }
            else
            {
                MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_QI() - unknown mBeta." );
            }
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

            // get the pressure
            real tPressure = tFIPressure->val()( 0 );

            // get the density property value
            real tDensity = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );

            // get the viscosity property value
            real tViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) )->val()( 0 );

            // get the maximum velocity property value
            real tUMax = mMasterProp( static_cast< uint >( Property_Type::VELOCITY_MAX ) )->val()( 0 );

            // get the diameter property value
            real tDiameter = mMasterProp( static_cast< uint >( Property_Type::DIAMETER ) )->val()( 0 );

            // compute strain-rate
            Matrix< DDRMat > tStrainRate = 0.5 * ( tFIVelocity->gradx( 1 ) + trans( tFIVelocity->gradx( 1 ) ) );

            // compute traction
            Matrix< DDRMat > tTraction = 2.0 * tViscosity * tStrainRate * mNormal - tPressure * mNormal;

            // compute deno = 0.5 rho U² D
            real tDeno = 0.5 * tDensity * std::pow( tUMax, 2.0 ) * tDiameter;

            if( mBeta == 1 ) // new drag coefficient
            {
                mSet->get_QI()( tQIIndex ) += aWStar * ( tTraction( 0 ) / tDeno );
            }
            else if( mBeta == -1 ) // new lift coefficient
            {
                mSet->get_QI()( tQIIndex ) += aWStar * ( tTraction( 1 ) / tDeno );
            }
            else
            {
                MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_QI() - unknown mBeta." );
            }
        }

        //------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


