/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Power_Dissipation_Bulk.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Power_Dissipation_Bulk.hpp"
#include "cl_FEM_CM_Fluid_Turbulence.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Power_Dissipation_Bulk::IQI_Power_Dissipation_Bulk()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::POWER_DISSIPATION_BULK;
            init_property("Select", IQI_Property_Type::SELECT);
            init_constitutive_model("Fluid", IQI_Constitutive_Type::FLUID);
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation_Bulk::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid = get_leader_constitutive_model(IQI_Constitutive_Type::FLUID);

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // get select property
            const std::shared_ptr< Property > & tPropSelect = get_leader_property(IQI_Property_Type::SELECT);

            // grab select property value
            real tSelect = 1.0;
            if( tPropSelect != nullptr )
            {
                tSelect = tPropSelect->val()( 0 );
            }

            // if select property value positive
            if ( tSelect > 0.0 )
            {
                // FIXME protect dof type
                // get velocity FI
                Field_Interpolator * tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

                // FIXME protect dof type
                // get the pressure FI
                Field_Interpolator* tFIPressure = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::P );

                // create identity matrix
                uint tSpaceDim = tFIVelocity->get_space_dim();
                Matrix< DDRMat > tI( tSpaceDim, 1, 1.0 );

                // evaluate pressure contribution to flux
                Matrix< DDRMat > tP( ( tSpaceDim - 1 ) * 3, 1, 0.0 );
                tP( { 0, tSpaceDim - 1 }, { 0, 0 } ) = tI * tFIPressure->val();

                // evaluate the QI
                aQI = tSelect * ( trans( tCMFluid->flux() + tP ) * tCMFluid->strain() + //
                        tPropDensity->val()( 0 )* tFIVelocity->val_trans() * tFIVelocity->val() );
            }
            else
            {
                aQI = 0.0;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation_Bulk::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // create matrix storage for IQI
            Matrix<DDRMat> tQI( 1, 1 );

            // compute QI
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation_Bulk::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid = get_leader_constitutive_model(IQI_Constitutive_Type::FLUID);

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity = tCMFluid->get_property( "Density" );

            // get select property
            const std::shared_ptr< Property > & tPropSelect = get_leader_property(IQI_Property_Type::SELECT);

            // grab select property value
            real tSelect = 1.0;
            if( tPropSelect != nullptr )
            {
                tSelect = tPropSelect->val()( 0 );
            }

            // if select property value positive
            if( tSelect > 0.0 )
            {
                // FIXME protect dof type
                // get velocity field interpolator
                Field_Interpolator * tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

                // FIXME protect dof type
                // get velocity field interpolator
                Field_Interpolator * tFIPressure = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::P );

                // get the number of leader dof type dependencies
                uint tNumDofDependencies = get_requested_leader_dof_types().size();

                // compute dQIdu for indirect dof dependencies
                for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
                {
                    // get the treated dof type
                    const Vector< MSI::Dof_Type > & tDofType = get_requested_leader_dof_types()( iDof );

                    // get leader index for residual dof type, indices for assembly
                    uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                    // if dof type is velocity
                    // FIXME protect dof type
                    if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex },
                                { 0, 0 } ) += aWStar * tSelect * //
                                2.0 * tPropDensity->val()( 0 ) * tFIVelocity->N_trans() * tFIVelocity->val();
                    }

                    // if dof type is pressure
                    // FIXME protect dof type
                    if ( tDofType( 0 ) == MSI::Dof_Type::P )
                    {
                        // get space dimension
                        uint tSpaceDim = tFIVelocity->get_space_dim();

                        // create identity matrix
                        Matrix< DDRMat > tI( tSpaceDim, 1, 1.0 );
                        Matrix< DDRMat > tII( ( tSpaceDim - 1 ) * 3, 1, 0.0 );
                        tII( { 0, tSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                        // add contribution to dQIdu
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex },
                                { 0, 0 } ) += aWStar * tSelect * //
                                trans( tII * tFIPressure->N() ) * tCMFluid->strain();
                    }

                    // if density property depends on dof type
                    if( tPropDensity->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex },
                                { 0, 0 } ) += aWStar * tSelect * //
                                tFIVelocity->val_trans() * tFIVelocity->val() * tPropDensity->dPropdDOF( tDofType );
                    }

                    // if fluid CM depends on dof type
                    if ( tCMFluid->check_dof_dependency( tDofType ) )
                    {
                        // create identity matrix
                        uint tSpaceDim = tFIVelocity->get_space_dim();
                        Matrix< DDRMat > tI( tSpaceDim, 1, 1.0 );

                        // evaluate pressure contribution to flux
                        Matrix< DDRMat > tP( ( tSpaceDim - 1 ) * 3, 1, 0.0 );
                        tP( { 0, tSpaceDim - 1 }, { 0, 0 } ) = tI * tFIPressure->val();

                        // add contribution to dQIdu
                        mSet->get_residual()( tQIIndex )(
                                { tLeaderDepStartIndex, tLeaderDepStopIndex },
                                { 0, 0 } ) += aWStar * tSelect * //
                                ( trans( tCMFluid->dFluxdDOF( tDofType ) ) * tCMFluid->strain()
                                        + trans( tCMFluid->dStraindDOF( tDofType ) ) * ( tCMFluid->flux() + tP ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation_Bulk::compute_dQIdu(
                Vector< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR( false, "IQI_Power_Dissipation_Bulk::compute_dQIdu - Not implemented.");
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

