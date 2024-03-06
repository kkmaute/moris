/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Property.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Property.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Property::IQI_Property()
        {
            init_property( "Property", IQI_Property_Type::PROPERTY );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Property::compute_QI( Matrix< DDRMat > &aQI )
        {
            // get index of property if defined; otherwise set to 0
            uint tTypeIndex = mIQITypeIndex != -1 ? mIQITypeIndex : 0;

            // evaluate the QI
            aQI = get_leader_property( IQI_Property_Type::PROPERTY )->val()( tTypeIndex );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Property::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the property
            std::shared_ptr< Property > const &tProperty = get_leader_property( IQI_Property_Type::PROPERTY );

            // check that pointer to property exists
            MORIS_ASSERT( tProperty,
                    "IQI_Property::compute_QI - property does not exist.\n" );

            // get index of property if defined; otherwise set to 0
            uint tTypeIndex = mIQITypeIndex != -1 ? mIQITypeIndex : 0;

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tProperty->val()( tTypeIndex ) );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Property::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the property
            std::shared_ptr< Property > const &tProperty = get_leader_property( IQI_Property_Type::PROPERTY );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();

            // compute dQIdu for indirect dof dependencies
            for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // if property depends on dof type
                if ( tProperty->check_dof_dependency( tDofType ) )
                {
                    // build selection matrix
                    uint tNumVecFieldComps = tProperty->val().numel();

                    Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                    // if no index is specified it will return the whole dqdu
                    if ( mIQITypeIndex == -1 )
                    {
                        tSelect.fill( 1.0 );
                    }
                    else
                    {
                        // if an index is specified it will only return dQ(index)d
                        tSelect( mIQITypeIndex, 0 ) = 1.0;
                    }

                    MORIS_ASSERT( tProperty->dPropdDOF( tDofType ).n_cols() == tLeaderDepStopIndex - tLeaderDepStartIndex + 1,
                            "IQI_Property::compute_dQIdu - Incorrect size of dof derivative of property" );

                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                            aWStar * trans( tProperty->dPropdDOF( tDofType ) ) * tSelect;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Property::compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu )
        {
            // get the property
            std::shared_ptr< Property > const &tProperty = get_leader_property( IQI_Property_Type::PROPERTY );

            // if property depends on dof type
            if ( tProperty->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = trans( tProperty->dPropdDOF( aDofType ) );
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
