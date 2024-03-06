/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_L2_Error_Analytic.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_L2_Error_Analytic::IQI_L2_Error_Analytic()
        {
            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::L2_ERROR_ANALYTIC;
            init_property("L2Check", IQI_Property_Type::L2_CHECK);
        }

        //------------------------------------------------------------------------------

        void
        IQI_L2_Error_Analytic::compute_QI( Matrix< DDRMat >& aQI )
        {
            // check that mQuantityDofType is set
            MORIS_ERROR( mQuantityDofType.size() > 0,
                    "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

            // get field interpolator
            Field_Interpolator* tFI = get_leader_fi_manager()->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            real tJumpNorm;

            Field_Interpolator* tFIField = nullptr;
            if ( get_field_type_list().size() != 0 )
            {
                tFIField = get_leader_fi_manager()->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

                tJumpNorm = norm( tFI->val() - tFIField->val() );
            }
            else
            {
                // get analytical solution property
               std::shared_ptr< Property > const &tPropL2Check = get_leader_property(IQI_Property_Type::L2_CHECK);
                tJumpNorm = norm( tFI->val() - tPropL2Check->val() );
            }

            // evaluate the QI
            aQI = tJumpNorm * tJumpNorm;
        }

        //------------------------------------------------------------------------------

        void
        IQI_L2_Error_Analytic::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // check that mQuantityDofType is set
            MORIS_ERROR( mQuantityDofType.size() > 0,
                    "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

            // get field interpolator
            Field_Interpolator* tFI = get_leader_fi_manager()->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // get analytical solution property
           std::shared_ptr< Property > const &tPropL2Check = get_leader_property(IQI_Property_Type::L2_CHECK);

            // compute jump
            real tJumpNorm = norm( tFI->val() - tPropL2Check->val() );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tJumpNorm * tJumpNorm );
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris

