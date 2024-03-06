/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Volume_Fraction.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Volume_Fraction.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Volume_Fraction::IQI_Volume_Fraction()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::VOLUME_FRACTION;
            init_stabilization_parameter( "Reciprocal_total_vol", IQI_Stabilization_Type::RECIPROCAL_TOTAL_VOLUME );
        }

        //------------------------------------------------------------------------------

        void IQI_Volume_Fraction::compute_QI( Matrix< DDRMat > &aQI )
        {
            // set the size for the QI
            aQI.set_size( 1, 1 );

            // evaluate the QI
            aQI( 0 ) = get_stabilization_parameter( IQI_Stabilization_Type::RECIPROCAL_TOTAL_VOLUME )->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void IQI_Volume_Fraction::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( get_stabilization_parameter( IQI_Stabilization_Type::RECIPROCAL_TOTAL_VOLUME )->val() );
        }

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
