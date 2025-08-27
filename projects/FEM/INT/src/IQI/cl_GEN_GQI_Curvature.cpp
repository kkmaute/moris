/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_GQI_Curvature.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_GEN_GQI_Curvature.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    GQI_Curvature::GQI_Curvature()
            : GQI()
    {
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::initialize()
    {
        if ( !mIsInitialized )
        {
            // set initialize flag to true
            mIsInitialized = true;
        }
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::compute_QI()
    {
        // initialize if needed
        this->initialize();

        // Have geometry compute GQI
        real tVal = mDesignVariableInterface->get_requested_GQI( mGeometryName, IQI_Type::CURVATURE );

        Matrix< DDRMat > tMat( 1, 1, tVal );

        this->evaluate_QI( tMat );
    }

    //------------------------------------------------------------------------------

    real
    GQI_Curvature::get_QI()
    {
        // initialize if needed
        this->initialize();

        Matrix< DDRMat > tMat( 1, 1 );

        this->evaluate_QI( tMat );

        return tMat( 0, 0 );
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::compute_QI( Matrix< DDRMat >& aQI )
    {
        // initialize if needed
        this->initialize();

        // evaluate QI
        this->evaluate_QI( aQI );
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::evaluate_QI( Matrix< DDRMat >& aMat )
    {
        MORIS_ERROR( false, "GQI_Curvature::evaluate_QI - Not implemented yet." );    // brendan todo

        aMat = { { 0.0 } };
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "GQI_Curvature::compute_dQIdu - Not implemented yet." );    // brendan todo
    }

    //------------------------------------------------------------------------------

    void
    GQI_Curvature::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "GQI_Curvature::compute_dQIdu - Not implemented yet." );    // brendan todo
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
