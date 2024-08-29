/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Measure.cpp
 *
 */

#include "cl_FEM_SP_Measure.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Model.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    SP_Measure::SP_Measure() {}

    //------------------------------------------------------------------------------

    void SP_Measure::set_cluster_measure_type_list(
            Vector< std::tuple<
                    fem::Measure_Type,
                    mtk::Primary_Void,
                    mtk::Leader_Follower > > &aClusterMeasureTuples,
            Vector< std::string >            &aClusterMeasureStrings )
    {
        // loop over provided cluster measure strings
        for ( uint iCMEA = 0; iCMEA < aClusterMeasureStrings.size(); iCMEA++ )
        {
            // if string is element size
            if ( aClusterMeasureStrings( iCMEA ) == "ElementSize" )
            {
                // set tuple for element size
                mElementSizeTuple = aClusterMeasureTuples( iCMEA );
            }
            else
            {
                // error unknown cluster measure string
                MORIS_ERROR( false,
                        "SP_Measure::set_cluster_measure_type_list - Unknown aClusterMeasureString : %s \n",
                        aClusterMeasureStrings( iCMEA ).c_str() );
            }
        }
    }

    //------------------------------------------------------------------------------

    Vector< std::tuple<
            fem::Measure_Type,
            mtk::Primary_Void,
            mtk::Leader_Follower > >
    SP_Measure::get_cluster_measure_tuple_list()
    {
        return { mElementSizeTuple };
    }

    //------------------------------------------------------------------------------

    void SP_Measure::eval_SP()
    {
        // compute stabilization parameter value
        mPPVal = mCluster->get_cluster_measure(
                                 std::get< 0 >( mElementSizeTuple ),
                                 std::get< 1 >( mElementSizeTuple ),
                                 std::get< 2 >( mElementSizeTuple ) )
                         ->val();
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
