/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Cluster_Measure.cpp
 *
 */

#include "cl_FEM_Cluster_Measure.hpp"
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Cluster_Measure::val()
        {
            MORIS_ERROR( mMEAEval, "Cluster_Measure::val -- Cluster measure was not initiated.");

            // return the cluster measure value
            return mMEAVal;
        }

        //------------------------------------------------------------------------------

        void Cluster_Measure::eval_cluster_measure()
        {
            switch( mMeasureType )
            {
                case fem::Measure_Type::CELL_MEASURE :
                {
                    // evaluate cell measure from the cluster
                    mMEAVal = mCluster->compute_cluster_cell_measure(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mMEAEval = true;

                    break;
                }
                case fem::Measure_Type::CELL_SIDE_MEASURE :
                {
                    // evaluate cell side measure from the cluster
                    mMEAVal = mCluster->compute_cluster_cell_side_measure(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mMEAEval = true;

                    break;
                }
                case fem::Measure_Type::CELL_LENGTH_MEASURE :
                {
                    // evaluate cell length measure from the cluster
                    mMEAVal = mCluster->compute_cluster_cell_length_measure(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mMEAEval = true;

                    break;
                }
                default:
                    MORIS_ERROR( false, " Cluster_Measure::eval_cluster_measure - No measure type specified. " );
            }
        }

        //------------------------------------------------------------------------------

        void Cluster_Measure::perturb_cluster_measure( moris::real aDeltaH )
        {
            // check that the cluste measure was initiated
            MORIS_ERROR( mMEAEval, "Cluster_Measure::perturb_cluster_measure -- Cluster measure was not initiated.");

            // update the cluster measure with perturbation
            mMEAVal += {{aDeltaH}};
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Cluster_Measure::dMEAdPDV()
        {
            MORIS_ERROR( mdMEAdPDVEval, "Cluster_Measure::dMEAdPDV -- Cluster measure derivatives were not evaluated.");

            // return the cluster measure derivatives value
            return mdMEAdPDV;
        }

        //------------------------------------------------------------------------------

        void Cluster_Measure::eval_cluster_measure_derivatives()
        {
            switch( mMeasureType )
            {
                case fem::Measure_Type::CELL_MEASURE :
                {
                    // evaluate cell measure from the cluster
                    mdMEAdPDV = mCluster->compute_cluster_cell_measure_derivative(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mdMEAdPDVEval = true;
                    break;
                }
                case fem::Measure_Type::CELL_SIDE_MEASURE :
                {
                    // evaluate cell side measure from the cluster
                    mdMEAdPDV = mCluster->compute_cluster_cell_side_measure_derivative(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mdMEAdPDVEval = true;
                    break;
                }
                case fem::Measure_Type::CELL_LENGTH_MEASURE :
                {
                    // evaluate cell length measure from the cluster
                    mdMEAdPDV = mCluster->compute_cluster_cell_length_measure_derivative(
                            mIsPrimary,
                            mIsMaster );

                    // set evaluation flag to true
                    mdMEAdPDVEval = true;
                    break;
                }
                default:
                    MORIS_ERROR( false, " Cluster_Measure::eval_cluster_measure_derivatives - No measure type specified. " );
            }
        }

        //------------------------------------------------------------------------------
    }/*end_fem_namespace */
}/*end_moris_namespace */

