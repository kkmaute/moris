/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Cluster_Measure.hpp
 *
 */

#pragma once

// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LNA/src
#include "cl_Matrix.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include <tuple>
#include <utility>

namespace moris::fem
{
    class Cluster;
    class Set;

    /**
     * Cluster Measure
     * This class provides cluster measure mainly for stabilization purpose.
     */
    class Cluster_Measure
    {
      public:
        using ClusterMeasureSpecification = std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower >;

        Cluster_Measure(
                Cluster_Measure::ClusterMeasureSpecification const &aClusterMeasureSpecification,
                fem::Cluster *const                                 aCluster )
                : mCluster( aCluster )
                , mMeasureType( std::get< 0 >( aClusterMeasureSpecification ) )
                , mIsPrimary( std::get< 1 >( aClusterMeasureSpecification ) )
                , mIsLeader( std::get< 2 >( aClusterMeasureSpecification ) ){};

        /**
         * trivial constructor for test purpose only
         */
        Cluster_Measure()
                : mMEAEval( true )
        {
            mMEAVal = { { 1.0 } };
        };

        /**
         * set cluster
         * @param[ in ] aCluster a fem cluster pointer
         */
        void set_cluster( fem::Cluster *const aCluster )
        {
            mCluster = aCluster;
            this->eval_cluster_measure();
        }

        /**
         * get measure type
         * @return fem::Measure_Type describing the cluster measure
         */
        fem::Measure_Type &get_measure_type()
        {
            return mMeasureType;
        }

        /**
         * get primary type
         * @return mtk::Primary_Void describing the cluster measure
         */
        mtk::Primary_Void &get_primary_type()
        {
            return mIsPrimary;
        }

        /**
         * get leader type
         * @return mtk::Leader_Follower describing the cluster measure
         */
        mtk::Leader_Follower &get_leader_type()
        {
            return mIsLeader;
        }

        /**
         * get the stabilization parameter value
         * @param[ out ] mMEAVal stabilization parameter value
         */
        const Matrix< DDRMat > &val();

        /**
         * evaluate cluster measure
         */
        void eval_cluster_measure();

        /**
         * perturb cluster measure
         * param[ in ] aDeltaH perturbation size
         */
        void perturb_cluster_measure( moris::real aDeltaH );

        /**
         * get the cluster measure derivative wrt pdv
         * @return mdMEAdPDV cluster measure derivative wrt leader dv
         */
        const Matrix< DDRMat > &dMEAdPDV();

        /**
         * evaluate the cluster measure derivatives wrt pdv
         */
        void eval_cluster_measure_derivatives();

      private:
        // cluster pointer
        fem::Cluster *mCluster;

        // enum storage
        fem::Measure_Type    mMeasureType = fem::Measure_Type::UNDEFINED;
        mtk::Primary_Void    mIsPrimary   = mtk::Primary_Void::UNDEFINED;
        mtk::Leader_Follower mIsLeader    = mtk::Leader_Follower::UNDEFINED;

        // storage
        Matrix< DDRMat > mMEAVal;
        Matrix< DDRMat > mdMEAdPDV;

        // flag for evaluation
        bool mMEAEval      = false;
        bool mdMEAdPDVEval = false;
    };
}    // namespace moris::fem
