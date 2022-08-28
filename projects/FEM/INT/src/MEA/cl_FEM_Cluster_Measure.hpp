/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Cluster_Measure.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CLUSTER_MEASURE_HPP_
#define SRC_FEM_CL_FEM_CLUSTER_MEASURE_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
//LNA/src
#include "cl_Matrix.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace fem
    {
        class Cluster;
        class Set;

        //------------------------------------------------------------------------------
        /**
         * Cluster Measure
         * This class provides cluster measure mainly for stabilization purpose.
         */
        class Cluster_Measure
        {
                //------------------------------------------------------------------------------
            protected :

                //------------------------------------------------------------------------------

                // cluster pointer
                fem::Cluster * mCluster = nullptr;

                // enum storage
                fem::Measure_Type mMeasureType = fem::Measure_Type::UNDEFINED;
                mtk::Primary_Void mIsPrimary = mtk::Primary_Void::UNDEFINED;
                mtk::Master_Slave mIsMaster = mtk::Master_Slave::UNDEFINED;

                // storage
                Matrix< DDRMat > mMEAVal;
                Matrix< DDRMat > mdMEAdPDV;

                // flag for evaluation
                bool mMEAEval     = false;
                bool mdMEAdPDVEval = false;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                Cluster_Measure(
                        fem::Measure_Type aMeasureType,
                        mtk::Primary_Void aIsPrimary,
                        mtk::Master_Slave aIsMaster,
                        fem::Cluster * aCluster )
                : mCluster( aCluster ),
                  mMeasureType( aMeasureType ),
                  mIsPrimary( aIsPrimary ),
                  mIsMaster( aIsMaster )
                {};

                /**
                 * trivial constructor for test purpose only
                 */
                Cluster_Measure()
                {
                    mMEAEval = true;
                    mMEAVal = {{1.0}};
                };

                //------------------------------------------------------------------------------
                /**
                 * virtual destructor
                 */
                virtual ~Cluster_Measure(){};

                //------------------------------------------------------------------------------
                /**
                 * set cluster
                 * @param[ in ] aCluster a fem cluster pointer
                 */
                void set_cluster( fem::Cluster * aCluster )
                {
                    // set a cluster
                    mCluster = aCluster;

                    // evaluate the cluster measure
                    this->eval_cluster_measure();
                }

                //------------------------------------------------------------------------------
                /**
                 * get measure type
                 * @return fem::Measure_Type describing the cluster measure
                 */
                fem::Measure_Type & get_measure_type()
                {
                    return  mMeasureType;
                }

                //------------------------------------------------------------------------------
                /**
                 * get primary type
                 * @return mtk::Primary_Void describing the cluster measure
                 */
                mtk::Primary_Void & get_primary_type()
                {
                    return  mIsPrimary;
                }

                //------------------------------------------------------------------------------
                /**
                 * get master type
                 * @return mtk::Master_Slave describing the cluster measure
                 */
                mtk::Master_Slave & get_master_type()
                {
                    return  mIsMaster;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the stabilization parameter value
                 * @param[ out ] mMEAVal stabilization parameter value
                 */
                const Matrix< DDRMat > & val();

                //------------------------------------------------------------------------------
                /**
                 * evaluate cluster measure
                 */
                void eval_cluster_measure();

                //------------------------------------------------------------------------------
                /**
                 * perturb cluster measure
                 * param[ in ] aDeltaH perturbation size
                 */
                void perturb_cluster_measure( moris::real aDeltaH );

                //------------------------------------------------------------------------------
                /**
                 * get the cluster measure derivative wrt pdv
                 * @return mdMEAdPDV cluster measure derivative wrt master dv
                 */
                const Matrix< DDRMat > & dMEAdPDV();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the cluster measure derivatives wrt pdv
                 */
                void eval_cluster_measure_derivatives();

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CLUSTER_MEASURE_HPP_ */

