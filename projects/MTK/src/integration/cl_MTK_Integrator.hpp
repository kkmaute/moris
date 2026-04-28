/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATOR_HPP_
#define SRC_MTK_CL_MTK_INTEGRATOR_HPP_

#include "moris_typedefs.hpp"               //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_MTK_Enums.hpp"                 //MTK/src
#include "cl_MTK_Integration_Rule.hpp"      //MTK/src
#include "cl_MTK_Integration_Coeffs.hpp"    //MTK/src
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    class Integrator
    {        
        // pointer to space rule, if specified
        std::unique_ptr< Integration_Coeffs_Base > mSpaceCoeffs;

        // pointer to time rule, if specified
        std::unique_ptr< Integration_Coeffs_Base > mTimeCoeffs;

        // number of points in space
        uint mNumOfSpacePoints;

        // number of points in time
        uint mNumOfTimePoints;

        // matrix with space points
        Matrix< DDRMat > mSpacePoints;

        // matrix with time points
        Matrix< DDRMat > mTimePoints;

        // matrix with space weights
        Matrix< DDRMat > mSpaceWeights;

        // matrix with time weights
        Matrix< DDRMat > mTimeWeights;

        // Function pointer for computing integration points and weights for each cluster. Currently moment fitting only works for bulk clusters, not for side clusters
        void ( Integrator::* m_compute_cluster_integration_points_and_weights )( const mtk::Cluster* ) = nullptr;

        // Function pointer for computing integration points and weights for each cluster for computing dRdP and dQIdP. Currently moment fitting only works for bulk clusters, not for side clusters
        void ( Integrator::*m_compute_cluster_integration_points_and_weights_perturbed )( const mtk::Cluster *, const moris_index, const Matrix< DDRMat > , const Matrix< DDRMat > ) = nullptr;

        // Moment fitting LHS inverse, so that moments need to be hit by this inverse to generate weights
        Matrix< DDRMat > mMomentFittingLHSinv;

        // Moment Fitting Quadrature points
        Matrix< DDRMat > mMomentFittingQuadPoints;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * constructs an integrator from an integration rule
         **/
        Integrator( const Integration_Rule &aIntegrationRule );

        //------------------------------------------------------------------------------

        /**
         * Sets function pointers to compute integration points and weights for each cluster. Currentoly moment fitting only works for bulk clusters, not for side clusters
         **/

        Integrator( const Integration_Rule &aIntegrationRule, const mtk::Set *aMeshSet, const bool aUseMomentFitting );


        //------------------------------------------------------------------------------
        /**
         * get the number of integration points
        **/

        uint get_number_of_points() const;

        //------------------------------------------------------------------------------
        /**
         * @brief get the integration points as a (d x n) matrix, where d is the dimension of the space
         * and n is the number of integration points (each column contains one point)
         * @param aIntegrationPoints
         */
        void get_points( Matrix< DDRMat > &aIntegrationPoints ) const;
       

        //------------------------------------------------------------------------------
        /**
         * @brief get the integration points as a (d x n) matrix, where d is the dimension of the space
         * and n is the number of integration points (each column contains one point)
         * @return
         */
        Matrix< DDRMat > get_points() const;

        //------------------------------------------------------------------------------
        /**
         * @brief get the integration points as a (d x n) matrix, where d is the dimension of the space
         * and n is the number of integration points (each column contains one point). The integration rule
         * need not be the one supplied to this class at instantiation.
         * @param aIntegrationPoints
         */
        void get_points( Matrix< DDRMat > &aIntegrationPoints, const Integration_Rule &aIntegrationRule ) const;

        Matrix< DDRMat > get_points( const Integration_Rule &aIntegrationRule ) const;

        //------------------------------------------------------------------------------
        /**
         * get the integration point weights
         **/
        void get_weights( Matrix< DDRMat > &aIntegrationWeights ) const;

        Matrix< DDRMat > get_weights() const;

        
        //------------------------------------------------------------------------------
        /**
         * get the integration point weights from any integrator, 
         * not necessarily the one supplied to this class at instantiation
         **/
        void get_weights( Matrix< DDRMat > &aIntegrationWeights, const Integration_Rule &aIntegrationRule ) const;

        Matrix< DDRMat > get_weights( const Integration_Rule &aIntegrationRule ) const;


        //------------------------------------------------------------------------------
        /**
         * Compute the integration points and weights for a bulk cluster using moment fitting
         **/

        void compute_bulk_cluster_integration_points_and_weights_moment_fitting( const Cluster* aCluster );
        //------------------------------------------------------------------------------
        /**
         * Compute the integration points and weights for a bulk cluster using moment fitting
         **/

        void compute_bulk_cluster_integration_points_and_weights_moment_fitting_perturbed( const Cluster *aCluster, const moris_index aNodeIndex, const Matrix< DDRMat > aPerturbation, const Matrix< DDRMat > aPhysicalPerturbation );
        //------------------------------------------------------------------------------

        /**
         * Compute the integration points and weights for a bulk cluster using standard quadrature
         **/

        void compute_bulk_cluster_integration_points_and_weights_standard( const Cluster *aCluster );
        
        //------------------------------------------------------------------------------

        /**
         * Compute the integration points and weights for a double side cluster using Standard scheme
         **/

        void compute_side_cluster_integration_points_and_weights_standard( const Cluster *aCluster );

        //------------------------------------------------------------------------------

        /**
         * Compute the integration points and weights for a side cluster using Standard scheme
         **/

        void compute_double_side_cluster_integration_points_and_weights_standard( const Cluster *aCluster );

        //------------------------------------------------------------------------------

        /** 
         * Compute the time integration points and weights
         */
        //void compute_time_integration_points_and_weights( const bool aIsTimeCon );

        //------------------------------------------------------------------------------

        /*
        *  Compute the integration points and weights for a cluster
        */
        void compute_cluster_integration_points_and_weights( const Cluster* aCluster )
        {
            ( this->*m_compute_cluster_integration_points_and_weights )( aCluster );
        }
        
        //------------------------------------------------------------------------------

        /*
         *  Compute the integration points and weights for a cluster when an IG node is perturbed for computing dRdP and dQIdP
         */
        void compute_cluster_integration_points_and_weights( const Cluster *aCluster, const moris_index aNodeIndex, const Matrix< DDRMat > aPerturbation, const Matrix< DDRMat > aPhysicalPerturbation )
        {
            ( this->*m_compute_cluster_integration_points_and_weights_perturbed )( aCluster, aNodeIndex, aPerturbation, aPhysicalPerturbation );
        }

        //------------------------------------------------------------------------------

        void restore_quadrature_weights_and_points( const Cluster *aCluster, Matrix< DDRMat > &aQuadraturePoints, Matrix< DDRMat > &aQuadratureWeights );

        //------------------------------------------------------------------------------

        /* 
        ** Determine the appropriate integration order for an unzipped IP cell based on the integration order assigned to its subphase cells. 
         *
        */
        
        mtk::Integration_Order get_ip_integration_order_from_cut_cell( const mtk::Geometry_Type &aGeometryType, const Integration_Rule &aIntegrationRule ) const;



        //------------------------------------------------------------------------------

    };    // class Integrator

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATOR_HPP_ */
