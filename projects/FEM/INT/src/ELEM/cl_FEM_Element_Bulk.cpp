/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Bulk.cpp
 *
 */

#include <iostream>
// FEM/INT/src
#include "cl_FEM_Element_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
// FEM/MSI/src
#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {
        //----------------------------------------------------------------------

        Element_Bulk::Element_Bulk(
                mtk::Cell const *  aCell,
                Set*               aSet,
                Cluster*           aCluster,
                moris::moris_index aCellIndexInCluster )
                : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {
        }

        //----------------------------------------------------------------------

        Element_Bulk::~Element_Bulk() {}

        //----------------------------------------------------------------------

        void
        Element_Bulk::init_ig_geometry_interpolator(
                Matrix< DDSMat >& aGeoLocalAssembly )
        {
            // get geometry interpolator for IG element
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mLeaderCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get the local cluster assembly indices
            if ( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices for IG element
                Matrix< IndexMat > tVertexIndices = mLeaderCell->get_vertex_inds();

                // get the requested geo pdv types
                moris::Cell< enum PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff( tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff( tIGParamTimeCoords );
        }

        //----------------------------------------------------------------------

        void
        Element_Bulk::init_ig_geometry_interpolator()
        {
            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mLeaderCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get geometry interpolator for IG element
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff( tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff( tIGParamTimeCoords );

            // set geometry interpolator of eigen vector FI
            if ( mSet->mNumEigenVectors > 0 )
            {
                Geometry_Interpolator* tIGGI =
                        mSet->get_field_interpolator_manager_eigen_vectors()->get_IG_geometry_interpolator();

                // set physical space and time coefficients for IG element GI
                tIGGI->set_space_coeff( tIGPhysSpaceCoords );
                tIGGI->set_time_coeff( tIGPhysTimeCoords );

                // set parametric space and time coefficients for IG element GI
                tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
                tIGGI->set_time_param_coeff( tIGParamTimeCoords );
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the current integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME: enforced nodal weak bcs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute residual at evaluation point
                    tReqIWG->compute_residual( tWStar );

                    // compute off-diagonal Jacobian for staggered solve
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_jacobian()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute Jacobian at evaluation point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_jacobian_and_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IWGs or IQIs
            if ( tNumIWGs == 0 && tNumIQIs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG = mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    if ( mSet->mEquationModel->get_is_forward_analysis() )
                    {
                        // compute residual at evaluation point
                        tReqIWG->compute_residual( tWStar );
                    }

                    // compute Jacobian at evaluation point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }

                // if sensitivity analysis and requested IQIs
                if ( ( !mSet->mEquationModel->get_is_forward_analysis() ) && ( tNumIQIs > 0 ) )
                {
                    // loop over the IQIs
                    for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // get requested IQI
                        const std::shared_ptr< IQI >& tReqIQI =
                                mSet->get_requested_IQIs()( iIQI );

                        // reset IQI
                        tReqIQI->reset_eval_flags();

                        // compute dQIdu at evaluation point
                        ( this->*m_compute_dQIdu )( tReqIQI, tWStar );
                    }
                }

                mSet->mFemModel->mBulkGaussPoints++;
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_dRdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_QI()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if ( tNumIQIs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // if eigen vectors
                if ( mSet->mNumEigenVectors > 0 )
                {
                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager_eigen_vectors()->set_space_time_from_local_IG_point( tLocalIntegPoint );
                }

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // compute QI at evaluation point
                    tReqIQI->compute_QI( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_dQIdp_explicit()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if ( tNumIQIs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IWG
                    tReqIQI->reset_eval_flags();

                    // compute dQIdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dQIdp )( tReqIQI, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_dRdp_and_dQIdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IWGs and IQIs
            if ( tNumIWGs == 0 && tNumIQIs == 0 )
            {
                return;
            }

            // get the vertices indices
            Matrix< IndexMat > tVertexIndices = mLeaderCell->get_vertex_inds();

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->    //
                        set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // compute dQIdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dQIdp )( tReqIQI, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_dQIdu()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if ( tNumIQIs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // compute dQIdu at evaluation point
                    ( this->*m_compute_dQIdu )( tReqIQI, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_quantity_of_interest_global( const uint aFemMeshIndex )
        {
            // get number of active local IQIs
            uint tNumLocalIQIs = mSet->get_number_of_requested_global_IQIs_for_visualization();

            // check that some IQIs need to be evaluated
            if ( tNumLocalIQIs == 0 )
            {
                return;
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // loop over IQI
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_global_IQIs_for_visualization()( iIQI );

                    // get IQI global index
                    moris_index tGlobalIndex =
                            mSet->get_requested_global_IQIs_global_indices_for_visualization()( iIQI );

                    // reset the requested IQI
                    tReqIQI->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIGlobal( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQIGlobal );

                    // assemble the global QI value on the set
                    ( *( mSet->mSetGlobalValues ) )( tGlobalIndex ) += tWStar * tQIGlobal( 0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_quantity_of_interest_elemental(
                const uint aFemMeshIndex,
                const bool aAverageOutput )
        {
            // get number of active local IQIs
            uint tNumLocalIQIs = mSet->get_number_of_requested_elemental_IQIs_for_visualization();

            // check that some IQIs need to be evaluated
            if ( tNumLocalIQIs == 0 )
            {
                return;
            }

            // get the VIS mesh index
            uint tVisMeshIndex = aFemMeshIndex - 1;

            // get the IG cell's index in the VIS mesh
            moris_index tCellIndex = mLeaderCell->get_index();
            moris_index tAssemblyIndex = mSet->mCellAssemblyMap( tVisMeshIndex )( tCellIndex );
            MORIS_ASSERT( tAssemblyIndex > -1,
                    "FEM::Element_Bulk::compute_quantity_of_interest_elemental() - "
                    "Element not part of VIS cell assembly map." );

            // FIXME: this seems inefficient and should be done elsewhere
            // set unused IQI values to not be NAN
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get IQI global index
                moris_index tGlobalIqiIndex =
                        mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                if ( ( *mSet->mSetElementalValues )( tAssemblyIndex, tGlobalIqiIndex ) == std::numeric_limits< real >::quiet_NaN() )
                {
                    ( *mSet->mSetElementalValues )( tAssemblyIndex, tGlobalIqiIndex ) = 0.0;
                }
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // initialize space - time volume
            real tSpaceTimeVolume = 0.0;

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // add contribution to space-time volume
                tSpaceTimeVolume += tWStar;

                // loop over IQI
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_elemental_IQIs_for_visualization()( iIQI );

                    // get IQI global index
                    moris_index tGlobalIqiIndex =
                            mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                    // reset the requested IQI
                    tReqIQI->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIElemental( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQIElemental );

                    // assemble the QI value on the set
                    ( *mSet->mSetElementalValues )( tAssemblyIndex, tGlobalIqiIndex ) += tWStar * tQIElemental( 0 );
                }
            }

            // loop over IQI and divide each elemental IQI by space-time volume if requested
            if ( aAverageOutput )
            {
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get IQI global index
                    moris_index tGlobalIqiIndex =
                            mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                    // divide by space-time volume
                    moris_index tLeaderCellIndex = mLeaderCell->get_index();
                    moris_index tAssemblyIndex = mSet->mCellAssemblyMap( tVisMeshIndex )( tLeaderCellIndex );
                    ( *mSet->mSetElementalValues )( tAssemblyIndex, tGlobalIqiIndex )  /= tSpaceTimeVolume;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Bulk::compute_quantity_of_interest_elemental(
                Matrix< DDRMat >& aValues,
                uint              aIQIIndex,
                real&             aSpaceTimeVolume )
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            Matrix< DDRMat > tValMat( 1, 1, 0.0 );

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat >& tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // add contribution to space-time volume
                aSpaceTimeVolume += tWStar;

                // get requested IQI
                const std::shared_ptr< IQI >& tReqIQI =
                        mSet->get_requested_field_IQIs()( aIQIIndex );

                // reset the requested IQI
                tReqIQI->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIElemental( 1, 1, 0.0 );
                tReqIQI->compute_QI( tQIElemental );

                tValMat( 0 ) += tWStar * tQIElemental( 0 );
            }

            // assemble the QI value on the set
            aValues( 0 ) += tValMat( 0 );
        }

        //------------------------------------------------------------------------------

        real
        Element_Bulk::compute_volume( mtk::Leader_Follower aIsLeader )
        {
            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator();

            // get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // get geometry interpolator
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // loop over integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                tIGGI->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += tIGGI->det_J() * mSet->get_integration_weights()( iGP );
            }

            // get time step
            real tTimeStep = tIGGI->get_time_step();

            // return the volume value
            return tVolume / tTimeStep;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
