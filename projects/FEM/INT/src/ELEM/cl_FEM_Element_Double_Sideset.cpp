/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Double_Sideset.cpp
 *
 */

#include <iostream>
// FEM/INT/src
#include "cl_FEM_Element_Double_Sideset.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_Side_Coordinate_Map.hpp"
// FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset(
                mtk::Cell const *  aLeaderIGCell,
                mtk::Cell const *  aFollowerIGCell,
                Set*               aSet,
                Cluster*           aCluster,
                moris::moris_index aCellIndexInCluster )
                : Element(
                        aLeaderIGCell,
                        aFollowerIGCell,
                        aSet,
                        aCluster,
                        aCellIndexInCluster )
        {
        }

        //------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset() {}

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::init_ig_geometry_interpolator(
                uint aLeaderSideOrdinal,
                uint aFollowerSideOrdinal )
        {
            // get leader IG geometry interpolator
            Geometry_Interpolator* tLeaderIGGI =
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                    get_IG_geometry_interpolator();

            // get follower IG geometry interpolator
            Geometry_Interpolator* tFollowerIGGI =
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                    get_IG_geometry_interpolator();

            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tLeaderIGPhysSpaceCoords =
                    mLeaderCell->get_cell_physical_coords_on_side_ordinal( aLeaderSideOrdinal );
            Matrix< DDRMat > tLeaderIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get follower physical space and time coordinates for IG element
            Matrix< DDRMat > tFollowerIGPhysSpaceCoords =
                    mFollowerCell->get_cell_physical_coords_on_side_ordinal( aFollowerSideOrdinal );
            Matrix< DDRMat > tFollowerIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tLeaderIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aLeaderSideOrdinal,
                            mtk::Leader_Follower::LEADER );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tLeaderIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get follower parametric space and time coordinates for IG element
            Matrix< DDRMat > tFollowerIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aFollowerSideOrdinal,
                            mtk::Leader_Follower::FOLLOWER );
            // FIXME not true if time is not linear
            Matrix< DDRMat > tFollowerIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // set physical space and time coefficients for leader IG element GI
            tLeaderIGGI->set_space_coeff( tLeaderIGPhysSpaceCoords );
            tLeaderIGGI->set_time_coeff( tLeaderIGPhysTimeCoords );

            // set physical space and time coefficients for follower IG element GI
            tFollowerIGGI->set_space_coeff( tFollowerIGPhysSpaceCoords );
            tFollowerIGGI->set_time_coeff( tFollowerIGPhysTimeCoords );

            // set parametric space and time coefficients for leader IG element GI
            tLeaderIGGI->set_space_param_coeff( tLeaderIGParamSpaceCoords );
            tLeaderIGGI->set_time_param_coeff( tLeaderIGParamTimeCoords );

            // set parametric space and time coefficients for follower IG element GI
            tFollowerIGGI->set_space_param_coeff( tFollowerIGParamSpaceCoords );
            tFollowerIGGI->set_time_param_coeff( tFollowerIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::init_ig_geometry_interpolator(
                uint              aLeaderSideOrdinal,
                uint              aFollowerSideOrdinal,
                Matrix< DDSMat >& aGeoLocalAssembly )
        {
            // get leader IG geometry interpolator
            Geometry_Interpolator* tLeaderIGGI =
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                    get_IG_geometry_interpolator();

            // get follower IG geometry interpolator
            Geometry_Interpolator* tFollowerIGGI =
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                    get_IG_geometry_interpolator();

            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tLeaderIGPhysSpaceCoords =
                    mLeaderCell->get_cell_physical_coords_on_side_ordinal( aLeaderSideOrdinal );

            Matrix< DDRMat > tLeaderIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get follower physical space and time coordinates for IG element
            Matrix< DDRMat > tFollowerIGPhysSpaceCoords =
                    mFollowerCell->get_cell_physical_coords_on_side_ordinal( aFollowerSideOrdinal );

            Matrix< DDRMat > tFollowerIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tLeaderIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aLeaderSideOrdinal,
                            mtk::Leader_Follower::LEADER );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tLeaderIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get follower parametric space and time coordinates for IG element
            Matrix< DDRMat > tFollowerIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell(
                            mCellIndexInCluster,
                            aFollowerSideOrdinal,
                            mtk::Leader_Follower::FOLLOWER );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tFollowerIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get the local cluster assembly indices
            if ( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices
                Matrix< IndexMat > tLeaderVertexIndices =
                        mLeaderCell->get_vertices_ind_on_side_ordinal( aLeaderSideOrdinal );

                Matrix< IndexMat > tFollowerVertexIndices =
                        mFollowerCell->get_vertices_ind_on_side_ordinal( aFollowerSideOrdinal );

                // get the requested geo pdv types
                moris::Cell< enum ge::PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tLeaderVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and time coefficients for leader IG element GI
            tLeaderIGGI->set_space_coeff( tLeaderIGPhysSpaceCoords );
            tLeaderIGGI->set_time_coeff( tLeaderIGPhysTimeCoords );

            // set physical space and time coefficients for follower IG element GI
            tFollowerIGGI->set_space_coeff( tFollowerIGPhysSpaceCoords );
            tFollowerIGGI->set_time_coeff( tFollowerIGPhysTimeCoords );

            // set parametric space and time coefficients for leader IG element GI
            tLeaderIGGI->set_space_param_coeff( tLeaderIGParamSpaceCoords );
            tLeaderIGGI->set_time_param_coeff( tLeaderIGParamTimeCoords );

            // set parametric space and time coefficients for follower IG element GI
            tFollowerIGGI->set_space_param_coeff( tFollowerIGParamSpaceCoords );
            tFollowerIGGI->set_time_param_coeff( tFollowerIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair(
                            mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );

            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the follower integration cell
                Matrix< DDRMat > tFollowerLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tFollowerNodeOrdOnSide,
                                tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute residual at integration point
                    tReqIWG->compute_residual( tWStar );

                    // compute Jacobian at integration point
                    // compute off-diagonal Jacobian for staggered solve
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_jacobian()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair(
                            mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );

            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the follower integration cell
                Matrix< DDRMat > tFollowerLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tFollowerNodeOrdOnSide,
                                tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute residual at integration point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_jacobian_and_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair( mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );
            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the follower integration cell
                Matrix< DDRMat > tFollowerLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tFollowerNodeOrdOnSide,
                                tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    if ( mSet->mEquationModel->get_is_forward_analysis() )
                    {
                        // compute residual at integration point
                        tReqIWG->compute_residual( tWStar );
                    }

                    // compute jacobian at integration point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }

                mSet->mFemModel->mDoubleSidedSideSetsGaussPoints++;
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_dRdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator(
                    tLeaderSideOrd,
                    tFollowerSideOrd,
                    tGeoLocalAssembly );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair(
                            mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );

            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // get the vertices indices
            moris::Cell< Matrix< IndexMat > > tVertexIndices( 2 );
            tVertexIndices( 0 ) =
                    mLeaderCell->get_vertices_ind_on_side_ordinal( tLeaderSideOrd );
            tVertexIndices( 1 ) =
                    mFollowerCell->get_vertices_ind_on_side_ordinal( tFollowerSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat > tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the follower integration cell
                const Matrix< DDRMat > tFollowerLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tFollowerNodeOrdOnSide,
                                tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute dRdpMat at evaluation point
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_dRdp_and_dQIdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if ( tNumIWGs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator(
                    tLeaderSideOrd,
                    tFollowerSideOrd,
                    tGeoLocalAssembly );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair(
                            mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );

            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // get the vertices indices
            moris::Cell< Matrix< IndexMat > > tVertexIndices( 2 );
            tVertexIndices( 0 ) =
                    mLeaderCell->get_vertices_ind_on_side_ordinal( tLeaderSideOrd );
            tVertexIndices( 1 ) =
                    mFollowerCell->get_vertices_ind_on_side_ordinal( tFollowerSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat > tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get local integration point for the follower integration cell
                const Matrix< DDRMat > tFollowerLocalIntegPoint =
                        side_coordinate_map(
                                mSet->get_IG_geometry_type(),
                                tFollowerNodeOrdOnSide,
                                tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                const real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                const real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh and set if for the IWG
                const Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IWGs
                for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG >& tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute dRdpMat at evaluation point
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }

                // FIXME: add part over IQIs
            }
        }

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_quantity_of_interest_global( const uint aFemMeshIndex )
        {
            // get number of active local IQIs
            uint tNumLocalIQIs = mSet->get_number_of_requested_global_IQIs_for_visualization();

            // check that some IQIs need to be evaluated
            if ( tNumLocalIQIs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get a leader follower vertex pair on the facet
            moris::mtk::Vertex const * tLeaderVertex   = mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 );
            moris::mtk::Vertex const * tFollowerVertex = mCluster->get_left_vertex_pair( tLeaderVertex );

            // get the ordinal of the follower vertex wrt the current facet
            moris_index tFollowerVertexOrdinalOnFacet = mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerVertex );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over integration points
            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get a copy of the local integration point for the follower integration cell
                Matrix< DDRMat > tFollowerLocalIntegPoint = side_coordinate_map(
                        mSet->get_IG_geometry_type(),
                        tFollowerVertexOrdinalOnFacet,
                        tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower FIs and GIs
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                        ->set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                        ->set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the requested IQIs on the element
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

                    // set the normal for the IWG
                    tReqIQI->set_normal( tNormal );

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIGlobal( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQIGlobal );

                    // assemble computed QI on the set
                    ( *( mSet->mSetGlobalValues ) )( tGlobalIndex ) += tWStar * tQIGlobal( 0 );

                }    // end for: each IQI on the current element

            }        // end for: each integration point

        }            // end function: Element_Double_Sideset::compute_quantity_of_interest_global()

        //------------------------------------------------------------------------------

        void
        Element_Double_Sideset::compute_quantity_of_interest_elemental(
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

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // get the leader IG cell's index in the VIS mesh
            moris_index tLeaderCellIndex          = mLeaderCell->get_index();
            moris_index tLeaderFacetIndexInVisSet = mSet->mFacetAssemblyMap( tVisMeshIndex )( tLeaderCellIndex, tLeaderSideOrd );
            MORIS_ASSERT( tLeaderFacetIndexInVisSet > -1,
                    "FEM::Element_Double_Sideset::compute_quantity_of_interest_elemental() - "
                    "Facet not part of VIS facet assembly map." );

            // FIXME: this seems inefficient and should be done elsewhere
            // set unused IQI values to not be NAN
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get IQI global index
                moris_index tGlobalIqiIndex =
                        mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                if ( ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) == std::numeric_limits< real >::quiet_NaN() )
                {
                    ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) = 0.0;
                }
            }

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // initialize space - time volume
            real tSpaceTimeVolume = 0.0;

            // get a leader follower vertex pair on the facet
            moris::mtk::Vertex const * tLeaderVertex   = mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 );
            moris::mtk::Vertex const * tFollowerVertex = mCluster->get_left_vertex_pair( tLeaderVertex );

            // get the ordinal of the follower vertex wrt the current facet
            moris_index tFollowerVertexOrdinalOnFacet = mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerVertex );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over integration points
            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get a copy of the local integration point for the follower integration cell
                Matrix< DDRMat > tFollowerLocalIntegPoint = side_coordinate_map(
                        mSet->get_IG_geometry_type(),
                        tFollowerVertexOrdinalOnFacet,
                        tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower FIs and GIs
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                        ->set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                        ->set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the requested IQIs on the element
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

                    // set the normal for the IWG
                    tReqIQI->set_normal( tNormal );

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIElemental( 1, 1, 0.0 );
                    tReqIQI->compute_QI( tQIElemental );

                    // assemble computed QI on the set
                    ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) += tWStar * tQIElemental( 0 );

                }    // end for: each IQI to be evaluated on cluster

            }        // end for: each integration point

            // loop over IQI and divide each elemental IQI by space-time volume if requested
            if ( aAverageOutput )
            {
                for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get IQI global index
                    moris_index tGlobalIqiIndex = mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                    // normalize by space-time volume
                    ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) /= tSpaceTimeVolume;
                }
            }
        }

        //------------------------------------------------------------------------------

        real
        Element_Double_Sideset::compute_volume( mtk::Leader_Follower aIsLeader )
        {
            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // initialize volume
            real tVolume = 0;

            // get geometry interpolator
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager( aIsLeader )->get_IG_geometry_interpolator();

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

        void
        Element_Double_Sideset::compute_QI()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if ( tNumIQIs == 0 )
            {
                return;
            }

            // get treated side ordinal on the leader and on the follower
            uint tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );
            uint tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( mCellIndexInCluster );

            // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
            this->init_ig_geometry_interpolator( tLeaderSideOrd, tFollowerSideOrd );

            // get first corresponding node from leader to follower
            moris::mtk::Vertex const * tFollowerNode =
                    mCluster->get_left_vertex_pair( mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 ) );
            moris_index tFollowerNodeOrdOnSide =
                    mCluster->get_right_vertex_ordinal_on_facet( mCellIndexInCluster, tFollowerNode );

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {

                // get local integration point for the leader integration cell
                const Matrix< DDRMat >& tLeaderLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // get copy of local integration point for the follower integration cell
                const Matrix< DDRMat > tFollowerLocalIntegPoint = side_coordinate_map(
                        mSet->get_IG_geometry_type(),
                        tFollowerNodeOrdOnSide,
                        tLeaderLocalIntegPoint );

                // set evaluation point for leader and follower interpolators
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                        set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                        set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

                // compute detJ of integration domain
                const real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // skip if detJ smaller than threshold
                if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                {
                    continue;
                }

                // compute integration point weight
                const real tWStar = mSet->get_integration_weights()( iGP ) * tDetJ;

                // get the normal from mesh
                const Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );

                // loop over the IQIs
                for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI >& tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // set the normal for the IQI
                    tReqIQI->set_normal( tNormal );

                    // compute QI at evaluation point
                    tReqIQI->compute_QI( tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
