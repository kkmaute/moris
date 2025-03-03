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
#include <limits>
#include <memory>
#include "cl_FEM_Element_Double_Sideset.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_FEM_Side_Coordinate_Map.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_assert.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"

namespace moris::fem
{
    Element_Double_Sideset::Element_Double_Sideset(
            mtk::Cell const *        aLeaderIGCell,
            mtk::Cell const *        aFollowerIGCell,
            Set*                     aSet,
            Cluster*                 aCluster,
            moris::moris_index const aCellIndexInCluster )
            : Element(
                      aLeaderIGCell,
                      aFollowerIGCell,
                      aSet,
                      aCluster,
                      aCellIndexInCluster )
    {
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::init_ig_geometry_interpolator() const
    {
        initialize_leader_follower_ig_interpolator( mtk::Leader_Follower::LEADER );
        initialize_leader_follower_ig_interpolator( mtk::Leader_Follower::FOLLOWER );
    }

    //------------------------------------------------------------------------------

    void Element_Double_Sideset::initialize_leader_follower_ig_interpolator( mtk::Leader_Follower const aLeaderFollowerType ) const
    {
        // get treated side ordinal on the leader and on the follower
        moris_index      tSideOrd;
        moris_index      tLocalCellIndex;
        const mtk::Cell* tCell;
        if ( aLeaderFollowerType == mtk::Leader_Follower::LEADER )
        {
            tCell           = mLeaderCell;
            tLocalCellIndex = get_leader_local_cell_index();
            tSideOrd        = mCluster->mLeaderListOfSideOrdinals( tLocalCellIndex );
        }
        else
        {
            tCell           = mFollowerCell;
            tLocalCellIndex = get_follower_local_cell_index();
            tSideOrd        = mCluster->mFollowerListOfSideOrdinals( tLocalCellIndex );
        }

        // get leader IG geometry interpolator
        Geometry_Interpolator* tIGInterpolator = mSet->get_field_interpolator_manager( aLeaderFollowerType )->get_IG_geometry_interpolator();

        // physical coefficients
        tIGInterpolator->set_space_coeff( tCell->get_cell_physical_coords_on_side_ordinal( tSideOrd ) );
        tIGInterpolator->set_time_coeff( mCluster->mInterpolationElement->get_time() );

        // parametric coefficients
        tIGInterpolator->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( tLocalCellIndex, tSideOrd, aLeaderFollowerType ) );
        tIGInterpolator->set_time_param_coeff( { { -1.0 }, { 1.0 } } );    // FIXME not true if time is not linear
    }

    //------------------------------------------------------------------------------

    Matrix< DDSMat > Element_Double_Sideset::get_local_cluster_assembly_indices( moris_index const aLeaderSideOrdinal, moris_index const aFollowerSideOrdinal ) const
    {
        Matrix< DDSMat > tGeoLocalAssembly;
        if ( mSet->get_geo_pdv_assembly_flag() )
        {
            // get the vertices indices
            Matrix< IndexMat > const tLeaderVertexIndices   = mLeaderCell->get_vertices_ind_on_side_ordinal( aLeaderSideOrdinal );
            Matrix< IndexMat > const tFollowerVertexIndices = mFollowerCell->get_vertices_ind_on_side_ordinal( aFollowerSideOrdinal );

            // get the requested geo pdv types
            Vector< enum gen::PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // get local assembly indices
            mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices( tLeaderVertexIndices, tGeoPdvType, tGeoLocalAssembly );
        }
        return tGeoLocalAssembly;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Element_Double_Sideset::get_follower_integration_point( uint const aGPIndex ) const
    {
        moris_index const   tLeaderSideOrd         = mCluster->mLeaderListOfSideOrdinals( get_leader_local_cell_index() );
        mtk::Vertex const * tLeaderVertex          = mLeaderCell->get_vertices_on_side_ordinal( tLeaderSideOrd )( 0 );
        mtk::Vertex const * tFollowerVertex        = mCluster->get_left_vertex_pair( tLeaderVertex );
        moris_index const   tFollowerNodeOrdOnSide = mCluster->get_right_vertex_ordinal_on_facet( get_follower_local_cell_index(), tFollowerVertex );

        return side_coordinate_map(
                mSet->get_IG_geometry_type(),
                tFollowerNodeOrdOnSide,
                this->get_leader_integration_point( aGPIndex ) );
    }

    //------------------------------------------------------------------------------

    moris_index Element_Double_Sideset::get_leader_local_cell_index() const
    {
        return mCellIndexInCluster;
    }

    //------------------------------------------------------------------------------

    moris_index Element_Double_Sideset::get_follower_local_cell_index() const
    {
        return mCellIndexInCluster;
    }

    //------------------------------------------------------------------------------

    void Element_Double_Sideset::initialize_quadrature_point( uint iGP )
    {
        // get local integration point for the leader and follower integration cell
        const Matrix< DDRMat >& tLeaderLocalIntegPoint   = get_leader_integration_point( iGP );
        const Matrix< DDRMat >& tFollowerLocalIntegPoint = get_follower_integration_point( iGP );

        // set evaluation point for leader and follower interpolators
        mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                ->set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

        mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                ->set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_residual()
    {
        uint const tNumIWGs = mSet->get_number_of_requested_IWGs();
        if ( tNumIWGs == 0 )
        {
            return;
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // loop over the integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // get the normal from mesh
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

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
                this->set_iwg_jacobian_strategy( tReqIWG );

                ( this->*m_compute_jacobian )( tReqIWG, tWStar );

                this->reset_iwg_jacobian_strategy();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_jacobian()
    {
        // get number of IWGs
        uint const tNumIWGs = mSet->get_number_of_requested_IWGs();

        // check for active IWGs
        if ( tNumIWGs == 0 )
        {
            return;
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // loop over the integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // get the normal from mesh and set if for the IWG
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get requested IWG
                std::shared_ptr< IWG > const & tReqIWG = mSet->get_requested_IWGs()( iIWG );

                // reset IWG
                tReqIWG->reset_eval_flags();

                // set the normal for the IWG
                tReqIWG->set_normal( tNormal );

                // compute residual at integration point
                this->set_iwg_jacobian_strategy( tReqIWG );

                ( this->*m_compute_jacobian )( tReqIWG, tWStar );

                this->reset_iwg_jacobian_strategy();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_jacobian_and_residual()
    {
        // get number of IWGs
        uint const tNumIWGs = mSet->get_number_of_requested_IWGs();

        // check for active IWGs
        if ( tNumIWGs == 0 )
        {
            return;
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // loop over the integration points
        uint const tNumIntegPoints = this->get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()
                                       ->get_IG_geometry_interpolator()
                                       ->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = this->get_integration_weight( iGP ) * tDetJ;

            // get normal at integration point
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get requested IWG
                const std::shared_ptr< IWG >& tReqIWG = mSet->get_requested_IWGs()( iIWG );

                // reset IWG
                tReqIWG->reset_eval_flags();

                // set the normal for the IWG
                tReqIWG->set_normal( tNormal );

                if ( mSet->mEquationModel->is_forward_analysis() )
                {
                    // compute residual at integration point
                    tReqIWG->compute_residual( tWStar );
                }

                // compute jacobian at integration point
                this->set_iwg_jacobian_strategy( tReqIWG );

                ( this->*m_compute_jacobian )( tReqIWG, tWStar );

                this->reset_iwg_jacobian_strategy();
            }

            mSet->mFemModel->mDoubleSidedSideSetsGaussPoints++;
        }
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_dRdp()
    {
        // get number of IWGs
        uint const tNumIWGs = mSet->get_number_of_requested_IWGs();

        // check for active IWGs
        if ( tNumIWGs == 0 )
        {
            return;
        }

        // get treated side ordinal on the leader and on the follower
        moris_index const tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( get_leader_local_cell_index() );
        moris_index const tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( get_follower_local_cell_index() );

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        Matrix< DDSMat > tGeoLocalAssembly = get_local_cluster_assembly_indices( tLeaderSideOrd, tFollowerSideOrd );

        // get the vertices indices
        Vector< Matrix< IndexMat > > tVertexIndices( 2 );
        tVertexIndices( 0 ) = mLeaderCell->get_vertices_ind_on_side_ordinal( tLeaderSideOrd );
        tVertexIndices( 1 ) = mFollowerCell->get_vertices_ind_on_side_ordinal( tFollowerSideOrd );

        // loop over integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // get the normal and set if for the IWG
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get requested IWG
                const std::shared_ptr< IWG >& tReqIWG = mSet->get_requested_IWGs()( iIWG );

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
        uint const tNumIWGs = mSet->get_number_of_requested_IWGs();

        // check for active IWGs
        if ( tNumIWGs == 0 )
        {
            return;
        }

        // get treated side ordinal on the leader and on the follower
        moris_index const tLeaderSideOrd   = mCluster->mLeaderListOfSideOrdinals( get_leader_local_cell_index() );
        moris_index const tFollowerSideOrd = mCluster->mFollowerListOfSideOrdinals( get_follower_local_cell_index() );

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        Matrix< DDSMat > tGeoLocalAssembly = get_local_cluster_assembly_indices( tLeaderSideOrd, tFollowerSideOrd );

        // get the vertices indices
        Vector< Matrix< IndexMat > > tVertexIndices( 2 );
        tVertexIndices( 0 ) = mLeaderCell->get_vertices_ind_on_side_ordinal( tLeaderSideOrd );
        tVertexIndices( 1 ) = mFollowerCell->get_vertices_ind_on_side_ordinal( tFollowerSideOrd );

        // loop over integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // get the normal and set if for the IWG
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get requested IWG
                const std::shared_ptr< IWG >& tReqIWG = mSet->get_requested_IWGs()( iIWG );

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
    Element_Double_Sideset::compute_quantity_of_interest_global( const uint /*aFemMeshIndex*/ )
    {
        // get number of active local IQIs
        uint const tNumLocalIQIs = mSet->get_number_of_requested_global_IQIs_for_visualization();

        // check that some IQIs need to be evaluated
        if ( tNumLocalIQIs == 0 )
        {
            return;
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // loop over integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // get normal at integration point
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the requested IQIs on the element
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get requested IQI
                const std::shared_ptr< IQI >& tReqIQI      = mSet->get_requested_global_IQIs_for_visualization()( iIQI );
                moris_index const             tGlobalIndex = mSet->get_requested_global_IQIs_global_indices_for_visualization()( iIQI );

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
        }    // end for: each integration point
    }    // end function: Element_Double_Sideset::compute_quantity_of_interest_global()

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_quantity_of_interest_elemental(
            const uint aFemMeshIndex,
            const bool aAverageOutput )
    {
        // get number of active local IQIs
        uint const tNumLocalIQIs = mSet->get_number_of_requested_elemental_IQIs_for_visualization();

        // check that some IQIs need to be evaluated
        if ( tNumLocalIQIs == 0 )
        {
            return;
        }

        // get the VIS mesh index
        uint const tVisMeshIndex = aFemMeshIndex - 1;

        // get the leader IG cell's index in the VIS mesh
        moris_index const tLeaderCellIndex          = mLeaderCell->get_index();
        moris_index const tLeaderSideOrd            = mCluster->mLeaderListOfSideOrdinals( get_leader_local_cell_index() );
        moris_index const tLeaderFacetIndexInVisSet = mSet->mFacetAssemblyMap( tVisMeshIndex )( tLeaderCellIndex, tLeaderSideOrd );

        MORIS_ASSERT( tLeaderFacetIndexInVisSet > -1,
                "FEM::Element_Double_Sideset::compute_quantity_of_interest_elemental() - "
                "Facet not part of VIS facet assembly map." );

        // FIXME: this seems inefficient and should be done elsewhere
        // set unused IQI values to not be NAN
        for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
        {
            // get IQI global index
            moris_index const tGlobalIqiIndex = mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

            if ( ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) == std::numeric_limits< real >::quiet_NaN() )
            {
                ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) = 0.0;
            }
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // initialize space - time volume
        real tSpaceTimeVolume = 0.0;

        // loop over integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            real const tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            real const tWStar = get_integration_weight( iGP ) * tDetJ;

            // add contribution to space-time volume
            tSpaceTimeVolume += tWStar;

            // get the normal from mesh
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the requested IQIs on the element
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get requested IQI
                const std::shared_ptr< IQI >& tReqIQI =
                        mSet->get_requested_elemental_IQIs_for_visualization()( iIQI );

                moris_index const tGlobalIqiIndex =
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
        }    // end for: each integration point

        // loop over IQI and divide each elemental IQI by space-time volume if requested
        if ( aAverageOutput )
        {
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get IQI global index
                moris_index const tGlobalIqiIndex = mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                // normalize by space-time volume
                ( *mSet->mSetElementalValues )( tLeaderFacetIndexInVisSet, tGlobalIqiIndex ) /= tSpaceTimeVolume;
            }
        }
    }

    //------------------------------------------------------------------------------

    real
    Element_Double_Sideset::compute_volume( mtk::Leader_Follower const aIsLeader )
    {
        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // get geometry interpolator
        Geometry_Interpolator* tIGGI = mSet->get_field_interpolator_manager( aIsLeader )->get_IG_geometry_interpolator();

        // loop over integration points
        real       tVolume           = 0;
        uint const tNumOfIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // set integration point for geometry interpolator
            tIGGI->set_space_time( get_leader_integration_point( iGP ) );

            // compute and add integration point contribution to volume
            tVolume += tIGGI->det_J() * get_integration_weight( iGP );
        }

        // return the volume value
        return tVolume / tIGGI->get_time_step();
    }

    //------------------------------------------------------------------------------

    void
    Element_Double_Sideset::compute_QI()
    {
        // get number of IQIs
        uint const tNumIQIs = mSet->get_number_of_requested_IQIs();

        // check for active IQIs
        if ( tNumIQIs == 0 )
        {
            return;
        }

        // set the leader/follower ig geometry interpolator physical/parametric space and time coefficients
        this->init_ig_geometry_interpolator();

        // loop over the integration points
        uint const tNumIntegPoints = get_number_of_integration_points();

        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            this->initialize_quadrature_point( iGP );

            // compute detJ of integration domain
            const real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

            // skip if detJ smaller than threshold
            if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
            {
                continue;
            }

            // compute integration point weight
            const real tWStar = get_integration_weight( iGP ) * tDetJ;

            // get the normal from mesh
            Matrix< DDRMat > const tNormal = get_leader_normal( iGP );

            // loop over the IQIs
            for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                // get requested IQI
                const std::shared_ptr< IQI >& tReqIQI = mSet->get_requested_IQIs()( iIQI );

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

}    // namespace moris::fem
