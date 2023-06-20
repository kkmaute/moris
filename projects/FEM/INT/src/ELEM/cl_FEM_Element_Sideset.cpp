/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Sideset.cpp
 *
 */

#include <iostream>
//FEM/INT/src
#include "cl_FEM_Element_Sideset.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Model.hpp"
// FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset(
                mtk::Cell const    * aCell,
                Set                * aSet,
                Cluster            * aCluster,
                moris::moris_index   aCellIndexInCluster)
                                : Element( aCell, aSet, aCluster, aCellIndexInCluster )
                                  {}

        //------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset(){}

        //------------------------------------------------------------------------------

        void Element_Sideset::init_ig_geometry_interpolator( uint aSideOrdinal )
        {
            // get geometry interpolator for IG element
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mLeaderCell->get_cell_physical_coords_on_side_ordinal( aSideOrdinal );
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster, aSideOrdinal );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::init_ig_geometry_interpolator(
                uint               aSideOrdinal,
                Matrix< DDSMat > & aGeoLocalAssembly )
        {
            // get geometry interpolator for IG element
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get leader physical space and time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords =
                    mLeaderCell->get_cell_physical_coords_on_side_ordinal( aSideOrdinal );
            Matrix< DDRMat > tIGPhysTimeCoords =
                    mCluster->mInterpolationElement->get_time();

            // get leader parametric space and time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster, aSideOrdinal );

            // FIXME not true if time is not linear
            Matrix< DDRMat > tIGParamTimeCoords = { { -1.0 }, { 1.0 } };

            // get the local cluster assembly indices
            if( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices for IG element
                Matrix< IndexMat > tVertexIndices =
                        mLeaderCell->get_vertices_ind_on_side_ordinal( aSideOrdinal );

                // get the requested geo pdv types
                moris::Cell < enum PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute residual at integration point
                    tReqIWG->compute_residual( tWStar );

                    // compute Jacobian at evaluation point
                    // compute off-diagonal Jacobian for staggered solve
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set BCs
                    tReqIWG->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute Jacobian at evaluation point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian_and_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IWGs or IQIs
            if (tNumIWGs == 0 && tNumIQIs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set BCs
                    tReqIWG->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    if( mSet->mEquationModel->get_is_forward_analysis() )
                    {
                        // compute residual at integration point
                        tReqIWG->compute_residual( tWStar );
                    }

                    // compute Jacobian at evaluation point
                    ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                }

                if( ( !mSet->mEquationModel->get_is_forward_analysis() ) && ( tNumIQIs > 0 ) )
                {
                    // loop over the IQIs
                    for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // get requested IQI
                        const std::shared_ptr< IQI > & tReqIQI =
                                mSet->get_requested_IQIs()( iIQI );

                        // reset IQI
                        tReqIQI->reset_eval_flags();

                        // set the normal for the IQI
                        tReqIQI->set_normal( tNormal );

                        // compute dQIdu at evaluation point
                        ( this->*m_compute_dQIdu )( tReqIQI, tWStar );
                    }
                }

                mSet->mFemModel->mSideSetsGaussPoints++;
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_dRdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tSideOrd, tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

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

        void Element_Sideset::compute_QI()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
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

        void Element_Sideset::compute_dQIdu()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIQI->set_normal( tNormal );

                    // compute dQIdu at evaluation point
                    ( this->*m_compute_dQIdu )( tReqIQI, tWStar );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_dQIdp_explicit()
        {
            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tSideOrd, tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IQI
                    tReqIQI->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIQI->set_normal( tNormal );

                    // compute dQIdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dQIdp )( tReqIQI, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_dRdp_and_dQIdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IWGs and IQIs
            if ( tNumIWGs == 0 && tNumIQIs == 0 )
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            Matrix< DDSMat > tGeoLocalAssembly;
            this->init_ig_geometry_interpolator( tSideOrd, tGeoLocalAssembly );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // get requested IWG
                    const std::shared_ptr< IWG > & tReqIWG =
                            mSet->get_requested_IWGs()( iIWG );

                    // reset IWG
                    tReqIWG->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    tReqIWG->set_nodal_weak_bcs(
                            mCluster->mInterpolationElement->get_weak_bcs() );

                    // set the normal for the IWG
                    tReqIWG->set_normal( tNormal );

                    // compute dRdp at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                }

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
                            mSet->get_requested_IQIs()( iIQI );

                    // reset IWG
                    tReqIQI->reset_eval_flags();

                    // set the normal for the IWG
                    tReqIQI->set_normal( tNormal );

                    // compute dQIdpMat at evaluation point
                    moris::Cell< Matrix< IndexMat > > tVertexIndices( 0 );
                    ( this->*m_compute_dQIdp )( tReqIQI, tWStar, tGeoLocalAssembly, tVertexIndices );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_quantity_of_interest_global( const uint aFemMeshIndex )
        {
            // get number of active local IQIs
            uint tNumLocalIQIs = mSet->get_number_of_requested_global_IQIs_for_visualization();

            // check that some IQIs need to be evaluated
            if ( tNumLocalIQIs == 0)
            {
                return;
            }

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over IQI
                for( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
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

                } // end for: each IQI on the current element

            } // end for: each integration point

        } // end function: Element_Sideset::compute_quantity_of_interest_global()

        //------------------------------------------------------------------------------

        void Element_Sideset::compute_quantity_of_interest_elemental(
                const uint aFemMeshIndex,
                const bool aAverageOutput )
        {
            // get number of active local IQIs
            uint tNumLocalIQIs = mSet->get_number_of_requested_elemental_IQIs_for_visualization();

            // check that some IQIs need to be evaluated
            if ( tNumLocalIQIs == 0)
            {
                return;
            }

            // get the VIS mesh index
            uint tVisMeshIndex = aFemMeshIndex - 1;

            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // get the leader IG cell's index in the VIS mesh
            moris_index tIgCellIndex = mLeaderCell->get_index();
            moris_index tFacetIndexInVisSet = mSet->mFacetAssemblyMap( tVisMeshIndex )( tIgCellIndex, tSideOrd );
            MORIS_ASSERT( tFacetIndexInVisSet > -1, 
                    "FEM::Element_Sideset::compute_quantity_of_interest_elemental() - "
                    "Facet not part of VIS facet assembly map." );

            // FIXME: this seems inefficient and should be done elsewhere
            // set unused IQI values to not be NAN
            for ( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
            {
                // get IQI global index
                moris_index tGlobalIqiIndex =
                        mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                if ( ( *mSet->mSetElementalValues )( tFacetIndexInVisSet, tGlobalIqiIndex ) == std::numeric_limits< real >::quiet_NaN() )
                {
                    ( *mSet->mSetElementalValues )( tFacetIndexInVisSet, tGlobalIqiIndex ) = 0.0;
                }
            }

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            // initialize space - time volume
            real tSpaceTimeVolume = 0.0;

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                const Matrix< DDRMat > & tLocalIntegPoint =
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

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeaderCell, tSideOrd );

                // loop over IQI
                for( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get requested IQI
                    const std::shared_ptr< IQI > & tReqIQI =
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
                    ( *mSet->mSetElementalValues )( tFacetIndexInVisSet, tGlobalIqiIndex ) += tWStar * tQIElemental( 0 );

                } // end for: each IQI to be evaluated on cluster

            } // end for: each integration point

            // loop over IQI and divide each elemental IQI by space-time volume if requested
            if ( aAverageOutput )
            {
                for( uint iIQI = 0; iIQI < tNumLocalIQIs; iIQI++ )
                {
                    // get IQI global index
                    moris_index tGlobalIqiIndex = mSet->get_requested_elemental_IQIs_global_indices_for_visualization()( iIQI );

                    // divide by space-time volume
                    ( *mSet->mSetElementalValues )( tFacetIndexInVisSet, tGlobalIqiIndex ) /= tSpaceTimeVolume;
                }
            }

        } // end function: Element_Sideset::compute_quantity_of_interest_elemental()

        //------------------------------------------------------------------------------

        real Element_Sideset::compute_volume( mtk::Leader_Follower aIsLeader )
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mLeaderListOfSideOrdinals( mCellIndexInCluster );

            // set physical and parametric space and time coefficients for IG element
            this->init_ig_geometry_interpolator( tSideOrd );

            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // get IG geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager( aIsLeader )->get_IG_geometry_interpolator();

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                tIGGI->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += tIGGI->det_J() * mSet->get_integration_weights()( iGP );
            }

            // get time step
            real tTimeStep = tIGGI->get_time_step();

            // return the volume value
            return tVolume/tTimeStep;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

