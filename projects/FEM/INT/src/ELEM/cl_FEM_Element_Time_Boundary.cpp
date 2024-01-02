/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Time_Boundary.cpp
 *
 */

#include <iostream>
//FEM/INT/src
#include "cl_FEM_Element_Time_Boundary.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
//FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Element_Time_Boundary::Element_Time_Boundary(
                mtk::Cell const  * aCell,
                Set              * aElementBlock,
                Cluster          * aCluster,
                moris::moris_index aCellIndexInCluster )
        : Element( aCell, aElementBlock, aCluster, aCellIndexInCluster )
        {}

        //------------------------------------------------------------------------------

        Element_Time_Boundary::~Element_Time_Boundary(){}

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::init_ig_geometry_interpolator(
                uint aTimeOrdinal )
        {
            // get param space time
            real tTimeParamCoeff = 2.0 * aTimeOrdinal - 1.0;

            // get geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get physical space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords = mLeaderCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords( 1, 1,
                    mCluster->mInterpolationElement->get_time()( aTimeOrdinal ) );

            // get leader parametric space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            Matrix< DDRMat > tIGParamTimeCoords( 1, 1, tTimeParamCoeff );

            // set physical space and current time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and current time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::init_ig_geometry_interpolator(
                uint               aTimeOrdinal,
                Matrix< DDSMat > & aGeoLocalAssembly )
        {
            // get param space time
            real tTimeParamCoeff = 2.0 * aTimeOrdinal - 1.0;

            // get geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get physical space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords = mLeaderCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords( 1, 1,
                    mCluster->mInterpolationElement->get_time()( aTimeOrdinal ) );

            // get leader parametric space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            Matrix< DDRMat > tIGParamTimeCoords( 1, 1, tTimeParamCoeff );

            // get the local cluster assembly indices
            if( mSet->get_geo_pdv_assembly_flag() )
            {
                // get the vertices indices for IG element
                Matrix< IndexMat > tVertexIndices = mLeaderCell->get_vertex_inds();

                // get the requested geo pdv types
                Vector < enum PDV_Type > tGeoPdvType;
                mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

                // get local assembly indices
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVertexIndices,
                        tGeoPdvType,
                        aGeoLocalAssembly );
            }

            // set physical space and current time coefficients for IG element GI
            tIGGI->set_space_coeff( tIGPhysSpaceCoords );
            tIGGI->set_time_coeff(  tIGPhysTimeCoords );

            // set parametric space and current time coefficients for IG element GI
            tIGGI->set_space_param_coeff( tIGParamSpaceCoords );
            tIGGI->set_time_param_coeff(  tIGParamTimeCoords );
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get number of integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // loop over integration points
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute detJ of integration domain
                    real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // skip if detJ smaller than threshold
                    if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                    {
                        continue;
                    }

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // get requested IWG
                        const std::shared_ptr< IWG > & tReqIWG =
                                mSet->get_requested_IWGs()( iIWG );

                        // reset IWG
                        tReqIWG->reset_eval_flags();

                        // compute Jacobian at evaluation point
                        tReqIWG->compute_residual( tWStar );

                        // compute Jacobian at evaluation point
                        // compute off-diagonal Jacobian for staggered solve
                        ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_jacobian()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get number of integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // loop over integration points
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get local integration point location
                    const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute detJ of integration domain
                    real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // skip if detJ smaller than threshold
                    if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                    {
                        continue;
                    }

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // get requested IWG
                        const std::shared_ptr< IWG > & tReqIWG =
                                mSet->get_requested_IWGs()( iIWG );

                        // reset IWG
                        tReqIWG->reset_eval_flags();

                        // compute Jacobian at evaluation point
                        ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_jacobian_and_residual()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get number of integration points integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // loop over integration points
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get local integration point location
                    const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute detJ of integration domain
                    real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // skip if detJ smaller than threshold
                    if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                    {
                        continue;
                    }

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // get requested IWG
                        const std::shared_ptr< IWG > & tReqIWG =
                                mSet->get_requested_IWGs()( iIWG );

                        // reset IWG
                        tReqIWG->reset_eval_flags();

                        if( mSet->mEquationModel->get_is_forward_analysis() )
                        {
                            // compute residual at evaluation point
                            tReqIWG->compute_residual( tWStar );
                        }

                        // compute Jacobian at evaluation point
                        ( this->*m_compute_jacobian )( tReqIWG, tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_dRdp()
        {
            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // check for active IWGs
            if (tNumIWGs == 0)
            {
                return;
            }

            // get number of integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                Matrix< DDSMat > tGeoLocalAssembly;
                this->init_ig_geometry_interpolator( iTimeBoundary, tGeoLocalAssembly );

                // loop over integration points
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    const Matrix< DDRMat > & tLocalIntegPoint =
                        mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute detJ of integration domain
                    real tDetJ = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // skip if detJ smaller than threshold
                    if ( tDetJ < Geometry_Interpolator::sDetJInvJacLowerLimit )
                    {
                        continue;
                    }

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // get requested IWG
                        const std::shared_ptr< IWG > & tReqIWG =
                                mSet->get_requested_IWGs()( iIWG );

                        // reset IWG
                        tReqIWG->reset_eval_flags();

                        // compute dRdp at evaluation point
                        Vector< Matrix< IndexMat > > tVertexIndices( 0 );
                        ( this->*m_compute_dRdp )( tReqIWG, tWStar, tGeoLocalAssembly, tVertexIndices );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_QI()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // get number of integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // loop over integration points
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
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IQIs
                    for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // get requested IQI
                        const std::shared_ptr< IQI > & tReqIQI =
                                mSet->get_requested_IQIs()( iIQI );

                        // reset IQI
                        tReqIQI->reset_eval_flags();

                        // compute QI at evaluation point
                        tReqIQI->compute_QI( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_dQIdu()
        {
            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // get number of integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // loop over integration points
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
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IQIs
                    for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // get requested IQI
                        const std::shared_ptr< IQI > & tReqIQI =
                                mSet->get_requested_IQIs()( iIQI );

                        // reset IQI
                        tReqIQI->reset_eval_flags();

                        // compute dQIdu at evaluation point
                        ( this->*m_compute_dQIdu )( tReqIQI, tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_dQIdp_explicit()
        {
            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // check for active IQIs
            if (tNumIQIs == 0)
            {
                return;
            }

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();

            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                Matrix< DDSMat > tGeoLocalAssembly;
                this->init_ig_geometry_interpolator( iTimeBoundary, tGeoLocalAssembly );

                // loop over integration points
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
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) * tDetJ;

                    // loop over the IQIs
                    for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                    {
                        // get requested IQI
                        const std::shared_ptr< IQI > & tReqIQI =
                                mSet->get_requested_IQIs()( iIQI );

                        // reset IQI
                        tReqIQI->reset_eval_flags();

                        // compute dQIdp at evaluation point
                        Vector< Matrix< IndexMat > > tVertexIndices( 0 );
                        ( this->*m_compute_dQIdp )( tReqIQI, tWStar, tGeoLocalAssembly, tVertexIndices );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

