#include <iostream>
#include "cl_FEM_Element_Time_Boundary.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

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
                uint                              aTimeOrdinal,
                moris::Cell< Matrix< DDSMat > > & aIsActiveDv )
        {
            // get param space time
            real tTimeParamCoeff = 2.0 * aTimeOrdinal - 1.0;

            // get geometry interpolator
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get physical space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGPhysSpaceCoords = mMasterCell->get_vertex_coords();
            Matrix< DDRMat > tIGPhysTimeCoords( 1, 1,
                    mCluster->mInterpolationElement->get_time()( aTimeOrdinal ) );

            // get master parametric space and current and previous time coordinates for IG element
            Matrix< DDRMat > tIGParamSpaceCoords =
                    mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster );
            Matrix< DDRMat > tIGParamTimeCoords( 1, 1, tTimeParamCoeff );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // determine if there are IG pdvs
            if ( tGeoPdvType.size() )
            {
                MORIS_ERROR( false, "init_ig_geometry_interpolator - pdv not handle so far on time boundary" );
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
            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                moris::Cell< Matrix< DDSMat > > tIsActiveDv;
                this->init_ig_geometry_interpolator( iTimeBoundary, tIsActiveDv );

                // get number of IWGs
                uint tNumIWGs = mSet->get_number_of_requested_IWGs();

                // loop over integration points
                uint tNumIntegPoints = mSet->get_number_of_integration_points();
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) *
                            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // reset IWG
                        mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                        // compute jacobian at evaluation point
                        mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                        // compute jacobian at evaluation point
                        // compute off-diagonal jacobian for staggered solve
                        mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_jacobian()
        {
            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                moris::Cell< Matrix< DDSMat > > tIsActiveDv;
                this->init_ig_geometry_interpolator( iTimeBoundary, tIsActiveDv );

                // get number of IWGs
                uint tNumIWGs = mSet->get_number_of_requested_IWGs();

                // loop over integration points
                uint tNumIntegPoints = mSet->get_number_of_integration_points();

                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get local integration point location
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) *
                            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // reset IWG
                        mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                        // compute jacobian at evaluation point
                        mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void Element_Time_Boundary::compute_jacobian_and_residual()
        {
            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set physical and parametric space and time coefficients for IG element
                moris::Cell< Matrix< DDSMat > > tIsActiveDv;
                this->init_ig_geometry_interpolator( iTimeBoundary, tIsActiveDv );

                // get number of IWGs
                uint tNumIWGs = mSet->get_number_of_requested_IWGs();

                // loop over integration points
                uint tNumIntegPoints = mSet->get_number_of_integration_points();
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get local integration point location
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()->
                            set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP ) *
                            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                    // loop over the IWGs
                    for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                    {
                        // reset IWG
                        mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                        // compute residual at evaluation point
                        mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                        // compute jacobian at evaluation point
                        mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
