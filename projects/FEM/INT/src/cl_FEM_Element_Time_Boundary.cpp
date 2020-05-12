#include <iostream>
#include "cl_FEM_Element_Time_Boundary.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        Element_Time_Boundary::Element_Time_Boundary( mtk::Cell const  * aCell,
                                                      Set              * aElementBlock,
                                                      Cluster          * aCluster,
                                                      moris::moris_index aCellIndexInCluster )
        : Element( aCell, aElementBlock, aCluster, aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------
        Element_Time_Boundary::~Element_Time_Boundary(){}

//------------------------------------------------------------------------------
        void Element_Time_Boundary::init_ig_geometry_interpolator( uint aTimeOrdinal )
        {
            // get param space time
            real tTimeParamCoeff = 2.0 * aTimeOrdinal - 1.0;

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_coeff( {{ mCluster->mInterpolationElement->get_time()( aTimeOrdinal ) }} );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_param_coeff( {{ tTimeParamCoeff }} );
        }

//------------------------------------------------------------------------------
        void Element_Time_Boundary::compute_residual()
        {
            // loop over time boundaries
            for ( uint iTimeBoundary = 0; iTimeBoundary < 2; iTimeBoundary++ )
            {
                // get param space time
                real tTimeParamCoeff = 2.0 * iTimeBoundary - 1.0;

                // set IG geometry interpolator
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // get number of IWGs
                uint tNumIWGs = mSet->get_number_of_requested_IWGs();

                // loop over integration points
                uint tNumIntegPoints = mSet->get_number_of_integration_points();
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()
                        ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP )
                                * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

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

                // set IG geometry interpolator
                this->init_ig_geometry_interpolator( iTimeBoundary );

                // get number of IWGs
                uint tNumIWGs = mSet->get_number_of_requested_IWGs();

                // loop over integration points
                uint tNumIntegPoints = mSet->get_number_of_integration_points();
                for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
                {
                    // get local integration point location
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // set evaluation point for interpolators (FIs and GIs)
                    mSet->get_field_interpolator_manager()
                        ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                    // compute integration point weight
                    real tWStar = tTimeParamCoeff * mSet->get_integration_weights()( iGP )
                                * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

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
            MORIS_ERROR( false, " Element_Time_Boundary::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
