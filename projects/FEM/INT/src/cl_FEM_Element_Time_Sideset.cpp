#include <iostream>
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Time_Sideset::Element_Time_Sideset( mtk::Cell const  * aCell,
                                                    Set              * aElementBlock,
                                                    Cluster          * aCluster,
                                                    moris::moris_index aCellIndexInCluster) : Element( aCell, aElementBlock, aCluster, aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------

        Element_Time_Sideset::~Element_Time_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_residual()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords() );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_param_coeff();

            // get number of field interpolator and properties
            uint tNumFI   = mSet->get_number_of_field_interpolators();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // get integration point location in the reference volume
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for IP geometry interpolator
                mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator()->set_space_time( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iFI = 0; iFI < tNumFI; iFI++ )
                {
                    mSet->mMasterFIManager->get_field_interpolators_for_type( mSet->mMasterDofTypes( iFI )( 0 ) )
                                          ->set_space_time( tGlobalIntegPoint );
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual(  tWStar );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );

                }
            }
//            // print residual for check
//            print( mCluster->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords() );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_param_coeff();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point location
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // get global integration point location
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for IP geometry interpolator
                mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator()->set_space_time( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                mSet->get_field_interpolator_manager()->set_space_time( tGlobalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
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

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Time_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
