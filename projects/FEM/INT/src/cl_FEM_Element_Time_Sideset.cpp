#include <iostream>
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src
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
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords() );
            mSet->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_IG_geometry_interpolator()->set_param_coeff();

            // get the treated time ordinal
            //moris_index tTreatedTimeOrdinal = mCluster->mListOfTimeOrdinals( mCellIndexInCluster );

//            // get the side phys and param coords
//            // fixme to be removed with new approach
//            mSet->get_IG_geometry_interpolator()->build_time_side_time_phys_coeff( tTreatedTimeOrdinal );
//            mSet->get_IG_geometry_interpolator()->build_time_side_time_param_coeff( tTreatedTimeOrdinal );

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                // get the index of the residual dof type for the ith IWG in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_field_interpolator() );

                //get number of integration points
                uint tNumOfIntegPoints = mSet->get_num_integration_points();

                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // get integration point location in the reference volume
                    Matrix< DDRMat > tGlobalIntegPoint;
                    mSet->get_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint,
                                                                                 tGlobalIntegPoint );

                    // set integration point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                    }

                    // compute integration point weight
//                    real tSurfDetJ = mSet->get_IG_geometry_interpolator()->time_surf_det_J( tLocalIntegPoint );
//                    real tWStar = mSet->get_integration_weights()( iGP ) * tSurfDetJ;
                    real tWStar = mSet->get_integration_weights()( iGP )
                                * mSet->get_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    tTreatedIWG->compute_residual( tResidual, tIWGInterpolators );

                    // get location of computed residual in global element residual
                    uint startDof = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                    uint stopDof  = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                    // add contribution to residual from evaluation point
                    mSet->mResidual( { startDof, stopDof }, { 0, 0 } )
                        = mSet->mResidual( { startDof, stopDof }, { 0, 0 } ) + tResidual * tWStar;
                }
            }
//            // print residual for check
//            print( mCluster->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords() );
            mSet->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_IG_geometry_interpolator()->set_param_coeff();

            // get the treated time ordinal
            //moris_index tTreatedTimeOrdinal = mCluster->mListOfTimeOrdinals( mCellIndexInCluster );

//            // get the side phys and param coords
//            // fixme to remove with new implementation of GI
//            mSet->get_IG_geometry_interpolator()->build_time_side_time_phys_coeff( tTreatedTimeOrdinal );
//            mSet->get_IG_geometry_interpolator()->build_time_side_time_param_coeff( tTreatedTimeOrdinal );

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                // get the index of the residual dof type for the ith IWG in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_field_interpolator() );

                //get number of integration points
                uint tNumOfIntegPoints = mSet->get_num_integration_points();

               // loop over the integration points
               for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
               {
                   // get local integration point location
                   Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                   // get global integration point location
//                   Matrix< DDRMat > tGlobalIntegPoint = mSet->get_IG_geometry_interpolator()->time_surf_val( tLocalIntegPoint );
                   Matrix< DDRMat > tGlobalIntegPoint;
                   mSet->get_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint,
                                                                                tGlobalIntegPoint );

                   // set integration point for field interpolator
                   for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                   {
                       tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                   }

                   // compute integration point weight
//                   real tWStar = mSet->get_integration_weights()( iGP )
//                               * mSet->get_IG_geometry_interpolator()->time_surf_det_J( tLocalIntegPoint );
                   real tWStar = mSet->get_integration_weights()( iGP )
                               * mSet->get_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                   // compute jacobian at evaluation point
                   moris::Cell< Matrix< DDRMat > > tJacobians;
                   tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                   // get location of computed jacobian in global element residual rows
                   uint startIDof = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                   uint stopIDof  = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                   // loop over the IWG active dof types
                   for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                   {
                       // get the index of the active dof type
                       uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                       // get location of computed jacobian in global element residual columns
                       uint startJDof = mSet->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 0 );
                       uint stopJDof  = mSet->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 1 );

                       // add contribution to jacobian from evaluation point
                       mSet->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                           = mSet->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                           + tWStar * tJacobians( iIWGFI );
                    }
                }
            }
//            // print residual for check
//            print( mCluster->mJacobian, " mJacobian " );
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
