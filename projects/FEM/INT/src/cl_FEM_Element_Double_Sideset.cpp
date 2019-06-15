#include <iostream>

#include "cl_FEM_Element_Double_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "fn_FEM_Rotation_Matrix.hpp"        //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset( mtk::Cell const  * aLeftIGCell,
                                                        mtk::Cell const  * aRightIGCell,
                                                        Set              * aSet,
                                                        Cluster          * aCluster,
                                                        moris::moris_index aCellIndexInCluster ) : Element( aLeftIGCell,
                                                                                                            aRightIGCell,
                                                                                                            aSet,
                                                                                                            aCluster,
                                                                                                            aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_residual()
        {
            // get treated side ordinal on the left and on the right
            uint tLeftSideOrd  = mCluster->mLeftListOfSideOrdinals( mCellIndexInCluster );
            uint tRightSideOrd = mCluster->mRightListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for left integration cell
            mSet->get_left_IG_geometry_interpolator()->set_space_coeff( mLeftCell->get_cell_physical_coords_on_side_ordinal( tLeftSideOrd ) );
            mSet->get_left_IG_geometry_interpolator()->set_time_coeff(  mCluster->mTime );

            // set the geometry interpolator physical space and time coefficients for right integration cell
            mSet->get_right_IG_geometry_interpolator()->set_space_coeff( mRightCell->get_cell_physical_coords_on_side_ordinal( tRightSideOrd ) );
            mSet->get_right_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for left integration cell
            mSet->get_left_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_left_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                            tLeftSideOrd ) );
            mSet->get_left_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for left integration cell
            mSet->get_right_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_right_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                              tRightSideOrd ) );
            mSet->get_right_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over the IWGs
             for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
             {
                 // get the treated IWG
                 IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                 // get the number of active Dof_type for the IWG
                 uint tNumOfIWGActiveDof = tTreatedIWG->get_active_dof_types().size();

                 // get the active field interpolators for the IWG for the left interpolation cell
                 Cell< Field_Interpolator* > tLeftIWGInterpolators  = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_left_field_interpolator() );
                 // get the active field interpolators for the IWG for the right interpolation cell
                 Cell< Field_Interpolator* > tRightIWGInterpolators = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_right_field_interpolator() );
                 //get number of integration points
                 uint tNumOfIntegPoints = mSet->get_num_integration_points();

                 for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                 {
                     // get local integration point for the left integration cell
                     Matrix< DDRMat > tLeftLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                     // get first corresponding node from left to right
                     //FIXME not sure it works, seems right
                     moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mLeftCell->get_vertices_on_side_ordinal( tLeftSideOrd )( 0 ) );
                     //std::cout<<tSlaveNode<<std::endl;

                     // get rotation matrix from left to right
                     Matrix< DDRMat> tR = rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode );

                     // get local integration point for the right integration cell
                     Matrix< DDRMat > tRightLocalIntegPoint = tLeftLocalIntegPoint;
                     tRightLocalIntegPoint({0,tRightLocalIntegPoint.numel()-2},{0,0})
                         = tR * tLeftLocalIntegPoint({0,tRightLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                     // get global integration point  for the left integration cell
                     Matrix< DDRMat > tLeftGlobalIntegPoint = mSet->get_left_IG_geometry_interpolator()->map_integration_point( tLeftLocalIntegPoint );
                     //print(tLeftGlobalIntegPoint,"tLeftGlobalIntegPoint");

                     // get global integration point for the right integration cell
                     Matrix< DDRMat > tRightGlobalIntegPoint = mSet->get_right_IG_geometry_interpolator()->map_integration_point( tRightLocalIntegPoint );
                     //print(tRightGlobalIntegPoint,"tRightGlobalIntegPoint");

                     // set integration point for left and right field interpolators
                     for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                     {
                         tLeftIWGInterpolators( iIWGFI )->set_space_time( tLeftGlobalIntegPoint );
                         tRightIWGInterpolators( iIWGFI )->set_space_time( tRightGlobalIntegPoint );
                     }

                     // compute the integration point weight // fixme both side?
                     real tWStar = mSet->get_integration_weights()( iGP )
                                 * mSet->get_left_IG_geometry_interpolator()->det_J( tLeftLocalIntegPoint );

                     // get the normal from mesh and set if for the IWG
                     Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeftCell, tLeftSideOrd, tLeftLocalIntegPoint );
                     tTreatedIWG->set_normal( tNormal );

                     // compute residual at integration point
                     Matrix< DDRMat > tResidual;
                     tTreatedIWG->compute_residual( tResidual, tLeftIWGInterpolators, tRightIWGInterpolators );

                     // fixme does it work as is?
                     // get the index of the residual dof type for the ith IWG in the list of element dof type
                     uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                     // get location of computed residual in global element residual
                     uint startDof = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                     uint stopDof  = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                     // add contribution to jacobian from evaluation point
                     mSet->mResidual( { startDof, stopDof }, { 0, 0 } )
                         = mSet->mResidual( { startDof, stopDof }, { 0, 0 } ) + tResidual * tWStar;
                 }
            }
//            // print residual for check
            print( mSet->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_jacobian()
        {
            // get treated side ordinal on the left and on the right
            uint tLeftSideOrd  = mCluster->mLeftListOfSideOrdinals( mCellIndexInCluster );
            uint tRightSideOrd = mCluster->mRightListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for left integration cell
            mSet->get_left_IG_geometry_interpolator()->set_space_coeff( mLeftCell->get_cell_physical_coords_on_side_ordinal( tLeftSideOrd ) );
            mSet->get_left_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator physical space and time coefficients for right integration cell
            mSet->get_right_IG_geometry_interpolator()->set_space_coeff( mRightCell->get_cell_physical_coords_on_side_ordinal( tRightSideOrd ) );
            mSet->get_right_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for left integration cell
            mSet->get_left_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_left_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                            tLeftSideOrd ) );
            mSet->get_left_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for left integration cell
            mSet->get_right_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_right_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                              tRightSideOrd ) );
            mSet->get_right_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over the IWGs
             for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
             {
                 // get the treated IWG
                 IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                 // get the index of the residual dof type for the ith IWG in the list of element dof type
                 uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                 // get the number of active Dof_type for the IWG
                 Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                 uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                 // get the active field interpolators for the IWG for the left interpolation cell
                 Cell< Field_Interpolator* > tLeftIWGInterpolators  = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_left_field_interpolator() );
                 // get the active field interpolators for the IWG for the right interpolation cell
                 Cell< Field_Interpolator* > tRightIWGInterpolators = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_right_field_interpolator() );
                 //get number of integration points
                 uint tNumOfIntegPoints = mSet->get_num_integration_points();

                 for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                 {
                     // get local integration point for the left integration cell
                     Matrix< DDRMat > tLeftLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                     // get first corresponding node from left to right
                     //FIXME not sure it works or provide the right thing
                     moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mLeftCell->get_vertices_on_side_ordinal( tLeftSideOrd )( 0 ) );

                     // get rotation matrix from left to right
                     Matrix< DDRMat> tR = rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode );

                     // get local integration point for the right integration cell
                     Matrix< DDRMat > tRightLocalIntegPoint = tLeftLocalIntegPoint;
                     tRightLocalIntegPoint({0,tLeftLocalIntegPoint.numel()-2},{0,0})
                         = tR * tLeftLocalIntegPoint({0,tRightLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                     // get global integration point  for the left integration cell
                     Matrix< DDRMat > tLeftGlobalIntegPoint = mSet->get_left_IG_geometry_interpolator()->map_integration_point( tLeftLocalIntegPoint );

                     // get global integration point for the right integration cell
                     Matrix< DDRMat > tRightGlobalIntegPoint = mSet->get_right_IG_geometry_interpolator()->map_integration_point( tRightLocalIntegPoint );

                     // set integration point for left and right field interpolators
                     for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                     {
                         tLeftIWGInterpolators( iIWGFI )->set_space_time( tLeftGlobalIntegPoint );
                         tRightIWGInterpolators( iIWGFI )->set_space_time( tRightGlobalIntegPoint );
                     }

                     // compute the integration point weight // fixme both side?
                     real tWStar = mSet->get_integration_weights()( iGP )
                                 * mSet->get_left_IG_geometry_interpolator()->det_J( tLeftLocalIntegPoint );

                     // get the normal from mesh and set if for the IWG
                     Matrix< DDRMat > tNormal = mCluster->get_side_normal( mLeftCell, tLeftSideOrd, tLeftLocalIntegPoint );
                     tTreatedIWG->set_normal( tNormal );

                     // compute residual at integration point
                     Cell< Matrix< DDRMat > > tJacobians;
                     tTreatedIWG->compute_jacobian( tJacobians, tLeftIWGInterpolators, tRightIWGInterpolators );
//                     print(tJacobians(0),"tJacobians");

//                     real tPerturbation = 1E-6;
//                     Cell< Matrix< DDRMat > > tJacobiansFD;
//                     tTreatedIWG->compute_jacobian_FD( tJacobiansFD,
//                                                       tLeftIWGInterpolators, tRightIWGInterpolators,
//                                                       tPerturbation );
//                     print(tJacobiansFD(0),"tJacobiansFD");

                     // fixme does it work as is?
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
//            // print jacobian for check
//            print( mSet->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Double_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
