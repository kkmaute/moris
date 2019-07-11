#include <iostream>

#include "cl_FEM_Element_Double_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "fn_FEM_Rotation_Matrix.hpp"        //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset( mtk::Cell const  * aMasterIGCell,
                                                        mtk::Cell const  * aSlaveIGCell,
                                                        Set              * aSet,
                                                        Cluster          * aCluster,
                                                        moris::moris_index aCellIndexInCluster ) : Element( aMasterIGCell,
                                                                                                            aSlaveIGCell,
                                                                                                            aSet,
                                                                                                            aCluster,
                                                                                                            aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_residual()
        {
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for master integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( tMasterSideOrd ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff(  mCluster->mTime );

            // set the geometry interpolator physical space and time coefficients for slave integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveCell->get_cell_physical_coords_on_side_ordinal( tSlaveSideOrd ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for master integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                                             tMasterSideOrd,
                                                                                                                                                             mtk::Master_Slave::MASTER ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for slave integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                                            tSlaveSideOrd,
                                                                                                                                                            mtk::Master_Slave::SLAVE ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over the IWGs
             for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
             {
                 // get the treated IWG
                 IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                 // get the number of active Dof_type for the IWG
                 uint tNumOfIWGActiveDof = tTreatedIWG->get_active_dof_types().size();

                 // get the active field interpolators for the IWG for the master interpolation cell
                 Cell< Field_Interpolator* > tMasterIWGInterpolators  = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                           mSet->get_field_interpolator( mtk::Master_Slave::MASTER ) );
                 // get the active field interpolators for the IWG for the slave interpolation cell
                 Cell< Field_Interpolator* > tSlaveIWGInterpolators = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_field_interpolator( mtk::Master_Slave::SLAVE ) );
                 //get number of integration points
                 uint tNumOfIntegPoints = mSet->get_num_integration_points();

                 for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                 {
                     // get local integration point for the master integration cell
                     Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                     // get first corresponding node from master to slave
                     //FIXME not sure it works, seems right
                     moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );
                     //std::cout<<tSlaveNode<<std::endl;

                     // get rotation matrix from left to right
                     Matrix< DDRMat> tR = rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode );

                     // get local integration point for the slave integration cell
                     Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                     tSlaveLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0})
                         = tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                     // get global integration point  for the master integration cell
                     Matrix< DDRMat > tMasterGlobalIntegPoint = mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->map_integration_point( tMasterLocalIntegPoint );
                     //print(tLeftGlobalIntegPoint,"tLeftGlobalIntegPoint");

                     // get global integration point for the slave integration cell
                     Matrix< DDRMat > tSlaveGlobalIntegPoint = mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->map_integration_point( tSlaveLocalIntegPoint );
                     //print(tRightGlobalIntegPoint,"tRightGlobalIntegPoint");

                     // set integration point for master and slave field interpolators
                     for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                     {
                         tMasterIWGInterpolators( iIWGFI )->set_space_time( tMasterGlobalIntegPoint );
                         tSlaveIWGInterpolators( iIWGFI )->set_space_time( tSlaveGlobalIntegPoint );
                     }

                     // compute the integration point weight // fixme both side?
                     real tWStar = mSet->get_integration_weights()( iGP )
                                 * mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->det_J( tMasterLocalIntegPoint );

                     // get the normal from mesh and set if for the IWG
                     Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd, tMasterLocalIntegPoint );
                     tTreatedIWG->set_normal( tNormal );

                     // compute residual at integration point
                     Matrix< DDRMat > tResidual;
                     tTreatedIWG->compute_residual( tResidual, tMasterIWGInterpolators, tSlaveIWGInterpolators );

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
            // get treated side ordinal on the master and on the slave
            uint tMasterSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );
            uint tSlaveSideOrd  = mCluster->mSlaveListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for master integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( tMasterSideOrd ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator physical space and time coefficients for slave integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveCell->get_cell_physical_coords_on_side_ordinal( tSlaveSideOrd ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for master integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                                             tMasterSideOrd,
                                                                                                                                                              mtk::Master_Slave::MASTER ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // set the geometry interpolator param space and time coefficients for slave integration cell
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                                            tSlaveSideOrd,
                                                                                                                                                            mtk::Master_Slave::SLAVE ) );
            mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

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

                 // get the active field interpolators for the IWG for the master interpolation cell
                 Cell< Field_Interpolator* > tMasterIWGInterpolators  = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                           mSet->get_field_interpolator( mtk::Master_Slave::MASTER ) );
                 // get the active field interpolators for the IWG for the slave interpolation cell
                 Cell< Field_Interpolator* > tSlaveIWGInterpolators = mSet->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                         mSet->get_field_interpolator( mtk::Master_Slave::SLAVE ) );
                 //get number of integration points
                 uint tNumOfIntegPoints = mSet->get_num_integration_points();

                 for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                 {
                     // get local integration point for the master integration cell
                     Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                     // get first corresponding node from master to slave
                     //FIXME not sure it works or provide the right thing
                     moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

                     // get rotation matrix from left to right
                     Matrix< DDRMat> tR = rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode );

                     // get local integration point for the slave integration cell
                     Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                     tSlaveLocalIntegPoint({0,tMasterLocalIntegPoint.numel()-2},{0,0})
                         = tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                     // get global integration point  for the master integration cell
                     Matrix< DDRMat > tMasterGlobalIntegPoint = mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->map_integration_point( tMasterLocalIntegPoint );

                     // get global integration point for the slave integration cell
                     Matrix< DDRMat > tSlaveGlobalIntegPoint = mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->map_integration_point( tSlaveLocalIntegPoint );

                     // set integration point for master and slave field interpolators
                     for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                     {
                         tMasterIWGInterpolators( iIWGFI )->set_space_time( tMasterGlobalIntegPoint );
                         tSlaveIWGInterpolators( iIWGFI )->set_space_time( tSlaveGlobalIntegPoint );
                     }

                     // compute the integration point weight // fixme both side?
                     real tWStar = mSet->get_integration_weights()( iGP )
                                 * mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->det_J( tMasterLocalIntegPoint );

                     // get the normal from mesh and set if for the IWG
                     Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd, tMasterLocalIntegPoint );
                     tTreatedIWG->set_normal( tNormal );

                     // compute residual at integration point
                     Cell< Matrix< DDRMat > > tJacobians;
                     tTreatedIWG->compute_jacobian( tJacobians, tMasterIWGInterpolators, tSlaveIWGInterpolators );
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
