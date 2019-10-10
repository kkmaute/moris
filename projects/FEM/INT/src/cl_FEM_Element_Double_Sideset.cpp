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

            // get first corresponding node from master to slave
            //FIXME not sure it works, seems right
            moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode, tR );

            // get number of field interpolator and properties for master and slave
            uint tMasterNumFI   = mSet->get_number_of_field_interpolators();
            uint tSlaveNumFI    = mSet->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );
            uint tMasterNumProp = mSet->get_number_of_properties();
            uint tSlaveNumProp  = mSet->get_number_of_properties( mtk::Master_Slave::SLAVE );
            uint tMasterNumCM   = mSet->get_number_of_constitutive_models();
            uint tSlaveNumCM    = mSet->get_number_of_constitutive_models( mtk::Master_Slave::SLAVE );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_IWGs();

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_time( tMasterLocalIntegPoint );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0})
                    = tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set the ith integration point in the IG param space for slave IG geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_time( tSlaveLocalIntegPoint );

                // get global integration point  for the master integration cell
                Matrix< DDRMat > tMasterGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->map_integration_point( tMasterGlobalIntegPoint );

                // get global integration point for the slave integration cell
                Matrix< DDRMat > tSlaveGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->map_integration_point( tSlaveGlobalIntegPoint );

                // set evaluation point for master IP geometry interpolator
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_time( tMasterGlobalIntegPoint );

                // set evaluation point for master IP geometry interpolator
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_time( tSlaveGlobalIntegPoint );

                // set evaluation point for master and slave field interpolator
                for ( uint iFI = 0; iFI < tMasterNumFI; iFI++ )
                {
                    mSet->get_field_interpolators( mtk::Master_Slave::MASTER )( iFI )->set_space_time( tMasterGlobalIntegPoint );
                }
                for ( uint iFI = 0; iFI < tSlaveNumFI; iFI++ )
                {
                    mSet->get_field_interpolators( mtk::Master_Slave::SLAVE )( iFI )->set_space_time( tSlaveGlobalIntegPoint );
                }

                // reset properties
                for ( uint iProp = 0; iProp < tMasterNumProp; iProp++ )
                {
                    mSet->get_properties()( iProp )->reset_eval_flags();
                }
                for ( uint iProp = 0; iProp < tSlaveNumProp; iProp++ )
                {
                    mSet->get_properties( mtk::Master_Slave::SLAVE )( iProp )->reset_eval_flags();
                }

                // reset constitutive models
                for ( uint iCM = 0; iCM < tMasterNumCM; iCM++ )
                {
                    mSet->get_constitutive_models()( iCM )->reset_eval_flags();
                }
                for ( uint iCM = 0; iCM < tSlaveNumCM; iCM++ )
                {
                    mSet->get_constitutive_models( mtk::Master_Slave::SLAVE )( iCM )->reset_eval_flags();
                }

                // compute the integration point weight // fixme both side?
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // set the normal for the IWG
                    mSet->get_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    moris::Cell< Matrix< DDRMat > > tResidual;
                    mSet->get_IWGs()( iIWG )->compute_residual( tResidual );

                    // add contribution to jacobian from evaluation point
                    mSet->mResidual( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 0 ), mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 1 ) },
                                     { 0, 0 } ) += tWStar * tResidual( 0 );
                    mSet->mResidual( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 1, 0 ), mSet->get_IWG_res_dof_assembly_map()( iIWG )( 1, 1 ) },
                                     { 0, 0 } ) += tWStar * tResidual( 1 );
                }
            }
//            // print residual for check
//            print( mSet->mResidual, " mResidual " );
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

            // get first corresponding node from master to slave
            //FIXME not sure it works or provide the right thing
            moris::moris_index tSlaveNode = mCluster->get_left_vertex_pair( mMasterCell->get_vertices_on_side_ordinal( tMasterSideOrd )( 0 ) );

            // get rotation matrix from left to right
            Matrix< DDRMat> tR;
            rotation_matrix( mSet->get_IG_geometry_type(), tSlaveNode, tR );

            // get number of field interpolator and properties for master and slave
            uint tMasterNumFI   = mSet->get_number_of_field_interpolators();
            uint tSlaveNumFI    = mSet->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );
            uint tMasterNumProp = mSet->get_number_of_properties();
            uint tSlaveNumProp  = mSet->get_number_of_properties( mtk::Master_Slave::SLAVE );
            uint tMasterNumCM   = mSet->get_number_of_constitutive_models();
            uint tSlaveNumCM    = mSet->get_number_of_constitutive_models( mtk::Master_Slave::SLAVE );

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_IWGs();

            // loop over the integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get local integration point for the master integration cell
                Matrix< DDRMat > tMasterLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_time( tMasterLocalIntegPoint );

                // get local integration point for the slave integration cell
                Matrix< DDRMat > tSlaveLocalIntegPoint = tMasterLocalIntegPoint;
                tSlaveLocalIntegPoint({0,tMasterLocalIntegPoint.numel()-2},{0,0})
                    = tR * tMasterLocalIntegPoint({0,tSlaveLocalIntegPoint.numel()-2},{0,0}); //fixme better way?

                // set the ith integration point in the IG param space for slave IG geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_time( tSlaveLocalIntegPoint );

                // get global integration point  for the master integration cell
                Matrix< DDRMat > tMasterGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->map_integration_point( tMasterGlobalIntegPoint );

                // get global integration point for the slave integration cell
                Matrix< DDRMat > tSlaveGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::SLAVE )->map_integration_point( tSlaveGlobalIntegPoint );

                // set evaluation point for master IP geometry interpolator
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_time( tMasterGlobalIntegPoint );

                // set evaluation point for master IP geometry interpolator
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_time( tSlaveGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iFI = 0; iFI < tMasterNumFI; iFI++ )
                {
                    mSet->get_field_interpolators()( iFI )->set_space_time( tMasterGlobalIntegPoint );
                }
                for ( uint iFI = 0; iFI < tSlaveNumFI; iFI++ )
                {
                    mSet->get_field_interpolators( mtk::Master_Slave::SLAVE  )( iFI )->set_space_time( tSlaveGlobalIntegPoint );
                }

                // reset properties
                for ( uint iProp = 0; iProp < tMasterNumProp; iProp++ )
                {
                    mSet->get_properties()( iProp )->reset_eval_flags();
                }
                for ( uint iProp = 0; iProp < tSlaveNumProp; iProp++ )
                {
                    mSet->get_properties( mtk::Master_Slave::SLAVE )( iProp )->reset_eval_flags();
                }

                // reset constitutive models
                for ( uint iCM = 0; iCM < tMasterNumCM; iCM++ )
                {
                    mSet->get_constitutive_models()( iCM )->reset_eval_flags();
                }
                for ( uint iCM = 0; iCM < tSlaveNumCM; iCM++ )
                {
                    mSet->get_constitutive_models( mtk::Master_Slave::SLAVE )( iCM )->reset_eval_flags();
                }

                // compute the integration point weight // fixme both side?
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->det_J();

                // get the normal from mesh and set if for the IWG
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tMasterSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // set the normal for the IWG
                    mSet->get_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    moris::Cell< moris::Cell< Matrix< DDRMat > > > tJacobians;
                    mSet->get_IWGs()( iIWG )->compute_jacobian( tJacobians );
//                    print( tJacobians(0), "tJacobians" );

//                    real tPerturbation = 1E-6;
//                    Cell< Matrix< DDRMat > > tJacobiansFD;
//                    mSet->get_IWGs()( iIWG )->compute_jacobian_FD( tJacobiansFD,
//                                                                   mSet->get_IWG_field_interpolators( mtk::Master_Slave::MASTER )( iIWG ),
//                                                                   mSet->get_IWG_field_interpolators( mtk::Master_Slave::SLAVE )( iIWG ),
//                                                                   tPerturbation );
//                    print( tJacobiansFD(0), "tJacobiansFD" );

                    // loop over the IWG active dof types
                    uint tNumIWGDof = mSet->get_IWGs()( iIWG )->get_global_dof_type_list().size()
                                    + mSet->get_IWGs()( iIWG )->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
                    for ( uint iIWGFI = 0; iIWGFI < tNumIWGDof; iIWGFI++ )
                    {
                        // add contribution to jacobian from evaluation point
                        mSet->mJacobian( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 0 ),      mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 1 ) },
                                         { mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 0 ), mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 1 ) } )
                                       += tWStar * tJacobians( 0 )( iIWGFI );

                        mSet->mJacobian( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 1, 0 ), mSet->get_IWG_res_dof_assembly_map()( iIWG )( 1, 1 ) },
                                         { mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 0 ), mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 1 ) } )
                                       += tWStar * tJacobians( 1 )( iIWGFI );
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
