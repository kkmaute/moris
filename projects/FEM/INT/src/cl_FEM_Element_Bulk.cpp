#include <iostream>

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

#include "cl_MSI_Design_Variable_Interface.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        Element_Bulk::Element_Bulk( mtk::Cell const    * aCell,
                                    Set                * aSet,
                                    Cluster            * aCluster,
                                    moris::moris_index   aCellIndexInCluster )
        : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------
        Element_Bulk::~Element_Bulk(){}

//------------------------------------------------------------------------------
        void Element_Bulk::compute_residual()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME: enforced nodal weak bcs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute residual at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute off-diagonal jacobian for staggered solve
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_jacobian()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff ( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} );//fixme

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()
                    ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, "Element_Bulk::compute_jacobian_and_residual - not implemented." );
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_dRdp()
        {
            // FIXME all this should be in the FI Manager
            // get the geometry pdv values
            Matrix< DDRMat > tPdvValues = mMasterCell->get_vertex_coords();

            // reshape the matrix tPdvValues into a cell of vectors tPdvValueList
            moris::Cell< Matrix< DDRMat > > tPdvValueList( 2 );
            tPdvValueList( 0 ) = tPdvValues.get_column( 0 );
            tPdvValueList( 1 ) = tPdvValues.get_column( 1 );

            // get the pdv values from the MSI/GEN interface
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            mSet->mDesignVariableInterface->get_ip_pdv_value( mMasterCell->get_vertex_inds(),
                                                              { PDV::X_COORDINATE, PDV::Y_COORDINATE },
                                                              tPdvValueList,
                                                              tIsActiveDv );//FIXME

            // reshape the cell of vectors tPdvValueList into a matrix tPdvValues
            mSet->mDesignVariableInterface->reshape_pdv_values( tPdvValueList,
                                                                tPdvValues );

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_coeff( tPdvValues );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_coeff ( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme
            // END FIXME all this should be in the FI Manager

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()
                    ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set nodal weak BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->mInterpolationElement->get_weak_bcs() );

                    // compute dRdpMat at evaluation point
                    real tPerturbation = 1E-6;
                    moris::Cell< Matrix< DDRMat > > tdRdpMatFD;
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_material( tWStar,
                                                                                  tPerturbation,
                                                                                  tdRdpMatFD );

                    // compute dRdpGeo at evaluation point
                    moris::Cell< Matrix< DDRMat > > tdRdpGeoFD;
                    mSet->get_requested_IWGs()( iIWG )->compute_dRdp_FD_geometry( tWStar,
                                                                                  tPerturbation,
                                                                                  tIsActiveDv,
                                                                                  tdRdpGeoFD );
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_QI()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // FIXME create a cell of IQI values
            moris::Cell< Matrix< DDRMat > > tQIValues( tNumIQIs );

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IQI
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute QI at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_QI( tWStar );
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_dQIdp()
        {
            // FIXME all this should be in the FI Manager
            // get the geometry pdv values
            Matrix< DDRMat > tPdvValues = mMasterCell->get_vertex_coords();

            // reshape the matrix tPdvValues into a cell of vectors tPdvValueList
            moris::Cell< Matrix< DDRMat > > tPdvValueList( 2 );
            tPdvValueList( 0 ) = tPdvValues.get_column( 0 );
            tPdvValueList( 1 ) = tPdvValues.get_column( 1 );

            // get the pdv values from the MSI/GEN interface
            moris::Cell< Matrix< DDSMat > > tIsActiveDv;
            mSet->mDesignVariableInterface->get_ip_pdv_value( mMasterCell->get_vertex_inds(),
                                                              { PDV::X_COORDINATE, PDV::Y_COORDINATE },
                                                              tPdvValueList, tIsActiveDv );//FIXME

            // reshape the cell of vectors tPdvValueList into a matrix tPdvValues
            mSet->mDesignVariableInterface->reshape_pdv_values( tPdvValueList,
                                                                tPdvValues );

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_coeff( tPdvValues );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_coeff ( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme
            // END FIXME all this should be in the FI Manager

            // get number of IWGs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()
                    ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute dQIdpMat at evaluation point
                    real tPerturbation = 1E-6;
                    Matrix< DDRMat > tdQIdpMatFD;
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_material( tWStar,
                                                                                   tPerturbation,
                                                                                   tdQIdpMatFD );

                    // compute dQIdpGeo at evaluation point
                    Matrix< DDRMat > tdQIdpGeoFD;
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdp_FD_geometry( tWStar,
                                                                                   tPerturbation,
                                                                                   tIsActiveDv,
                                                                                   tdQIdpGeoFD );
                    // fill list of dQIdp
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_dQIdu()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff ( mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} );//fixme

            // get number of IQIs
            uint tNumIQIs = mSet->get_number_of_requested_IQIs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()
                    ->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // loop over the IQIs
                for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
                {
                    // reset IWG
                    mSet->get_requested_IQIs()( iIQI )->reset_eval_flags();

                    // compute dQIdu at evaluation point
                    mSet->get_requested_IQIs()( iIQI )->compute_dQIdu( tWStar );
                }
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_quantity_of_interest_global( const uint aMeshIndex,
                                                                enum  vis::Output_Type aOutputType )
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                *( mSet->mSetGlobalValues ) += tQIValue( 0 ) * tWStar;
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_quantity_of_interest_nodal( const uint aMeshIndex,
                                                               enum vis::Output_Type aOutputType )
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ) );
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // get the vertices
            moris::Cell< mtk::Vertex * > tVertices = mMasterCell->get_vertex_pointers();

            // loop over the vertices
            uint tNumNodes = tVertices.size();
            for( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
            {
                // get the ith vertex coordinates in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint = mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ).get_row( iVertex );
                tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) = mCluster->mInterpolationElement->get_time()( 0 );
                tGlobalIntegPoint = trans( tGlobalIntegPoint );

                // set vertex coordinates for field interpolator
                mSet->get_field_interpolator_manager()->set_space_time( tGlobalIntegPoint );

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                // FIXME add up on shared node and divide or overwrite
                (*mSet->mSetNodalValues)( tVertices( iVertex )->get_index(), 0 ) += tQIValue( 0 );

                mSet->mSetNodalCounter( tVertices( iVertex )->get_index(), 0 ) += 1;
            }
        }

//------------------------------------------------------------------------------
        void Element_Bulk::compute_quantity_of_interest_elemental( const uint aMeshIndex,
                                                                   enum vis::Output_Type aOutputType )
        {
            // set the IG geometry interpolator physical space and time coefficients
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the IP geometry interpolator param space and time coefficients
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ) );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

                // reset the requested IQI
                mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQIValue;
                mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                // FIXME assemble on the set here or inside the compute QI?
                ( *mSet->mSetElementalValues )( mSet->mCellAssemblyMap( aMeshIndex )( mMasterCell->get_index() ), 0 )
                += tQIValue( 0 ) * tWStar / tNumIntegPoints;
            }
        }

//------------------------------------------------------------------------------
        real Element_Bulk::compute_volume( mtk::Master_Slave aIsMaster )
        {
            // set the IG geometry interpolator physical space and time coefficients
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_time_coeff(  mCluster->mInterpolationElement->get_time() );

            // set the IP geometry interpolator param space and time coefficients
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
                ->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster ) );
            mSet->get_field_interpolator_manager()
                ->get_IG_geometry_interpolator()
              ->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                mSet->get_field_interpolator_manager()
                    ->get_IG_geometry_interpolator()
                    ->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J()
                         * mSet->get_integration_weights()( iGP );
            }

            // return the volume value
            return tVolume;
        }

//------------------------------------------------------------------------------

//        real Element::compute_integration_error(
//                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
//        {
//            // create field interpolation rule
//            Interpolation_Rule tFieldInterpolationRule( this->get_geometry_type(),
//                                                        Interpolation_Type::LAGRANGE,
//                                                        this->get_interpolation_order() );
//            // <- add second type in order
//            //    to interpolate in space and time
//
//            // create geometry interpolation rule
//            Interpolation_Rule tGeometryInterpolationRule( this->get_geometry_type(),
//                                                           Interpolation_Type::LAGRANGE,
//                                                           this->get_interpolation_order() );
//
//            // create integration rule
//            Integration_Rule tIntegration_Rule( this->get_geometry_type(),
//                                                Integration_Type::GAUSS,
//                                                this->get_auto_integration_order() );
//
//            // set number of fields
//            uint tNumberOfFields = 1;
//
//            // create interpolator
//            Interpolator tInterpolator( this->get_node_coords(),
//                                        tNumberOfFields,
//                                        tFieldInterpolationRule,
//                                        tGeometryInterpolationRule,
//                                        tIntegration_Rule );
//
//            // get number of points
//            auto tNumberOfIntegrationPoints
//                = tInterpolator.get_number_of_integration_points();
//
//            real aError = 0.0;
//
//            mIWG->create_matrices( &tInterpolator );
//
//            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
//            {
//                 //evaluate shape function at given integration point
//                aError += mIWG->compute_integration_error( mPdofValues,
//                                                           aFunction,
//                                                           k )
//                        * tInterpolator.get_det_J( k )
//                        * tInterpolator.get_integration_weight( k );
//            }
//
//            //std::cout << "Element error " << aError << std::endl;
//            mIWG->delete_matrices();
//
//            return aError;
//        }

//------------------------------------------------------------------------------

//        real Element::compute_element_average_of_scalar_field()
//        {
//
//            // create field interpolation rule
//            Interpolation_Rule tFieldInterpolationRule( this->get_geometry_type(),
//                                                        Interpolation_Type::LAGRANGE,
//                                                        this->get_interpolation_order() );
//            // <- add second type in order
//            //    to interpolate in space and time
//
//            // create geometry interpolation rule
//            Interpolation_Rule tGeometryInterpolationRule( this->get_geometry_type(),
//                                                           Interpolation_Type::LAGRANGE,
//                                                           mtk::Interpolation_Order::LINEAR );
//
//            // create integration rule
//            Integration_Rule tIntegration_Rule( this->get_geometry_type(),
//                                                Integration_Type::GAUSS,
//                                                this->get_auto_integration_order() );
//
//            // set number of fields
//            uint tNumberOfFields = 1;
//
//            // create interpolator
//            Interpolator tInterpolator( this->get_node_coords(),
//                                        tNumberOfFields,
//                                        tFieldInterpolationRule,
//                                        tGeometryInterpolationRule,
//                                        tIntegration_Rule );
//
//            // get number of points
//            auto tNumberOfIntegrationPoints
//                = tInterpolator.get_number_of_integration_points();
//
//            mIWG->create_matrices( &tInterpolator );
//
//            real aValue  = 0.0;Cell< Field_Interpolator* > tFieldInterpolators
//            real tWeight = 0.0;
//
//            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
//            {
//                real tScale = tInterpolator.get_integration_weight( k )
//                            * tInterpolator.get_det_J( k );
//
//                aValue += mIWG->interpolate_scalar_at_point( mNodalWeakBCs, k )
//                        * tScale;
//
//                tWeight += tScale;
//
//            }
//
//            // close IWG object
//            mIWG->delete_matrices();
//
//            return aValue / tWeight;
//
//        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
