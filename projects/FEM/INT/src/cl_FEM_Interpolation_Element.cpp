/*
 * cl_FEM_Interpolation_Element.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Element.hpp"                    //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"      //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_MSI_Design_Variable_Interface.hpp"   //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                   //FEM/INT/src
#include "cl_FEM_Set.hpp"                   //FEM/INT/src
#include "cl_FEM_Model.hpp"                   //FEM/INT/src

#include "cl_SOL_Dist_Vector.hpp"

#include "fn_isfinite.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Interpolation_Element::Interpolation_Element( const Element_Type                aElementType,
                                                  const Cell< const mtk::Cell * > & aInterpolationCell,
                                                        moris::Cell< Node_Base* > & aNodes,
                                                        Set                       * aSet )
        : MSI::Equation_Object( aSet ),
          mSet( aSet ),
          mElementType( aElementType )
        {
            // fill the master interpolation cell
            mMasterInterpolationCell = aInterpolationCell( 0 );

            // get vertices from cell
            moris::Cell< mtk::Vertex* > tVertices = mMasterInterpolationCell->get_vertex_pointers();

            // get number of vertices from cell
            uint tNumOfVertices = tVertices.size();

            // assign node object
            mNodeObj.resize( 1 );
            mNodeObj( 0 ).resize( tNumOfVertices, nullptr );

            // fill master node objects
            for( uint iVertex = 0; iVertex < tNumOfVertices; iVertex++)
            {
                mNodeObj( 0 )( iVertex ) = aNodes( tVertices( iVertex )->get_index() );
            }

//            // set size of Weak BCs
//            mNodalWeakBCs.set_size( tNumOfVertices, 1 );

            // switch on the element type
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // fill the slave interpolation cell
                mSlaveInterpolationCell  = aInterpolationCell( 1 );

                // get vertices from cell
                moris::Cell< mtk::Vertex* > tSlaveVertices  = mSlaveInterpolationCell->get_vertex_pointers();

                // get number of vertices from cell
                uint tNumOfSlaveVertices = tSlaveVertices.size();

                // assign node object
                mNodeObj.resize( 2 );
                mNodeObj( 1 ).resize( tNumOfSlaveVertices , nullptr );

                // fill slave node objects
                for( uint iVertex = 0; iVertex < tNumOfSlaveVertices; iVertex++)
                {
                    mNodeObj( 1 )( iVertex ) = aNodes( tSlaveVertices( iVertex )->get_index() );
                }
            }
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::set_cluster( std::shared_ptr< fem::Cluster > aCluster,
                                                 const uint                      aMeshIndex )
        {
            // if mesh index is 0 (i.e., forward analysis mesh, IG mesh)
            if( aMeshIndex == 0 )
            {
                // fem cluster with index 0 should be set only once and shall not be changed
                MORIS_ASSERT( !( mFemCluster.size() >= 1 ),
                              "Interpolation_Element::set_cluster() - first fem cluster is already set");
            }

            // get max size for fem cluster list
            sint tSize = std::max( ( sint )mFemCluster.size(), ( sint )aMeshIndex + 1 );

            // resize fem cluster list
            mFemCluster.resize( tSize );

            // add the fem cluster to the list
            mFemCluster( aMeshIndex ) = aCluster;
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::set_field_interpolators_coefficients()
        {
            // dof field interpolators------------------------------------------

            // get number of master dof types
             uint tMasterNumDofTypes = mSet->get_dof_type_list().size();

             // loop on the master dof types
             for( uint iDOF = 0; iDOF < tMasterNumDofTypes; iDOF++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_dof_type_list()( iDOF );

                 // get the pdof values for the ith dof type group
                 Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
                 this->get_my_pdof_values( mPdofValues, tDofTypeGroup, tCoeff_Original );

                 // reshape tCoeffs into the order the cluster expects them
                 Matrix< DDRMat > tCoeff;
                 this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                 // set field interpolator coefficients
                 mSet->get_field_interpolator_manager()
                     ->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );

                 if( mSet->get_time_continuity() )
                 {
                     // get the pdof values for the ith dof type group
                     Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
                     this->get_my_pdof_values( mPreviousPdofValues, tDofTypeGroup, tCoeff_Original );

                     // reshape tCoeffs into the order the cluster expects them
                     Matrix< DDRMat > tCoeff;
                     this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                     // set field interpolator coefficients
                     mSet->get_field_interpolator_manager_previous_time()
                         ->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
                 }
             }

             // get number of slave dof types
             uint tSlaveNumDofTypes = mSet->get_dof_type_list( mtk::Master_Slave::SLAVE ).size();

             // loop on the slave dof types
             for( uint iDOF = 0; iDOF < tSlaveNumDofTypes; iDOF++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup
                 = mSet->get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF );

                 // get the pdof values for the ith dof type group
                 Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
                 this->get_my_pdof_values( mPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

                 // reshape tCoeffs into the order the cluster expects them
                 Matrix< DDRMat > tCoeff;
                 this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );

                 // set the field coefficients
                 mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )
                     ->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
             }

             // dv field interpolators------------------------------------------

             // get number of master dv types
             uint tMasterNumDvTypes = mSet->get_dv_type_list().size();

             // loop on the master dv types
             for( uint iDv = 0; iDv < tMasterNumDvTypes; iDv++ )
             {
                 // get the dv type group
                 moris::Cell< PDV_Type > tDvTypeGroup
                 = mSet->get_dv_type_list()( iDv );

                // get the pdv values for the ith dv type group
                Cell< Matrix< DDRMat > > tCoeff_Original;
                mSet->mDesignVariableInterface->get_ip_pdv_value( mMasterInterpolationCell->get_vertex_inds(),
                                                                  tDvTypeGroup,
                                                                  tCoeff_Original );

                // reshape tCoeffs into the order the FI expects them
                Matrix< DDRMat > tCoeff;
                mSet->mDesignVariableInterface
                    ->reshape_pdv_values( tCoeff_Original, tCoeff );

                // set field interpolator coefficients
                mSet->get_field_interpolator_manager()
                    ->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
             }

             // get number of slave dv types
             uint tSlaveNumDvTypes
             = mSet->get_dv_type_list( mtk::Master_Slave::SLAVE ).size();

             // loop on the slave dv types
             for( uint iDv = 0; iDv < tSlaveNumDvTypes; iDv++ )
             {
                 // get the dv type group
                 moris::Cell< PDV_Type > tDvTypeGroup
                 = mSet->get_dv_type_list( mtk::Master_Slave::SLAVE )( iDv );

                 // get the pdv values for the ith dv type group
                 Cell< Matrix< DDRMat > > tCoeff_Original;
                 mSet->mDesignVariableInterface->get_ip_pdv_value( mSlaveInterpolationCell->get_vertex_inds(),
                                                                   tDvTypeGroup,
                                                                   tCoeff_Original );

                 // reshape tCoeffs into the order the FI expects them
                 Matrix< DDRMat > tCoeff;
                 mSet->mDesignVariableInterface
                     ->reshape_pdv_values( tCoeff_Original, tCoeff );

                 // set the field coefficients
                 mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )
                     ->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
             }

             // geometry interpolators------------------------------------------
             // set the IP geometry interpolator physical space and time coefficients for the master
             mSet->get_field_interpolator_manager()
                 ->get_IP_geometry_interpolator()
                 ->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
             mSet->get_field_interpolator_manager()
                 ->get_IP_geometry_interpolator()
                 ->set_time_coeff( this->get_time() );

             // if double sideset
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave
                 mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )
                     ->get_IP_geometry_interpolator()
                     ->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )
                     ->get_IP_geometry_interpolator()
                     ->set_time_coeff( this->get_time() );
             }

             // if time sideset
             if( mElementType == fem::Element_Type::TIME_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the previous
                 mSet->get_field_interpolator_manager_previous_time( mtk::Master_Slave::MASTER )
                     ->get_IP_geometry_interpolator()
                     ->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
                 mSet->get_field_interpolator_manager_previous_time( mtk::Master_Slave::MASTER )
                     ->get_IP_geometry_interpolator()
                     ->set_time_coeff( this->get_previous_time() );
             }
         }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_jacobian()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // init the jacobian
            mSet->initialize_mJacobian();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_jacobian();

//            // check that jacobian is finite
//            // FIXME
//            bool tIsFinite = true;
//            uint tNumRows = mSet->get_jacobian().n_rows();
//            uint tNumCols = mSet->get_jacobian().n_cols();
//            for( uint iRow = 0; iRow < tNumRows; iRow++ )
//            {
//                for( uint iCol = 0; iCol < tNumCols; iCol++ )
//                {
//                    Matrix< DDRMat> tValue = {{mSet->get_jacobian()( iRow, iCol )}};
//                    if( !isfinite( tValue ) )
//                    {
//                        tIsFinite = false;
//                    }
//                }
//            }
//            if( !tIsFinite )
//            {
//                std::cout<<"Interpolation_Element::compute_jacobian - non finite values in jacobian."<<std::endl;
//            }
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_residual()
        {
            //Fixme do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // init the residual
            mSet->initialize_mResidual();

            // init the jacobian
            mSet->initialize_mJacobian();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            if( mSet->mEquationModel->get_is_forward_analysis() )
            {
                // FIXME should not be like this
                mSet->set_IWG_field_interpolator_managers();

                // set cluster for stabilization parameter
                mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

                // ask cluster to compute residual
                mFemCluster( 0 )->compute_residual();
            }
            else
            {
                // FIXME should not be like this
                mSet->set_IQI_field_interpolator_managers();

                // set cluster for stabilization parameter
                mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

                // ask cluster to compute jacobian
                mFemCluster( 0 )->compute_dQIdu();
            }
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_jacobian_and_residual()
         {
             //Fixme do this only once
             this->compute_my_pdof_values();

             // if time continuity set
             if ( mSet->get_time_continuity() )
             {
                 // compute pdof values for previous time step
                 // FIXME do this only once
                 this->compute_previous_pdof_values();
             }

             // init the jacobian
             mSet->initialize_mJacobian();

             // init the residual
             mSet->initialize_mResidual();

             // set the field interpolators coefficients
             this->set_field_interpolators_coefficients();

             // FIXME should not be like this
             mSet->set_IWG_field_interpolator_managers();

             // set cluster for stabilization parameter
             mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

             MORIS_ERROR( false, "Interpolation_Element::compute_jacobian_and_residual(), function not tested and works only non staggered");

             // ask cluster to compute jacobian and residual
             mFemCluster( 0 )->compute_jacobian_and_residual();
         }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_dRdp()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // init geo pdv assembly map
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // init dRdp
            mSet->initialize_mdRdpMat();
            mSet->initialize_mdRdpGeo( mFemCluster( 0 ) );

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IWG_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dRdp();
        }

//------------------------------------------------------------------------------

        void Interpolation_Element::compute_dQIdp_explicit()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // init geo pdv assembly map
            mSet->create_geo_pdv_assembly_map( mFemCluster( 0 ) );

            // init dRdp
            mSet->initialize_mdQIdpMat();
            mSet->initialize_mdQIdpGeo( mFemCluster( 0 ) );

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dQIdp_explicit();

            //----------------------------------------------------------------------------------------

            for( uint Ik = 0; Ik < mSet->mdQIdp( 0 ).size(); Ik++ )
            {
                Cell< enum PDV_Type > tRequestedIPDvTypes;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_requested_dv_types( tRequestedIPDvTypes );

                // get vertices from cell
                Matrix< IndexMat > tVerticesInds = mMasterInterpolationCell->get_vertex_inds();

                //FIXME add Slave

                moris::Cell< moris::Matrix< IdMat > > tTypeListOfLocalToGlobalIds;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_dv_ids_for_type_and_ind( tVerticesInds,
                                                              tRequestedIPDvTypes,
                                                              tTypeListOfLocalToGlobalIds );   //FIXME add type and nodei inds

                moris::uint tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                moris::Matrix< IdMat > tLocalToGlobalIds( tCounter, 1, moris::gNoIndex );

                tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tLocalToGlobalIds( { tCounter, tTypeListOfLocalToGlobalIds( Ii ).numel() -1 },{ 0, 0 } )
                            = tTypeListOfLocalToGlobalIds( Ii ).matrix_data();

                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                mEquationSet->get_equation_model()
                            ->get_explicit_dQidu()
                            ->sum_into_global_values( tLocalToGlobalIds,
                                                      mSet->mdQIdp( 0 )( Ik ),
                                                      Ik );
            }

            for( uint Ik = 0; Ik < mSet->mdQIdp( 1 ).size(); Ik++ )
            {
                Cell< enum PDV_Type > tRequestedIGDvTypes;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_requested_dv_types( tRequestedIGDvTypes );

                moris::Cell< moris::Matrix< IdMat > > tTypeListOfLocalToGlobalIds;

                // get vertices from cell
                Matrix< IndexMat > tVerticesInds = mFemCluster( 0 )->get_mesh_cluster()
                                                                   ->get_vertex_indices_in_cluster();

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ig_dv_ids_for_type_and_ind( tVerticesInds,
                                                              tRequestedIGDvTypes,
                                                              tTypeListOfLocalToGlobalIds );   //FIXME add type and nodei inds

                moris::uint tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                moris::Matrix< IdMat > tLocalToGlobalIds( tCounter, 1, moris::gNoIndex );

                tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tLocalToGlobalIds( { tCounter, tTypeListOfLocalToGlobalIds( Ii ).numel() -1 },{ 0, 0 } )
                            = tTypeListOfLocalToGlobalIds( Ii ).matrix_data();

                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                mEquationSet->get_equation_model()
                            ->get_implicit_dQidu()
                            ->sum_into_global_values( tLocalToGlobalIds,
                                                      mSet->mdQIdp( 1 )( Ik ),
                                                      Ik );
            }
        }

//-------------------------------------------------------------------------------------------------

        void Interpolation_Element::compute_dQIdp()
        {
            this->compute_dRdp();

            moris::Cell< Matrix< DDRMat > > & tdRdp = mEquationSet->get_drdp();

            this->compute_my_adjoint_values();

    //        moris::Matrix< DDRMat > tMyAdjointValues( mAdjointPdofValues.numel(), 1, 0.0 );
    //
    //        // get number of master dof types
    //        uint tMasterNumDofTypes = mEquationSet->get_dof_type_list().size();
    //
    //        // loop on the master dof types
    //        for( uint iDOF = 0; iDOF < tMasterNumDofTypes; iDOF++ )
    //        {
    //            // get the ith dof type group
    //            moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_dof_type_list()( iDOF );
    //
    //            // get the pdof values for the ith dof type group. Outer cell are multi-vector entries
    //            Cell< Cell< Matrix< DDRMat > > > tCoeff_Original;
    //            this->get_my_pdof_values( mAdjointPdofValues, tDofTypeGroup, tCoeff_Original );
    //
    //        //FIXME reshape correctly
    ////            // reshape tCoeffs into the order the cluster expects them
    ////            Matrix< DDRMat > tCoeff;
    ////            this->reshape_pdof_values( tCoeff_Original( 0 ), tCoeff );
    //        }

            for( uint Ik = 0; Ik < mAdjointPdofValues.size(); Ik++ )
            {
                moris::Matrix< DDRMat > tLocalIPdQiDp = trans( mAdjointPdofValues( Ik ) ) * tdRdp( 0 );

                Cell< enum PDV_Type > tRequestedIPDvTypes;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_requested_dv_types( tRequestedIPDvTypes );

                // get vertices from cell
                Matrix< IndexMat > tVerticesInds = mMasterInterpolationCell->get_vertex_inds();

                //FIXME add Slave

                moris::Cell< moris::Matrix< IdMat > > tTypeListOfLocalToGlobalIds;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_dv_ids_for_type_and_ind( tVerticesInds,
                                                              tRequestedIPDvTypes,
                                                              tTypeListOfLocalToGlobalIds );   //FIXME add type and nodei inds

                moris::uint tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                moris::Matrix< IdMat > tLocalToGlobalIds( tCounter, 1, moris::gNoIndex );

                tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tLocalToGlobalIds( { tCounter, tTypeListOfLocalToGlobalIds( Ii ).numel() -1 },{ 0, 0 } )
                            = tTypeListOfLocalToGlobalIds( Ii ).matrix_data();

                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                mEquationSet->get_equation_model()->get_implicit_dQidu()->sum_into_global_values( tLocalToGlobalIds,
                                                                                       tLocalIPdQiDp,
                                                                                       Ik );
            }

            for( uint Ik = 0; Ik < mAdjointPdofValues.size(); Ik++ )
            {
                moris::Matrix< DDRMat > tLocalIGdQiDp = trans( mAdjointPdofValues( Ik ) ) * tdRdp( 1 );

                Cell< enum PDV_Type > tRequestedIGDvTypes;

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ip_requested_dv_types( tRequestedIGDvTypes );

                moris::Cell< moris::Matrix< IdMat > > tTypeListOfLocalToGlobalIds;

                // get vertices from cell
                Matrix< IndexMat > tVerticesInds = mFemCluster( 0 )->get_mesh_cluster()
                                                                   ->get_vertex_indices_in_cluster();

                mEquationSet->get_equation_model()
                            ->get_design_variable_interface()
                            ->get_ig_dv_ids_for_type_and_ind( tVerticesInds,
                                                              tRequestedIGDvTypes,
                                                              tTypeListOfLocalToGlobalIds );   //FIXME add type and nodei inds

                moris::uint tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                moris::Matrix< IdMat > tLocalToGlobalIds( tCounter, 1, moris::gNoIndex );

                tCounter = 0;

                for( uint Ii = 0; Ii < tTypeListOfLocalToGlobalIds.size(); Ii++ )
                {
                    tLocalToGlobalIds( { tCounter, tTypeListOfLocalToGlobalIds( Ii ).numel() -1 },{ 0, 0 } )
                            = tTypeListOfLocalToGlobalIds( Ii ).matrix_data();

                    tCounter += tTypeListOfLocalToGlobalIds( Ii ).numel();
                }

                mEquationSet->get_equation_model()->get_implicit_dQidu()->sum_into_global_values( tLocalToGlobalIds,
                                                                                       tLocalIGdQiDp,
                                                                                       Ik );
            }
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_dQIdu()
        {
            // compute pdof values
            //FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute jacobian
            mFemCluster( 0 )->compute_dQIdu();
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_QI()
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            mSet->initialize_mQI();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IQI_field_interpolator_managers();

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            // ask cluster to compute quantity of interest
            mFemCluster( 0 )->compute_QI();
        }

//------------------------------------------------------------------------------
        void
        Interpolation_Element::compute_quantity_of_interest( const uint            aMeshIndex,
                                                             enum vis::Output_Type aOutputType,
                                                             enum vis::Field_Type  aFieldType )
        {
            // compute pdof values
            // FIXME do this only once
            this->compute_my_pdof_values();

            // if time continuity set
            if ( mSet->get_time_continuity() )
            {
                // compute pdof values for previous time step
                // FIXME do this only once
                this->compute_previous_pdof_values();
            }

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->get_IQI_for_vis( aOutputType )
                ->set_field_interpolator_manager( mSet->get_field_interpolator_manager() );

            // set cluster for stabilization parameter
            mSet->set_IQI_cluster_for_stabilization_parameters( mFemCluster( 0 ).get() );

            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                mSet->get_IQI_for_vis( aOutputType )
                    ->set_field_interpolator_manager( mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE ),
                                                      mtk::Master_Slave::SLAVE );
            }

            if( aFieldType == vis::Field_Type::NODAL )
            {
                // get the master vertices indices on the mesh cluster
                moris::Cell< moris_index > tVertexIndices
                = mFemCluster( aMeshIndex )->get_vertex_indices_in_cluster();

                // get the master vertices local coordinates on the mesh cluster
                moris::Matrix<moris::DDRMat> tVertexLocalCoords
                = mFemCluster( aMeshIndex )->get_vertices_local_coordinates_wrt_interp_cell();

                // get number of vertices on the treated mesh cluster
                uint tNumNodes = tVertexIndices.size();

                // loop over the vertices on the treated mesh cluster
                for( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
                {
                    // get the ith vertex coordinates in the IP param space
                    Matrix< DDRMat > tGlobalIntegPoint = tVertexLocalCoords.get_row( iVertex );
                    tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                    tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) = -1.0;
                    tGlobalIntegPoint = trans( tGlobalIntegPoint );

                    // set vertex coordinates for field interpolator
                    mSet->get_field_interpolator_manager()->set_space_time( tGlobalIntegPoint );

                    // reset the requested IQI
                    mSet->get_IQI_for_vis( aOutputType )->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIValue;
                    mSet->get_IQI_for_vis( aOutputType )->compute_QI( tQIValue );

                    // fill in the nodal set values
                    ( * mSet->mSetNodalValues )( tVertexIndices( iVertex ), 0 )
                             = tQIValue( 0 );
                }
            }
            else
            {
                // ask cluster to compute quantity of interest
                mFemCluster( aMeshIndex )->compute_quantity_of_interest( aMeshIndex, aOutputType, aFieldType );
            }
        }

//------------------------------------------------------------------------------
        real Interpolation_Element::compute_volume()
        {
            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // ask cluster to compute colume
            return mFemCluster( 0 )->compute_volume();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
