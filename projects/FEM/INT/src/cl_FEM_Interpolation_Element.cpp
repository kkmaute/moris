/*
 * cl_FEM_Interpolation_Element.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Element.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"   //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src

#include "cl_FEM_Cluster.hpp"                   //FEM/INT/src

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

            // set size of Weak BCs
            mNodalWeakBCs.set_size( tNumOfVertices, 1 );

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
                 Cell< Matrix< DDRMat > > tCoeff_Original;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original );

                 // reshape tCoeffs into the order the cluster expects them
                 Matrix< DDRMat > tCoeff;
                 this->reshape_pdof_values( tCoeff_Original, tCoeff );

                 // set field interpolator coefficients
                 mSet->mMasterFIManager->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
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
                 Cell< Matrix< DDRMat > > tCoeff_Original;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

                 // reshape tCoeffs into the order the cluster expects them
                 Matrix< DDRMat > tCoeff;
                 this->reshape_pdof_values( tCoeff_Original, tCoeff );

                 // set the field coefficients
                 mSet->mSlaveFIManager->set_coeff_for_type( tDofTypeGroup( 0 ), tCoeff );
             }

             // dv field interpolators------------------------------------------

             // get number of master dv types
             uint tMasterNumDvTypes = mSet->get_dv_type_list().size();

             // loop on the master dv types
             for( uint iDv = 0; iDv < tMasterNumDvTypes; iDv++ )
             {
                 // get the dv type group
                 moris::Cell< MSI::Dv_Type > tDvTypeGroup
                 = mSet->get_dv_type_list()( iDv );

                 // FIXME get the pdv values for the dv type group
                 Matrix< DDRMat > tCoeff;
                 //this->get_my_pdv_values( tDvTypeGroup, tCoeff );

                 // set field interpolator coefficients
                 mSet->mMasterFIManager->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
             }

             // get number of slave dv types
             uint tSlaveNumDvTypes
             = mSet->get_dv_type_list( mtk::Master_Slave::SLAVE ).size();

             // loop on the slave dv types
             for( uint iDv = 0; iDv < tSlaveNumDvTypes; iDv++ )
             {
                 // get the dv type group
                 moris::Cell< MSI::Dv_Type > tDvTypeGroup
                 = mSet->get_dv_type_list( mtk::Master_Slave::SLAVE )( iDv );

                 // FIXME get the pdv values for the dv type group
                 Matrix< DDRMat > tCoeff;
                 //this->get_my_pdv_values( tDvTypeGroup, tCoeff, mtk::Master_Slave::SLAVE );

                 // set the field coefficients
                 mSet->mSlaveFIManager->set_coeff_for_type( tDvTypeGroup( 0 ), tCoeff );
             }

             // geometry interpolators------------------------------------------

             // set the IP geometry interpolator physical space and time coefficients for the master
             mSet->mMasterFIManager->get_IP_geometry_interpolator()->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
             mSet->mMasterFIManager->get_IP_geometry_interpolator()->set_time_coeff( this->mTime );

             // if double sideset
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave
                 mSet->mSlaveFIManager->get_IP_geometry_interpolator()->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->mSlaveFIManager->get_IP_geometry_interpolator()->set_time_coeff( this->mTime );
             }
         }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_jacobian()
        {
             //Fixme do this only once
             this->compute_my_pdof_values();

             // init the jacobian
             mSet->initialize_mJacobian();

             // set the field interpolators coefficients
             this->set_field_interpolators_coefficients();

             // FIXME should not be like this
             mSet->set_IWG_field_interpolator_managers();

             // ask cluster to compute jacobian
             mFemCluster( 0 )->compute_jacobian();
         }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_residual()
        {
            //Fixme do this only once
            this->compute_my_pdof_values();

            // init the residual
            mSet->initialize_mResidual();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();

            // ask cluster to compute residual
            mFemCluster( 0 )->compute_residual();
        }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_jacobian_and_residual()
         {
             //Fixme do this only once
             this->compute_my_pdof_values();

             // init the jacobian
             mSet->initialize_mJacobian();

             // init the residual
             mSet->initialize_mResidual();

             // set the field interpolators coefficients
             this->set_field_interpolators_coefficients();

             // FIXME should not be like this
             mSet->set_IWG_field_interpolator_managers();

             //
             MORIS_ERROR( false, "Interpolation_Element::compute_jacobian_and_residual(), function not tested and works only non staggered");

             // ask cluster to compute jacobian and residual
             mFemCluster( 0 )->compute_jacobian_and_residual();
         }

//------------------------------------------------------------------------------
        void Interpolation_Element::compute_quantity_of_interest( const uint            aMeshIndex,
                                                                  enum vis::Output_Type aOutputType,
                                                                  enum vis::Field_Type  aFieldType )
        {
            // FIXME do this only once
            this->compute_my_pdof_values();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->get_requested_IQI( aOutputType )
                ->set_field_interpolator_manager( mSet->get_field_interpolator_manager() );

            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                mSet->get_requested_IQI( aOutputType )
                    ->set_field_interpolator_manager( mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE ) );
            }

            if( aFieldType == vis::Field_Type::NODAL_IP )
            {
                // get the vertices on the treated mesh cluster
                // FIXME also slave?
                const moris::Cell< const moris::mtk::Vertex * > tVertices
                = mFemCluster( aMeshIndex )->get_mesh_cluster()
                                           ->get_vertices_in_cluster();

                moris::Matrix<moris::DDRMat> tVertexLocalCoords
                = mFemCluster( aMeshIndex )->get_mesh_cluster()
                                           ->get_vertices_local_coordinates_wrt_interp_cell();

                // get number of vertices on the treated mesh cluster
                uint tNumNodes = tVertices.size();

                // loop over the vertices on the treated mesh cluster
                for( uint iVertex = 0; iVertex < tNumNodes; iVertex++ )
                {
                    // get the ith vertex coordinates in the IP param space
                    Matrix< DDRMat > tGlobalIntegPoint = tVertexLocalCoords.get_row( iVertex );
                    tGlobalIntegPoint.resize( 1, tGlobalIntegPoint.numel() + 1 );
                    tGlobalIntegPoint( tGlobalIntegPoint.numel() - 1 ) = this->get_time()( 0 );
                    tGlobalIntegPoint = trans( tGlobalIntegPoint );

                    // set vertex coordinates for IP geometry interpolator
                    mSet->get_field_interpolator_manager()
                        ->get_IP_geometry_interpolator()
                        ->set_space_time( tGlobalIntegPoint );

                    // set vertex coordinates for field interpolator
                    mSet->mMasterFIManager->set_space_time( tGlobalIntegPoint );

                    // reset the requested IQI
                    mSet->get_requested_IQI( aOutputType )->reset_eval_flags();

                    // compute quantity of interest at evaluation point
                    Matrix< DDRMat > tQIValue;
                    mSet->get_requested_IQI( aOutputType )->compute_QI( tQIValue );

                    // fill in the nodal set values
                    ( * mSet->mSetNodalValues )( tVertices( iVertex )->get_index(), 0 )
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
            // ask cluster to compute colume
            return mFemCluster( 0 )->compute_volume();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
