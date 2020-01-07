/*
 * cl_FEM_Cluster.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Element.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"   //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        Cluster::Cluster( const Element_Type                aElementType,
                          const mtk::Cluster              * aMeshCluster,
                                moris::Cell< Node_Base* > & aNodes,
                                Set                       * aSet )
        : MSI::Equation_Object( aSet ),
          mSet( aSet ),
          mElementType( aElementType )
        {
            // fill the cell cluster pointer
            mMeshCluster = aMeshCluster;

            // fill the master interpolation cell
            mMasterInterpolationCell = & aMeshCluster->get_interpolation_cell();

            // fill the master integration cells
            mMasterIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

            // get the number of IWGs //FIXME
            mNumOfIWGs = mSet->get_number_of_IWGs();

            // switch on the element type
            switch ( mElementType )
            {
                case ( fem::Element_Type::BULK ):
                {
                    // select the element nodes from aNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tVertices = mMasterInterpolationCell->get_vertex_pointers();

                    // get number of nodes from cell
                    uint tNumOfNodes = tVertices.size();

                    // assign node object
                    mNodeObj.resize( 1 );
                    mNodeObj( 0 ).resize( tNumOfNodes, nullptr );

                    // fill node objects
                    for( uint i = 0; i < tNumOfNodes; i++)
                    {
                        mNodeObj( 0 )( i ) = aNodes( tVertices( i )->get_index() );
                    }

                    // set size of Weak BCs
                    mNodalWeakBCs.set_size( tNumOfNodes, 1 );

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik);
                    }
                    break;
                }
                case ( fem::Element_Type::SIDESET ):
                {
                    // set the side ordinals for the IG cells in the cluster
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals();

                    // select the element nodes from aNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tVertices = mMasterInterpolationCell->get_vertex_pointers();

                    // get number of nodes from cell
                    uint tNumOfNodes = tVertices.size();

//                    std::cout<< aNodes.size() <<" num nodes"<<std::endl;

                    // assign node object
                    mNodeObj.resize( 1 );
                    mNodeObj( 0 ).resize( tNumOfNodes, nullptr );

                    // fill node objects
                    for( uint i = 0; i < tNumOfNodes; i++)
                    {
//                        std::cout<< tVertices( i )->get_index() <<" vertex index"<<std::endl;
                        mNodeObj( 0 )( i ) = aNodes( tVertices( i )->get_index() );
                    }

                    // set size of Weak BCs
                    mNodalWeakBCs.set_size( tNumOfNodes, 1 );

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik );
                    }
                    break;
                }
                case( fem::Element_Type::DOUBLE_SIDESET ):
                {
                    // fill the slave interpolation cell
                    mSlaveInterpolationCell  = & aMeshCluster->get_interpolation_cell( mtk::Master_Slave::SLAVE );

                    // fill the slave integration cells
                    mSlaveIntegrationCells  = aMeshCluster->get_primary_cells_in_cluster( mtk::Master_Slave::SLAVE );

                    // set the side ordinals for the master and slave IG cells
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::MASTER );
                    mSlaveListOfSideOrdinals  = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::SLAVE );

                    // select the element nodes from aIPNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tMasterVertices = mMasterInterpolationCell->get_vertex_pointers();
                    moris::Cell< mtk::Vertex* > tSlaveVertices  = mSlaveInterpolationCell->get_vertex_pointers();

                    // get number of nodes from cell
                    uint tNumOfNodesMaster = tMasterVertices.size();
                    uint tNumOfNodesSlave = tSlaveVertices.size();

                    // assign node object
                    mNodeObj.resize( 2 );
                    mNodeObj( 0 ).resize( tNumOfNodesMaster, nullptr );
                    mNodeObj( 1 ).resize( tNumOfNodesSlave , nullptr );

                    // fill node objects
                    for( uint Ik = 0; Ik < tNumOfNodesMaster; Ik++)
                    {
                        mNodeObj( 0 )( Ik ) = aNodes( tMasterVertices( Ik )->get_index() );
                    }
                    for( uint Ik = 0; Ik < tNumOfNodesSlave; Ik++)
                    {
                        mNodeObj( 1 )( Ik ) = aNodes( tSlaveVertices( Ik )->get_index() );
                    }

                    // set size of Weak BCs
                    mNodalWeakBCs.set_size( tNumOfNodesMaster, 1 );             // FIXME  replace this

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSlaveIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "Cluster::Cluster - No element type specified" );
                    break;
            }
        }

//------------------------------------------------------------------------------
        Cluster::~Cluster()
        {
            for( auto tElements : mElements )
            {
                delete tElements;
            }
            mElements.clear();

            for( auto tElements : mVisElements )
            {
                delete tElements;
            }
            mVisElements.clear();
        }

//------------------------------------------------------------------------------
        moris::Matrix< moris::DDRMat > Cluster::get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                               moris::moris_index aSideOrdinal,
                                                                                               mtk::Master_Slave  aIsMaster,
                                                                                               bool               aIsVis )
        {
            // check that side cluster
            MORIS_ASSERT( ( mElementType == fem::Element_Type::DOUBLE_SIDESET ) || ( mElementType == fem::Element_Type::SIDESET ),
                          "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - not a side or double side cluster.");

            // is trivial master or slave?
            bool tIsTrivial = mSet->mIsTrivialMaster;
            if ( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                tIsTrivial = mSet->mIsTrivialSlave;
            }

            // if trivial cluster IP cell = IG cell
            if( tIsTrivial )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_IP_geometry_interpolator( aIsMaster )->extract_space_side_space_param_coeff( aSideOrdinal,
                                                                                                              mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                if ( !aIsVis )
                {
                    // check that the mesh cluster was set
                    MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

                    // get the side param coords from the side cluster
                    return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, aIsMaster );
                }
                else
                {
                    // check that the mesh cluster was set
                    MORIS_ASSERT( mVisMeshCluster != NULL, "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

                    // get the side param coords from the side cluster
                    return mVisMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, aIsMaster );
                }
            }
        }

//------------------------------------------------------------------------------
        moris::Matrix< moris::DDRMat > Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellIndexInCluster,
                                                                                                       bool aIsVis )
        {
            // check that bulk cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::BULK, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a bulk cluster.");

            // if trivial cluster IP cell = IG cell
            if( mSet->mIsTrivialMaster )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_IP_geometry_interpolator()->extract_space_param_coeff( mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                if ( !aIsVis )
                {
                    // check that the mesh cluster was set
                    MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

                    // get the side param coords from the side cluster
                    return mMeshCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellIndexInCluster );
                }
                else
                {
                    // check that the mesh cluster was set
                    MORIS_ASSERT( mVisMeshCluster != NULL, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

                    // get the side param coords from the side cluster
                    return mVisMeshCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellIndexInCluster );
                }
            }
        }

//------------------------------------------------------------------------------
        Matrix< DDRMat > Cluster::get_side_normal( const mtk::Cell   * aCell,
                                                  moris::moris_index   aSideOrdinal )
        {
            // init normal
            Matrix < DDRMat > tNormal;

            // if interpolation cell is linear
//            if( mSet->get_IG_space_interpolation_order() == mtk::Interpolation_Order::LINEAR )
//            {
//                // get normal from the mesh
//                tNormal = aCell->compute_outward_side_normal( aSideOrdinal );
//            }
//            // if integration cell is higher order
//            else
//            {
                // get normal from the integration cell geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->get_normal( tNormal );
//            }

            return tNormal;
        }

//------------------------------------------------------------------------------
        moris::mtk::Vertex const * Cluster::get_left_vertex_pair( moris::mtk::Vertex const * aLeftVertex )
        {
            // check that a double sided cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                          "Cluster::get_left_vertex_pair - not a double side cluster.");

            // get the paired vertex on the right
            return mMeshCluster->get_left_vertex_pair( aLeftVertex );
        }

//------------------------------------------------------------------------------
        moris::moris_index Cluster::get_right_vertex_ordinal_on_facet( moris_index aCellIndexInCluster,
                                                                     moris::mtk::Vertex const * aVertex )
        {
            // check that a double sided cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                          "Cluster::get_left_vertex_pair - not a double side cluster.");

            // return the index of the paired vertex on the right
            return mMeshCluster->get_right_vertex_ord_on_facet(aCellIndexInCluster, aVertex);
        }

//------------------------------------------------------------------------------
        void Cluster::set_field_interpolators_coefficients()
         {
             // get number of master dof types
             uint tMasterNumDofTypes = mSet->mMasterDofTypes.size();

             // loop on the master dof types
             for( uint iDOF = 0; iDOF < tMasterNumDofTypes; iDOF++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_dof_type_list()( iDOF );

                 // get the pdof values for the ith dof type group
                 Cell< Matrix< DDRMat > > tCoeff_Original;
                 Matrix< DDRMat > tCoeff;

                 this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original );

                 // reshape tCoeffs into the order the cluster expects them
                 this->reshape_pdof_values( tCoeff_Original, tCoeff );

                 // set field interpolator coefficients
                 mSet->mMasterFIManager->get_field_interpolators_for_type( tDofTypeGroup( 0 ) )
                                       ->set_coeff( tCoeff );
             }

             // get number of slave dof types
             uint tSlaveNumDofTypes = mSet->mSlaveDofTypes.size();

             // loop on the slave dof types
             for( uint iDOF = 0; iDOF < tSlaveNumDofTypes; iDOF++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF );

                 // get the pdof values for the ith dof type group
                 Cell< Matrix< DDRMat > > tCoeff_Original;
                 Matrix< DDRMat > tCoeff;

                 this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

                 // reshape tCoeffs into the order the cluster expects them
                 this->reshape_pdof_values( tCoeff_Original, tCoeff );

                 // set the field coefficients
                 mSet->mSlaveFIManager->get_field_interpolators_for_type( tDofTypeGroup( 0 ) )
                                      ->set_coeff( tCoeff );
             }
         }

//------------------------------------------------------------------------------
        void Cluster::compute_jacobian()
         {

            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
             }

             //Fixme do this only once
             this->compute_my_pdof_values();

             // init the jacobian //fixme still ok?
             mSet->initialize_mJacobian();

             // set the field interpolators coefficients
             this->set_field_interpolators_coefficients();

             // FIXME should not be like this
             mSet->set_IWG_field_interpolator_managers();
             mSet->set_IWG_geometry_interpolators();

             // loop over the elements
             for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
             {
                 // compute the jacobian for the element
                 mElements( iElem )->compute_jacobian();
             }
//             print( mSet->mJacobian,"Jacobian");
         }

//------------------------------------------------------------------------------
        void Cluster::compute_residual()
        {
            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
            }

            //Fixme do this only once
            this->compute_my_pdof_values();

            // init the residual
            mSet->initialize_mResidual();

            // set the field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // FIXME should not be like this
            mSet->set_IWG_field_interpolator_managers();
            mSet->set_IWG_geometry_interpolators();

            // loop over the elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the residual for the element
                mElements( iElem )->compute_residual();
            }
//            print( mSet->mResidual,"Residual");
        }

//------------------------------------------------------------------------------
        void Cluster::compute_jacobian_and_residual()
         {
            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
             }

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
             mSet->set_IWG_geometry_interpolators();

             // loop over the elements
             for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
             {
                 // compute the jacobian for the element
                 mElements( iElem )->compute_jacobian();

                 // compute the residual for the element
                 mElements( iElem )->compute_residual();
             }
         }

//------------------------------------------------------------------------------
        //void Cluster::compute_quantity_of_interest( fem::QI_Compute_Type aQIComputeType )
        void Cluster::compute_quantity_of_interest( enum vis::Output_Type aOutputType,
                                                    enum vis::Field_Type  aFieldType )
        {
            // set the IP geometry interpolator physical space and time coefficients for the master interpolation cell
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_space_coeff( mMasterInterpolationCell->get_vertex_coords() );
            mSet->get_IP_geometry_interpolator( mtk::Master_Slave::MASTER )->set_time_coeff( this->mTime );

            // if double side cluster
             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_space_coeff( mSlaveInterpolationCell->get_vertex_coords() );
                 mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE )->set_time_coeff( this->mTime );
             }

             // FIXME do this only once
             this->compute_my_pdof_values();

             // set the field interpolators coefficients
             this->set_field_interpolators_coefficients();

             // FIXME should not be like this
             mSet->get_requested_IQI( aOutputType )
                 ->set_field_interpolator_manager( mSet->get_field_interpolator_manager() );
             mSet->get_requested_IQI( aOutputType )
                 ->set_geometry_interpolator( mSet->get_IP_geometry_interpolator() );

             if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
             {
                 // set the IP geometry interpolator physical space and time coefficients for the slave interpolation cell
                 mSet->get_requested_IQI( aOutputType )
                     ->set_field_interpolator_manager( mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE ) );
                 mSet->get_requested_IQI( aOutputType )
                     ->set_geometry_interpolator( mSet->get_IP_geometry_interpolator( mtk::Master_Slave::SLAVE) );
             }

             // FIXME choose between visualization and integration
             // flag to choose between mVisElements and mElements

             // loop over the elements
             for ( uint iElem = 0; iElem < mVisElements.size(); iElem++ )
             {
                 // compute the quantity of interest for the element
                 mVisElements( iElem )->compute_quantity_of_interest( aOutputType, aFieldType );
             }
         }

//------------------------------------------------------------------------------
        real Cluster::compute_volume()
        {
            // set cluster volume
            real tClusterVolume = 0;

            // loop over the elements in cluster
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the volume for the element
                tClusterVolume += mElements( iElem )->compute_volume();
            }
            // return cluster volume value
            return tClusterVolume;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
