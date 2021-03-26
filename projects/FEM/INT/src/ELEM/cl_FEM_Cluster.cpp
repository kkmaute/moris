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
// FEM/MSI/src
#include "cl_MSI_Equation_Model.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Cluster::Cluster(
                const Element_Type         aElementType,
                const mtk::Cluster        * aMeshCluster,
                Set                       * aSet,
                MSI::Equation_Object      * aEquationObject )
        : mInterpolationElement( aEquationObject ),
          mSet( aSet ),
          mElementType( aElementType )
        {
            // fill the cell cluster pointer
            mMeshCluster = aMeshCluster;

            // fill the master integration cells
            mMasterIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

            // get number of subelements (IG cells)
            uint tNumMasterIGCells = mMasterIntegrationCells.size();

            // element factory
            fem::Element_Factory tElementFactory;

            // set size for the number of IG cells
            mElements.resize( tNumMasterIGCells, nullptr );

            // switch on the element type
            switch ( mElementType )
            {
                case fem::Element_Type::BULK :
                case fem::Element_Type::TIME_SIDESET :
                {
                    // loop over the IG cells
                    for ( uint iIGCell = 0; iIGCell < tNumMasterIGCells; iIGCell++ )
                    {
                        // create an element
                        mElements( iIGCell ) = tElementFactory.create_element(
                                mElementType,
                                mMasterIntegrationCells( iIGCell ),
                                mSet,
                                this,
                                iIGCell );
                    }
                    break;
                }
                case fem::Element_Type::SIDESET :
                {
                    // set the side ordinals for the IG cells in the cluster
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals();

                    // loop over the IG cells
                    for ( uint iIGCell = 0; iIGCell < tNumMasterIGCells; iIGCell++ )
                    {
                        // create an element
                        mElements( iIGCell ) = tElementFactory.create_element(
                                mElementType,
                                mMasterIntegrationCells( iIGCell ),
                                mSet,
                                this,
                                iIGCell );
                    }
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET :
                {
                    // fill the slave integration cells
                    mSlaveIntegrationCells  = aMeshCluster->get_primary_cells_in_cluster( mtk::Master_Slave::SLAVE );

                    // set the side ordinals for the master and slave IG cells
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::MASTER );
                    mSlaveListOfSideOrdinals  = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::SLAVE );

                    // loop over the IG cells
                    for( moris::uint iIGCell = 0; iIGCell < tNumMasterIGCells; iIGCell++)
                    {
                        // create element
                        mElements( iIGCell ) = tElementFactory.create_element(
                                mElementType,
                                mMasterIntegrationCells( iIGCell ),
                                mSlaveIntegrationCells( iIGCell ),
                                mSet,
                                this,
                                iIGCell );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "Cluster::Cluster - No element type specified" );
            }

            // get cluster measure map from set
            mClusterMEAMap = mSet->get_cluster_measure_map();

            // get cluster measure tuples from set
            moris::Cell< std::tuple<
            fem::Measure_Type,
            mtk::Primary_Void,
            mtk::Master_Slave > > tClusterMEATuples = mSet->get_cluster_measure_tuples();

            // build the cluster measures from tuples
            mClusterMEA.resize( tClusterMEATuples.size(), nullptr );
            for( uint iCMEA = 0; iCMEA < tClusterMEATuples.size(); iCMEA++ )
            {
                mClusterMEA( iCMEA ) = std::make_shared< Cluster_Measure >(
                        std::get<0>( tClusterMEATuples( iCMEA ) ),
                        std::get<1>( tClusterMEATuples( iCMEA ) ),
                        std::get<2>( tClusterMEATuples( iCMEA ) ),
                        this );
            }
        }

        //------------------------------------------------------------------------------

        Cluster::Cluster()
        {
            // FIXME could only collect from SP
            // create default cluster measure
            mClusterMEA.resize( 4, nullptr );
            mClusterMEA( 0 ) = std::make_shared< Cluster_Measure >();
            mClusterMEA( 1 ) = std::make_shared< Cluster_Measure >();
            mClusterMEA( 2 ) = std::make_shared< Cluster_Measure >();
            mClusterMEA( 3 ) = std::make_shared< Cluster_Measure >();

            // FIXME could only collect from SP
            // fill the cluster measure access map
            mClusterMEAMap[ std::make_tuple(
                    fem::Measure_Type::CELL_SIDE_MEASURE,
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER ) ] = 0;
            mClusterMEAMap[ std::make_tuple(
                    fem::Measure_Type::CELL_MEASURE,
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER ) ] = 1;
            mClusterMEAMap[ std::make_tuple(
                    fem::Measure_Type::CELL_MEASURE,
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::SLAVE ) ] = 2;
            mClusterMEAMap[ std::make_tuple(
                    fem::Measure_Type::CELL_LENGTH_MEASURE,
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER ) ] = 3;
        }

        //------------------------------------------------------------------------------

        Cluster::~Cluster()
        {
            for( auto tElements : mElements )
            {
                delete tElements;
            }
            mElements.clear();
        }

        //------------------------------------------------------------------------------

        Matrix< IndexMat > & Cluster::get_side_ordinal_info(
                mtk::Master_Slave aIsMaster )
        {
            switch( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                    return mMasterListOfSideOrdinals;

                case mtk::Master_Slave::SLAVE :
                    return mSlaveListOfSideOrdinals;

                default:
                    MORIS_ERROR( false, "Cluster::get_side_ordinal_info - can only be master or slave." );
                    return mMasterListOfSideOrdinals;
            }
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > Cluster::get_vertices_local_coordinates_wrt_interp_cell()
        {
            // check that side cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::BULK,
                    "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a bulk cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_vertices_local_coordinates_wrt_interp_cell - empty cluster.");

                // get the vertices param coords from the cluster
                return mMeshCluster->get_vertices_local_coordinates_wrt_interp_cell();
        }

        //------------------------------------------------------------------------------

        void Cluster::get_vertex_indices_in_cluster_for_sensitivity(
                moris::Matrix< moris::IndexMat > & aVerticesIndices )
        {
            // if mesh cluster is not trivial
            if( !mMeshCluster->is_trivial( mtk::Master_Slave::MASTER ) )
            {
                // get vertices indices on cluster
                aVerticesIndices = mMeshCluster->get_vertex_indices_in_cluster();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::get_vertex_indices_in_cluster_for_visualization(
                moris::Matrix< moris::IndexMat > & aVerticesIndices )
        {
            // check that bulk cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::BULK,
                    "Cluster::get_vertex_indices_in_cluster - not a bulk cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_vertex_indices_in_cluster - empty cluster.");

            // fill vertex indices matrix
            aVerticesIndices = mMeshCluster->get_vertex_indices_in_cluster();
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index aCellIndexInCluster,
                moris::moris_index aSideOrdinal,
                mtk::Master_Slave  aIsMaster )
        {
            // check that side cluster
            MORIS_ASSERT( ( mElementType == fem::Element_Type::DOUBLE_SIDESET ) || ( mElementType == fem::Element_Type::SIDESET ),
                    "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - not a side or double side cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

            // get the side param coords from the side cluster
            return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, aIsMaster );
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index aPrimaryCellIndexInCluster )
        {
            // check that bulk cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::BULK || mElementType == fem::Element_Type::TIME_SIDESET,
                    "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a bulk cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");

            // get the side param coords from the side cluster
            return mMeshCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellIndexInCluster );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Cluster::get_side_normal(
                const mtk::Cell    * aCell,
                moris::moris_index   aSideOrdinal )
        {
            // init normal
            Matrix < DDRMat > tNormal;

            //            // if interpolation cell is linear
            //            if( mSet->get_IG_space_interpolation_order() == mtk::Interpolation_Order::LINEAR )
            //            {
            //                // get normal from the mesh
            //                tNormal = aCell->compute_outward_side_normal( aSideOrdinal );
            //            }
            //            // if integration cell is higher order
            //            else
            //            {
            // get normal from the integration cell geometry interpolator
            mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->get_normal( tNormal );
            //            }

            return tNormal;
        }

        //------------------------------------------------------------------------------

        moris::mtk::Vertex const * Cluster::get_left_vertex_pair(
                moris::mtk::Vertex const * aLeftVertex )
        {
            // check that a double sided cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                    "Cluster::get_left_vertex_pair - not a double side cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_left_vertex_pair - empty cluster.");

            // get the paired vertex on the right
            return mMeshCluster->get_master_vertex_pair( aLeftVertex );
        }

        //------------------------------------------------------------------------------

        moris::moris_index Cluster::get_right_vertex_ordinal_on_facet(
                moris_index                aCellIndexInCluster,
                moris::mtk::Vertex const * aVertex )
        {
            // check that a double sided cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                    "Cluster::get_left_vertex_pair - not a double side cluster.");

            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::get_right_vertex_ordinal_on_facet - empty cluster.");

            // return the index of the paired vertex on the right
            return mMeshCluster->get_slave_vertex_ord_on_facet(aCellIndexInCluster, aVertex);
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_jacobian()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the jacobian for the IG element
                mElements( iElem )->compute_jacobian();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_residual()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the residual for the IG element
                mElements( iElem )->compute_residual();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_jacobian_and_residual()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the jacobian and residual for the IG element
                mElements( iElem )->compute_jacobian_and_residual();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_dRdp()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // reset cluster measures derivatives
            this->reset_cluster_measure_derivatives();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute dRdp for the IG element
                mElements( iElem )->compute_dRdp();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_dQIdp_explicit()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // reset cluster measures derivatives
            this->reset_cluster_measure_derivatives();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute dQIdp for the IG element
                mElements( iElem )->compute_dQIdp_explicit();
            }
        }
        
        //------------------------------------------------------------------------------

        void Cluster::compute_dRdp_and_dQIdp()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // reset cluster measures derivatives
            this->reset_cluster_measure_derivatives();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute dRdp and dQIdp for the IG element
                mElements( iElem )->compute_dRdp_and_dQIdp();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_dQIdu()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the dQIdu for the IG element
                mElements( iElem )->compute_dQIdu();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_QI()
        {
            // reset cluster measures
            this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the quantities of interest for the IG element
                mElements( iElem )->compute_QI();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::compute_quantity_of_interest(
                const uint                         aMeshIndex,
                enum vis::Field_Type               aFieldType )
        {
            // FIXME
            // cannot do it here cause vis mesh
            // reset cluster measures
            //this->reset_cluster_measure();

            // loop over the IG elements
            for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
            {
                // compute the quantity of interest for the IG element
                mElements( iElem )->compute_quantity_of_interest( aMeshIndex, aFieldType );
            }
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< Cluster_Measure > & Cluster::get_cluster_measure(
                fem::Measure_Type aMeasureType,
                mtk::Primary_Void aIsPrimary,
                mtk::Master_Slave aIsMaster )
        {
            // init cluster index
            uint tClusterMEAIndex = UINT_MAX;

            // check if the cluster measure exists in map
            if( mClusterMEAMap.find( std::make_tuple(
                    aMeasureType,
                    aIsPrimary,
                    aIsMaster ) ) != mClusterMEAMap.end() )
            {
                // add the mesh set name map
                tClusterMEAIndex = mClusterMEAMap[ std::make_tuple(
                        aMeasureType,
                        aIsPrimary,
                        aIsMaster ) ];
            }

            MORIS_ERROR( tClusterMEAIndex != UINT_MAX, "Cluster measure not found!" );

            return mClusterMEA( tClusterMEAIndex );
        }

        //------------------------------------------------------------------------------

        void Cluster::reset_cluster_measure()
        {
            // loop over cluster measures
            for( uint iCMEA = 0; iCMEA < mClusterMEA.size(); iCMEA++ )
            {
                // evaluate each cluster measure
                mClusterMEA( iCMEA )->eval_cluster_measure();
            }
        }

        //------------------------------------------------------------------------------

        void Cluster::reset_cluster_measure_derivatives()
        {
            // loop over cluster measures
            for( uint iCMEA = 0; iCMEA < mClusterMEA.size(); iCMEA++ )
            {
                // evaluate each cluster measure
                mClusterMEA( iCMEA )->eval_cluster_measure_derivatives();
            }
        }

        //------------------------------------------------------------------------------

        moris::real Cluster::compute_cluster_cell_measure(
                const mtk::Primary_Void aPrimaryOrVoid ,
                const mtk::Master_Slave aIsMaster ) const
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_measure - empty cluster.");

            return mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsMaster );
        }

        //------------------------------------------------------------------------------

        moris::Matrix< DDRMat > Cluster::compute_cluster_cell_measure_derivative(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster )
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_measure_derivative - empty cluster.");

            // get vertices indices for sensitivity
            moris::Matrix< moris::IndexMat > tVerticesIndices;
            this->get_vertex_indices_in_cluster_for_sensitivity( tVerticesIndices );

            // get the vertex pointers in cluster
            moris::Cell<moris::mtk::Vertex const *> tVertices =
                    mMeshCluster->get_vertices_in_cluster( aIsMaster );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // get local assembly indices
            Matrix< DDSMat > tGeoLocalAssembly;
            if( mSet->get_geo_pdv_assembly_flag() )
            {
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVerticesIndices,
                        tGeoPdvType,
                        tGeoLocalAssembly );
            }

            // get number of pdv on cluster
            uint tNumPDV = mSet->get_geo_pdv_assembly_vector().numel();

            // fill matrix with derivatives
            Matrix< DDRMat > tDerivatives( 1, tNumPDV, 0.0 );

            // create container for the ig cells in cluster
            moris::Cell<moris::mtk::Cell const *> const* tCells = nullptr;

            // get the primary ig cells in cluster
            if(aPrimaryOrVoid == mtk::Primary_Void::PRIMARY)
            {
                tCells = &mMeshCluster->get_primary_cells_in_cluster();
            }
            // get the void ig cells in cluster
            else
            {
                tCells = &mMeshCluster->get_void_cells_in_cluster();
            }

            // loop over the vertices in cluster
            for( uint iClusterNode = 0; iClusterNode < tVertices.size(); iClusterNode++ )
            {
                // get the node coordinates
                Matrix< DDRMat > tPerturbedNodeCoords =
                        tVertices( iClusterNode )->get_coords();

                // loop over the ig cells in cluster
                for( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
                {
                    // get the cell coordinates
                    Matrix< DDRMat > tCellCoords = (*iC)->get_vertex_coords();

                    // get number of nodes in this cell
                    uint tNumNodesInCell = tCellCoords.n_rows(); // number of nodes in this cell

                    // FIXME could be done with node index
                    // check if this cell in cluster is affected by the perturbed node
                    // flag true if cell is affected by perturbed cluster node
                    bool tIsAffected = false;

                    // init cell local node index
                    uint tLocalVertexID = UINT_MAX;

                    // loop over the nodes of the cell
                    for( uint iCellNode = 0; iCellNode < tNumNodesInCell; iCellNode++ )
                    {
                        // check if perturbed cluster node affects this cell by using the distance between two nodes
                        tIsAffected = tIsAffected || ( moris::norm( tPerturbedNodeCoords - tCellCoords.get_row( iCellNode ) ) < 1e-12 );

                        // if the cell is affected by perturbed cluster node
                        if( tIsAffected == true )
                        {
                            // get cell local node index
                            tLocalVertexID = iCellNode;

                            // break the loop as node correspondence was found
                            break;
                        }
                    }

                    // if the cell is affected by perturbed cluster node
                    if( tIsAffected == true )
                    {
                        // loop over the space directions
                        for( uint iSpace = 0; iSpace < tGeoLocalAssembly.n_cols(); iSpace++ )
                        {
                            // get the pdv assembly index
                            sint tPDVAssemblyIndex = tGeoLocalAssembly( iClusterNode, iSpace );

                            // if pdv assembly index is set
                            if( tPDVAssemblyIndex != -1 )
                            {
                                // add contribution from the cell to the cluster measure
                                tDerivatives( tPDVAssemblyIndex ) +=
                                        (*iC)->compute_cell_measure_deriv( tLocalVertexID, iSpace );
                            }
                        }
                    }
                }
            }
            return tDerivatives;
        }

        //------------------------------------------------------------------------------

        moris::real Cluster::compute_cluster_cell_side_measure(
                const mtk::Primary_Void aPrimaryOrVoid ,
                const mtk::Master_Slave aIsMaster ) const
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_side_measure - empty cluster.");

            return mMeshCluster->compute_cluster_cell_side_measure(aPrimaryOrVoid,aIsMaster);
        }

        //------------------------------------------------------------------------------

        moris::Matrix< DDRMat > Cluster::compute_cluster_cell_side_measure_derivative(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster )
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_side_measure_derivative - empty cluster.");

            // get vertices indices for sensitivity
            moris::Matrix< moris::IndexMat > tVerticesIndices;
            this->get_vertex_indices_in_cluster_for_sensitivity( tVerticesIndices );

            // get the vertex pointers in cluster
            moris::Cell<moris::mtk::Vertex const *> tVertices =
                    mMeshCluster->get_vertices_in_cluster( aIsMaster );

            // get the requested geo pdv types
            moris::Cell < enum PDV_Type > tGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

            // get local assembly indices
            Matrix< DDSMat > tGeoLocalAssembly;
            if( mSet->get_geo_pdv_assembly_flag() )
            {
                mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                        tVerticesIndices,
                        tGeoPdvType,
                        tGeoLocalAssembly );
            }

            // get number of pdv on cluster
            uint tNumPDV = mSet->get_geo_pdv_assembly_vector().numel();

            // fill matrix with derivatives
            Matrix< DDRMat > tDerivatives( 1, tNumPDV, 0.0 );

            //get the primary ig cells in cluster
            moris::Cell<mtk::Cell const *> const & tCells =
                    mMeshCluster->get_primary_cells_in_cluster();

            // get the ig cell side ordinals
            moris::Matrix<IndexMat> tSideOrds =
                    mMeshCluster->get_cell_side_ordinals( aIsMaster );

            // loop over the vertices in cluster
            for( uint iClusterNode = 0; iClusterNode < tVertices.size(); iClusterNode++ )
            {
                // get the node coordinates
                Matrix< DDRMat > tPerturbedNodeCoords =
                        tVertices( iClusterNode )->get_coords();

                // loop over the ig cells in cluster
                for( moris::uint iC = 0 ; iC < tCells.size(); iC++ )
                {
                    // get the cell coordinates
                    Matrix< DDRMat > tCellCoords =
                            tCells( iC )->get_cell_physical_coords_on_side_ordinal(tSideOrds(iC));

                    // check if this cell in cluster is affected by the perturbed node
                    uint tNumNodesInCell = tCellCoords.n_rows(); // number of nodes in this cell

                    // FIXME could be done with node index
                    // check if this cell in cluster is affected by the perturbed node
                    // flag true if cell is affected by perturbed cluster node
                    bool tIsAffected = false;

                    // init cell local node index
                    uint tLocalVertexID = UINT_MAX;

                    // loop over the nodes of the cell
                    for( uint iCellNode = 0; iCellNode < tNumNodesInCell; iCellNode++ )
                    {
                        // check if perturbed cluster node affects this cell by using the distance between two nodes
                        tIsAffected = tIsAffected || ( moris::norm( tPerturbedNodeCoords - tCellCoords.get_row( iCellNode ) ) < 1e-12 );

                        // if the cell is affected by perturbed cluster node
                        if( tIsAffected == true )
                        {
                            // get cell local node index
                            tLocalVertexID = iCellNode;

                            // break the loop as node correspondence was found
                            break;
                        }
                    }

                    // if the cell is affected by perturbed cluster node
                    if( tIsAffected == true )
                    {
                        // loop over the space directions
                        for( uint iSpace = 0; iSpace < tGeoLocalAssembly.n_cols(); iSpace++ )
                        {
                            // get the pdv assembly index
                            sint tPDVAssemblyIndex = tGeoLocalAssembly( iClusterNode, iSpace );

                            // if pdv assembly index is set
                            if( tPDVAssemblyIndex != -1 )
                            {
                                // add contribution from the cell to the cluster side measure
                                tDerivatives( tPDVAssemblyIndex ) +=
                                        tCells( iC )->compute_cell_side_measure_deriv(
                                                tSideOrds( iC ),
                                                tLocalVertexID,
                                                iSpace );
                            }
                        }
                    }
                }
            }
            return tDerivatives;
        }

        //------------------------------------------------------------------------------

        moris::real Cluster::compute_cluster_cell_length_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_length_measure - empty cluster.");

            // init volume of the custer
            real tVolume = 0.0;

            // get spatial dimension from IP geometry interpolator
            uint tSpaceDim = mSet->get_field_interpolator_manager( aIsMaster )->
                    get_IP_geometry_interpolator()->
                    get_number_of_space_dimensions();

            // switch on set type
            fem::Element_Type tElementType = mSet->get_element_type();
            switch( tElementType )
            {
                case fem::Element_Type::BULK :
                case fem::Element_Type::TIME_SIDESET :
                case fem::Element_Type::TIME_BOUNDARY :
                {
                    tVolume = mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsMaster );
                    break;
                }
                case fem::Element_Type::SIDESET :
                case fem::Element_Type::DOUBLE_SIDESET :
                {
                    tVolume = mMeshCluster->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsMaster );
                    tSpaceDim = tSpaceDim - 1;
                    break;
                }
                default:
                    MORIS_ERROR( false, "Undefined set type" );
            }

            // compute element size
            switch ( tSpaceDim )
            {
                case 3 :
                {
                    // compute length from volume
                    return std::pow( 6.0 * tVolume / M_PI, 1.0 / 3.0 );
                }
                case 2 :
                {
                    // compute length from volume
                    return std::pow( 4.0 * tVolume / M_PI, 1.0 / 2.0 );
                }
                case 1 :
                {
                    return tVolume;
                }

                default:
                    MORIS_ERROR( false, "Cluster::compute_cluster_cell_length_measure - space dimension can only be 1, 2, or 3. ");
                    return 0.0;
            }
        }

        //------------------------------------------------------------------------------

        moris::Matrix< DDRMat > Cluster::compute_cluster_cell_length_measure_derivative(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster )
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mMeshCluster != NULL,
                    "Cluster::compute_cluster_cell_length_measure_derivative - empty cluster.");

            // get spatial dimension from IP geometry interpolator
            uint tSpaceDim = mSet->get_field_interpolator_manager( aIsMaster )->
                    get_IP_geometry_interpolator()->
                    get_number_of_space_dimensions();

            // init volume of the custer
            real tVolume = 0.0;

            // init derivatives
            Matrix< DDRMat > tDerivatives;

            // switch on set type
            fem::Element_Type tElementType = mSet->get_element_type();
            switch( tElementType )
            {
                case fem::Element_Type::BULK :
                case fem::Element_Type::TIME_SIDESET :
                case fem::Element_Type::TIME_BOUNDARY :
                {
                    tVolume = mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsMaster );
                    tDerivatives = this->compute_cluster_cell_measure_derivative( aPrimaryOrVoid, aIsMaster );
                    break;
                }
                case fem::Element_Type::SIDESET :
                case fem::Element_Type::DOUBLE_SIDESET :
                {
                    tVolume = mMeshCluster->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsMaster );
                    tDerivatives = this->compute_cluster_cell_side_measure_derivative( aPrimaryOrVoid, aIsMaster );
                    tSpaceDim = tSpaceDim - 1;
                    break;
                }
                default:
                    MORIS_ERROR( false, "Undefined set type" );
            }

            // compute element size
            switch ( tSpaceDim )
            {
                case 3 :
                {
                    // compute derivative of length from volume
                    real tdLengthdVolume = 2.0 * std::pow( 6.0 * tVolume / M_PI, -2.0 / 3.0 ) / M_PI;
                    return tdLengthdVolume * tDerivatives;
                }
                case 2 :
                {
                    // compute derivative of length from volume
                    real tdLengthdVolume = 2.0 * std::pow( 4.0 * tVolume / M_PI, -1.0 / 2.0 ) / M_PI;
                    return tdLengthdVolume * tDerivatives;
                }
                case 1 :
                {
                    return tDerivatives;
                }

                default:
                    MORIS_ERROR( false, "Cluster::compute_cluster_cell_length_measure_derivative - space dimension can only be 1, 2, or 3. ");
                    return tDerivatives;
            }
        }

        //------------------------------------------------------------------------------

        moris::real Cluster::compute_ip_cell_length_measure(
                const mtk::Master_Slave aIsMaster ) const
        {
            // check that the mesh cluster was set
            MORIS_ASSERT( mInterpolationElement != NULL,
                    "Cluster::compute_ip_cell_length_measure - empty interpolation element.");

            // compute the volume of the master IP cell
            real tVolume = reinterpret_cast< fem::Interpolation_Element* >( mInterpolationElement )->
                    get_ip_cell( mtk::Master_Slave::MASTER )->compute_cell_measure();

            // get spatial dimension from IP geometry interpolator
            uint tSpaceDim = mSet->get_field_interpolator_manager( aIsMaster )->
                    get_IP_geometry_interpolator()->
                    get_number_of_space_dimensions();

            // compute element size
            switch ( tSpaceDim )
            {
                case 3 :
                {
                    // compute length from volume
                    return std::pow( 6.0 * tVolume / M_PI, 1.0 / 3.0 );
                }
                case 2 :
                {
                    // compute length from volume
                    return std::pow( 4.0 * tVolume / M_PI, 1.0 / 2.0 );
                }
                case 1 :
                {
                    return tVolume;
                }

                default:
                    MORIS_ERROR( false, "Cluster::compute_ip_cell_length_measure - space dimension can only be 1, 2, or 3. ");
                    return 0.0;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
