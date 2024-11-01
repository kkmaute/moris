/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Cluster.cpp
 *
 */

#include <iostream>

#include "cl_FEM_Cluster_Measure.hpp"
#include "cl_FEM_Element.hpp"                       //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src

#include "cl_MSI_Equation_Model.hpp"

#include "cl_MTK_PointPairs.hpp"
#include "fn_norm.hpp"
#include "fn_sort.hpp"
#include "fn_sum.hpp"
#include <algorithm>
#include <cl_FEM_Element_Nonconformal_Sideset.hpp>
#include <cl_MTK_Nonconformal_Side_Cluster.hpp>
#include <memory>

namespace moris::fem
{

    //------------------------------------------------------------------------------

    Cluster::Cluster(
            const Element_Type    aElementType,
            const mtk::Cluster   *aMeshCluster,
            Set                  *aSet,
            MSI::Equation_Object *aEquationObject,
            bool                  aIsVisCluster )
            : mInterpolationElement( aEquationObject )
            , mSet( aSet )
            , mElementType( aElementType )
            , mIsVisCluster( aIsVisCluster )
    {
        // fill the cell cluster pointer
        mMeshCluster = aMeshCluster;

        // fill the leader integration cells
        mLeaderIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

        // get number of subelements (IG cells)
        uint tNumLeaderIGCells = mLeaderIntegrationCells.size();

        // element factory
        fem::Element_Factory tElementFactory;

        // set size for the number of IG cells
        mElements.resize( tNumLeaderIGCells, nullptr );

        // switch on the element type
        switch ( mElementType )
        {
            case fem::Element_Type::BULK:
            case fem::Element_Type::TIME_SIDESET:
            {
                // loop over the IG cells
                for ( uint iIGCell = 0; iIGCell < tNumLeaderIGCells; iIGCell++ )
                {
                    // create an element
                    mElements( iIGCell ) = tElementFactory.create_single_sided_element(
                            mElementType,
                            mLeaderIntegrationCells( iIGCell ),
                            mSet,
                            this,
                            iIGCell );
                }
                break;
            }
            case fem::Element_Type::SIDESET:
            {
                // set the side ordinals for the IG cells in the cluster
                mLeaderListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals();

                // loop over the IG cells
                for ( uint iIGCell = 0; iIGCell < tNumLeaderIGCells; iIGCell++ )
                {
                    // create an element
                    mElements( iIGCell ) = tElementFactory.create_single_sided_element(
                            mElementType,
                            mLeaderIntegrationCells( iIGCell ),
                            mSet,
                            this,
                            iIGCell );
                }
                break;
            }
            case fem::Element_Type::DOUBLE_SIDESET:
            {
                // fill the follower integration cells
                mFollowerIntegrationCells = aMeshCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::FOLLOWER );

                // set the side ordinals for the leader and follower IG cells
                mLeaderListOfSideOrdinals   = aMeshCluster->get_cell_side_ordinals( mtk::Leader_Follower::LEADER );
                mFollowerListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals( mtk::Leader_Follower::FOLLOWER );

                // loop over the IG cells
                for ( moris::uint iIGCell = 0; iIGCell < tNumLeaderIGCells; iIGCell++ )
                {
                    // create element
                    mElements( iIGCell ) = tElementFactory.create_double_sided_element(
                            mElementType,
                            mLeaderIntegrationCells( iIGCell ),
                            mFollowerIntegrationCells( iIGCell ),
                            mSet,
                            this,
                            iIGCell );
                }
                break;
            }
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                const auto *tNonconformalSideCluster = dynamic_cast< mtk::Nonconformal_Side_Cluster const * >( aMeshCluster );
                auto const  tIntegrationPointPairs   = tNonconformalSideCluster->get_integration_point_pairs();
                auto const  tNodalPointPairs         = tNonconformalSideCluster->get_nodal_point_pairs();

                // the same cells might be used multiple times (i.e. they are paired with multiple leader cells)
                mtk::Nonconformal_Side_Cluster::NonconformalCellPairing const &tNonconformalPairing = tNonconformalSideCluster->get_nonconformal_cell_pairing();

                // the nonconformal side set provides different cell-pairs (i.e. one follower-cell might be paired with multiple leader-cells)
                mElements.resize( tNonconformalPairing.size(), nullptr );

                // fill the follower integration cells
                mLeaderIntegrationCells   = tNonconformalSideCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                mFollowerIntegrationCells = tNonconformalSideCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::FOLLOWER );

                // set the side ordinals for the leader and follower IG cells
                mLeaderListOfSideOrdinals   = tNonconformalSideCluster->get_cell_side_ordinals( mtk::Leader_Follower::LEADER );
                mFollowerListOfSideOrdinals = tNonconformalSideCluster->get_cell_side_ordinals( mtk::Leader_Follower::FOLLOWER );

                // loop over the nonconformal pairings
                size_t iElementIndex = 0;
                for ( auto const &tNonconformalPair : tNonconformalPairing )
                {
                    auto [ tLeaderCellLocalIndex, tFollowerCellLocalIndex ] = tNonconformalPair.first;
                    auto [ tIntegrationPairIndex, tNodalPairIndex ]         = tNonconformalPair.second;

                    mElements( iElementIndex ) =
                            new Element_Nonconformal_Sideset(
                                    mLeaderIntegrationCells( tLeaderCellLocalIndex ),
                                    mFollowerIntegrationCells( tFollowerCellLocalIndex ),
                                    aSet,
                                    this,
                                    tLeaderCellLocalIndex,
                                    tFollowerCellLocalIndex,
                                    tIntegrationPairIndex == -1 ? mtk::IntegrationPointPairs() : tIntegrationPointPairs( tIntegrationPairIndex ),
                                    tNodalPairIndex == -1 ? mtk::NodalPointPairs() : tNodalPointPairs( tNodalPairIndex ) );
                    iElementIndex++;
                }
                break;
            }
            default:
                MORIS_ERROR( false, "Cluster::Cluster - No element type specified" );
        }

        // determine IG cells to be ignored in residual computation
        // FIXME: should be used to ignore cell at time of construction
        this->determine_elements_for_residual_and_iqi_computation();

        // get cluster measure map from set
        std::set< Cluster_Measure::ClusterMeasureSpecification > const tClusterSpecifications =
                mSet->get_cluster_measure_specifications();

        // build the cluster measures from tuples
        for ( auto const &tClusterSpec : tClusterSpecifications )
        {
            mClusterMeasures[ tClusterSpec ] = std::make_shared< Cluster_Measure >( tClusterSpec, this );
        }
    }

    //------------------------------------------------------------------------------

    Cluster::Cluster()
    {
        Vector< Cluster_Measure::ClusterMeasureSpecification > tClusterSpecifications;

        // create default cluster measure
        tClusterSpecifications.push_back( std::make_tuple(
                fem::Measure_Type::CELL_SIDE_MEASURE,
                mtk::Primary_Void::PRIMARY,
                mtk::Leader_Follower::LEADER ) );
        tClusterSpecifications.push_back( std::make_tuple(
                fem::Measure_Type::CELL_MEASURE,
                mtk::Primary_Void::PRIMARY,
                mtk::Leader_Follower::LEADER ) );
        tClusterSpecifications.push_back( std::make_tuple(
                fem::Measure_Type::CELL_MEASURE,
                mtk::Primary_Void::PRIMARY,
                mtk::Leader_Follower::FOLLOWER ) );
        tClusterSpecifications.push_back( std::make_tuple(
                fem::Measure_Type::CELL_LENGTH_MEASURE,
                mtk::Primary_Void::PRIMARY,
                mtk::Leader_Follower::LEADER ) );

        // build the cluster measures from tuples
        for ( auto const &tClusterSpec : tClusterSpecifications )
        {
            mClusterMeasures[ tClusterSpec ] = std::make_shared< Cluster_Measure >();
        }
    }

    //------------------------------------------------------------------------------

    Cluster::~Cluster()
    {
        for ( auto tElements : mElements )
        {
            delete tElements;
        }
        mElements.clear();
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat > &
    Cluster::get_side_ordinal_info(
            mtk::Leader_Follower aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
                return mLeaderListOfSideOrdinals;

            case mtk::Leader_Follower::FOLLOWER:
                return mFollowerListOfSideOrdinals;

            default:
                MORIS_ERROR( false, "Cluster::get_side_ordinal_info - can only be leader or follower." );
                return mLeaderListOfSideOrdinals;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< moris::DDRMat >
    Cluster::get_vertices_local_coordinates_wrt_interp_cell( mtk::Leader_Follower aLeaderFollower )
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "FEM::Cluster::get_vertices_local_coordinates_wrt_interp_cell() - Empty cluster." );

        // get the vertices param coords from the cluster
        return mMeshCluster->get_vertices_local_coordinates_wrt_interp_cell( aLeaderFollower );
    }

    //------------------------------------------------------------------------------

    void
    Cluster::get_vertex_indices_in_cluster_for_sensitivity(
            Matrix< IndexMat > &aVerticesIndices )
    {
        // if mesh cluster is not trivial
        if ( !mMeshCluster->is_trivial( mtk::Leader_Follower::LEADER ) )
        {
            // for bulk and single-sided elements, get the leader vertices, for double-sided both leader and follower
            mtk::Leader_Follower tLeaderFollower = mtk::Leader_Follower::LEADER;

            if ( this->get_element_type() == Element_Type::DOUBLE_SIDESET || this->get_element_type() == Element_Type::NONCONFORMAL_SIDESET )
            {
                tLeaderFollower = mtk::Leader_Follower::UNDEFINED;
            }

            // get vertices indices on cluster
            aVerticesIndices = mMeshCluster->get_vertex_indices_in_cluster( tLeaderFollower );
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::get_vertex_indices_in_cluster_for_visualization(
            Matrix< IndexMat >  &aVerticesIndices,
            mtk::Leader_Follower aLeaderFollower )
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr, "FEM::Cluster::get_vertex_indices_in_cluster() - Cluster is empty." );

        // fill vertex indices matrix
        aVerticesIndices = mMeshCluster->get_vertex_indices_in_cluster( aLeaderFollower );
    }

    //------------------------------------------------------------------------------

    Matrix< moris::DDRMat >
    Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index   aCellIndexInCluster,
            moris::moris_index   aSideOrdinal,
            mtk::Leader_Follower aIsLeader )
    {
        // check that side cluster
        MORIS_ASSERT(
                ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                        || ( mElementType == fem::Element_Type::NONCONFORMAL_SIDESET )
                        || ( mElementType == fem::Element_Type::SIDESET ),
                "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - not a side or double side cluster." );

        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - empty cluster." );

        // get the side param coords from the side cluster
        return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, aIsLeader );
    }

    //------------------------------------------------------------------------------

    Matrix< moris::DDRMat >
    Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index aPrimaryCellIndexInCluster )
    {
        // check that bulk cluster
        MORIS_ASSERT( mElementType == fem::Element_Type::BULK || mElementType == fem::Element_Type::TIME_SIDESET,
                "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a bulk cluster." );

        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - empty cluster." );

        // get the side param coords from the side cluster
        return mMeshCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellIndexInCluster );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::get_side_normal(
            const mtk::Cell   *aCell,
            moris::moris_index aSideOrdinal )
    {
        // init normal
        Matrix< DDRMat > tNormal;

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

    moris::mtk::Vertex const *
    Cluster::get_left_vertex_pair(
            moris::mtk::Vertex const *aLeftVertex )
    {
        // check that a double sided cluster
        MORIS_ASSERT(
                ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                        || ( mElementType == fem::Element_Type::NONCONFORMAL_SIDESET ),
                "Cluster::get_left_vertex_pair - not a double side cluster." );

        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::get_left_vertex_pair - empty cluster." );

        // get the paired vertex on the right
        return mMeshCluster->get_leader_vertex_pair( aLeftVertex );
    }

    //------------------------------------------------------------------------------

    moris::moris_index
    Cluster::get_right_vertex_ordinal_on_facet(
            moris_index               aCellIndexInCluster,
            moris::mtk::Vertex const *aVertex )
    {
        // check that a double sided cluster
        MORIS_ASSERT( ( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                              || ( mElementType == fem::Element_Type::NONCONFORMAL_SIDESET ),
                "Cluster::get_left_vertex_pair - not a double side cluster." );

        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::get_right_vertex_ordinal_on_facet - empty cluster." );

        // return the index of the paired vertex on the right
        return mMeshCluster->get_follower_vertex_ord_on_facet( aCellIndexInCluster, aVertex );
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_jacobian()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute jacobian
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the jacobian for the IG element
                mElements( iElem )->compute_jacobian();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_residual()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute residual
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the residual for the IG element
                mElements( iElem )->compute_residual();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_jacobian_and_residual()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute residual and jacobian for this element
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the jacobian and residual for the IG element
                mElements( iElem )->compute_jacobian_and_residual();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_dRdp()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // reset cluster measures derivatives
        this->reset_cluster_measure_derivatives();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute derivative of residual
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute dRdp for the IG element
                mElements( iElem )->compute_dRdp();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_dQIdp_explicit()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // reset cluster measures derivatives
        this->reset_cluster_measure_derivatives();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute derivative of QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute dQIdp for the IG element
                mElements( iElem )->compute_dQIdp_explicit();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_dRdp_and_dQIdp()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // reset cluster measures derivatives
        this->reset_cluster_measure_derivatives();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute derivative of residual and QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute dRdp and dQIdp for the IG element
                mElements( iElem )->compute_dRdp_and_dQIdp();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_dQIdu()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute derivative of QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the dQIdu for the IG element
                mElements( iElem )->compute_dQIdu();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_QI()
    {
        // reset cluster measures
        this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the quantities of interest for the IG element
                mElements( iElem )->compute_QI();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_quantity_of_interest(
            const uint           aFemMeshIndex,
            enum vis::Field_Type aFieldType )
    {
        // FIXME
        // cannot do it here cause vis mesh
        // reset cluster measures
        // this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the quantity of interest for the IG element
                mElements( iElem )->compute_quantity_of_interest( aFemMeshIndex, aFieldType );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::compute_quantity_of_interest(
            Matrix< DDRMat >      &aValues,
            mtk::Field_Entity_Type aFieldType,
            uint                   aIQIIndex,
            real                  &aSpaceTimeVolume )
    {
        // FIXME
        // cannot do it here cause vis mesh
        // reset cluster measures
        // this->reset_cluster_measure();

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // check whether to compute QI
            if ( mComputeResidualAndIQI( iElem ) )
            {
                // compute the quantity of interest for the IG element
                mElements( iElem )->compute_quantity_of_interest( aValues, aFieldType, aIQIIndex, aSpaceTimeVolume );
            }
        }
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< Cluster_Measure >
    Cluster::get_cluster_measure(
            fem::Measure_Type    aMeasureType,
            mtk::Primary_Void    aIsPrimary,
            mtk::Leader_Follower aIsLeader )
    {
        auto tClusterMeasure = mClusterMeasures.find( std::make_tuple( aMeasureType, aIsPrimary, aIsLeader ) );

        MORIS_ASSERT( tClusterMeasure != mClusterMeasures.end(),
                "Cluster::get_cluster_measure - Cluster measure not found!" );

        return tClusterMeasure->second;
    }

    //------------------------------------------------------------------------------

    void
    Cluster::reset_cluster_measure()
    {
        for ( auto &[ tClusterSpecification, tClusterMeasure ] : mClusterMeasures )
        {
            tClusterMeasure->eval_cluster_measure();
        }
    }

    //------------------------------------------------------------------------------

    void
    Cluster::reset_cluster_measure_derivatives()
    {
        for ( auto &[ tClusterSpecification, tClusterMeasure ] : mClusterMeasures )
        {
            tClusterMeasure->eval_cluster_measure_derivatives();
        }
    }

    //------------------------------------------------------------------------------

    moris::real
    Cluster::compute_cluster_cell_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_measure - empty cluster." );

        return mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_cluster_cell_measure_derivative(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader )
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_measure_derivative - empty cluster." );

        // check that this function is not called on the VIS clusters, only makes sense on actual FEM clusters
        MORIS_ASSERT( !this->is_VIS_cluster(), "FEM cluster is actually VIS cluster." );

        // get adv ids on cluster
        const Vector< sint > &tAdvIds = mSet->get_ig_adv_ids();

        // get number of pdv on cluster
        uint tNumADV = tAdvIds.size();

        // fill matrix with derivatives
        Matrix< DDRMat > tDerivatives( 1, tNumADV, 0.0 );

        // if pdv defined on cluster
        if ( tNumADV == 0 )
        {
            return tDerivatives;
        }

        // get the vertex pointers in cluster
        Vector< moris::mtk::Vertex const * > tVertices =
                mMeshCluster->get_vertices_in_cluster( aIsLeader );

        // get the vertex indices in cluster
        Matrix< IndexMat > tVertexIndices = mMeshCluster->get_vertex_indices_in_cluster( aIsLeader );

        // build design extraction operator for all nodes in cluster
        mSet->create_geo_adv_assembly_data( tVertexIndices, aIsLeader );

        // get design extraction operator
        const Vector< Matrix< DDRMat > > &tGeoWeights = mSet->get_adv_geo_weights( aIsLeader );

        //        // get the requested geo pdv types
        //        Vector< enum gen::PDV_Type > tGeoPdvType;
        //        mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

        //        // get number of pdv on cluster
        //        uint tNumPDV = mSet->get_geo_pdv_assembly_vector().numel();
        //
        //        // fill matrix with derivatives
        //        Matrix< DDRMat > tDerivatives( 1, tNumPDV, 0.0 );

        //        // if pdv defined on cluster
        //        if ( tNumPDV > 0 )
        //        {
        //            // loop over the vertices in cluster

        for ( uint iClusterNode = 0; iClusterNode < tVertices.size(); iClusterNode++ )
        {
            // check that adv influence this node
            if ( tGeoWeights( iClusterNode ).n_rows() == 0 )
            {
                continue;
            }

            //
            //            // get local assembly indices
            //            Matrix< DDSMat >     tGeoLocalAssembly;
            //            MSI::Equation_Model *tEquationModel    = mSet->get_equation_model();
            //            Matrix< IndexMat >   tClusterNodeIndex = { { tVertices( iClusterNode )->get_index() } };
            //            tEquationModel->get_integration_xyz_pdv_assembly_indices(
            //                    tClusterNodeIndex,
            //                    tGeoPdvType,
            //                    tGeoLocalAssembly );

            //            // if cluster node is not associated with pdv
            //            if ( sum( tGeoLocalAssembly ) == -2 || sum( tGeoLocalAssembly ) == -3 )
            //            {
            //                continue;
            //            }

            // get the node coordinates
            const Matrix< DDRMat > &tPerturbedNodeCoords = tVertices( iClusterNode )->get_coords();

            // loop over the space directions
            for ( uint iSpace = 0; iSpace < tPerturbedNodeCoords.n_cols(); iSpace++ )
            {
                //                // get the pdv assembly index
                //                moris_index tPDVAssemblyIndex = tGeoLocalAssembly( iSpace );
                //
                //                // sanity check
                //                MORIS_ASSERT( tPDVAssemblyIndex != MORIS_INDEX_MAX,
                //                        "FEM::Cluster::compute_cluster_cell_measure_derivative() - "
                //                        "tPDVAssemblyIndex is MORIS_INDEX_MAX which likely means that mXYZLocalAssemblyIndices"
                //                        "in the FEM-Model has not been set for this node and PDV type." );

                Matrix< DDRMat > dIGNodeCorddAdv = tGeoWeights( iClusterNode ).get_row( iSpace );

                //                // if pdv assembly index is set
                //                if ( tPDVAssemblyIndex != -1 )
                if ( norm( dIGNodeCorddAdv ) > MORIS_REAL_EPS )
                {
                    //                    // check
                    //                    MORIS_ASSERT( tPDVAssemblyIndex < (moris_index)tDerivatives.numel(),
                    //                            "FEM::Cluster::compute_cluster_cell_measure_derivative() - "
                    //                            "tPDVAssemblyIndex is out of bounds. Assembly index is %i, but number of PDVs only is: %zu",
                    //                            tPDVAssemblyIndex,
                    //                            tDerivatives.numel() );

                    // add contribution from the cell to the cluster measure
                    real tDeriv = mMeshCluster->compute_cluster_cell_measure_derivative(
                            tPerturbedNodeCoords,
                            iSpace,
                            aPrimaryOrVoid,
                            aIsLeader );

                    tDerivatives += tDeriv * dIGNodeCorddAdv;
                }
            }
            //        }
        }
        return tDerivatives;
    }

    //------------------------------------------------------------------------------

    moris::real
    Cluster::compute_cluster_cell_side_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_side_measure - empty cluster." );

        return mMeshCluster->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsLeader );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_cluster_cell_side_measure_derivative(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader )
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_side_measure_derivative - empty cluster." );

        // get adv ids on cluster
        const Vector< sint > &tAdvIds = mSet->get_ig_adv_ids();

        // get number of pdv on cluster
        uint tNumADV = tAdvIds.size();

        // fill matrix with derivatives
        Matrix< DDRMat > tDerivatives( 1, tNumADV, 0.0 );

        // if pdv defined on cluster
        if ( tNumADV == 0 )
        {
            return tDerivatives;
        }

        // get the vertex pointers in cluster
        Vector< moris::mtk::Vertex const * > tVertices =
                mMeshCluster->get_vertices_in_cluster( aIsLeader );

        // get the vertex indices in cluster
        Matrix< IndexMat > tVertexIndices = mMeshCluster->get_vertex_indices_in_cluster( aIsLeader );

        //        // get the requested geo pdv types
        //        Vector< enum gen::PDV_Type >
        //                tGeoPdvType;
        //        mSet->get_ig_unique_dv_types_for_set( tGeoPdvType );

        // build design extraction operator for all nodes in cluster
        mSet->create_geo_adv_assembly_data( tVertexIndices );

        // get design extraction operator
        const Vector< Matrix< DDRMat > > &tGeoWeights = mSet->get_adv_geo_weights();

        // loop over the vertices in cluster
        for ( uint iClusterNode = 0; iClusterNode < tVertices.size(); iClusterNode++ )
        {
            // check that adv influence this node
            if ( tGeoWeights( iClusterNode ).n_rows() == 0 )
            {
                continue;
            }

            //            // get local assembly indices
            //            Matrix< DDSMat >     tGeoLocalAssembly;
            //            MSI::Equation_Model *tEquationModel    = mSet->get_equation_model();
            //            Matrix< IndexMat >   tClusterNodeIndex = { { tVertices( iClusterNode )->get_index() } };
            //            tEquationModel->get_integration_xyz_pdv_assembly_indices(
            //                    tClusterNodeIndex,
            //                    tGeoPdvType,
            //                    tGeoLocalAssembly );
            //
            //            // if cluster node is not associated with pdv
            //            if ( sum( tGeoLocalAssembly ) == -2 || sum( tGeoLocalAssembly ) == -3 )
            //            {
            //                continue;
            //            }

            // get the node coordinates
            const Matrix< DDRMat > &tPerturbedNodeCoords = tVertices( iClusterNode )->get_coords();

            // loop over the space directions
            for ( uint iSpace = 0; iSpace < tPerturbedNodeCoords.n_cols(); iSpace++ )
            {
                //                // get the pdv assembly index
                //                sint tPDVAssemblyIndex = tGeoLocalAssembly( iSpace );
                //
                //                // sanity check
                //                MORIS_ASSERT( tPDVAssemblyIndex != MORIS_INDEX_MAX,
                //                        "FEM::Cluster::compute_cluster_cell_side_measure_derivative() - "
                //                        "tPDVAssemblyIndex is MORIS_INDEX_MAX which likely means that mXYZLocalAssemblyIndices"
                //                        "in the FEM-Model has not been set for this node and PDV type." );

                Matrix< DDRMat > dIGNodeCorddAdv = tGeoWeights( iClusterNode ).get_row( iSpace );
                // print( dIGNodeCorddAdv, "dIGNodeCorddAdv" );

                //                // if pdv assembly index is set
                //                if ( tPDVAssemblyIndex != -1 )
                if ( norm( dIGNodeCorddAdv ) > MORIS_REAL_EPS )
                {
                    // add contribution from the cell to the cluster side measure
                    real tDeriv = mMeshCluster->compute_cluster_cell_side_measure_derivative(
                            tPerturbedNodeCoords,
                            iSpace,
                            aPrimaryOrVoid,
                            aIsLeader );

                    tDerivatives += tDeriv * dIGNodeCorddAdv;
                }
            }
        }

        // print( tDerivatives, "tDerivatives" );
        return tDerivatives;
    }

    //------------------------------------------------------------------------------

    moris::real
    Cluster::compute_cluster_cell_length_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_length_measure - empty cluster." );

        // initialize volume of the cluster
        real tVolume = 0.0;

        // get spatial dimension from IP geometry interpolator
        uint tSpaceDim = mSet->get_field_interpolator_manager( aIsLeader )->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        // switch on set type
        fem::Element_Type tElementType = mSet->get_element_type();
        switch ( tElementType )
        {
            case fem::Element_Type::BULK:
            case fem::Element_Type::TIME_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
            {
                tVolume = mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
                break;
            }
            case fem::Element_Type::SIDESET:
            case fem::Element_Type::DOUBLE_SIDESET:
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                tVolume   = mMeshCluster->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsLeader );
                tSpaceDim = tSpaceDim - 1;
                break;
            }
            default:
                MORIS_ERROR( false,
                        "Cluster::compute_cluster_cell_length_measure - Undefined element type" );
        }

        // compute element size
        switch ( tSpaceDim )
        {
            case 3:
            {
                // compute length from volume
                return std::pow( 6.0 * tVolume / M_PI, 1.0 / 3.0 );
            }
            case 2:
            {
                // compute length from volume
                return std::pow( 4.0 * tVolume / M_PI, 1.0 / 2.0 );
            }
            case 1:
            {
                return tVolume;
            }

            default:
                MORIS_ERROR( false,
                        "Cluster::compute_cluster_cell_length_measure - space dimension can only be 1, 2, or 3. " );
                return 0.0;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_cluster_cell_length_measure_derivative(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader )
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mMeshCluster != nullptr,
                "Cluster::compute_cluster_cell_length_measure_derivative - empty cluster." );

        // get spatial dimension from IP geometry interpolator
        uint tSpaceDim = mSet->get_field_interpolator_manager( aIsLeader )->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        // initialize volume of the cluster
        real tVolume = 0.0;

        // init derivatives
        Matrix< DDRMat > tDerivatives;

        // switch on set type
        fem::Element_Type tElementType = mSet->get_element_type();
        switch ( tElementType )
        {
            case fem::Element_Type::BULK:
            case fem::Element_Type::TIME_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
            {
                tVolume      = mMeshCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
                tDerivatives = this->compute_cluster_cell_measure_derivative( aPrimaryOrVoid, aIsLeader );
                break;
            }
            case fem::Element_Type::SIDESET:
            case fem::Element_Type::DOUBLE_SIDESET:
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                tVolume      = mMeshCluster->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsLeader );
                tDerivatives = this->compute_cluster_cell_side_measure_derivative( aPrimaryOrVoid, aIsLeader );
                tSpaceDim    = tSpaceDim - 1;
                break;
            }
            default:
                MORIS_ERROR( false,
                        "Cluster::compute_cluster_cell_length_measure_derivative - Undefined element type" );
        }

        // compute element size
        switch ( tSpaceDim )
        {
            case 3:
            {
                // compute derivative of length from volume
                real tdLengthdVolume = 2.0 * std::pow( 6.0 * tVolume / M_PI, -2.0 / 3.0 ) / M_PI;
                return tdLengthdVolume * tDerivatives;
            }
            case 2:
            {
                // compute derivative of length from volume
                real tdLengthdVolume = 2.0 * std::pow( 4.0 * tVolume / M_PI, -1.0 / 2.0 ) / M_PI;
                return tdLengthdVolume * tDerivatives;
            }
            case 1:
            {
                return tDerivatives;
            }

            default:
                MORIS_ERROR( false,
                        "Cluster::compute_cluster_cell_length_measure_derivative - space dimension can only be 1, 2, or 3. " );
                return tDerivatives;
        }
    }

    //------------------------------------------------------------------------------

    moris::real
    Cluster::compute_ip_cell_length_measure(
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the mesh cluster was set
        MORIS_ASSERT( mInterpolationElement != nullptr,
                "Cluster::compute_ip_cell_length_measure - empty interpolation element." );

        // compute the volume of the leader IP cell
        real tVolume = reinterpret_cast< fem::Interpolation_Element * >( mInterpolationElement )->get_ip_cell( mtk::Leader_Follower::LEADER )->compute_cell_measure();

        // get spatial dimension from IP geometry interpolator
        uint tSpaceDim = mSet->get_field_interpolator_manager( aIsLeader )->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        // compute element size
        switch ( tSpaceDim )
        {
            case 3:
            {
                // compute length from volume
                return std::pow( 6.0 * tVolume / M_PI, 1.0 / 3.0 );
            }
            case 2:
            {
                // compute length from volume
                return std::pow( 4.0 * tVolume / M_PI, 1.0 / 2.0 );
            }
            case 1:
            {
                return tVolume;
            }

            default:
                MORIS_ERROR( false,
                        "Cluster::compute_ip_cell_length_measure - space dimension can only be 1, 2, or 3. " );
                return 0.0;
        }
    }

    //------------------------------------------------------------------------------

    real Cluster::compute_volume()
    {
        // new way to compute cluster volume

        const mtk::Primary_Void    tPrimaryOrVoid = mtk::Primary_Void::PRIMARY;
        const mtk::Leader_Follower tIsLeader      = mtk::Leader_Follower::LEADER;

        real tClusterVolume = 0.0;

        // switch on set type
        fem::Element_Type tElementType = mSet->get_element_type();

        switch ( tElementType )
        {
            case fem::Element_Type::BULK:
            case fem::Element_Type::TIME_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
            {
                tClusterVolume = mMeshCluster->compute_cluster_cell_measure( tPrimaryOrVoid, tIsLeader );
                break;
            }
            case fem::Element_Type::SIDESET:
            case fem::Element_Type::DOUBLE_SIDESET:
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                tClusterVolume = mMeshCluster->compute_cluster_cell_side_measure( tPrimaryOrVoid, tIsLeader );
                break;
            }
            default:
                MORIS_ERROR( false, "Cluster::compute_volume - Undefined element type" );
        }

        // check for zero or negative volume
        MORIS_ASSERT( tClusterVolume > MORIS_REAL_MIN,
                "Cluster::compute_volume - volume of cluster is smaller / equal zero: %15.9e\n",
                tClusterVolume );

        return tClusterVolume;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_element_volumes()
    {
        // new way to compute cluster volume

        const mtk::Primary_Void    tPrimaryOrVoid = mtk::Primary_Void::PRIMARY;
        const mtk::Leader_Follower tIsLeader      = mtk::Leader_Follower::LEADER;

        // switch on set type
        fem::Element_Type tElementType = mSet->get_element_type();

        switch ( tElementType )
        {
            case fem::Element_Type::BULK:
            case fem::Element_Type::TIME_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
            {
                return mMeshCluster->compute_cluster_ig_cell_measures( tPrimaryOrVoid, tIsLeader );
                break;
            }
            case fem::Element_Type::SIDESET:
            case fem::Element_Type::DOUBLE_SIDESET:
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                return mMeshCluster->compute_cluster_ig_cell_side_measures( tPrimaryOrVoid, tIsLeader );
                break;
            }
            default:
                MORIS_ERROR( false, "Cluster::compute_element_volumes - Undefined element type" );
        }

        return { { 0 } };
    }

    //------------------------------------------------------------------------------

    real Cluster::compute_volume_in_fem()
    {
        // compute cluster volume by numerical integration in FEM

        // initialize cluster volume
        real tClusterVolume = 0;

        // loop over the IG elements
        for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
        {
            // add volume contribution for the IG element
            tClusterVolume += mElements( iElem )->compute_volume();
        }

        // check for zero or negative volume
        MORIS_ASSERT( tClusterVolume > MORIS_REAL_MIN,
                "Cluster::compute_volume - volume of cluster is smaller / equal zero: %15.9e\n",
                tClusterVolume );

        // check for consistency between old and new way to compute cluster volume
        MORIS_ASSERT( std::abs( tClusterVolume - this->compute_volume() ) < std::max( 1e-8 * tClusterVolume, 10.0 * MORIS_REAL_EPS ),
                "Cluster::compute_volume - inconsistent volume computation: %15.9e vs %15.9e\n",
                tClusterVolume,
                this->compute_volume() );

        // return cluster volume value
        return tClusterVolume;
    }

    //------------------------------------------------------------------------------

    std::map< mtk::Cell const *, Vector< moris_index > > Cluster::get_cell_to_element_map( mtk::Leader_Follower aLeaderFollowerType ) const
    {
        std::map< mtk::Cell const *, Vector< moris_index > > tUniqueElements;

        for ( size_t iElement = 0; iElement < mElements.size(); ++iElement )
        {
            fem::Element    *tElement = mElements( iElement );
            mtk::Cell const *tCell    = tElement->get_mtk_cell( aLeaderFollowerType );
            tUniqueElements[ tCell ].push_back( iElement );
        }
        return tUniqueElements;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_relative_volume()
    {
        // get (unique) cells in cluster
        Vector< mtk::Cell const * > const tUniqueCells = mMeshCluster->get_primary_cells_in_cluster();

        uint const tNumberOfCells = tUniqueCells.size();

        // initialize cluster volume
        real tClusterVolume = 0.0;

        // compute volumes/areas of IG cells
        Matrix< DDRMat > tRelativeVolume = this->compute_element_volumes();

        // check for correct number of IG cells
        MORIS_ERROR( tRelativeVolume.numel() == tNumberOfCells,
                "Cluster::compute_relative_volume - inconsistent number of IG cells.\n" );

        // loop over the IG cells and drop zero elements
        for ( uint iCell = 0; iCell < tNumberOfCells; iCell++ )
        {
            // if it smaller than zero; set it to zero
            tRelativeVolume( iCell ) = std::max( tRelativeVolume( iCell ), 0.0 );

            // add volume contribution for the IG element
            tClusterVolume += tRelativeVolume( iCell );
        }

        // check for consistent cluster volume computation
        MORIS_ASSERT( std::abs( tClusterVolume - this->compute_volume() ) < std::max( 1e-8 * tClusterVolume, 10.0 * MORIS_REAL_EPS ),
                "Cluster::compute_relative_volume - inconsistent volume computation: %15.9e vs %15.9e\n",
                tClusterVolume,
                this->compute_volume() );

        // compute relative volumes of each element in cluster; if total volume is close to zero
        // fill vector with negative one
        if ( tClusterVolume > MORIS_REAL_MIN )
        {
            tRelativeVolume = 1.0 / tClusterVolume * tRelativeVolume;

            // check that summation of volumes has been done correctly
            MORIS_ASSERT( std::abs( 1.0 - sum( tRelativeVolume ) ) < 10. * MORIS_REAL_EPS,
                    "Cluster::compute_relative_volume - relative volumes do not sum up to one: %15.9e.\n",
                    sum( tRelativeVolume ) );
        }
        else
        {
            tRelativeVolume.fill( -1.0 );
        }

        // return relative volumes
        return tRelativeVolume;
    }

    //------------------------------------------------------------------------------

    void Cluster::determine_elements_for_residual_and_iqi_computation()
    {
        // get (unique) cells in cluster
        // we need to use cells instead of elements because elements can be repeated (e.g. in a nonconformal setting)
        Vector< mtk::Cell const * > const tCells = mMeshCluster->get_primary_cells_in_cluster();

        uint const tNumberOfCells = tCells.size();

        // initialize flags for computing residuals and IQIs (default: on)
        // this is done per element, not per cell
        mComputeResidualAndIQI.resize( mElements.size(), true );

        // skip remainder if there is only one element
        if ( tNumberOfCells == 1 )
        {
            return;
        }

        // determine whether IG element should be used for computation of residual and and jacobian
        Matrix< DDRMat > tRelativeCellVolume = this->compute_relative_volume();

        // check for degenerated element, i.e. all components of tRelativeElementVolume are negative
        if ( tRelativeCellVolume( 0 ) < 0.0 )    // if first element is negative, all elements are negative
        {
            for ( size_t iElement = 0; iElement < mElements.size(); ++iElement )
            {
                mComputeResidualAndIQI( iElement ) = false;
            }
            return;
        }

        // get drop tolerance
        real const tCellDropTolerance = this->compute_volume_drop_threshold( tRelativeCellVolume, mVolumeError );

        // check that drop tolerance is positive
        MORIS_ASSERT( tCellDropTolerance > -MORIS_REAL_MIN,
                "Cluster::determine_elements_for_residual_and_iqi_computation - %s",
                "drop tolerance is negative.\n" );

        std::map< mtk::Cell const *, Vector< int > > tCellToElementMap = this->get_cell_to_element_map();

        // loop over cells and deactivate all element that use this cell if necessary
        for ( size_t iCell = 0; iCell < tCells.size(); ++iCell )
        {
            if ( tRelativeCellVolume( iCell ) < tCellDropTolerance )
            {
                // if cell volume is smaller than threshold
                // set flag for computing residual and IQI to false
                // for all elements that use this cell
                for ( auto const &iElement : tCellToElementMap[ tCells( iCell ) ] )
                {
                    mComputeResidualAndIQI( iElement ) = false;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    real Cluster::compute_volume_drop_threshold(
            const Matrix< DDRMat > &tRelativeCellVolume,
            const real             &tVolumeError )
    {
        // create copy of vector of relative cell volumes
        Matrix< DDRMat > tSortedVolumes = tRelativeCellVolume;

        // sort volumes
        sort( tRelativeCellVolume, tSortedVolumes, "descend", 0 );

        // initialize sum of relative volumes
        real tSum = 0;

        // find threshold
        for ( uint iCell = 0; iCell < tRelativeCellVolume.numel(); ++iCell )
        {
            // check if error criterion is satisfied and if so
            // return current relative cell volume as threshold
            if ( 1.0 - tSum < tVolumeError )
            {
                return tSortedVolumes( iCell );
            }

            // add to sum of relative element volumes
            tSum += tSortedVolumes( iCell );
        }

        // if all elements are needed return 0
        return 0.0;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
