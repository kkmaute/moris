/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster.cpp
 *
 */

#include "cl_MTK_Side_Cluster.hpp"
#include "moris_typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell.hpp"

#include "cl_MTK_Cluster.hpp"

#include "fn_norm.hpp"

namespace moris::mtk
{

    // ----------------------------------------------------------------------------------

    Side_Cluster::Side_Cluster() {}

    // ----------------------------------------------------------------------------------

    Vector< mtk::Cell const * > const &
    Side_Cluster::get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return this->get_cells_in_side_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::real tVolume = 0.0;

        Vector< moris::mtk::Cell const * > const *tCells = nullptr;

        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
        {
            tCells = &this->get_primary_cells_in_cluster();
        }
        else
        {
            tCells = &this->get_void_cells_in_cluster();
        }

        for ( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
        {
            tVolume += ( *iC )->compute_cell_measure();
        }

        return tVolume;
    }

    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Side_Cluster::compute_cluster_ig_cell_measures(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        Vector< moris::mtk::Cell const * > const *tCells = nullptr;

        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
        {
            tCells = &this->get_primary_cells_in_cluster();
        }
        else
        {
            tCells = &this->get_void_cells_in_cluster();
        }

        Matrix< DDRMat > tMeasureVec( tCells->size(), 1 );

        for ( uint iC = 0; iC < tCells->size(); iC++ )
        {
            tMeasureVec( iC ) = ( *tCells )( iC )->compute_cell_measure();
        }

        return tMeasureVec;
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_measure_derivative(
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::real tDerivative = 0.0;

        Vector< moris::mtk::Cell const * > const *tCells = nullptr;

        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
        {
            tCells = &this->get_primary_cells_in_cluster();
        }
        else
        {
            tCells = &this->get_void_cells_in_cluster();
        }

        // loop over the ig cells in cluster
        for ( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
        {
            // get the cell coordinates
            Matrix< DDRMat > tCellCoords =
                    ( *iC )->get_vertex_coords();
            ;

            // check if this cell in cluster is affected by the perturbed node
            uint tNumNodesInCell = tCellCoords.n_rows();    // number of nodes in this cell

            // FIXME could be done with node index
            // check if this cell in cluster is affected by the perturbed node
            // flag true if cell is affected by perturbed cluster node
            bool tIsAffected = false;

            // init cell local node index
            uint tLocalVertexID = UINT_MAX;

            // loop over the nodes of the cell
            for ( uint iCellNode = 0; iCellNode < tNumNodesInCell; iCellNode++ )
            {
                // check if perturbed cluster node affects this cell by using the distance between two nodes
                tIsAffected = tIsAffected || ( moris::norm( aPerturbedVertexCoords - tCellCoords.get_row( iCellNode ) ) < 1e-12 );

                // if the cell is affected by perturbed cluster node
                if ( tIsAffected == true )
                {
                    // get cell local node index
                    tLocalVertexID = iCellNode;

                    // break the loop as node correspondence was found
                    break;
                }
            }

            // if the cell is affected by perturbed cluster node
            if ( tIsAffected == true )
            {
                // add contribution from the cell to the cluster side measure
                tDerivative += ( *iC )->compute_cell_measure_deriv(
                        tLocalVertexID,
                        aDirection );
            }
        }

        return tDerivative;
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_side_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        MORIS_ASSERT( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY,
                "Side cluster only operates on primary cells." );

        moris::real tMeasure = 0.0;

        Vector< mtk::Cell const * > const &tCells = this->get_primary_cells_in_cluster();

        moris::Matrix< IndexMat > tSideOrds = this->get_cell_side_ordinals( aIsLeader );

        for ( moris::uint iC = 0; iC < tCells.size(); iC++ )
        {
            tMeasure += tCells( iC )->compute_cell_side_measure( tSideOrds( iC ) );
        }

        return tMeasure;
    }

    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Side_Cluster::compute_cluster_ig_cell_side_measures(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        MORIS_ASSERT( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY,
                "Side cluster only operates on primary cells." );

        Vector< mtk::Cell const * > const &tCells = this->get_primary_cells_in_cluster();

        Matrix< DDRMat > tMeasureVec( tCells.size(), 1 );

        moris::Matrix< IndexMat > tSideOrds = this->get_cell_side_ordinals( aIsLeader );

        for ( moris::uint iC = 0; iC < tCells.size(); iC++ )
        {
            tMeasureVec( iC ) = tCells( iC )->compute_cell_side_measure( tSideOrds( iC ) );
        }

        return tMeasureVec;
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_side_measure_derivative(
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::real tDerivative = 0.0;

        Vector< mtk::Cell const * > const &tCells = this->get_primary_cells_in_cluster();

        moris::Matrix< IndexMat > tSideOrds = this->get_cell_side_ordinals( aIsLeader );

        // loop over the ig cells in cluster
        for ( moris::uint iC = 0; iC < tCells.size(); iC++ )
        {
            // get the cell coordinates
            Matrix< DDRMat > tCellCoords =
                    tCells( iC )->get_cell_physical_coords_on_side_ordinal( tSideOrds( iC ) );

            // check if this cell in cluster is affected by the perturbed node
            uint tNumNodesInCell = tCellCoords.n_rows();    // number of nodes in this cell

            // FIXME could be done with node index
            // check if this cell in cluster is affected by the perturbed node
            // flag true if cell is affected by perturbed cluster node
            bool tIsAffected = false;

            // init cell local node index
            uint tLocalVertexID = UINT_MAX;

            // loop over the nodes of the cell
            for ( uint iCellNode = 0; iCellNode < tNumNodesInCell; iCellNode++ )
            {
                // check if perturbed cluster node affects this cell by using the distance between two nodes
                tIsAffected = tIsAffected || ( moris::norm( aPerturbedVertexCoords - tCellCoords.get_row( iCellNode ) ) < 1e-12 );

                // if the cell is affected by perturbed cluster node
                if ( tIsAffected == true )
                {
                    // get cell local node index
                    tLocalVertexID = iCellNode;

                    // break the loop as node correspondence was found
                    break;
                }
            }

            // if the cell is affected by perturbed cluster node
            if ( tIsAffected == true )
            {
                // add contribution from the cell to the cluster side measure
                tDerivative += tCells( iC )->compute_cell_side_measure_deriv(
                        tSideOrds( iC ),
                        tLocalVertexID,
                        aDirection );
            }
        }

        return tDerivative;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Side_Cluster::get_cell_indices_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        Vector< moris::mtk::Cell const * > const &tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix< moris::IndexMat > tCellIndices( 1, tNumCells );

        // get cell indices and store
        for ( moris::uint i = 0; i < tNumCells; i++ )
        {
            tCellIndices( i ) = tCells( i )->get_index();
        }

        return tCellIndices;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Side_Cluster::get_primary_cell_indices_in_cluster() const
    {
        return this->get_cell_indices_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index
    Side_Cluster::get_interpolation_cell_index() const
    {
        return this->get_interpolation_cell().get_index();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Side_Cluster::get_vertex_indices_in_cluster( mtk::Leader_Follower aLeaderFollower ) const
    {
        // number of cells in cluster
        moris::uint tNumVertices = this->get_num_vertices_in_cluster();

        // cell access
        Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_in_cluster();

        // initialize output
        moris::Matrix< moris::IndexMat > tVertexIndices( 1, tNumVertices );

        // get cell indices and store
        for ( moris::uint i = 0; i < tNumVertices; i++ )
        {
            tVertexIndices( i ) = tVertices( i )->get_index();
        }

        return tVertexIndices;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Side_Cluster::get_cell_ids_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        Vector< moris::mtk::Cell const * > const &tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix< moris::IdMat > tCellIds( 1, tNumCells );

        // get cell indices and store
        for ( moris::uint i = 0; i < tNumCells; i++ )
        {
            tCellIds( i ) = tCells( i )->get_id();
        }

        return tCellIds;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Side_Cluster::get_primary_cell_ids_in_cluster() const
    {
        return this->get_cell_ids_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Side_Cluster::get_vertex_ids_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumVertices = this->get_num_vertices_in_cluster();

        // cell access
        Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_in_cluster();

        // initialize output
        moris::Matrix< moris::IndexMat > tVertexIds( 1, tNumVertices );

        // get cell indices and store
        for ( moris::uint i = 0; i < tNumVertices; i++ )
        {
            tVertexIds( i ) = tVertices( i )->get_id();
        }

        return tVertexIds;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Side_Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index         aClusterLocalIndex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        MORIS_ASSERT( aClusterLocalIndex < (moris_index)this->get_num_sides_in_cluster(),
                "Integration Cell Cluster index out of bounds" );

        // get side ordinal of interest
        moris_index tSideOrdinal = this->get_cell_side_ordinal( aClusterLocalIndex );

        // get the integration cell of interest
        moris::mtk::Cell const *tIntegrationCell =
                this->get_cells_in_side_cluster()( aClusterLocalIndex );

        // get the vertex pointers on the side
        Vector< moris::mtk::Vertex const * > tVerticesOnSide =
                tIntegrationCell->get_vertices_on_side_ordinal( tSideOrdinal );

        // allocate output (nnode x dim_xsi)
        moris::Matrix< moris::DDRMat > tVertexParamCoords(
                tVerticesOnSide.size(),
                this->get_dim_of_param_coord() );

        // iterate through vertices and collect local coordinates
        for ( moris::uint i = 0; i < tVerticesOnSide.size(); i++ )
        {
            tVertexParamCoords.get_row( i ) =
                    this->get_vertex_local_coordinate_wrt_interp_cell( tVerticesOnSide( i ) ).get_row( 0 );
        }

        return tVertexParamCoords;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Side_Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index aPrimaryCellClusterIndex ) const
    {
        return this->get_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellClusterIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_sides_in_cluster() const
    {
        return this->get_cells_in_side_cluster().size();
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_primary_cells() const
    {
        return this->get_num_sides_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_vertices_in_cluster() const
    {
        return this->get_vertices_in_cluster().size();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Side_Cluster::get_vertex_coords_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumVertices = this->get_num_vertices_in_cluster();

        // cell access
        Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_in_cluster();

        // initialize output
        moris::Matrix< moris::DDRMat > tVertexCoords( tNumVertices, this->get_dim_of_param_coord() );

        // get cell indices and store
        for ( moris::uint i = 0; i < tNumVertices; i++ )
        {
            tVertexCoords.get_row( i ) = tVertices( i )->get_coords().get_row( 0 );
        }

        return tVertexCoords;
    }

    // ----------------------------------------------------------------------------------
}    // namespace moris::mtk
