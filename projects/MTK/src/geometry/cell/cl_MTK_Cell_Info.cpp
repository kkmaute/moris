/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info.cpp
 *
 */

#include "cl_MTK_Cell_Info.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Interpolation_Function.hpp"
#include "fn_trans.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include <cmath>

#include <memory>

namespace moris::mtk
{
    // ---------------------------------------------------------------------------------

    bool
    Cell_Info::normalize_shape_vector( Matrix< DDRMat > &aVector )
    {
        // hypot avoids squaring dimensional components, which can underflow or
        // overflow before normalization. No absolute length cutoff is imposed.
        real tLength = 0.0;

        for ( uint i = 0; i < aVector.numel(); ++i )
        {
            tLength = std::hypot( tLength, aVector( i ) );
        }

        if ( !std::isfinite( tLength ) || tLength == 0.0 )
        {
            return false;
        }

        for ( uint i = 0; i < aVector.numel(); ++i )
        {
            aVector( i ) /= tLength;
        }

        return true;
    }

    // ---------------------------------------------------------------------------------

    bool
    Cell_Info::shape_vectors_parallel( const Matrix< DDRMat > &aFirst,
            const Matrix< DDRMat >                            &aSecond )
    {
        // Inputs are unit vectors; the tolerance is dimensionless.
        constexpr real tTolerance = 1.0e-8;

        if ( aFirst.numel() == 2 )
        {
            return std::abs( aFirst( 0 ) * aSecond( 1 ) - aFirst( 1 ) * aSecond( 0 ) ) <= tTolerance;
        }

        Matrix< DDRMat > tCross = cross( aFirst, aSecond );

        return std::hypot( tCross( 0 ), tCross( 1 ), tCross( 2 ) ) <= tTolerance;
    }

    // ---------------------------------------------------------------------------------

    CellShape
    Cell_Info::compute_quad_cell_shape( const Cell *aCell ) const
    {
        constexpr real             tTolerance = 1.0e-8;
        const Vector< Vertex * >   tVertices  = aCell->get_vertex_pointers();
        Vector< Matrix< DDRMat > > tEdges( 4 );

        for ( uint iEdge = 0; iEdge < 4; ++iEdge )
        {
            // Edge maps put endpoints first, followed by higher-order nodes.
            const Matrix< IndexMat > tNodes  = this->get_node_to_edge_map( iEdge );
            const Matrix< DDRMat >   tOrigin = tVertices( tNodes( 0 ) )->get_coords();
            tEdges( iEdge )                  = tVertices( tNodes( 1 ) )->get_coords() - tOrigin;

            if ( !normalize_shape_vector( tEdges( iEdge ) ) )
            {
                return CellShape::GENERAL;
            }

            for ( uint iNode = 2; iNode < tNodes.numel(); ++iNode )
            {
                Matrix< DDRMat > tDirection = tVertices( tNodes( iNode ) )->get_coords() - tOrigin;

                if ( !normalize_shape_vector( tDirection ) || !shape_vectors_parallel( tEdges( iEdge ), tDirection ) )
                {
                    return CellShape::GENERAL;
                }
            }
        }

        if ( !shape_vectors_parallel( tEdges( 0 ), tEdges( 2 ) ) || !shape_vectors_parallel( tEdges( 1 ), tEdges( 3 ) ) )
        {
            return CellShape::STRAIGHT;
        }

        // Preserve QUAD4's existing orientation convention.
        if ( std::abs( dot( tEdges( 0 ), tEdges( 1 ) ) ) > tTolerance
                || std::abs( tEdges( 0 )( 1 ) ) > tTolerance
                || ( this->get_cell_topology() != CellTopology::QUAD4 && tEdges( 0 )( 0 ) < 0.0 ) )
        {
            return CellShape::PARALLEL;
        }

        return CellShape::RECTANGULAR;
    }

    // ---------------------------------------------------------------------------------

    CellShape
    Cell_Info::compute_hex_cell_shape( const Cell *aCell ) const
    {
        constexpr real             tTolerance = 1.0e-8;
        const Vector< Vertex * >   tVertices  = aCell->get_vertex_pointers();
        Vector< Matrix< DDRMat > > tNormals( 6 );

        for ( uint iFace = 0; iFace < 6; ++iFace )
        {
            // Face maps put the four corners first. Normalize the two edge
            // directions before crossing them, then normalize the face normal.
            const Matrix< IndexMat > tNodes  = this->get_node_to_face_map( iFace );
            const Matrix< DDRMat >   tOrigin = tVertices( tNodes( 0 ) )->get_coords();
            Matrix< DDRMat >         tEdge0  = tVertices( tNodes( 1 ) )->get_coords() - tOrigin;
            Matrix< DDRMat >         tEdge1  = tVertices( tNodes( 2 ) )->get_coords() - tVertices( tNodes( 1 ) )->get_coords();

            if ( !normalize_shape_vector( tEdge0 ) || !normalize_shape_vector( tEdge1 ) )
            {
                return CellShape::GENERAL;
            }

            tNormals( iFace ) = cross( tEdge0, tEdge1 );

            if ( !normalize_shape_vector( tNormals( iFace ) ) )
            {
                return CellShape::GENERAL;
            }

            // Every remaining face node must lie in the corner plane. Unit
            // displacement vectors keep this test independent of cell size.
            for ( uint iNode = 3; iNode < tNodes.numel(); ++iNode )
            {
                Matrix< DDRMat > tDirection = tVertices( tNodes( iNode ) )->get_coords() - tOrigin;

                if ( !normalize_shape_vector( tDirection ) || std::abs( dot( tNormals( iFace ), tDirection ) ) > tTolerance )
                {
                    return CellShape::GENERAL;
                }
            }
        }

        if ( !shape_vectors_parallel( tNormals( 0 ), tNormals( 2 ) )
                || !shape_vectors_parallel( tNormals( 1 ), tNormals( 3 ) )
                || !shape_vectors_parallel( tNormals( 4 ), tNormals( 5 ) ) )
        {
            return CellShape::STRAIGHT;
        }

        Matrix< DDRMat > tEdge0 = tVertices( 1 )->get_coords() - tVertices( 0 )->get_coords();
        Matrix< DDRMat > tEdge1 = tVertices( 2 )->get_coords() - tVertices( 1 )->get_coords();

        if ( !normalize_shape_vector( tEdge0 ) || !normalize_shape_vector( tEdge1 ) )
        {
            return CellShape::GENERAL;
        }

        if ( std::abs( dot( tNormals( 0 ), tNormals( 1 ) ) ) > tTolerance
                || std::abs( dot( tNormals( 0 ), tNormals( 4 ) ) ) > tTolerance
                || std::abs( dot( tNormals( 1 ), tNormals( 4 ) ) ) > tTolerance
                || std::abs( tEdge0( 1 ) ) > tTolerance || std::abs( tEdge0( 2 ) ) > tTolerance
                || std::abs( tEdge1( 0 ) ) > tTolerance || std::abs( tEdge1( 2 ) ) > tTolerance
                || tEdge0( 0 ) < 0.0 || tEdge1( 1 ) < 0.0 )
        {
            return CellShape::PARALLEL;
        }

        return CellShape::RECTANGULAR;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info::get_geometric_node_to_facet_map() const
    {
        MORIS_ERROR( 0, "ERROR" );
        return moris::Matrix< moris::IndexMat >();
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info::get_geometric_node_to_facet_map( moris::uint aSideOrdinal ) const
    {
        MORIS_ERROR( 0, "ERROR" );
        return moris::Matrix< moris::IndexMat >();
    }

    // ---------------------------------------------------------------------------------

    moris::uint
    Cell_Info::get_adjacent_side_ordinal( moris::uint aSideOrdinal ) const
    {
        MORIS_ERROR( 0, "get_adjacent_side_ordinal only makes sense for hex/quad type cells" );
        return MORIS_UINT_MAX;
    }

    // ---------------------------------------------------------------------------------
    Matrix< DDRMat >
    Cell_Info::get_vertex_loc_coord( moris_index const &aVertexOrdinal ) const
    {
        MORIS_ERROR( 0, "get_vertex_loc_coord not implemented for given cell info type" );
        return Matrix< DDRMat >( 0, 0 );
    }

    Vector< moris_index >
    Cell_Info::get_vertex_path_to_entity_rank_and_ordinal(
            moris_index aVertexOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
    {
        MORIS_ERROR( 0, "get_vertex_path_to_entity_rank_and_ordinal not implemented." );
        return Vector< moris_index >( 0, 0 );
    }

    Vector< moris_index >
    Cell_Info::get_edge_path_to_entity_rank_and_ordinal(
            moris_index aVertexOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
    {
        MORIS_ERROR( 0, "get_edge_path_to_entity_rank_and_ordinal not implemented." );
        return Vector< moris_index >( 0, 0 );
    }

    bool
    Cell_Info::is_entity_connected_to_facet(
            moris_index aEntityOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
    {
        MORIS_ERROR( 0, "get_edge_path_to_entity_rank_and_ordinal not implemented." );
        return false;
    }

    // ---------------------------------------------------------------------------------
    moris_index
    Cell_Info::get_shared_vertex_ordinal_between_edges(
            moris_index aEdgeOrdinal1,
            moris_index aEdgeOrdinal2 ) const
    {
        MORIS_ERROR( 0, "No default implementation." );
        return MORIS_INDEX_MAX;
    }

    // ---------------------------------------------------------------------------------

    void
    Cell_Info::get_loc_coord_on_side_ordinal(
            moris::uint const &aSideOrdinal,
            Matrix< DDRMat >  &aXi ) const
    {
        // number of vertices on facet side
        uint tNumVertsPerFacet = this->get_num_verts_per_facet();

        // dimension of local coordinates
        uint tDimXi = this->get_loc_coord_dim();

        // allocate space
        aXi.resize( tNumVertsPerFacet, tDimXi );

        // get the vertex ordinals on facet
        Matrix< IndexMat > tVerticesOnSide = this->get_node_to_facet_map( aSideOrdinal );

        MORIS_ASSERT( tNumVertsPerFacet == tVerticesOnSide.numel(),
                "Dimension mismatch between get_num_verts_per_facet() and get_node_to_facet_map()" );

        // iterate through vertices and get their local coordinate
        for ( moris::uint i = 0; i < tNumVertsPerFacet; i++ )
        {
            aXi.get_row( i ) = this->get_vertex_loc_coord( tVerticesOnSide( i ) ).get_row( 0 );
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Cell_Info::get_loc_coords_of_cell( Matrix< DDRMat > &aXi ) const
    {
        // number of vertices
        uint tNumVerts = this->get_num_verts();

        // dimension of local coordinates
        uint tDimXi = this->get_loc_coord_dim();

        // allocate space
        aXi.resize( tNumVerts, tDimXi );

        // iterate through vertices and get their local coordinate
        for ( moris::uint i = 0; i < tNumVerts; i++ )
        {
            aXi.get_row( i ) = this->get_vertex_loc_coord( (moris_index)i ).get_row( 0 );
        }
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info::get_edge_to_face_map() const
    {
        MORIS_ERROR( 0, "get edge to face map not implemented for this CELL_INFO, this CELL_INFO currently used in XTK only for tets/tris" );
        return moris::Matrix< moris::IndexMat >( 0, 0 );
    }

    // ---------------------------------------------------------------------------------

    void
    Cell_Info::eval_N( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                  &aNXi ) const
    {
        MORIS_ERROR( 0, "eval_N not implemented for this type of cell info" );
    }

    // ---------------------------------------------------------------------------------

    moris::real
    Cell_Info::compute_cell_size( moris::mtk::Cell const *aCell ) const
    {
        // Linear simplices have a valid geometric shortcut for every admissible shape.
        const CellTopology tTopology = this->get_cell_topology();

        if ( tTopology == CellTopology::TRI3 || tTopology == CellTopology::TET4 )
        {
            return this->compute_cell_size_special( aCell );
        }

        // RECTANGULAR includes the axis alignment required by the box formulas.
        if ( aCell->get_cell_shape() == CellShape::RECTANGULAR )
        {
            return this->compute_cell_size_special( aCell );
        }

        return this->compute_cell_size_general( aCell );
    }

    // ---------------------------------------------------------------------------------

    moris::real
    Cell_Info::compute_cell_size_general( moris::mtk::Cell const *aCell ) const
    {
        // pulling geometry
        Geometry_Type       tGeoType = aCell->get_geometry_type();
        Interpolation_Order tIPOrder = aCell->get_interpolation_order();

        // creating interpolation rule
        Interpolation_Rule tInterpolationRule( tGeoType,
                Interpolation_Type::LAGRANGE,
                tIPOrder,
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR );

        // create a space interpolator
        Space_Interpolator tSpaceInterpolator( tInterpolationRule );

        // get the coefficients xHat and XiHat
        Matrix< DDRMat > tXHat = aCell->get_vertex_coords();

        // Reuse the interpolator's own reference interpolation function.
        tSpaceInterpolator.set_space_coeff( tXHat );
        tSpaceInterpolator.set_param_coeff();

        // create a integration rule
        Integration_Order tIGOrder = aCell->get_integration_order();

        Integration_Rule tIntegrationRule( tGeoType,
                Integration_Type::GAUSS,
                tIGOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create space rule
        std::unique_ptr< Integration_Coeffs_Base > tSpaceCoeffs = tIntegrationRule.create_space_coeffs();

        // get number of integration points, integration points and weights
        Matrix< DDRMat > tIntegPoints;
        tSpaceCoeffs->get_points( tIntegPoints );
        Matrix< DDRMat > tIntegWeights;
        tSpaceCoeffs->get_weights( tIntegWeights );
        uint tNumIntegPoints = tIntegWeights.numel();

        // init the surface of the integration mesh
        real tCellSize = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tSpaceInterpolator.set_space( tIntegPointI );

            // add to cell size
            tCellSize += tSpaceInterpolator.space_det_J() * tIntegWeights( iGP );
        }

        return tCellSize;
    }

    moris::real
    Cell_Info::compute_cell_size_straight( moris::mtk::Cell const *aCell ) const
    {
        MORIS_ERROR( false, "Cell_Info::compute_cell_size_straight not implemented for this cell" );
        return 0.0;
    }

    moris::real
    Cell_Info::compute_cell_size_deriv(
            moris::mtk::Cell const *aCell,
            uint                    aLocalVertexID,
            uint                    aDirection ) const
    {
        MORIS_ERROR( false,
                "Cell_Info::compute_cell_size_deriv not implemented for this cell. Could you use "
                "Cell_Info::compute_cell_size_deriv_general for any shape" );
        return 0.0;
    }

    moris::real
    Cell_Info::compute_cell_size_deriv_general(
            moris::mtk::Cell const *aCell,
            uint                    aLocalVertexID,
            uint                    aDirection ) const
    {
        // eventually this will be the non-virtual portion of "compute_cell_size_deriv"
        // pulling geometry
        Geometry_Type       tGeoType = aCell->get_geometry_type();
        Interpolation_Order tIPOrder = aCell->get_interpolation_order();

        // creating interpolation rule
        Interpolation_Rule tInterpolationRule(
                tGeoType,
                Interpolation_Type::LAGRANGE,
                tIPOrder,
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR );

        // create a space interpolator
        Space_Interpolator tSpaceInterpolator( tInterpolationRule );

        // get the coefficients xHat and XiHat
        Matrix< DDRMat > tXHat = aCell->get_vertex_coords();

        // Reuse the interpolator's own reference interpolation function.
        tSpaceInterpolator.set_space_coeff( tXHat );
        tSpaceInterpolator.set_param_coeff();

        // create a integration rule
        Integration_Order tIGOrder = aCell->get_integration_order();

        Integration_Rule tIntegrationRule(
                tGeoType,
                Integration_Type::GAUSS,
                tIGOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create space rule
        std::unique_ptr< Integration_Coeffs_Base > tSpaceCoeffs = tIntegrationRule.create_space_coeffs();

        // get number of integration points, integration points and weights
        Matrix< DDRMat > tIntegPoints;
        tSpaceCoeffs->get_points( tIntegPoints );
        Matrix< DDRMat > tIntegWeights;
        tSpaceCoeffs->get_weights( tIntegWeights );
        uint tNumIntegPoints = tIntegWeights.numel();

        // init the surface of the integration mesh
        real tCellSizeDeriv = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tSpaceInterpolator.set_space( tIntegPointI );

            // add to area
            tCellSizeDeriv +=
                    tSpaceInterpolator.space_det_J_deriv( aLocalVertexID, aDirection ) * tIntegWeights( iGP );
        }

        return tCellSizeDeriv;
    }

    // ---------------------------------------------------------------------------------

    void
    Cell_Info::compute_cell_size_derivatives_general(
            const Cell             *aCell,
            Matrix< DDRMat >       &aDerivatives,
            const Matrix< DDSMat > &aRequested ) const
    {
        const Matrix< DDRMat > tCoordinates = aCell->get_vertex_coords();
        MORIS_ERROR( aRequested.numel() == 0 || ( aRequested.n_rows() == tCoordinates.n_rows() && aRequested.n_cols() == tCoordinates.n_cols() ),
                "Cell_Info::compute_cell_size_derivatives_general - Invalid request mask." );

        aDerivatives.set_size( tCoordinates.n_rows(), tCoordinates.n_cols(), 0.0 );

        if ( aRequested.numel() > 0 )
        {
            bool tAnyRequested = false;

            for ( uint i = 0; i < aRequested.numel(); ++i )
            {
                tAnyRequested = tAnyRequested || aRequested.data()[ i ] != 0;
            }

            if ( !tAnyRequested )
            {
                return;
            }
        }

        // Reuse one interpolator and quadrature rule for all requested entries.
        Interpolation_Rule tRule( aCell->get_geometry_type(), Interpolation_Type::LAGRANGE, aCell->get_interpolation_order(), Interpolation_Type::LAGRANGE, Interpolation_Order::LINEAR );
        Space_Interpolator tInterpolator( tRule );
        tInterpolator.set_space_coeff( tCoordinates );
        tInterpolator.set_param_coeff();

        Integration_Rule tIntegrationRule( aCell->get_geometry_type(), Integration_Type::GAUSS, aCell->get_integration_order(), Integration_Type::GAUSS, Integration_Order::BAR_1 );
        auto             tCoefficients = tIntegrationRule.create_space_coeffs();
        Matrix< DDRMat > tPoints;
        Matrix< DDRMat > tWeights;
        tCoefficients->get_points( tPoints );
        tCoefficients->get_weights( tWeights );

        Matrix< DDRMat > tPointDerivatives;

        for ( uint iPoint = 0; iPoint < tWeights.numel(); ++iPoint )
        {
            tInterpolator.set_space( tPoints.get_column( iPoint ) );
            tInterpolator.space_det_J_derivatives( tPointDerivatives, aRequested );
            aDerivatives.matrix_data() += tWeights( iPoint ) * tPointDerivatives.matrix_data();
        }
    }

    // ---------------------------------------------------------------------------------

    moris::real
    Cell_Info::compute_cell_side_size_deriv(
            moris::mtk::Cell const *aCell,
            moris_index const      &aSideOrd,
            uint                    aLocalVertexID,
            uint                    aDirection ) const
    {
        MORIS_ERROR( false, "Cell_Info::compute_cell_side_size_deriv not implemented for this cell." );
        return 0.0;
    }

}    // namespace moris::mtk
