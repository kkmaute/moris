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

namespace moris
{
    namespace mtk
    {
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
            Matrix< DDRMat >             tXHat           = aCell->get_vertex_coords();
            Interpolation_Function_Base *tInterpFunction = tInterpolationRule.create_space_interpolation_function();
            Matrix< DDRMat >             tXiHat;
            tInterpFunction->get_param_coords( tXiHat );
            tXiHat = trans( tXiHat );

            // setting space coefficients
            tSpaceInterpolator.set_space_coeff( tXHat );
            tSpaceInterpolator.set_space_param_coeff( tXiHat );

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
            Matrix< DDRMat >             tXHat           = aCell->get_vertex_coords();
            Interpolation_Function_Base *tInterpFunction = tInterpolationRule.create_space_interpolation_function();
            Matrix< DDRMat >             tXiHat;
            tInterpFunction->get_param_coords( tXiHat );
            tXiHat = trans( tXiHat );

            // setting space coefficients
            tSpaceInterpolator.set_space_coeff( tXHat );
            tSpaceInterpolator.set_space_param_coeff( tXiHat );

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

    }    // namespace mtk
}    // namespace moris
