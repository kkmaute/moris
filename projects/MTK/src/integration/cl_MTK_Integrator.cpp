/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator.cpp
 *
 */

#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Cell_Info.hpp"

// LINALG/src
#include "op_times.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"
#include "fn_cross.hpp"
#include "fn_inv.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Integrator::Integrator( const Integration_Rule &aIntegrationRule )
            : mSpaceCoeffs( aIntegrationRule.create_space_coeffs() )
            , mTimeCoeffs( aIntegrationRule.create_time_coeffs() )
            , mNumOfSpacePoints( mSpaceCoeffs->get_number_of_points() )
            , mNumOfTimePoints( mTimeCoeffs->get_number_of_points() )
    {
        // matrix with space points
        mSpaceCoeffs->get_points( mSpacePoints );

        // matrix with time points
        mTimeCoeffs->get_points( mTimePoints );

        // matrix with space weights
        mSpaceCoeffs->get_weights( mSpaceWeights );

        // matrix with time weights
        mTimeCoeffs->get_weights( mTimeWeights );
    }

    //------------------------------------------------------------------------------

    Integrator::Integrator( const Integration_Rule &aIntegrationRule, const mtk::Set* aMeshSet, const bool aUseMomentFitting )
            : mSpaceCoeffs( aIntegrationRule.create_space_coeffs() )
            , mTimeCoeffs( aIntegrationRule.create_time_coeffs() )
            , mNumOfSpacePoints( mSpaceCoeffs->get_number_of_points() )
            , mNumOfTimePoints( mTimeCoeffs->get_number_of_points() )
    {
        // Set function pointers to compute integration points and weights depending on type of set and whether moment fitting is used,
        if ( aMeshSet->get_set_type() == mtk::SetType::BULK )
        {
            m_compute_cluster_integration_points_and_weights =
                    aUseMomentFitting ?
                            &Integrator::compute_bulk_cluster_integration_points_and_weights_moment_fitting :
                            &Integrator::compute_bulk_cluster_integration_points_and_weights_standard;
        }
        else if ( aMeshSet->get_set_type() == mtk::SetType::SIDESET )
        {
            // No moment fitting implementation currently implemented.
            m_compute_cluster_integration_points_and_weights =
                    &Integrator::compute_side_cluster_integration_points_and_weights_standard;
        }
        else if ( aMeshSet->get_set_type() == mtk::SetType::DOUBLE_SIDED_SIDESET )
        {
            // No moment fitting implementation currently implemented.
            m_compute_cluster_integration_points_and_weights =
                    &Integrator::compute_double_side_cluster_integration_points_and_weights_standard;
        }
        else
        {
            MORIS_ERROR( false, "Integrator constructor: mesh set has undefined geometry type" );
        }

        // matrix with space points
        mSpaceCoeffs->get_points( mSpacePoints );

        // matrix with time points
        mTimeCoeffs->get_points( mTimePoints );

        // matrix with space weights
        mSpaceCoeffs->get_weights( mSpaceWeights );

        // matrix with time weights
        mTimeCoeffs->get_weights( mTimeWeights );


        if ( aMeshSet->get_set_type() == mtk::SetType::BULK && aUseMomentFitting )
        {
            // If using moment fitting compute the value of the moment fitting LHS beforehand
            // Define interpolation rule to get moment fitting polynomials
            mtk::Interpolation_Rule tIPInterpolationRule( aMeshSet->get_interpolation_cell_geometry_type(), mtk::Interpolation_Type::LAGRANGE, aMeshSet->get_interpolation_cell_interpolation_order(),
                                                         mtk::Geometry_Type::LINE, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );

            // Create interpolation function to evaluate moment fitting polynomials at quadrature points
            mtk::Interpolation_Function_Base *tIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            // Get integration rule for moment fitting (cut phases will have triangle/tet interpolation rule, moment fitting integration rule alweays defined over the IP cell)
            mtk::Integration_Rule tMomentFittingIntegrationRule( aMeshSet->get_interpolation_cell_geometry_type(), mtk::Integration_Type::GAUSS,
                                                                this->get_ip_integration_order_from_cut_cell( aMeshSet->get_interpolation_cell_geometry_type(), aIntegrationRule ), 
                                                                mtk::Geometry_Type::LINE, mtk::Integration_Type::GAUSS, aIntegrationRule.get_time_integration_order() );

            // Get quadrature points for moment fitting
            mMomentFittingQuadPoints = this->get_points( tMomentFittingIntegrationRule );

            // Get the number of vertices on each IP cell
            uint tNumVertices = aMeshSet->get_clusters_on_set()( 0 )->get_interpolation_cell().get_cell_info()->get_num_verts();

            // Initialize the moment fitting LHS matrix
            Matrix< DDRMat > tMomentFittingLHS( tNumVertices , mMomentFittingQuadPoints.n_cols(), 0.0 );

            for ( uint iQuadPointIndex = 0; iQuadPointIndex < tMomentFittingLHS.n_cols(); iQuadPointIndex++ )
            {
                // Declare matrix for basis function values
                Matrix< DDRMat > tN;

                // Get quad point
                Matrix< DDRMat > tXi = mMomentFittingQuadPoints.get_column( iQuadPointIndex );

                // Get value of basis functions at quad point
                tIPInterp->eval_N( tXi, tN );

                // Place it in LHS
                tMomentFittingLHS.set_column( iQuadPointIndex, trans( tN ) );
            }

            mMomentFittingLHSinv = inv( tMomentFittingLHS );
        }
            
    }
    //------------------------------------------------------------------------------
    
    uint Integrator::get_number_of_points() const
    {
        return mNumOfSpacePoints * mNumOfTimePoints;
    }

    //------------------------------------------------------------------------------

    void Integrator::get_points( Matrix< DDRMat > &aIntegrationPoints ) const
    {
        // get number of dimensions in space
        uint tNumOfSpaceDim = mSpaceCoeffs->get_number_of_dimensions();

        // set output matrix size for space time
        aIntegrationPoints.set_size( tNumOfSpaceDim + 1,
                mNumOfSpacePoints * mNumOfTimePoints );

        Matrix< DDRMat > tOnes( 1, mNumOfSpacePoints, 1.0 );

        // loop over time
        uint startCol, stopCol;
        for ( uint k = 0; k < mNumOfTimePoints; ++k )
        {
            // indices for columns
            startCol = k * mNumOfSpacePoints;
            stopCol  = ( k + 1 ) * mNumOfSpacePoints - 1;

            // fill in the space points coordinates
            aIntegrationPoints( { 0, tNumOfSpaceDim - 1 }, { startCol, stopCol } ) =
                    mSpacePoints.matrix_data();

            // fill in the time point coordinates
            aIntegrationPoints( { tNumOfSpaceDim, tNumOfSpaceDim }, { startCol, stopCol } ) =
                    mTimePoints( k ) * tOnes;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Integrator::get_points() const
    {
        Matrix< DDRMat > tPoints;
        get_points( tPoints );
        return tPoints;
    }

    //------------------------------------------------------------------------------

    void Integrator::get_weights( Matrix< DDRMat > &aIntegrationWeights ) const
    {
        // get weights
        aIntegrationWeights = trans(
                vectorize( trans( mSpaceWeights ) * mTimeWeights ) );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Integrator::get_weights() const
    {
        Matrix< DDRMat > tWeights;
        get_weights( tWeights );
        return tWeights;
    }

    //------------------------------------------------------------------------------

    void Integrator::get_points( Matrix< DDRMat > &aIntegrationPoints, const Integration_Rule &aIntegrationRule ) const
    {
        // get number of dimensions in space
        uint tNumOfSpaceDim = aIntegrationRule.create_space_coeffs()->get_number_of_dimensions();

        // Get number of points in space from the integration rule
        uint tNumOfSpacePoints = aIntegrationRule.create_space_coeffs()->get_number_of_points();

        // Get number of points in time from the integration rule
        uint tNumOfTimePoints = aIntegrationRule.create_time_coeffs()->get_number_of_points();

        // Get the space points from the integration rule
        Matrix< DDRMat > tSpacePoints;
        aIntegrationRule.create_space_coeffs()->get_points( tSpacePoints );

        // Get the time points from the integration rule
        Matrix< DDRMat > tTimePoints;
        aIntegrationRule.create_time_coeffs()->get_points( tTimePoints );   

        // set output matrix size for space time
        aIntegrationPoints.set_size( tNumOfSpaceDim + 1,
                tNumOfSpacePoints * tNumOfTimePoints );

        Matrix< DDRMat > tOnes( 1, tNumOfSpacePoints, 1.0 );

        // loop over time
        uint startCol, stopCol;
        for ( uint k = 0; k < tNumOfTimePoints; ++k )
        {
            // indices for columns
            startCol = k * tNumOfSpacePoints;
            stopCol  = ( k + 1 ) * tNumOfSpacePoints - 1;

            // fill in the space points coordinates
            aIntegrationPoints( { 0, tNumOfSpaceDim - 1 }, { startCol, stopCol } ) =
                    tSpacePoints.matrix_data();

            // fill in the time point coordinates
            aIntegrationPoints( { tNumOfSpaceDim, tNumOfSpaceDim }, { startCol, stopCol } ) =
                    tTimePoints( k ) * tOnes;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Integrator::get_points( const Integration_Rule &aIntegrationRule ) const
    {
        Matrix< DDRMat > tPoints;
        get_points( tPoints, aIntegrationRule );
        return tPoints;
    }

    //------------------------------------------------------------------------------

    void Integrator::get_weights( Matrix< DDRMat > &aIntegrationWeights, const Integration_Rule &aIntegrationRule ) const
    {
        // get space weights from integration rule
        Matrix< DDRMat > tSpaceWeights;
        aIntegrationRule.create_space_coeffs()->get_weights( tSpaceWeights );
        
        // get time weights from integration rule
        Matrix< DDRMat > tTimeWeights;
        aIntegrationRule.create_time_coeffs()->get_weights( tTimeWeights );

        // get weights
        aIntegrationWeights = trans(
                vectorize( trans( tSpaceWeights ) * tTimeWeights ) );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Integrator::get_weights( const Integration_Rule &aIntegrationRule ) const
    {
        Matrix< DDRMat > tWeights;
        get_weights( tWeights, aIntegrationRule );
        return tWeights;
    }
    //------------------------------------------------------------------------------

    void Integrator::compute_bulk_cluster_integration_points_and_weights_standard( const Cluster *aCluster )
    {
        // Drop const and static cast into tCellCluster because we know func is for cell clusters and because we need to modify it
        Cell_Cluster *tCellCluster = static_cast< Cell_Cluster * >( const_cast< Cluster * >( aCluster ) );

        // if trivial, just get quadrature pts and weights for BG element
        // TODO : Add serendipity element support for trivial and non trivial cell clusters
        if ( tCellCluster->is_trivial() )
        {
            // assign quadrature weights and points
            tCellCluster->set_quadrature_points( this->get_points() );
            tCellCluster->set_quadrature_weights( this->get_weights() );
        }
        else
        {
            // Get IP element
            moris::mtk::Cell const &tInterpolationCell = tCellCluster->get_interpolation_cell();

            // Create IP cell interpolation rule
            mtk::Interpolation_Rule tIPInterpolationRule( tInterpolationCell.get_geometry_type(), mtk::Interpolation_Type::LAGRANGE, tInterpolationCell.get_interpolation_order(), mtk::Geometry_Type::LINE, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );

            // Get IG elements
            Vector< moris::mtk::Cell const * > const &tIGCells = tCellCluster->get_primary_cells_in_cluster();

            // Set the size of the quadrature points and weights matrices
            Matrix < DDRMat > tQuadraturePoints;
            Matrix < DDRMat > tQuadratureWeights;
            tQuadraturePoints.reshape( this->get_points().n_rows(), tIGCells.size() * this->get_weights().n_cols() );
            tQuadratureWeights.reshape( this->get_weights().n_rows(), tIGCells.size() * this->get_weights().n_cols() );

            // Loop over IG Elements
            for ( uint iPrimaryCell = 0; iPrimaryCell < tIGCells.size(); iPrimaryCell++ )
            {
                // Get cell info from tIGCell (required by space interpolator)
                const mtk::Cell_Info *tCellInfo = tIGCells( iPrimaryCell )->get_cell_info();

                // Create Interpolation Rule
                mtk::Interpolation_Rule tIGInterpolationRule( tIGCells( iPrimaryCell )->get_geometry_type(), mtk::Interpolation_Type::LAGRANGE, tIGCells( iPrimaryCell )->get_interpolation_order(), mtk::Geometry_Type::LINE, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );

                // Create space interpolator
                mtk::Space_Interpolator tIGSpaceInterpolator = mtk::Space_Interpolator(
                        tIGInterpolationRule,
                        tIPInterpolationRule,
                        tCellInfo->compute_cell_shape( tIGCells( iPrimaryCell ) ),
                        false );    // since the cell cluster is not a sideset this is false

                // Get parametric points from cluster and set them in space interpolator
                Matrix< DDRMat > tCellLocalCoords = tCellCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( iPrimaryCell );
                tIGSpaceInterpolator.set_space_param_coeff( tCellLocalCoords );

                // Set coeffs for the jacobian
                tIGSpaceInterpolator.set_space_coeff( tCellLocalCoords );

                // Loop over all quadrature points in local IG element space to IP element local space
                for ( uint iQuadPoint = 0; iQuadPoint < this->get_weights().numel(); iQuadPoint++ )
                {
                    // Set local coord point
                    tIGSpaceInterpolator.set_space_time( this->get_points().get_column( iQuadPoint ) );

                    // get mapped quadrature point (from IG element parent element space to IP element parent element space)
                    Matrix< DDRMat > tMappedQuadraturePoint = tIGSpaceInterpolator.map_integration_point();

                    // Place in quadrature point vector
                    tQuadraturePoints.set_column( iPrimaryCell * this->get_weights().numel() + iQuadPoint, ( tMappedQuadraturePoint ) );

                    // Get determinant
                    real tDetJ = tIGSpaceInterpolator.space_det_J();

                    // Get modified weight
                    real tWStar = this->get_weights()( iQuadPoint ) * tDetJ;

                    // Set quadrature point as multiplied with det_J
                    tQuadratureWeights( iPrimaryCell * this->get_weights().numel() + iQuadPoint ) = tWStar;
                }
            }

            // Set Quadrature weights and points inside the cluster
            tCellCluster->set_quadrature_points( tQuadraturePoints );
            tCellCluster->set_quadrature_weights( tQuadratureWeights );
        }
        
    }

    //------------------------------------------------------------------------------

    void Integrator::compute_bulk_cluster_integration_points_and_weights_moment_fitting( const Cluster *aCluster )
    {
        // Drop const and static cast into tCellCluster because we know func is for cell clusters and because we need to modify it
        Cell_Cluster *tCellCluster = static_cast< Cell_Cluster * >( const_cast< Cluster * >( aCluster ) );

        // If cell cluster is trivial simply set points and weights from the FEM Set Integration rule
        if ( tCellCluster->is_trivial() )
        {
            // assign quadrature weights and points
            tCellCluster->set_quadrature_points( this->get_points() );
            tCellCluster->set_quadrature_weights( this->get_weights() );
            return;
        }
        
        // Get the interpolation order and spatial dimension from the cluster
        Interpolation_Order tInterpOrder = tCellCluster->get_interpolation_cell().get_interpolation_order();
        uint tDim = mSpaceCoeffs->get_number_of_dimensions();

        uint tOrder = tInterpOrder == Interpolation_Order::LINEAR ? 1 : tInterpOrder == Interpolation_Order::QUADRATIC ? 2 : tInterpOrder == Interpolation_Order::CUBIC ? 3 : 0;

        // Determine the number of moments
        uint tNmoments = std::pow( tOrder + 1, tDim );

        // tMomentFittingLHS.reshape( tNmoments , tNmoments );

        // Generate LHS
        /*for (uint iQuadPointIndex = 0; iQuadPointIndex < tMomentFittingLHS.n_cols() ; iQuadPointIndex++)
        {
            // Declare matrix for basis function values
            Matrix< DDRMat > tN;

            // Get quad point
            Matrix< DDRMat > tXi = mQuadraturePoints.get_column( iQuadPointIndex );

            // Get value of basis functions at quad point
            mIPInterp->eval_N( tXi , tN );

            // Place it in LHS
            tMomentFittingLHS.set_column( iQuadPointIndex , trans( tN ) );

        }*/

        // create RHS vector
        Matrix< DDRMat > tMomentFittingRHS;

        // Allocate memory for RHS
        tMomentFittingRHS.reshape( tNmoments, 1 );

        // Obtain boundary facet element ordinals from cluster
        const Matrix< DDRMat >& tBoundaryFacetElementOrdinals = tCellCluster->get_boundary_facet_element_ordinals();

        // Compute the RHS
        for ( uint iFacetIndex = 0; iFacetIndex < tBoundaryFacetElementOrdinals.n_rows(); iFacetIndex++ )
        {
            // Get ID of cell which the facet belongs to
            moris_id tCellId = tBoundaryFacetElementOrdinals( iFacetIndex, 0 );

            // Get cell ordinal for the facet
            moris_index tElementOrdinal = tBoundaryFacetElementOrdinals( iFacetIndex, 1 );

            // Preallocate storage for facet nodal coordinates
            Matrix< DDRMat > tFacetVertexCoordinates( tDim, tDim, 0.0 );

            // Find the element with the given ID in the list of primary cells in the cluster
            const Vector< moris::mtk::Cell const * > &tPrimaryCells = tCellCluster->get_primary_cells_in_cluster();
            const moris::mtk::Cell *tCell = nullptr;
            for ( uint iCell = 0; iCell < tPrimaryCells.size(); iCell++ )
            {
                if ( tPrimaryCells( iCell )->get_id() == tCellId )
                {
                    tCell = tPrimaryCells( iCell );
                    break;
                }   
            }

            // Get node-to-facet map for the cell 
            moris::Matrix< moris::IndexMat > tNodeToFacetMap = tCell->get_cell_info_sp()->get_node_to_facet_map();
            
            // Declare vector to store vertex pointers
            Vector< mtk::Vertex* > tVertexPointers;
            tVertexPointers.resize( tNodeToFacetMap.n_cols() );

            for ( uint iNodeIndex = 0; iNodeIndex < tNodeToFacetMap.n_cols(); iNodeIndex++ )
            {
                // Get the vertex pointer corresponding to vertex index in node-to-facet map
                tVertexPointers( iNodeIndex ) = tCell->get_vertex_pointers()( tNodeToFacetMap( tElementOrdinal, iNodeIndex ) );
                
                // Get the coordinates of the vertex with respect to the interpolation cell 
                tFacetVertexCoordinates.set_row( iNodeIndex, tCellCluster->get_vertex_local_coordinate_wrt_interp_cell( tVertexPointers( iNodeIndex ) ) );
            }

            // Get subphase facet coordinates
            //Matrix< DDRMat > tFacetCoords = mFacetVertexCoordinates( iFacetIndex );

            // Get facet normal
            Matrix< DDRMat > tFacetNormal = tCell->compute_outward_side_normal( tElementOrdinal );

            // Get individual facet normal coordinates
            real tNx = tFacetNormal( 0, 0 );
            real tNy = tFacetNormal( 1, 0 );
            real tNz = 0.0;
            if ( tDim == 3 )
            {
                tNz = tFacetNormal( 2, 0 );
            }

            // Determine facet geometry type
            const mtk::Geometry_Type tGeometryType = ( tDim == 2 ) ? mtk::Geometry_Type::LINE : ( tDim == 3 ) ? mtk::Geometry_Type::TRI
                                                                                                              : mtk::Geometry_Type::UNDEFINED;
            // Determine facet integration order
            const mtk::Integration_Order tIntegrationOrder = ( tDim == 2 ) ? mtk::Integration_Order::BAR_6 : ( tDim == 3 ) ? mtk::Integration_Order::TRI_12
                                                                                                                           : mtk::Integration_Order::BAR_1;

            // Construct the rule using the const variables
            mtk::Integration_Rule tIntObjFacet(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );


            // Get quadrature points for facet
            Matrix< DDRMat > tIntPointsFacet;
            this->get_points( tIntPointsFacet, tIntObjFacet );

            // Get quadrature weights for facet
            Matrix< DDRMat > tIntWeightsFacet;
            this->get_weights( tIntWeightsFacet, tIntObjFacet );

            // Get geometric jacobian
            real tD = 0.0;

            if ( tDim == 2 )
            {
                tD = std::sqrt( std::pow( tFacetVertexCoordinates( 1, 0 ) - tFacetVertexCoordinates( 0, 0 ), 2 ) + std::pow( tFacetVertexCoordinates( 1, 1 ) - tFacetVertexCoordinates( 0, 1 ), 2 ) );
            }
            if ( tDim == 3 && tNz == 0 )
            {
                Matrix< DDRMat > tZeroMoments;
                tZeroMoments.set_size( tNmoments, 1, 0.0 );
                tMomentFittingRHS += tZeroMoments;
            }
            else
            {
                // Loop over facet quadrature points to compute the integral
                for ( uint iQuadPtIndex = 0; iQuadPtIndex < tIntWeightsFacet.numel(); iQuadPtIndex++ )
                {
                    // Get quad point
                    Matrix< DDRMat > tQuadPoint = tIntPointsFacet.get_column( iQuadPtIndex );

                    // Get quad weight
                    real tQuadWeight = tIntWeightsFacet( iQuadPtIndex );

                    real tXm = 0.0;
                    real tYm = 0.0;
                    real tZm = 0.0;
                    
                    // Get all x,y and (if applicable) z coordinates of facet vertices in separate vectors for ease of interpolation
                    Matrix< DDRMat > tXvector;
                    Matrix< DDRMat > tYvector;
                    Matrix< DDRMat > tZvector;

                    if ( tDim == 3 )
                    {
                        tXvector = { { tFacetVertexCoordinates( 0, 0 ) }, { tFacetVertexCoordinates( 1, 0 ) }, { tFacetVertexCoordinates( 2, 0 ) } };
                        tYvector = { { tFacetVertexCoordinates( 0, 1 ) }, { tFacetVertexCoordinates( 1, 1 ) }, { tFacetVertexCoordinates( 2, 1 ) } };
                        tZvector = { { tFacetVertexCoordinates( 0, 2 ) }, { tFacetVertexCoordinates( 1, 2 ) }, { tFacetVertexCoordinates( 2, 2 ) } };
                    }

                    if ( tDim == 2 )
                    {
                        // Get the mapped quad point value x
                        tXm = 0.5 * ( 1.0 - tQuadPoint( 0 ) ) * tFacetVertexCoordinates( 0, 0 ) + 0.5 * ( 1.0 + tQuadPoint( 0 ) ) * tFacetVertexCoordinates( 1, 0 );

                        // Get the mapped quad point value x
                        tYm = 0.5 * ( 1.0 - tQuadPoint( 0 ) ) * tFacetVertexCoordinates( 0, 1 ) + 0.5 * ( 1.0 + tQuadPoint( 0 ) ) * tFacetVertexCoordinates( 1, 1 );
                    }
                    else if ( tDim == 3 )
                    {
                        // Define interpolation rule for facet
                        mtk::Interpolation_Rule           tIGInterpolationRule( mtk::Geometry_Type::TRI, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Geometry_Type::LINE, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
                        mtk::Interpolation_Function_Base *tIGInterp = tIGInterpolationRule.create_space_interpolation_function();

                        // Determine the integration order to compute the antiderivative based on polynomial order of mesh
                        mtk::Integration_Order tAntiDerivIntOrder = mtk::Integration_Order::BAR_1;

                        if ( tOrder == 1 )
                        {
                            tAntiDerivIntOrder = mtk::Integration_Order::BAR_6;
                        }
                        if ( tOrder == 2 )
                        {
                            tAntiDerivIntOrder = mtk::Integration_Order::BAR_16;
                        }
                        if ( tOrder == 3 )
                        {
                            tAntiDerivIntOrder = mtk::Integration_Order::BAR_32;
                        }

                        // Define interpolation rule for getting the moment expressions
                        mtk::Interpolation_Rule           tMomentInterpolationRule( mtk::Geometry_Type::HEX, mtk::Interpolation_Type::LAGRANGE, tInterpOrder, mtk::Geometry_Type::LINE, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
                        mtk::Interpolation_Function_Base *tMomentInterp = tMomentInterpolationRule.create_space_interpolation_function();

                        // Declare matrix for basis functions
                        Matrix< DDRMat > tNmapping;

                        // Get the mapped quad point value x
                        tIGInterp->eval_N( tQuadPoint, tNmapping );

                        // get mapped quad point value x
                        Matrix< DDRMat > tXmat = tNmapping * tXvector;
                        tXm                    = tXmat( 0 );

                        // get mapped quad point value y
                        Matrix< DDRMat > tYmat = tNmapping * tYvector;
                        tYm                    = tYmat( 0 );

                        // get mapped quad point value z
                        Matrix< DDRMat > tZmat = tNmapping * tZvector;
                        tZm                    = tZmat( 0 );

                        // Allocate matrix for storing mapping
                        Matrix< DDRMat > tNximapping;

                        // Compute geometric Jacobian
                        tIGInterp->eval_dNdXi( tQuadPoint, tNximapping );

                        Matrix< DDRMat > tRxi_x = tNximapping.get_row( 0 ) * tXvector;
                        Matrix< DDRMat > tRxi_y = tNximapping.get_row( 0 ) * tYvector;
                        Matrix< DDRMat > tRxi_z = tNximapping.get_row( 0 ) * tZvector;


                        Matrix< DDRMat > tReta_x = tNximapping.get_row( 1 ) * tXvector;
                        Matrix< DDRMat > tReta_y = tNximapping.get_row( 1 ) * tYvector;
                        Matrix< DDRMat > tReta_z = tNximapping.get_row( 1 ) * tZvector;

                        Matrix< DDRMat > tRxi = { { tRxi_x( 0 ) }, { tRxi_y( 0 ) }, { tRxi_z( 0 ) } };

                        Matrix< DDRMat > tReta = { { tReta_x( 0 ) }, { tReta_y( 0 ) }, { tReta_z( 0 ) } };

                        Matrix< DDRMat > tCrossProdRxiReta = cross( tRxi, tReta );

                        tD = std::sqrt( std::pow( tCrossProdRxiReta( 0 ), 2 ) + std::pow( tCrossProdRxiReta( 1 ), 2 ) + std::pow( tCrossProdRxiReta( 2 ), 2 ) );

                        // Define 1D quadrature for evaluating the antiderivative
                        mtk::Integration_Rule tIntObjAntiDeriv( mtk::Geometry_Type::LINE, mtk::Integration_Type::GAUSS, tAntiDerivIntOrder, mtk::Geometry_Type::LINE, mtk::Integration_Type::GAUSS, mtk::Integration_Order::BAR_1 );
                        
                        // Get 1D quadrature points and weights
                        Matrix< DDRMat > tIntPointsAntiDeriv;
                        Matrix< DDRMat > tIntWeightsAntiDeriv;

                        this->get_weights( tIntWeightsAntiDeriv, tIntObjAntiDeriv );
                        this->get_points( tIntPointsAntiDeriv, tIntObjAntiDeriv );

                        for ( uint iAntiDerivQuadPoint = 0; iAntiDerivQuadPoint < tIntWeightsAntiDeriv.numel(); iAntiDerivQuadPoint++ )
                        {
                            // Get 1D quadrature point
                            Matrix< DDRMat > tAntiDerivQuadPoint = tIntPointsAntiDeriv.get_column( iAntiDerivQuadPoint );

                            // Develop mapping for x,y,z
                            real tSigma_z = 0.5 * ( 1.0 + tAntiDerivQuadPoint( 0 ) ) * ( tZm );

                            // Get final point coordinates for the moment functions
                            Matrix< DDRMat > tAntiDerivativeEvalPoint_z = { { tXm }, { tYm }, { tSigma_z } };

                            // Evaluate Moment fitting polynomial at this point
                            Matrix< DDRMat > tMomentValue_z;

                            tMomentInterp->eval_N( tAntiDerivativeEvalPoint_z, tMomentValue_z );

                            // Geometric jacobian will be 0.5 * tZm
                            // Compute value of moment
                            tMomentFittingRHS += 0.25 * trans( tMomentValue_z ) * tD * 0.5 * tZm * tNz * tQuadWeight * tIntWeightsAntiDeriv( iAntiDerivQuadPoint );
                        }
                    }


                    // Evaluate the moments
                    if ( tOrder == 1 && tDim == 2 )
                    {
                        tMomentFittingRHS( 0, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm - 0.5 * ( tXm * tXm ) ) * ( 1.0 - tYm ) + 0.25 * tNy * ( tYm - 0.5 * ( tYm * tYm ) ) * ( 1.0 - tXm ) );

                        tMomentFittingRHS( 1, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm + 0.5 * ( tXm * tXm ) ) * ( 1.0 - tYm ) + 0.25 * tNy * ( tYm - 0.5 * ( tYm * tYm ) ) * ( 1.0 + tXm ) );

                        tMomentFittingRHS( 2, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm + 0.5 * ( tXm * tXm ) ) * ( 1.0 + tYm ) + 0.25 * tNy * ( tYm + 0.5 * ( tYm * tYm ) ) * ( 1.0 + tXm ) );

                        tMomentFittingRHS( 3, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm - 0.5 * ( tXm * tXm ) ) * ( 1.0 + tYm ) + 0.25 * tNy * ( tYm + 0.5 * ( tYm * tYm ) ) * ( 1.0 - tXm ) );
                    }
                    else if ( tOrder == 2 && tDim == 2 )
                    {
                        tMomentFittingRHS( 0, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( -0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm ) * ( -tYm ) + 0.25 * tNy * ( -0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tYm, 3 ) ) * ( 1.0 - tXm ) * ( -tXm ) );

                        tMomentFittingRHS( 1, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( 0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm ) * ( -tYm ) + 0.25 * tNy * ( -0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tYm, 3 ) ) * ( 1.0 + tXm ) * ( tXm ) );

                        tMomentFittingRHS( 2, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( 0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 + tYm ) * ( tYm ) + 0.25 * tNy * ( 0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tYm, 3 ) ) * ( 1.0 + tXm ) * ( tXm ) );

                        tMomentFittingRHS( 3, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( -0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 + tYm ) * ( tYm ) + 0.25 * tNy * ( 0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tYm, 3 ) ) * ( 1.0 - tXm ) * ( -tXm ) );

                        tMomentFittingRHS( 4, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( 1.0 - tYm ) * ( -tYm ) * ( tXm - ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) + 0.5 * tNy * ( 1.0 - tXm * tXm ) * ( -0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tYm, 3 ) ) );

                        tMomentFittingRHS( 5, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( 1.0 - tYm * tYm ) * ( 0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) + 0.5 * tNy * ( tXm + 1.0 ) * ( tXm ) * ( tYm - ( 1.0 / 3.0 ) * ( std::pow( tYm, 3 ) ) ) );

                        tMomentFittingRHS( 6, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( tXm - ( 1.0 / 3.0 ) * ( std::pow( tXm, 3 ) ) ) * ( 1.0 + tYm ) * ( tYm ) + 0.5 * tNy * ( 1.0 - tXm * tXm ) * ( 0.5 * std::pow( tYm, 2 ) + ( 1.0 / 3.0 ) * ( std::pow( tYm, 3 ) ) ) );

                        tMomentFittingRHS( 7, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( -0.5 * std::pow( tXm, 2 ) + ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm * tYm ) + 0.5 * tNy * ( 1.0 - tXm ) * ( -tXm ) * ( tYm - ( 1.0 / 3.0 ) * ( std::pow( tYm, 3 ) ) ) );

                        tMomentFittingRHS( 8, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 1.0 * tNx * ( tXm - ( 1.0 / 3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm * tYm ) + 1.0 * tNy * ( tYm - ( 1.0 / 3.0 ) * ( std::pow( tYm, 3 ) ) ) * ( 1.0 - tXm * tXm ) );
                    }
                    else if ( tOrder == 3 && tDim == 2 )
                    {
                        real tt2 = std::pow( tXm, 2 );
                        real tt3 = std::pow( tXm, 3 );
                        real tt4 = std::pow( tYm, 2 );
                        real tt5 = std::pow( tYm, 3 );

                        tMomentFittingRHS( 0, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) - tt5 * ( 9.0 / 5.12e+2 ) + tYm / 5.12e+2 - 1.0 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 - tt5 * 7.91015625e-2 + tYm * 8.7890625e-3 - 8.7890625e-3 ) - tXm * ( tt4 * ( 9.0 / 2.56e+2 ) - tt5 * ( 9.0 / 2.56e+2 ) + tYm / 2.56e+2 - 1.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 0, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) - tt3 * ( 9.0 / 5.12e+2 ) + tXm / 5.12e+2 - 1.0 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 - tt3 * 7.91015625e-2 + tXm * 8.7890625e-3 - 8.7890625e-3 ) - tYm * ( tt2 * ( 9.0 / 2.56e+2 ) - tt3 * ( 9.0 / 2.56e+2 ) + tXm / 2.56e+2 - 1.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 1, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) - tt5 * ( 9.0 / 5.12e+2 ) + tYm / 5.12e+2 - 1.0 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 - tt5 * 7.91015625e-2 + tYm * 8.7890625e-3 - 8.7890625e-3 ) - tXm * ( tt4 * ( 9.0 / 2.56e+2 ) - tt5 * ( 9.0 / 2.56e+2 ) + tYm / 2.56e+2 - 1.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 1, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) + tt3 * ( 9.0 / 5.12e+2 ) - tXm / 5.12e+2 - 1.0 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 + tt3 * 7.91015625e-2 - tXm * 8.7890625e-3 - 8.7890625e-3 ) - tYm * ( tt2 * ( 9.0 / 2.56e+2 ) + tt3 * ( 9.0 / 2.56e+2 ) - tXm / 2.56e+2 - 1.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 2, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) + tt5 * ( 9.0 / 5.12e+2 ) - tYm / 5.12e+2 - 1.0 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 + tt5 * 7.91015625e-2 - tYm * 8.7890625e-3 - 8.7890625e-3 ) - tXm * ( tt4 * ( 9.0 / 2.56e+2 ) + tt5 * ( 9.0 / 2.56e+2 ) - tYm / 2.56e+2 - 1.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 2, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) + tt3 * ( 9.0 / 5.12e+2 ) - tXm / 5.12e+2 - 1.0 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 + tt3 * 7.91015625e-2 - tXm * 8.7890625e-3 - 8.7890625e-3 ) - tYm * ( tt2 * ( 9.0 / 2.56e+2 ) + tt3 * ( 9.0 / 2.56e+2 ) - tXm / 2.56e+2 - 1.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 3, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) + tt5 * ( 9.0 / 5.12e+2 ) - tYm / 5.12e+2 - 1.0 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 + tt5 * 7.91015625e-2 - tYm * 8.7890625e-3 - 8.7890625e-3 ) - tXm * ( tt4 * ( 9.0 / 2.56e+2 ) + tt5 * ( 9.0 / 2.56e+2 ) - tYm / 2.56e+2 - 1.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 3, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) - tt3 * ( 9.0 / 5.12e+2 ) + tXm / 5.12e+2 - 1.0 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 - tt3 * 7.91015625e-2 + tXm * 8.7890625e-3 - 8.7890625e-3 ) - tYm * ( tt2 * ( 9.0 / 2.56e+2 ) - tt3 * ( 9.0 / 2.56e+2 ) + tXm / 2.56e+2 - 1.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 4, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) - tt5 * ( 2.43e+2 / 5.12e+2 ) + tYm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 - tt5 * 2.373046875e-1 + tYm * 2.63671875e-2 - 2.63671875e-2 ) + tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 4, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) - tt3 * ( 2.7e+1 / 5.12e+2 ) + tXm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 - tt3 * 2.373046875e-1 + tXm * 2.373046875e-1 - 7.91015625e-2 ) + tYm * ( tt2 * ( 9.0 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 5, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) - tt5 * ( 2.43e+2 / 5.12e+2 ) + tYm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 - tt5 * 2.373046875e-1 + tYm * 2.63671875e-2 - 2.63671875e-2 ) + tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 5, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) + tt3 * ( 2.7e+1 / 5.12e+2 ) - tXm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 + tt3 * 2.373046875e-1 - tXm * 2.373046875e-1 - 7.91015625e-2 ) + tYm * ( tt2 * ( 9.0 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 6, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) - tt5 * ( 2.7e+1 / 5.12e+2 ) + tYm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 - tt5 * 2.373046875e-1 + tYm * 2.373046875e-1 - 7.91015625e-2 ) + tXm * ( tt4 * ( 9.0 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 6, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) + tt3 * ( 2.43e+2 / 5.12e+2 ) - tXm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 + tt3 * 2.373046875e-1 - tXm * 2.63671875e-2 - 2.63671875e-2 ) + tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 7, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) + tt5 * ( 2.7e+1 / 5.12e+2 ) - tYm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 + tt5 * 2.373046875e-1 - tYm * 2.373046875e-1 - 7.91015625e-2 ) + tXm * ( tt4 * ( 9.0 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 7, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) + tt3 * ( 2.43e+2 / 5.12e+2 ) - tXm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 + tt3 * 2.373046875e-1 - tXm * 2.63671875e-2 - 2.63671875e-2 ) + tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 8, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) + tt5 * ( 2.43e+2 / 5.12e+2 ) - tYm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 + tt5 * 2.373046875e-1 - tYm * 2.63671875e-2 - 2.63671875e-2 ) + tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 8, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) + tt3 * ( 2.7e+1 / 5.12e+2 ) - tXm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 + tt3 * 2.373046875e-1 - tXm * 2.373046875e-1 - 7.91015625e-2 ) + tYm * ( tt2 * ( 9.0 / 2.56e+2 ) + tt3 * ( 2.7e+1 / 2.56e+2 ) - tXm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 9, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) + tt5 * ( 2.43e+2 / 5.12e+2 ) - tYm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 + tt5 * 2.373046875e-1 - tYm * 2.63671875e-2 - 2.63671875e-2 ) + tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 9, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 9.0 / 5.12e+2 ) - tt3 * ( 2.7e+1 / 5.12e+2 ) + tXm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 7.91015625e-2 - tt3 * 2.373046875e-1 + tXm * 2.373046875e-1 - 7.91015625e-2 ) + tYm * ( tt2 * ( 9.0 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 10, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) + tt5 * ( 2.7e+1 / 5.12e+2 ) - tYm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 + tt5 * 2.373046875e-1 - tYm * 2.373046875e-1 - 7.91015625e-2 ) + tXm * ( tt4 * ( 9.0 / 2.56e+2 ) + tt5 * ( 2.7e+1 / 2.56e+2 ) - tYm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 10, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) - tt3 * ( 2.43e+2 / 5.12e+2 ) + tXm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 - tt3 * 2.373046875e-1 + tXm * 2.63671875e-2 - 2.63671875e-2 ) + tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 11, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( -std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 9.0 / 5.12e+2 ) - tt5 * ( 2.7e+1 / 5.12e+2 ) + tYm * ( 2.7e+1 / 5.12e+2 ) - 9.0 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 7.91015625e-2 - tt5 * 2.373046875e-1 + tYm * 2.373046875e-1 - 7.91015625e-2 ) + tXm * ( tt4 * ( 9.0 / 2.56e+2 ) - tt5 * ( 2.7e+1 / 2.56e+2 ) + tYm * ( 2.7e+1 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );
                        tMomentFittingRHS( 11, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( -std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 2.7e+1 / 2.56e+2 ) + tXm * ( 3.0 / 2.56e+2 ) - 3.0 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) - tt3 * ( 2.43e+2 / 5.12e+2 ) + tXm * ( 2.7e+1 / 5.12e+2 ) - 2.7e+1 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 - tt3 * 2.373046875e-1 + tXm * 2.63671875e-2 - 2.63671875e-2 ) + tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 9.0 / 2.56e+2 ) - 9.0 / 2.56e+2 ) );

                        tMomentFittingRHS( 12, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) - tt5 * ( 7.29e+2 / 5.12e+2 ) + tYm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 - tt5 * 7.119140625e-1 + tYm * 7.119140625e-1 - 2.373046875e-1 ) - tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) - tt5 * ( 2.43e+2 / 2.56e+2 ) + tYm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );
                        tMomentFittingRHS( 12, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) - tt3 * ( 7.29e+2 / 5.12e+2 ) + tXm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 - tt3 * 7.119140625e-1 + tXm * 7.119140625e-1 - 2.373046875e-1 ) - tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) - tt3 * ( 2.43e+2 / 2.56e+2 ) + tXm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );

                        tMomentFittingRHS( 13, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) - tt5 * ( 8.1e+1 / 2.56e+2 ) + tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) - tt5 * ( 7.29e+2 / 5.12e+2 ) + tYm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 - tt5 * 7.119140625e-1 + tYm * 7.119140625e-1 - 2.373046875e-1 ) - tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) - tt5 * ( 2.43e+2 / 2.56e+2 ) + tYm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );
                        tMomentFittingRHS( 13, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) + tt3 * ( 7.29e+2 / 5.12e+2 ) - tXm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) - std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 + tt3 * 7.119140625e-1 - tXm * 7.119140625e-1 - 2.373046875e-1 ) - tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) + tt3 * ( 2.43e+2 / 2.56e+2 ) - tXm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );

                        tMomentFittingRHS( 14, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) + tt5 * ( 7.29e+2 / 5.12e+2 ) - tYm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) + std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 + tt5 * 7.119140625e-1 - tYm * 7.119140625e-1 - 2.373046875e-1 ) - tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) + tt5 * ( 2.43e+2 / 2.56e+2 ) - tYm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );
                        tMomentFittingRHS( 14, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) + tt3 * ( 8.1e+1 / 2.56e+2 ) - tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) + tt3 * ( 7.29e+2 / 5.12e+2 ) - tXm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 + tt3 * 7.119140625e-1 - tXm * 7.119140625e-1 - 2.373046875e-1 ) - tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) + tt3 * ( 2.43e+2 / 2.56e+2 ) - tXm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );

                        tMomentFittingRHS( 15, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx * ( std::pow( tXm, 3 ) * ( tt4 * ( 2.7e+1 / 2.56e+2 ) + tt5 * ( 8.1e+1 / 2.56e+2 ) - tYm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) + std::pow( tXm, 2 ) * ( tt4 * ( 2.43e+2 / 5.12e+2 ) + tt5 * ( 7.29e+2 / 5.12e+2 ) - tYm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) - std::pow( tXm, 4 ) * ( tt4 * 2.373046875e-1 + tt5 * 7.119140625e-1 - tYm * 7.119140625e-1 - 2.373046875e-1 ) - tXm * ( tt4 * ( 8.1e+1 / 2.56e+2 ) + tt5 * ( 2.43e+2 / 2.56e+2 ) - tYm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );
                        tMomentFittingRHS( 15, 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy * ( std::pow( tYm, 3 ) * ( tt2 * ( 2.7e+1 / 2.56e+2 ) - tt3 * ( 8.1e+1 / 2.56e+2 ) + tXm * ( 8.1e+1 / 2.56e+2 ) - 2.7e+1 / 2.56e+2 ) - std::pow( tYm, 2 ) * ( tt2 * ( 2.43e+2 / 5.12e+2 ) - tt3 * ( 7.29e+2 / 5.12e+2 ) + tXm * ( 7.29e+2 / 5.12e+2 ) - 2.43e+2 / 5.12e+2 ) + std::pow( tYm, 4 ) * ( tt2 * 2.373046875e-1 - tt3 * 7.119140625e-1 + tXm * 7.119140625e-1 - 2.373046875e-1 ) - tYm * ( tt2 * ( 8.1e+1 / 2.56e+2 ) - tt3 * ( 2.43e+2 / 2.56e+2 ) + tXm * ( 2.43e+2 / 2.56e+2 ) - 8.1e+1 / 2.56e+2 ) );
                    }
                }


                /*else if ( tOrder == 1 && aDim == 3 )
                {
                    Matrix< DDRMat > tMomentFittingRHSTest;
                    tMomentFittingRHSTest.reshape( 8 , 1 );
                    tMomentFittingRHSTest( 0 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm - 1),2)*(tYm - 1.0)*(tZm - 1.0)) / 16.0 + tNy * -((tXm - 1.0)*std::pow((tYm - 1.0),2)*(tZm - 1.0))/16.0 + tNz * -((tXm - 1.0)*(tYm - 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 1 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm + 1.0),2)*(tYm - 1.0)*(tZm - 1.0))/16.0 + tNy * ((tXm + 1.0)*std::pow((tYm - 1.0),2)*(tZm - 1.0))/16.0 + tNz * ((tXm + 1.0)*(tYm - 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 2 , 0 ) +=  0.5 *( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm + 1.0),2)*(tYm + 1.0)*(tZm - 1.0))/16.0 + tNy * -((tXm + 1.0)*std::pow((tYm + 1.0),2)*(tZm - 1.0))/16.0 + tNz * -((tXm + 1.0)*(tYm + 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 3 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm - 1.0),2)*(tYm + 1.0)*(tZm - 1.0))/16.0 + tNy * ((tXm - 1.0)*std::pow((tYm + 1.0),2)*(tZm - 1.0))/16.0 + tNz * ((tXm - 1.0)*(tYm + 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 4 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm - 1.0),2)*(tYm - 1.0)*(tZm + 1.0))/16.0 + tNy * ((tXm - 1.0)*std::pow((tYm - 1.0),2)*(tZm + 1.0))/16.0 + tNz * ((tXm - 1.0)*(tYm - 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 5 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm + 1.0),2)*(tYm - 1.0)*(tZm + 1.0))/16.0 + tNy * -((tXm + 1.0)*std::pow((tYm - 1.0),2)*(tZm + 1.0))/16.0 + tNz * -((tXm + 1.0)*(tYm - 1.0)*std::pow((tZm + 1.0),2))/16.0 );

                    tMomentFittingRHSTest( 6 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm + 1.0),2)*(tYm + 1.0)*(tZm + 1.0))/16.0 + tNy * ((tXm + 1.0)*std::pow((tYm + 1.0),2)*(tZm + 1.0))/16.0 + tNz * ((tXm + 1.0)*(tYm + 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;

                    tMomentFittingRHSTest( 7 , 0 ) +=  0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx* -(std::pow((tXm - 1),2)*(tYm + 1.0)*(tZm + 1.0))/16.0 + tNy * -((tXm - 1.0)*std::pow((tYm + 1.0),2)*(tZm + 1.0))/16.0 + tNz * -((tXm - 1.0)*(tYm + 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;







                }*/
            }


            // Get coordinate transform


            /*real tu1  = tFacetCoords( 0 , 0 );
            real tv1  = tFacetCoords( 0 , 1 );

            real tu2  = tFacetCoords( 1 , 0 );
            real tv2  = tFacetCoords( 1 , 1 );

            real tNx  = tFacetNormal( 0 , 0 );
            real tNy  = tFacetNormal( 1 , 0 );

            real tD   = std::sqrt(std::pow( tu2-tu1 , 2) + std::pow( tv2-tv1 , 2 ));

            real t21 = tu1*tu1;
            real t31 = -tu2;
            real t41 = -tv2;
            real t51 = tv1-1.0;
            real t61 = tu1/4.0;
            real t71 = tu2/4.0;
            real t81 = t31+tu1;
            real t91 = t41+tv1;
            real t101 = -t61;
            real t121 = t21/8.0;
            real t111 = t81*t81;
            //real t131 = -t121;
            real t141 = t61*t81;
            //real t151 = t61+t131;
            real t161 = t71+t101+t141;

            tMomentFittingRHS( 0 , 0 ) += tD*0.5*tNx*( t51*(t101+t121) - (t91*(t101+t121))/2.0 + (t51*t111)/2.4e+1 - (t91*t111)/3.2e+1 - (t51*t161)/2.0 + (t91*t161)/3.0);

            real t22 = tv1*tv1;
            real t32 = -tv1;
            real t42 = -tv2;
            real t52 = tu1/4.0;
            real t62 = tu2/4.0;
            real t72 = t42+tv1;
            real t82 = -t62;
            real t92 = t22/2.0;
            real t142 = t52-(1.0/4.0);
            real t102 = t72*t72;
            //real t112 = -t92;
            real t122 = t72*tv1;
            real t152 = t52+t82;
            //real t132 = t112+tv1;
            real t162 = t32+t122+tv2;

            tMomentFittingRHS( 0 , 0 ) += tD*0.5*tNy*(t142*(t32+t92) - (t152*(t32+t92))/2.0 + (t102*t142)/6.0 - (t102*t152)/8.0 - (t142*t162)/2.0 + (t152*t162)/3.0);

            real t23 = tu1*tu1;
            real t33 = -tu2;
            real t43 = -tv2;
            real t53 = tv1-1.0;
            real t63 = tu1/4.0;
            real t73 = tu2/4.0;
            real t83 = t33+tu1;
            real t93 = t43+tv1;
            real t103 = -t73;
            real t123 = t23/8.0;
            real t113 = t83*t83;
            real t133 = t63*t83;
            real t143 = t63+t123;
            real t153 = t63+t103+t133;

            tMomentFittingRHS( 1 , 0 ) += tD*0.5*tNx*(t53*t113*(-1.0/2.4e+1) - t53*t143 + (t53*t153)/2.0 + (t93*t113)/(3.2e+1) + (t93*t143)/2.0 - (t93*t153)/3.0);

            real t24 = tv1*tv1;
            real t34 = -tv1;
            real t44 = -tv2;
            real t54 = tu1/4.0;
            real t64 = tu2/4.0;
            real t74 = t44+tv1;
            real t84 = -t64;
            real t94 = t24/2.0;
            real t134 = t54+1.0/4.0;
            real t104 = t74*t74;
            //real t114 = -t94;
            real t124 = t74*tv1;
            real t154 = t54+t84;
            //real t144 = t114+tv1;
            real t164 = t34+t124+tv2;

            tMomentFittingRHS( 1 , 0 ) += tD*0.5*tNy*((-t134)*(t34+t94) + (t154*(t34+t94))/2.0 -(t104*t134)/6.0 + (t104*t154)/8.0 + (t134*t164)/2.0 - (t154*t164)/3.0);

            real t25 = tu1*tu1;
            real t35 = tv1+1.0;
            real t45 = -tu2;
            real t55 = -tv2;
            real t65 = tu1/4.0;
            real t75 = tu2/4.0;
            real t85 = t45+tu1;
            real t95 = t55+tv1;
            real t105 = -t75;
            real t125 = t25/8.0;
            real t115 = t85*t85;
            real t135 = t65*t85;
            real t145 = t65+t125;
            real t155 = t65+t105+t135;

            tMomentFittingRHS( 2 , 0 ) += tD*0.5*tNx*( (t35*t115)/2.4e+1 + t35*t145 - (t35*t155)/2.0 - (t95*t115)/3.2e+1 - (t95*t145)/2.0 + (t95*t155)/3.0);

            real t26 = tv1*tv1;
            real t36 = -tv2;
            real t46 = tu1/4.0;
            real t56 = tu2/4.0;
            real t66 = t36+tv1;
            real t76 = -t56;
            real t86 = t26/2.0;
            real t126 = t46+1.0/4.0;
            real t96 = t66*t66;
            real t106 = t66*tv1;
            real t116 = t86+tv1;
            real t136 = t46+t76;
            real t146 = t66+t106;

            tMomentFittingRHS( 2 , 0 ) += tD*0.5*tNy*((t96*t126)/6.0 - (t96*t136)/8.0 + t116*t126 - (t116*t136)/2.0 - (t126*t146)/2.0 + (t136*t146)/3.0);

            real t27 = tu1*tu1;
            real t37 = tv1+1.0;
            real t47 = -tu2;
            real t57 = -tv2;
            real t67 = tu1/4.0;
            real t77 = tu2/4.0;
            real t87 = t47+tu1;
            real t97 = t57+tv1;
            real t107 = -t67;
            real t127 = t27/8.0;
            real t117 = t87*t87;
            //real t137 = -t127;
            real t147 = t67*t87;
            //real t157 = t67+t137;
            real t167 = t77+t107+t147;

            tMomentFittingRHS( 3 , 0 ) += tD*0.5*tNx*( -t37*(t107+t127) + (t97*(t107+t127))/2.0 - (t37*t117)/(2.4e+1) + (t37*t167)/2.0 + (t97*t117)/(3.2e+1) - (t97*t167)/3.0);

            real t28 = tv1*tv1;
            real t38 = -tv2;
            real t48 = tu1/4.0;
            real t58 = tu2/4.0;
            real t68 = t38+tv1;
            real t78 = -t58;
            real t88 = t28/2.0;
            real t128 = t48-(1.0/4.0);
            real t98 = t68*t68;
            real t108 = t68*tv1;
            real t118 = t88+tv1;
            real t138 = t48+t78;
            real t148 = t68+t108;

            tMomentFittingRHS( 3 , 0 ) += tD*0.5*tNy*(t98*t128*(-1.0/6.0) + (t98*t138)/8.0 - t118*t128 + (t118*t138)/2.0 + (t128*t148)/2.0- (t138*t148)/3.0);*/


            /*real ta11 = tx0 - tx1;
            real ta12 = ty0 - ty1;
            real tb11  = 1.0 - tx0;
            real tb12  = 1.0 - ty0;

            real ta21 = tx1 - tx0;
            real ta22 = ty0 - ty1;
            real tb21  = 1.0 + tx0;
            real tb22  = 1.0 - ty0;

            real ta31 = tx1 - tx0;
            real ta32 = ty1 - ty0;
            real tb31  = 1.0 + tx0;
            real tb32  = 1.0 + ty0;

            real ta41 = tx0 - tx1;
            real ta42 = ty1 - ty0;
            real tb41  = 1.0 - tx0;
            real tb42  = 1.0 + ty0;

            // Compute the RHS
            tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) - ta11 * tNx * (1.0/4.0) * ( ( ta12/ 4*ta11 )*( std::pow( (ta11 + tb11) , 4 ) -  std::pow( tb11 , 4 ) ) + ( ta12* tb11 / (3 * ta11) )*( std::pow( (ta11 + tb11) , 3 ) -  std::pow( tb11 , 3 ) ) + (tb12 / 3)*( std::pow( (ta11 + tb11) , 3 ) -  std::pow( tb11 , 3 ) ) );

            tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) - ta12 * tNy * (1.0/4.0) * ( ( ta11/ 4*ta12 )*( std::pow( (ta12 + tb12) , 4 ) -  std::pow( tb12 , 4 ) ) + ( ta11* tb12 / (3 * ta12) )*( std::pow( (ta12 + tb12) , 3 ) -  std::pow( tb11 , 3 ) ) + (tb12 / 3)*( std::pow( (ta12 + tb12) , 3 ) -  std::pow( tb12 , 3 ) ) );


            tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) + ta21 * tNx * (1.0/4.0) * ( ( ta22/ 4*ta21 )*( std::pow( (ta21 + tb21) , 4 ) -  std::pow( tb21 , 4 ) ) + ( ta22* tb21 / (3 * ta21) )*( std::pow( (ta21 + tb21) , 3 ) -  std::pow( tb21 , 3 ) ) + (tb22 / 3)*( std::pow( (ta21 + tb21) , 3 ) -  std::pow( tb21 , 3 ) ) );

            tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) - ta22 * tNy * (1.0/4.0) * ( ( ta21/ 4*ta22 )*( std::pow( (ta22 + tb22) , 4 ) -  std::pow( tb22 , 4 ) ) + ( ta21* tb22 / (3 * ta22) )*( std::pow( (ta22 + tb22) , 3 ) -  std::pow( tb21 , 3 ) ) + (tb22 / 3)*( std::pow( (ta22 + tb22) , 3 ) -  std::pow( tb22 , 3 ) ) );


            tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + ta31 * tNx * (1.0/4.0) * ( ( ta32/ 4*ta31 )*( std::pow( (ta31 + tb31) , 4 ) -  std::pow( tb31 , 4 ) ) + ( ta32* tb31 / (3 * ta31) )*( std::pow( (ta31 + tb31) , 3 ) -  std::pow( tb31 , 3 ) ) + (tb32 / 3)*( std::pow( (ta31 + tb31) , 3 ) -  std::pow( tb31 , 3 ) ) );

            tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + ta32 * tNy * (1.0/4.0) * ( ( ta31/ 4*ta32 )*( std::pow( (ta32 + tb32) , 4 ) -  std::pow( tb32 , 4 ) ) + ( ta31* tb32 / (3 * ta32) )*( std::pow( (ta32 + tb32) , 3 ) -  std::pow( tb31 , 3 ) ) + (tb32 / 3)*( std::pow( (ta32 + tb32) , 3 ) -  std::pow( tb32 , 3 ) ) );


            tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) - ta41 * tNx * (1.0/4.0) * ( ( ta42/ 4*ta41 )*( std::pow( (ta41 + tb41) , 4 ) -  std::pow( tb41 , 4 ) ) + ( ta42* tb41 / (3 * ta41) )*( std::pow( (ta41 + tb41) , 3 ) -  std::pow( tb41 , 3 ) ) + (tb42 / 3)*( std::pow( (ta41 + tb41) , 3 ) -  std::pow( tb41 , 3 ) ) );

            tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) + ta42 * tNy * (1.0/4.0) * ( ( ta41/ 4*ta42 )*( std::pow( (ta42 + tb42) , 4 ) -  std::pow( tb42 , 4 ) ) + ( ta41* tb42 / (3 * ta42) )*( std::pow( (ta42 + tb42) , 3 ) -  std::pow( tb41 , 3 ) ) + (tb42 / 3)*( std::pow( (ta42 + tb42) , 3 ) -  std::pow( tb42 , 3 ) ) );*/

            // Compute the RHS

            // tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5 * ( (std::pow( (tx0 - tx1) , 2)*( ty0 - 1.0 ))/24.0 - (std::pow( (tx0 - tx1) , 2)*( ty0 - ty1 ))/32.0 + ((( 1 / 8.0 ) * std::pow( tx0 , 2 ) -  tx0 / 4.0 )*( ty0 - ty1 ))/2.0 )*/

            // tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5*tNx*(( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + 0.5 * (( 1.0 / 8.0 )* (- std::pow( tx0 , 2 )) + ( 1.0 / 4.0 )*( tx0 ) )*( ty0 - ty1 ) + (1.0 / 3.0)*( ty0 - ty1 )*( tx1/4.0 - tx0/4.0 + (tx0/4.0) * (tx0 - tx1) ) - ( -std::pow( tx0 , 2 )/8.0 + tx0/4.0 )*(ty0 - 1) - 0.5*(ty0 - 1.0)*(tx1/4.0 - tx0/4.0 + (tx0/4.0)*(tx0 - tx1) ));

            // tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5*tNy*( 0.25*0.5*(tx0 - tx1)*(0.5 * std::pow( ty0 , 2 ) + ty0) - (1.0/32.0)*(tx0-tx1)*(std::pow( (ty0-ty1), 2)) + (1.0/12.0)*(tx0 - tx1)*(ty1 - ty0 + ty0*(ty0 - ty1)) + (1.0/24.0)*(tx0 - 1.0)*(std::pow( (ty0-ty1), 2)) - ( -0.5*std::pow( ty0, 2 ) + ty0 )*0.25*( tx0 - 1.0 ) - 0.5*0.25*( tx0 - 1 )*( ty1- ty0 + ty0*(ty0 - ty1) ) );

            // tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) + 0.5*tNx*( ( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + (( 1.0/8.0 )*std::pow( tx0 , 2) - tx0/4.0 )*0.5*( ty0 -ty1 ) - (1.0/3.0)*(ty0 - ty1)*( tx0/4.0 - tx1/4.0 + (tx0/4.0)*( tx0 - tx1 )) - (std::pow(tx0,2)/8.0 + tx0/4.0)*(ty0 - 1.0) + 0.5*(ty0 - 1.0)*(tx0/4.0 - tx1/4.0 + (tx0/4.0)*(tx0 - tx1) ) );

            // tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) + 0.5*tNy*((1.0/32.0)*(tx0 - tx1)*(std::pow( (ty0 - ty1), 2)) - (1.0/8.0)*(tx0 - tx1)*( 0.5*std::pow( ty0 , 2 ) + ty0 ) - (1.0/12.0)*(tx0 - tx1)*(ty1 - ty0 + ty0*(ty0 - ty1)) -  (1.0/24.0)*(std::pow(ty0 - ty1, 2))*( tx0 + 1.0) + 0.25*(-std::pow(ty0,2) + ty0 ) + 0.25*0.5*(tx0 + 1)*(ty1 - ty0 + ty0*(ty0 - ty1)));

            // tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + 0.5*tNx*( ( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 + 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + 0.5*((1.0/8.0)*std::pow(tx0,2)+ tx0/4.0 )*(ty0 - ty1) + (1.0/3.0)*(ty0 -ty1)*0.25*(tx0 - tx1 + tx0*(tx0 - tx1)) + 0.25*(0.5*std::pow(tx0,2) + tx0)*(ty0 + 1.0) - 0.5*0.25*(ty0 + 1.0)*(tx0 - tx1 + tx0*(tx0 - tx1)));

            // tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + 0.5*tNy*( (1.0/12.0)*(tx0 - tx1)*(ty0 - ty1 + ty0*(ty0 - ty1)) - 0.5*0.25*(tx0 - tx1)*(0.5*std::pow(ty0,2)+ty0) - (1.0/32.0)*(tx0 - tx1)*(std::pow(ty0 - ty1,2)) + (1.0/24.0)*(tx0 + 1.0)*(std::pow(ty0 - ty1, 2)) + 0.25*(0.5*std::pow(ty0,2) + ty0)*(tx0 + 1.0) - 0.5*0.25*(tx0 + 1.0)*(ty0 - ty1 +ty0*(ty0 - ty1)));

            // tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) + 0.5*tNx*( -( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 + 1.0 ) + ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) - 0.25*(-0.5*std::pow(tx0 , 2) + tx0 )*( ty0 - ty1 ) - (1.0/12.0)*(ty0 - 1)*( tx1 - tx0 + tx0*(tx0 - tx1)) + 0.25*(ty0 + 1.0)*( -0.5*std::pow(tx0 , 2) + tx0 ) + (1.0/8.0)*(ty0 + 1)*( tx1 - tx0 + tx0*(tx0 - tx1)) );

            // tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) + 0.5*tNy*( (1.0/32.0)*(tx0 -tx1)*(std::pow(ty0 - ty1, 2)) + 0.25*0.5*(tx0 - tx1)*( 0.5*std::pow(ty0, 2) + ty0 ) - (1.0/12.0)*(tx0 - tx1)*(ty0 - ty1 + ty0*(ty0 - ty1)) - (1.0/24.0)*(tx0 - 1.0)*(std::pow(ty0 - ty1, 2)) + 0.25*(tx0 - 1)*(0.25*std::pow(ty0,2)) + 0.5*0.25*(tx0 - 1)*( ty0 - ty1 + ty0*(ty0 - ty1)));

            // fprintf(stdout, "Facet Index %d \n", iFacetIndex );
        }

        // Solve the system
        // mQuadratureWeights = {{0.0}, {0.0}, {0.0}, {0.0}};


        Matrix< DDRMat > tQuadratureWeights = ( mMomentFittingLHSinv * tMomentFittingRHS );

        // Reshape so as to be compatible with format required downstream
        tQuadratureWeights.reshape( 1, tNmoments );

        // add in the points and weights inside the cluster
        tCellCluster->set_quadrature_points( mMomentFittingQuadPoints );
        tCellCluster->set_quadrature_weights( tQuadratureWeights );
    }

    //------------------------------------------------------------------------------

    void Integrator::compute_side_cluster_integration_points_and_weights_standard( const Cluster *aCluster )
    {
        MORIS_ERROR( false, "Standard integration for side clusters not yet implemented" );
    }

    //------------------------------------------------------------------------------

    void Integrator::compute_double_side_cluster_integration_points_and_weights_standard( const Cluster *aCluster )
    {
        MORIS_ERROR( false, "Moment fitting for side clusters not yet implemented" );
    }

    //------------------------------------------------------------------------------

    mtk::Integration_Order Integrator::get_ip_integration_order_from_cut_cell( const mtk::Geometry_Type &aGeometryType, const Integration_Rule &aIntegrationRule ) const
    {
        // Implementation for determining integration order from cut cell
        mtk::Integration_Order tIntegrationOrder = aIntegrationRule.get_space_integration_order();
        switch ( tIntegrationOrder )
        {
            case mtk::Integration_Order::TRI_7:
                return mtk::Integration_Order::QUAD_2x2;
            
            case mtk::Integration_Order::TRI_12:
                return mtk::Integration_Order::QUAD_3x3;
            
            case mtk::Integration_Order::TRI_25:
                return mtk::Integration_Order::QUAD_4x4;

            case mtk::Integration_Order::QUAD_2x2:
                return mtk::Integration_Order::QUAD_2x2;

            case mtk::Integration_Order::QUAD_3x3:
                return mtk::Integration_Order::QUAD_3x3;

            case mtk::Integration_Order::QUAD_4x4:
                return mtk::Integration_Order::QUAD_4x4;

            case mtk::Integration_Order::TET_11:
                return mtk::Integration_Order::HEX_2x2x2;

            case mtk::Integration_Order::TET_35:
                return mtk::Integration_Order::HEX_3x3x3;

            case mtk::Integration_Order::TET_56:
                return mtk::Integration_Order::HEX_4x4x4;

            case mtk::Integration_Order::HEX_2x2x2:
                return mtk::Integration_Order::HEX_2x2x2;

            case mtk::Integration_Order::HEX_3x3x3:
                return mtk::Integration_Order::HEX_3x3x3;

            case mtk::Integration_Order::HEX_4x4x4:
                return mtk::Integration_Order::HEX_4x4x4;

            default:
                return mtk::Integration_Order::UNDEFINED;
        }
    }

}    // namespace moris::mtk



