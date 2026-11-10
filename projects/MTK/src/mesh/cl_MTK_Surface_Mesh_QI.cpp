/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh_QI.cpp
 *
 */
#include "cl_MTK_Surface_Mesh_QI.hpp"
#include "cl_Library_IO.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "op_elemwise_mult.hpp"

namespace moris::mtk
{
    Surface_Mesh_QI::Surface_Mesh_QI( const Surface_Mesh& aSurfaceMesh )
            : mSurfaceMesh( aSurfaceMesh )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_raycast_wedge_area( const Vector< moris_index >& aSourceFacets,
            Matrix< DDRMat >&                                                      aVertexSensitivities ) const
    {
        MORIS_ERROR( mSurfaceMesh.get_spatial_dimension() == 2,
                "Surface_Mesh_QI::compute_raycast_wedge_area - Requires a 2D surface mesh." );
        MORIS_ERROR( aSourceFacets.size() > 0,
                "Surface_Mesh_QI::compute_raycast_wedge_area - The source loop is empty." );

        // Collect each source vertex once: neighboring source facets share a ray and its hit.
        // Keep the mesh's local vertex indices for assembling the final sensitivity rows.
        const uint     tNumVertices = mSurfaceMesh.get_number_of_vertices();
        Vector< bool > tIsSource( tNumVertices, false );
        Vector< uint > tSourceVertices;
        for ( moris_index iF : aSourceFacets )
        {
            MORIS_ERROR( iF >= 0 && static_cast< uint >( iF ) < mSurfaceMesh.get_number_of_facets(),
                    "Surface_Mesh_QI::compute_raycast_wedge_area - Source facet index out of bounds." );
            for ( moris_index iV : mSurfaceMesh.get_facets_vertex_indices( iF ) )
            {
                if ( !tIsSource( iV ) )
                {
                    tIsSource( iV ) = true;
                    tSourceVertices.push_back( iV );
                }
            }
        }

        // Reverse the outward length-weighted normals to cast inward toward the target.
        // Subtract from +0.0 to keep zero components positive for ArborX ray-box tests.
        const Matrix< DDRMat >     tCoordinates = mSurfaceMesh.get_all_vertex_coordinates();
        const Matrix< DDRMat >     tNormals     = 0.0 - mSurfaceMesh.compute_vertex_normals();
        Matrix< DDRMat >           tOrigins( 2, tSourceVertices.size() );
        Vector< Matrix< DDRMat > > tDirections( tSourceVertices.size() );
        for ( uint i = 0; i < tSourceVertices.size(); i++ )
        {
            tOrigins.set_column( i, tCoordinates.get_column( tSourceVertices( i ) ) );
            tDirections( i ) = tNormals.get_column( tSourceVertices( i ) );
        }

        // The batch is indexed by source vertex, then ray direction (one direction here).
        // Store hit data by mesh vertex so either incident source facet can reuse it.
        Vector< Vector< bool > > tWarnings;
        auto                     tRayIntersections = mSurfaceMesh.cast_batch_of_rays( tOrigins, tDirections, tWarnings );
        Vector< real >           tDistances( tNumVertices, 0.0 );
        Vector< uint >           tHitFacets( tNumVertices, MORIS_UINT_MAX );
        Matrix< DDRMat >         tIntersections( 2, tNumVertices, 0.0 );
        for ( uint i = 0; i < tSourceVertices.size(); i++ )
        {
            uint iV = tSourceVertices( i );
            for ( const auto& tHit : tRayIntersections( i )( 0 ) )
            {
                // Hits are sorted by distance. Skip the source vertex's self-intersection
                // and retain the closest hit strictly ahead of it.
                if ( tHit.second > 1e-12 )
                {
                    tDistances( iV ) = tHit.second;
                    tHitFacets( iV ) = tHit.first;
                    break;
                }
            }

            MORIS_ERROR( tHitFacets( iV ) != MORIS_UINT_MAX, "Surface_Mesh_QI::compute_raycast_wedge_area - No nontrivial intersection for source vertex %u.", iV );

            // The hit facet must belong entirely to the fixed target boundary; hitting the
            // moving source loop would violate the sensitivity assumptions below.
            for ( moris_index iHitVertex : mSurfaceMesh.get_facets_vertex_indices( tHitFacets( iV ) ) )
            {
                MORIS_ERROR( !tIsSource( iHitVertex ), "Surface_Mesh_QI::compute_raycast_wedge_area - Source normals must point toward the fixed target boundary." );
            }
            tIntersections.set_column( iV, tCoordinates.get_column( iV ) + tDistances( iV ) * tNormals.get_column( iV ) );
        }

        // Separate direct source-coordinate derivatives from derivatives through ray hits.
        // Both are row gradients; target vertex rows remain zero because the target is fixed.
        aVertexSensitivities.set_size( tNumVertices, 2, 0.0 );
        Matrix< DDRMat > tAreaDIntersection( tNumVertices, 2, 0.0 );
        real             tArea = 0.0;
        for ( moris_index iF : aSourceFacets )
        {
            const auto& tVertices = mSurfaceMesh.get_facets_vertex_indices( iF );
            // Quadrilateral order: source vertex 0, source vertex 1, ray hit 1, ray hit 0.
            Matrix< DDRMat > tQuad( 2, 4 );
            tQuad.set_column( 0, tCoordinates.get_column( tVertices( 0 ) ) );
            tQuad.set_column( 1, tCoordinates.get_column( tVertices( 1 ) ) );
            tQuad.set_column( 2, tIntersections.get_column( tVertices( 1 ) ) );
            tQuad.set_column( 3, tIntersections.get_column( tVertices( 0 ) ) );
            // Shoelace formula: 2*A_signed = sum_i (x_i*y_(i+1) - x_(i+1)*y_i).
            // Sum absolute wedge areas so reversing the quadrilateral order preserves area.
            real tTwiceSignedArea = 0.0;
            for ( uint i = 0; i < 4; i++ )
            {
                uint tNext = ( i + 1 ) % 4;
                tTwiceSignedArea += tQuad( 0, i ) * tQuad( 1, tNext ) - tQuad( 0, tNext ) * tQuad( 1, i );
            }
            MORIS_ERROR( std::abs( tTwiceSignedArea ) > MORIS_REAL_EPS,
                    "Surface_Mesh_QI::compute_raycast_wedge_area - A zero-area wedge has no unique absolute-area derivative." );
            tArea += 0.5 * std::abs( tTwiceSignedArea );
            // With a fixed orientation, d|A_signed| = sign(A_signed)*dA_signed.
            real tSign = tTwiceSignedArea > 0.0 ? 1.0 : -1.0;
            for ( uint i = 0; i < 4; i++ )
            {
                uint             tPrevious = ( i + 3 ) % 4;
                uint             tNext     = ( i + 1 ) % 4;
                Matrix< DDRMat > tGradient = { { 0.5 * tSign * ( tQuad( 1, tNext ) - tQuad( 1, tPrevious ) ),
                        0.5 * tSign * ( tQuad( 0, tPrevious ) - tQuad( 0, tNext ) ) } };
                // Corners 0 and 1 move directly with the source. Corners 2 and 3
                // are hits for source vertices 1 and 0, respectively, and need the chain rule.
                if ( i < 2 )
                {
                    uint iV = tVertices( i );
                    aVertexSensitivities.set_row( iV, aVertexSensitivities.get_row( iV ) + tGradient );
                }
                else
                {
                    uint iV = tVertices( 3 - i );
                    tAreaDIntersection.set_row( iV, tAreaDIntersection.get_row( iV ) + tGradient );
                }
            }
        }

        // For the hit q = o + t*n, the distance t changes with both origin o and direction n:
        // dq = (I + n*dt/do)*do + (t*I + n*dt/dn)*dn.
        // Apply the accumulated dArea/dq once per ray, including both neighboring wedges.
        // The target facet is fixed; these derivatives assume its identity stays unchanged.
        Matrix< DDRMat > tIdentity = eye( 2, 2 );
        for ( uint iV : tSourceVertices )
        {
            Matrix< DDRMat > tOrigin      = tCoordinates.get_column( iV );
            Matrix< DDRMat > tNormal      = tNormals.get_column( iV );
            Matrix< DDRMat > tDHitDOrigin = tIdentity + tNormal * mSurfaceMesh.compute_draycast_dorigin( tOrigin, tNormal, tHitFacets( iV ) );
            Matrix< DDRMat > tDHitDNormal = tDistances( iV ) * tIdentity + tNormal * mSurfaceMesh.compute_draycast_ddirection( tOrigin, tNormal, tHitFacets( iV ) );
            Matrix< DDRMat > tDADNormal   = tAreaDIntersection.get_row( iV ) * tDHitDNormal;
            aVertexSensitivities.set_row( iV, aVertexSensitivities.get_row( iV ) + tAreaDIntersection.get_row( iV ) * tDHitDOrigin );

            // The unit normal is formed from sum_f(length_f * unit_normal_f).
            // Its derivative includes facet-normal motion, changing length weights, and
            // normalization. Distribute the direction term to its incident source vertices.
            Vector< moris_index > tDependencies = mSurfaceMesh.get_vertex_neighbors( iV );
            tDependencies.push_back( iV );
            for ( moris_index iDependency : tDependencies )
            {
                if ( tIsSource( iDependency ) )
                {
                    // The ray direction is the negative mesh normal, so its Jacobian
                    // must carry the same sign as the direction used in the raycast.
                    Matrix< DDRMat > tDNormalDVertex = 0.0 - mSurfaceMesh.compute_dvertex_normal_dvertex( iV, iDependency );
                    aVertexSensitivities.set_row( iDependency, aVertexSensitivities.get_row( iDependency ) + tDADNormal * tDNormalDVertex );
                }
            }
        }
        return tArea;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_QI( QI_Type aType, const Parameter_List& aParameters, const std::shared_ptr< Library_IO >& aLibrary )
    {
        // A failed evaluation must not leave callbacks or sensitivities from a previous QI active.
        mQIType              = QI_Type::END_ENUM;
        mSensitivityFunction = nullptr;
        real tValue          = 0.0;

        switch ( aType )
        {
            case QI_Type::USER_DEFINED:
            {
                const std::string& tFunctionName    = aParameters.get< std::string >( "function_name" );
                const std::string& tSensitivityName = aParameters.get< std::string >( "sensitivity_function_name" );
                MORIS_ERROR( !tFunctionName.empty() && !tSensitivityName.empty(),
                        "Surface_Mesh_QI::compute_QI - USER_DEFINED requires function_name and sensitivity_function_name." );
                MORIS_ERROR( aLibrary != nullptr,
                        "Surface_Mesh_QI::compute_QI - USER_DEFINED requires an input library." );

                auto tFunction    = aLibrary->load_function< Surface_Mesh_QI_Function >( tFunctionName );
                auto tSensitivity = aLibrary->load_function< Surface_Mesh_QI_Sensitivity_Function >( tSensitivityName );
                return this->compute_QI( tFunction, tSensitivity );
            }
            case QI_Type::VOLUME:
                tValue = this->compute_volume();
                break;

            case QI_Type::RAYCAST_SHAPE_DIAMETER:
            case QI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER:
            case QI_Type::SHORTEST_DISTANCE_SHAPE_DIAMETER:
            {
                auto tAgglom    = create_agglomeration_parameters( aParameters, aLibrary );
                real tConeAngle = aParameters.get< real >( "cone_angle" );

                if ( aType == QI_Type::RAYCAST_SHAPE_DIAMETER )
                {
                    tValue = this->compute_global_shape_diameter_raycast( *tAgglom, tConeAngle, aParameters.get< moris_index >( "number_of_polar_rays" ), aParameters.get< moris_index >( "number_of_azimuth_rays" ) );
                }
                else if ( aType == QI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER )
                {
                    tValue = this->compute_global_shape_diameter_inscribed_circle( *tAgglom, tConeAngle, aParameters.get< real >( "minimum_relative_chord_length" ), aParameters.get< moris_index >( "number_of_samples" ) );
                }
                else
                {
                    tValue = this->compute_global_shape_diameter_shortest_distance( *tAgglom, tConeAngle, aParameters.get< moris_index >( "number_of_samples" ) );
                }
                break;
            }
            default:
                MORIS_ERROR( false, "Surface_Mesh_QI::compute_QI - Unsupported QI type." );
                return 0.0;
        }

        mQIType = aType;
        return tValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_QI( Surface_Mesh_QI_Function aFunction,
            Surface_Mesh_QI_Sensitivity_Function               aSensitivityFunction )
    {
        mQIType              = QI_Type::END_ENUM;
        mSensitivityFunction = nullptr;
        MORIS_ERROR( aFunction != nullptr && aSensitivityFunction != nullptr,
                "Surface_Mesh_QI::compute_QI - Both user-defined callbacks are required." );
        real tValue          = aFunction( mSurfaceMesh );
        mQIType              = QI_Type::USER_DEFINED;
        mSensitivityFunction = aSensitivityFunction;
        return tValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_QI::compute_dQI_dvertex( uint aVertexIndex ) const
    {
        MORIS_ERROR( aVertexIndex < mSurfaceMesh.get_number_of_vertices(),
                "Surface_Mesh_QI::compute_dQI_dvertex - Vertex index out of bounds." );

        Matrix< DDRMat > tSensitivity;
        switch ( mQIType )
        {
            case QI_Type::USER_DEFINED:
                tSensitivity = mSensitivityFunction( mSurfaceMesh, aVertexIndex );
                break;
            case QI_Type::VOLUME:
                tSensitivity = this->compute_dvolume_dvertex( aVertexIndex );
                break;
            case QI_Type::RAYCAST_SHAPE_DIAMETER:
            case QI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER:
            case QI_Type::SHORTEST_DISTANCE_SHAPE_DIAMETER:
                tSensitivity = this->compute_ddiameter_dvertex( aVertexIndex );
                break;
            default:
                MORIS_ERROR( false, "Surface_Mesh_QI::compute_dQI_dvertex - Compute the QI before requesting sensitivities." );
        }

        MORIS_ERROR( tSensitivity.n_rows() == 1 && tSensitivity.n_cols() == mSurfaceMesh.get_spatial_dimension(),
                "Surface_Mesh_QI::compute_dQI_dvertex - Sensitivity must be a 1 x spatial dimension row matrix." );
        return tSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_QI::compute_QI_sensitivities( const Vector< Vector< moris_index > >& aVertexPDVIDs,
            sol::Dist_Vector*                                                              aSensitivities,
            uint                                                                           aRequestIndex ) const
    {
        MORIS_ERROR( aSensitivities != nullptr, "Surface_Mesh_QI::compute_QI_sensitivities - Missing sensitivity vector." );
        MORIS_ERROR( aVertexPDVIDs.size() == mSurfaceMesh.get_spatial_dimension(),
                "Surface_Mesh_QI::compute_QI_sensitivities - PDV map must have one entry per spatial dimension." );
        for ( const auto& tIDs : aVertexPDVIDs )
        {
            MORIS_ERROR( tIDs.size() == mSurfaceMesh.get_number_of_vertices(),
                    "Surface_Mesh_QI::compute_QI_sensitivities - PDV map must contain every surface vertex." );
        }

        for ( uint iV = 0; iV < mSurfaceMesh.get_number_of_vertices(); iV++ )
        {
            bool tHasPDV = false;
            for ( const auto& tIDs : aVertexPDVIDs )
            {
                tHasPDV = tHasPDV || tIDs( iV ) != -1;
            }
            if ( !tHasPDV )
            {
                continue;
            }

            Matrix< DDRMat > tSensitivity = this->compute_dQI_dvertex( iV );
            for ( uint iDim = 0; iDim < aVertexPDVIDs.size(); iDim++ )
            {
                moris_index tPDVID = aVertexPDVIDs( iDim )( iV );
                if ( tPDVID != -1 )
                {
                    ( *aSensitivities )( tPDVID, aRequestIndex ) += tSensitivity( iDim );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::sech( real aValue )
    {
        return 1.0 / std::cosh( aValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::tanh_clip( real aValue, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex )
    {
        return std::pow( aValue - aAgglom.get_ref( aCoords, aIndex ), aAgglom.get_exp() ) * ( 1.0 - std::tanh( 10.0 * ( aValue - aAgglom.get_ref( aCoords, aIndex ) ) ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::dtanh_clip( real aValue, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex )
    {
        return aAgglom.get_exp() * std::pow( aValue - aAgglom.get_ref( aCoords, aIndex ), aAgglom.get_exp() - 1.0 ) * ( 1.0 - std::tanh( 10.0 * ( aValue - aAgglom.get_ref( aCoords, aIndex ) ) ) )
             - 10.0 * std::pow( aValue - aAgglom.get_ref( aCoords, aIndex ), aAgglom.get_exp() ) * Surface_Mesh_QI::sech( 10.0 * ( aValue - aAgglom.get_ref( aCoords, aIndex ) ) ) * Surface_Mesh_QI::sech( 10.0 * ( aValue - aAgglom.get_ref( aCoords, aIndex ) ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::max_clip( real aValue, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex )
    {
        return std::pow( std::max( 1.0 - ( ( aValue - aAgglom.get_shift() ) / aAgglom.get_ref( aCoords, aIndex ) ), 0.0 ), aAgglom.get_exp() );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::dmax_clip( real aValue, const Agglomeration_Parameters& aAgglom, const Matrix< DDRMat >& aCoords, uint aIndex )
    {
        return ( ( 1.0 - ( ( aValue - aAgglom.get_shift() ) / aAgglom.get_ref( aCoords, aIndex ) ) ) > 0.0 ) ? aAgglom.get_exp() * std::pow( 1.0 - ( ( aValue - aAgglom.get_shift() ) / aAgglom.get_ref( aCoords, aIndex ) ), aAgglom.get_exp() - 1.0 ) / -aAgglom.get_ref( aCoords, aIndex ) : 0.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::dnormal_dsigma( real aValue, real aStdDev )
    {
        return ( aValue * aValue - aStdDev * aStdDev ) * std::exp( -0.5 * std::pow( aValue / aStdDev, 2.0 ) ) / ( aStdDev * aStdDev * aStdDev * aStdDev * std::sqrt( 2.0 * M_PI ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_QI::write_shape_diameters_to_file( const std::string& aFilePath ) const
    {
        std::ofstream tFile( aFilePath );
        tFile << std::fixed << std::setprecision( 8 );

        for ( uint iF = 0; iF < mSurfaceMesh.get_number_of_facets(); iF++ )
        {
            tFile << mShapeDiameters( iF ) << "\n";
        }
        tFile.close();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_volume() const
    {
        real tVolume = 0.0;

        // Apply shoelace formula for 2D surface mesh BRENDAN NOT GENERALIZABLE
        if ( mSurfaceMesh.get_spatial_dimension() == 2 )
        {
            for ( uint iFacet = 0; iFacet < mSurfaceMesh.get_number_of_facets(); iFacet++ )
            {
                Matrix< DDRMat > tCoords = mSurfaceMesh.get_all_vertex_coordinates_of_facet( iFacet );

                tVolume += 0.5 * ( tCoords( 0, 0 ) * tCoords( 1, 1 ) - tCoords( 0, 1 ) * tCoords( 1, 0 ) );
            }
        }
        else
        {
            MORIS_ERROR( false, "Surface Mesh volume computation only implemented for 2D (lines) meshes" );
        }

        return tVolume;
    }

    //--------------------------------------------------------------------------------------------------------------

    Shape_Diameter_Data Surface_Mesh_QI::cast_shape_diameter_ray_cones( const real aConeAngle, const uint aNumPolarRays, uint aNumAzimuthRays ) const
    {
        // Get the current vertex coordinates for the entire surface mesh
        Matrix< DDRMat > tRayOrigins = mSurfaceMesh.compute_facet_centroids();

        // Initialize the output struct - get the direction weights from the ray cones
        Shape_Diameter_Data tDiameterData( aConeAngle, aNumPolarRays, aNumAzimuthRays, mSurfaceMesh.get_all_facet_normals(), aNumPolarRays * aNumAzimuthRays );

        // Batch process all rays for all vertices
        Vector< Vector< bool > >                tWarnings;
        Vector< Vector< Intersection_Vector > > tIntersections = mSurfaceMesh.cast_batch_of_rays( tRayOrigins, tDiameterData.get_directions(), tWarnings );

        // Get just the nearest intersection for each ray
        tDiameterData.store_nearest_nontrivial_intersections( tIntersections );

        return tDiameterData;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real > Surface_Mesh_QI::compute_raycast_shape_diameter(
            const Agglomeration_Parameters& aAgglom,
            real                            aConeAngle,
            uint                            aNumPolarRays,
            uint                            aNumAzimuthRays )
    {
        bool tRayDropping = false;    // brendan delete, temporary for debugging
        real tAlphaPower  = 4.0;      // brendan delete, temporary for debugging

        // Get the minimum distance for every ray cone at every facet
        Shape_Diameter_Data tDiameterData = this->cast_shape_diameter_ray_cones( aConeAngle, aNumPolarRays, aNumAzimuthRays );

        // Get the number of facets in the mesh
        uint tNumFacets = mSurfaceMesh.get_number_of_facets();

        // Get the number of rays per facet
        uint tNumRays = tDiameterData.get_num_rays_per_cone();

        // Check that we got a distance for every facet and ray
        MORIS_ASSERT( tDiameterData.get_distances().size() == tNumFacets, "Inconsistent number of facets in the shape diameter distances." );

        // Initialize return variable
        Vector< real > tShapeDiameterValues( tNumFacets, MORIS_REAL_MAX );

        // Reset the derivative
        mdShapeDiameterdVertex.set_size( mSurfaceMesh.get_number_of_vertices(), mSurfaceMesh.get_spatial_dimension(), 0.0 );

        // Centroids needed for raycast sensitivities
        Matrix< DDRMat > tCentroids = mSurfaceMesh.compute_facet_centroids();

        // Get the facet measures for sensitivity calculations
        Vector< real > tFacetMeasures = mSurfaceMesh.compute_facet_measure();

        // Loop through every facet in this surface mesh
        for ( uint iF = 0; iF < tNumFacets; iF++ )
        {
            MORIS_ASSERT( tDiameterData.get_distances()( iF ).size() == tNumRays, "Inconsistent number of rays in cone for facet %d (Expected %d, Got %ld)", iF, tNumRays, tDiameterData.get_distances()( iF ).size() );

            // Get the ray distances for this facet (one distance per ray)
            const Intersection_Vector& tRaysOnFacet = tDiameterData.get_distances()( iF );

            Intersection_Vector tSortedRaysOnFacet = tRaysOnFacet;    // FIXME annoying copy

            // Get the vertex indices for the origin facet, these vertices need their sensitivities updated
            const Vector< moris_index >& tFacetVertices          = mSurfaceMesh.get_facets_vertex_indices( iF );
            Matrix< DDRMat >             tFacetVertexCoordinates = mSurfaceMesh.get_all_vertex_coordinates_of_facet( iF );

            // Compute the median ray distance for this facet
            std::sort( tSortedRaysOnFacet.begin(), tSortedRaysOnFacet.end(), []( const std::pair< moris_index, real >& aR1, const std::pair< moris_index, real >& aR2 ) { return aR1.second < aR2.second; } );
            real tMedian;
            if ( tNumRays % 2 == 0 )
            {
                tMedian = ( tSortedRaysOnFacet( tNumRays / 2 - 1 ).second + tSortedRaysOnFacet( tNumRays / 2 ).second ) / 2;
            }
            else
            {
                tMedian = tSortedRaysOnFacet( tNumRays / 2 ).second;
            }

            // Compute the mean ray distance for this facet
            real tMean = 0.0;
            for ( const auto& tDistance : tRaysOnFacet )
            {
                tMean += tDistance.second;
            }
            tMean = tMean / tNumRays;

            // Compute the standard deviation ray distances for this facet
            real tStdDeviation = 0.0;
            for ( const auto& tDistance : tRaysOnFacet )
            {
                tStdDeviation += ( tDistance.second - tMean ) * ( tDistance.second - tMean );    // variance
            }
            tStdDeviation = std::sqrt( tStdDeviation / (real)tNumRays );    // take sqrt for std deviation

            // Compute the weighted sum for entries within one standard deviation from the median
            real tDiamWeightedSum   = 0.0;
            real tSumOfValidWeights = 0.0;

            // Compute the sum of weights. Need to do this first to have ready for sensitivities
            for ( uint iR = 0; iR < tNumRays; iR++ )
            {
                if ( !tRayDropping or std::abs( tRaysOnFacet( iR ).second - tMedian ) <= tStdDeviation )
                {
                    // Compute the angle between the ray direction and the intersected facets normal
                    Matrix< DDRMat > tIntersectedFacetNormal = mSurfaceMesh.get_facet_normal( tRaysOnFacet( iR ).first );

                    real tAngle = Surface_Mesh::compute_dot_absolute( tDiameterData.get_ray_direction( iF, iR ), tIntersectedFacetNormal );    // Angle factor raised to a power to increase differentiation between rays

                    real tAugmentedWeight = tDiameterData.get_ray_weight( iR ) * std::pow( tAngle, tAlphaPower );    // Weight augmented by angle factor

                    // Sum the weights and the weighted shape diameter values
                    tSumOfValidWeights += tAugmentedWeight;
                    tDiamWeightedSum += tAugmentedWeight * this->max_clip( tDiameterData.get_distances()( iF )( iR ).second, aAgglom, tFacetVertexCoordinates, iF );
                }
            }

            MORIS_ERROR( tSumOfValidWeights > 0.0, "Sum of valid weights is zero for facet %d. Check ray directions and intersections.", iF );

            // Divide by the sum of weights to get the weighted average
            tShapeDiameterValues( iF ) = tDiamWeightedSum / tSumOfValidWeights;

            // ---------------------------------------------------------------------------
            // STEP 2: Compute the sensitivities of this facet, loop back through the rays
            // ---------------------------------------------------------------------------

            // Each ray intersection has sensitivities wrt to 2*dim vertices (dim vertices that create the origin facet and dim vertices that create the hit facet).
            // Note that if a ray hits the neighboring facet, the origin facet vertices and the intersected facet vertices can share some vertices
            // We compute the derivative of the shape diameter to the vertices and store the accumulation of sensitivity for all vertices
            // For each vertex ddiameter/dv = dr/dx * weight / sum_of_weights

            for ( uint iR = 0; iR < tNumRays; iR++ )
            {
                if ( !tRayDropping or std::abs( tRaysOnFacet( iR ).second - tMedian ) <= tStdDeviation )
                {
                    // Augment the theta weight by the angle between the ray and the intersected facet normal
                    Matrix< DDRMat > tIntersectedFacetNormal = mSurfaceMesh.get_facet_normal( tRaysOnFacet( iR ).first );
                    real             tAngle                  = Surface_Mesh::compute_dot_absolute( tDiameterData.get_ray_direction( iF, iR ), tIntersectedFacetNormal );
                    real             tAugmentedWeight        = tDiameterData.get_ray_weight( iR ) * std::pow( tAngle, tAlphaPower );    // Weight augmented by angle factor
                    real             tDistance               = tDiameterData.get_distance( iF, iR ).second;
                    real             tAgglomeratedDistance   = this->max_clip( tDistance, aAgglom, tFacetVertexCoordinates, iF );

                    // Compute the sensitivity factor components
                    real tAgglomerationTerm = this->dmax_clip( tDistance, aAgglom, tFacetVertexCoordinates, iF );    // derivative of tanh clipping function
                    real tWeightTerm        = tAugmentedWeight / tSumOfValidWeights;                                 // Weighted contribution of this ray
                    real tSensitivityFactor = tFacetMeasures( iF ) * tWeightTerm * tAgglomerationTerm;               // Full sensitivity factor for this ray

                    // Origin facets vertices
                    for ( uint iV : tFacetVertices )
                    {
                        Matrix< DDRMat > tRotation = -tDiameterData.get_ray_rotation_matrix( iR );

                        // Raycast distance is a function of the origin and the normal, which are both influenced by these vertices
                        Matrix< DDRMat > tCentroidSens = mSurfaceMesh.compute_draycast_dorigin( tCentroids.get_column( iF ), tDiameterData.get_ray_direction( iF, iR ), tRaysOnFacet( iR ).first )
                                                       * mSurfaceMesh.compute_dfacet_centroid_dvertex( iF, iV, true );
                        Matrix< DDRMat > tNormalSens   = mSurfaceMesh.compute_draycast_ddirection( tCentroids.get_column( iF ), tDiameterData.get_ray_direction( iF, iR ), tRaysOnFacet( iR ).first )
                                                       * tRotation * mSurfaceMesh.compute_dfacet_normal_dvertex( iF, iV, true );

                        // Weight sensitivity is a function of the normal only for these vertices (rotation since the ray direction is the rotated origin facet normal)
                        Matrix< DDRMat > tWeightSens = tFacetMeasures( iF ) * tAlphaPower * std::pow( tAngle, tAlphaPower - 1 ) * mSurfaceMesh.compute_ddot_absolute( tDiameterData.get_ray_direction( iF, iR ), tIntersectedFacetNormal )
                                                     * tRotation * mSurfaceMesh.compute_dfacet_normal_dvertex( iF, iV, true ) / tSumOfValidWeights * tDiameterData.get_ray_weight( iR );

                        // Full sensitivity for this vertex
                        Matrix< DDRMat > tVertexSensitivity = tSensitivityFactor * ( tCentroidSens + tNormalSens ) + tWeightSens * ( tAgglomeratedDistance - tShapeDiameterValues( iF ) );

                        // Accumulate the sensitivity for this vertex
                        mdShapeDiameterdVertex.set_row( iV, mdShapeDiameterdVertex.get_row( iV ) + tVertexSensitivity );
                    }

                    // Intersected facet vertices
                    // Get the vertex indices for the intersected facet, these vertices need their sensitivities updated
                    const Vector< moris_index >& tIntersectedFacetVertices = mSurfaceMesh.get_facets_vertex_indices( tRaysOnFacet( iR ).first );
                    Matrix< DDRMat >             tVertexSens               = mSurfaceMesh.compute_draycast_dvertices( tCentroids.get_column( iF ), tDiameterData.get_ray_direction( iF, iR ), tRaysOnFacet( iR ).first );
                    for ( uint iV = 0; iV < tIntersectedFacetVertices.size(); iV++ )
                    {
                        // Weight sensitivity is a function of the normal only for these vertices (no rotation here since we are directly measuring against the intersected normal)
                        Matrix< DDRMat > tWeightSens = tFacetMeasures( iF ) * tAlphaPower * std::pow( tAngle, tAlphaPower - 1 ) * mSurfaceMesh.compute_ddot_absolute( tIntersectedFacetNormal, tDiameterData.get_ray_direction( iF, iR ) )
                                                     * mSurfaceMesh.compute_dfacet_normal_dvertex( tRaysOnFacet( iR ).first, tIntersectedFacetVertices( iV ), true ) / tSumOfValidWeights * tDiameterData.get_ray_weight( iR );

                        Matrix< DDRMat > tVertexSensitivity = tSensitivityFactor * tVertexSens.get_row( iV ) + tWeightSens * ( tAgglomeratedDistance - tShapeDiameterValues( iF ) );

                        // Raycast distance is a function of the facets vertices explicitly for these vertices. Computed for all vertices at once above
                        mdShapeDiameterdVertex.set_row( tIntersectedFacetVertices( iV ), mdShapeDiameterdVertex.get_row( tIntersectedFacetVertices( iV ) ) + tVertexSensitivity );
                    }
                }
            }
        }

        return tShapeDiameterValues;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real > Surface_Mesh_QI::compute_inscribed_circle_shape_diameter( const Agglomeration_Parameters& aAgglom, real aAngle, real aRelativeChord, uint aNumCircles )
    {
        MORIS_ERROR( mSurfaceMesh.get_spatial_dimension() == 2, "Inscribed circle shape diameter is only implemented for 2D surface meshes." );

        // Initialize output vector
        Vector< real > tAverageVertexDiameter( mSurfaceMesh.get_number_of_vertices(), std::numeric_limits< real >::quiet_NaN() );
        Vector< real > tFacetShapeDiameterValues( mSurfaceMesh.get_number_of_facets(), std::numeric_limits< real >::quiet_NaN() );

        // Get some quantities of the surface mesh, used for sensitivities and other checks
        Matrix< DDRMat > tVertexNormals = 0.0 - mSurfaceMesh.compute_vertex_normals();
        Vector< real >   tFacetMeasures = mSurfaceMesh.compute_facet_measure();

        // Minimum dot product of chord and normal for an inscribed circle to be valid (360 because we want the half cone angle)
        real tMinDot = std::cos( aAngle * M_PI / 360.0 );

        // Reset the derivative
        mdShapeDiameterdVertex.set_size( mSurfaceMesh.get_number_of_vertices(), mSurfaceMesh.get_spatial_dimension(), 0.0 );

        // Loop over vertices for circles to originate from
        for ( uint iVO = 0; iVO < mSurfaceMesh.get_number_of_vertices(); iVO++ )
        {
            // Create vector for distance values for this vertex
            Vector< real > tCircleDistances( aNumCircles, MORIS_REAL_MAX );    // To compute average of closest aNumDist distances
            Vector< uint > tHitVertices( aNumCircles, MORIS_INDEX_MAX );       // needed for sensitivities
            uint           tNumValidCircles = 0;                               /// In case there are not aNumDist valid distances within the cone

            // Get the facets connected to this vertex, needed for sensitivities
            Vector< moris_index > tFacetsConnectedToOriginVertex = mSurfaceMesh.get_vertexs_facet_indices( iVO );
            Matrix< DDRMat >      tVertexCoordinates             = mSurfaceMesh.get_vertex_coordinates( iVO );

            // Loop over vertices to check inscribed circle against
            for ( uint iVC = 0; iVC < mSurfaceMesh.get_number_of_vertices(); iVC++ )
            {
                // Get the vector from the origin vertex to the check vertex
                Matrix< DDRMat > tChord     = mSurfaceMesh.get_vertex_coordinates( iVC ) - mSurfaceMesh.get_vertex_coordinates( iVO );
                real             tChordNorm = norm( tChord );

                // Check that the chord is within the angle threshold of the vertex normal and that we aren't checking the same vertex
                if ( dot( tChord, tVertexNormals.get_column( iVO ) ) > tMinDot * tChordNorm )
                {
                    // Compute the radius of the inscribed circle between these two points
                    real tDiameter = dot( tChord, tChord ) / dot( tChord, tVertexNormals.get_column( iVO ) );

                    if ( tDiameter > 0 && tChordNorm / tDiameter > aRelativeChord )
                    {
                        // Check if the valid diameter should replace a minimum
                        for ( uint iMin = 0; iMin < aNumCircles; iMin++ )
                        {
                            if ( tDiameter < tCircleDistances( iMin ) )
                            {
                                // Shift the larger values down the list
                                for ( uint iMaxBack = aNumCircles - 1; iMaxBack > iMin; iMaxBack-- )
                                {
                                    tCircleDistances( iMaxBack ) = tCircleDistances( iMaxBack - 1 );
                                    tHitVertices( iMaxBack )     = tHitVertices( iMaxBack - 1 );
                                }

                                // Insert the new value
                                tCircleDistances( iMin ) = tDiameter;
                                tHitVertices( iMin )     = iVC;
                                tNumValidCircles++;
                                break;
                            }
                        }    // end loop over minima
                    }    // end check for positive diameter
                }    // end check for valid inscribed circle
            }    // end loop over check vertices

            tNumValidCircles = std::min( tNumValidCircles, aNumCircles );

            // Loop over valid circles, take agglomerated mean, and compute sensitivities
            real tMeanViolation = 0.0;
            for ( uint iC = 0; iC < tNumValidCircles; iC++ )
            {
                tMeanViolation += this->max_clip( tCircleDistances( iC ), aAgglom, tVertexCoordinates, iVO );

                // Sensitivity calculations for this circle, each circle has sensitivity to the origin vertex and the hit vertex
                uint iVC = tHitVertices( iC );

                Matrix< DDRMat > tChord = mSurfaceMesh.get_vertex_coordinates( iVC ) - mSurfaceMesh.get_vertex_coordinates( iVO );

                // Agglomeration derivative and other scaling terms
                real tSensitivityFactor = 0.5 * this->dmax_clip( tCircleDistances( iC ), aAgglom, tVertexCoordinates, iVO ) / (real)tNumValidCircles;

                // Derivative of radius wrt chord is given by the quotient rule, with the numerator and denominator as follows
                real tRadiusNumerator   = dot( tChord, tChord );
                real tRadiusDenominator = dot( tChord, tVertexNormals.get_column( iVO ) );

                Matrix< DDRMat > tDradiusDchord  = trans( ( 2.0 * tChord * tRadiusDenominator - tRadiusNumerator * tVertexNormals.get_column( iVO ) ) / ( tRadiusDenominator * tRadiusDenominator ) );
                Matrix< DDRMat > tDradiusDnormal = -tRadiusNumerator / ( tRadiusDenominator * tRadiusDenominator ) * trans( tChord );

                // Derivative of chord wrt origin vertex is -I, derivative of chord wrt hit vertex is I.
                // Loop over facets connected to the origin vertex and add sensitivities
                real tLengthSum = tFacetMeasures( tFacetsConnectedToOriginVertex( 0 ) ) + tFacetMeasures( tFacetsConnectedToOriginVertex( 1 ) );
                for ( uint iF : tFacetsConnectedToOriginVertex )
                {
                    // Add chord sensitivity, chord is only dependent on my vertex
                    mdShapeDiameterdVertex.set_row( iVO, mdShapeDiameterdVertex.get_row( iVO ) - tSensitivityFactor * tFacetMeasures( iF ) * tDradiusDchord );

                    // Hit vertex sensitivity is only dependent on the chord. If the hit vertex is a neighbor and thus can influence the normal, this is accounted for above
                    mdShapeDiameterdVertex.set_row( iVC, mdShapeDiameterdVertex.get_row( iVC ) + tSensitivityFactor * tFacetMeasures( iF ) * tDradiusDchord );

                    // Loop over vertices on this facet to add normal sensitivity
                    for ( uint iV : mSurfaceMesh.get_facets_vertex_indices( iF ) )
                    {
                        mdShapeDiameterdVertex.set_row( iV, mdShapeDiameterdVertex.get_row( iV ) - tSensitivityFactor * tLengthSum * tDradiusDnormal * mSurfaceMesh.compute_dvertex_normal_dvertex( iVO, iV ) );
                    }
                }
            }

            // Set the shape diameter value for this vertex
            tAverageVertexDiameter( iVO ) = tMeanViolation / std::max( 1.0, (real)tNumValidCircles );
        }

        // Compute the value on the facets by taking the mean of the vertices and return
        return mSurfaceMesh.average_vertex_quantity_over_facets( tAverageVertexDiameter );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real > Surface_Mesh_QI::compute_shortest_distance_shape_diameter( const Agglomeration_Parameters& aAgglom, real aAngle, uint aNumDist )
    {
        MORIS_ERROR( mSurfaceMesh.get_spatial_dimension() == 2, "Closest point diameter is only implemented for 2D surface meshes." );

        // Initialize output vector
        Vector< real > tAverageVertexDiameter( mSurfaceMesh.get_number_of_vertices(), std::numeric_limits< real >::quiet_NaN() );
        Vector< real > tFacetShapeDiameterValues( mSurfaceMesh.get_number_of_facets(), std::numeric_limits< real >::quiet_NaN() );

        // Get some quantities of the surface mesh, used for sensitivities and other checks
        Matrix< DDRMat > tVertexNormals = 0.0 - mSurfaceMesh.compute_vertex_normals();     // Inward normals with positive zero components
        Matrix< DDRMat > tCenters       = mSurfaceMesh.compute_facet_centroids();     // Used for bounding circle checks
        Vector< real >   tFacetMeasures = mSurfaceMesh.compute_facet_measure();       // Used for sensitivityes and bounding circle checks
        Vector< real >   tMeasure2      = Vector< real >( tFacetMeasures.size() );    // Used for bounding circle
        std::transform( tFacetMeasures.cbegin(), tFacetMeasures.cend(), tMeasure2.begin(), []( real aVal ) { return aVal * aVal; } );

        // Minimum dot product of chord and normal for an inscribed circle to be valid (360 because we want the half cone angle)
        real tMinDot = std::cos( aAngle * M_PI / 360.0 );

        // Reset the derivative
        mdShapeDiameterdVertex.set_size( mSurfaceMesh.get_number_of_vertices(), mSurfaceMesh.get_spatial_dimension(), 0.0 );

        // Loop over vertices for distance to originate from
        for ( uint iV = 0; iV < mSurfaceMesh.get_number_of_vertices(); iV++ )
        {
            Matrix< DDRMat > tMyCoords = mSurfaceMesh.get_vertex_coordinates( iV );

            // Create vector for distance values for this vertex
            Vector< real >             tDistances( aNumDist, MORIS_REAL_MAX );                                                                   // To compute average of closest aNumDist distances
            Vector< uint >             tHitFacets( aNumDist, MORIS_INDEX_MAX );                                                                  // needed for sensitivities
            Vector< real >             tTMins( aNumDist, MORIS_REAL_MAX );                                                                       // Stored local distance along hit vertex, used for sensitivities
            Vector< Matrix< DDRMat > > tClosestDirs( aNumDist, Matrix< DDRMat >( mSurfaceMesh.get_spatial_dimension(), 1, MORIS_REAL_MAX ) );    // Stored direction from vertex to closest point on facet, used for sensitivities
            uint                       tNumValidDistances = 0;                                                                                   // In case there are not aNumDist valid distances within the cone
            real                       tBottomDist        = tDistances( aNumDist - 1 );                                                          // Initialize upper bound on minima

            // Get the facets connected to this vertex, needed for sensitivities
            Vector< moris_index > tFacetsConnectedToOriginVertex = mSurfaceMesh.get_vertexs_facet_indices( iV );
            Matrix< DDRMat >      tVertexCoordinates             = mSurfaceMesh.get_vertex_coordinates( iV );

            // Loop over facets to get distance to
            for ( uint iF = 0; iF < mSurfaceMesh.get_number_of_facets(); iF++ )
            {
                // Get the vertices of this facet and check that my vertex is not part of it, as the distance would be in plane with the facet
                const Vector< moris_index >& tVerticesOnFacet = mSurfaceMesh.get_facets_vertex_indices( iF );
                if ( std::any_of( tVerticesOnFacet.cbegin(), tVerticesOnFacet.cend(), [ iV ]( moris_index aVertex ) { return aVertex == (moris_index)iV; } ) )
                {
                    continue;
                }

                // Compute bounding circle of this facet and see if its possible that this radius would be a minima
                real tDistCenter2 = dot( tCenters.get_column( iF ) - tMyCoords, tCenters.get_column( iF ) - tMyCoords );
                real tBound       = tBottomDist + 0.5 * tFacetMeasures( iF );
                if ( tDistCenter2 > tBound * tBound )
                {
                    continue;
                }

                // Compute edge vector
                Matrix< DDRMat > tV1   = mSurfaceMesh.get_vertex_coordinates( tVerticesOnFacet( 0 ) );
                Matrix< DDRMat > tEdge = mSurfaceMesh.get_vertex_coordinates( tVerticesOnFacet( 1 ) ) - tV1;

                // Compute where along this facet segment the closest point is and the distance to that point
                real tT                    = dot( tMyCoords - mSurfaceMesh.get_vertex_coordinates( tVerticesOnFacet( 0 ) ), tEdge ) / tMeasure2( iF );
                tT                         = std::max( 0.0, std::min( 1.0, tT ) );    // Clamp to facet segment
                Matrix< DDRMat > tQ        = tV1 + tT * tEdge;                        // Closest point on facet to my vertex
                Matrix< DDRMat > tDir      = tQ - tMyCoords;                          // Vector from my vertex to closest point
                real             tDistance = norm( tDir );                            // Squared distance from my vertex to closest point on facet

                // Check that the distance is within the cone and in the direction of the normal
                if ( tDistance > 0 and dot( tDir, tVertexNormals.get_column( iV ) ) > tMinDot * tDistance )
                {
                    // Check if the valid diameter should replace a minimum
                    for ( uint iMin = 0; iMin < aNumDist; iMin++ )
                    {
                        if ( tDistance < tDistances( iMin ) )
                        {
                            // Shift the larger values down the list
                            for ( uint iMaxBack = aNumDist - 1; iMaxBack > iMin; iMaxBack-- )
                            {
                                tDistances( iMaxBack )   = tDistances( iMaxBack - 1 );
                                tHitFacets( iMaxBack )   = tHitFacets( iMaxBack - 1 );
                                tTMins( iMaxBack )       = tTMins( iMaxBack - 1 );
                                tClosestDirs( iMaxBack ) = tClosestDirs( iMaxBack - 1 );
                            }

                            // Insert the new value
                            tDistances( iMin )   = tDistance;
                            tHitFacets( iMin )   = iF;
                            tTMins( iMin )       = tT;
                            tClosestDirs( iMin ) = tDir / tDistance;    // Store the normalized direction for sensitivities
                            tNumValidDistances++;
                            tBottomDist = tDistances( aNumDist - 1 );    // Update the upper bound on minima
                            break;
                        }
                    }    // end loop over minima
                }    // end check for valid distance
            }    // end loop over facets

            tNumValidDistances = std::min( tNumValidDistances, aNumDist );

            // Loop over valid circles, take agglomerated mean, and compute sensitivities
            real tMeanViolation = 0.0;
            for ( uint iD = 0; iD < tNumValidDistances; iD++ )
            {
                // Get the vertices of the facet that we hit
                const Vector< moris_index >& tVerticesOnFacet = mSurfaceMesh.get_facets_vertex_indices( tHitFacets( iD ) );

                tMeanViolation += this->max_clip( tDistances( iD ), aAgglom, tVertexCoordinates, iV );

                // Agglomeration derivative and other scaling terms
                real tSensitivityFactor = 0.5 * this->dmax_clip( tDistances( iD ), aAgglom, tVertexCoordinates, iV ) / (real)tNumValidDistances;

                // Derivative of distance wrt my vertex is given by the normalized direction vector
                real tLengthSum = tFacetMeasures( tFacetsConnectedToOriginVertex( 0 ) ) + tFacetMeasures( tFacetsConnectedToOriginVertex( 1 ) );
                mdShapeDiameterdVertex.set_row( iV, mdShapeDiameterdVertex.get_row( iV ) - tSensitivityFactor * tLengthSum * trans( tClosestDirs( iD ) ) );

                // Derivative of distance wrt to hit facet is also only dependent on the normalized direction vector and the tmin value
                for ( uint iVC = 0; iVC < tVerticesOnFacet.size(); iVC++ )
                {
                    real tVertexFactor;
                    if ( iVC == 0 )
                    {
                        tVertexFactor = ( 1.0 - tTMins( iD ) );
                    }
                    else
                    {
                        tVertexFactor = tTMins( iD );
                    }

                    mdShapeDiameterdVertex.set_row( tVerticesOnFacet( iVC ), mdShapeDiameterdVertex.get_row( tVerticesOnFacet( iVC ) ) + tSensitivityFactor * tVertexFactor * tLengthSum * trans( tClosestDirs( iD ) ) );
                }
            }    // end loop over sensitivities and valid distances

            // Set the shape diameter value for this vertex
            tAverageVertexDiameter( iV ) = tMeanViolation / std::max( 1.0, (real)tNumValidDistances );

        }    // end loop over compute vertex distances

        // Compute the value on the facets by taking the mean of the vertices and return
        return mSurfaceMesh.average_vertex_quantity_over_facets( tAverageVertexDiameter );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::integrate_shape_diameter_over_surface()
    {
        // Compute the area of every facet in the surface mesh
        Vector< real > tFacetMeasure = mSurfaceMesh.compute_facet_measure();

        // Need to integrate the shape diameter over the surface mesh, loop over facets
        mIntegratedShapeDiameter = 0.0;
        for ( uint iF = 0; iF < mSurfaceMesh.get_number_of_facets(); iF++ )
        {
            mIntegratedShapeDiameter += mShapeDiameters( iF ) * tFacetMeasure( iF );
        }

        // Store the agglomeration before we take the power(need it for sensitivities later)
        return mIntegratedShapeDiameter;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_global_shape_diameter_raycast(
            const Agglomeration_Parameters& aAgglom,
            real                            aConeAngle,
            uint                            aNumPolarRays,
            uint                            aNumAzimuthRays )
    {
        mShapeDiameters = this->compute_raycast_shape_diameter( aAgglom, aConeAngle, aNumPolarRays, aNumAzimuthRays );
        return this->integrate_shape_diameter_over_surface();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_global_shape_diameter_inscribed_circle(
            const Agglomeration_Parameters& aAgglom,
            real                            aAngle,
            real                            aRelativeChord,
            uint                            aNumCircles )
    {
        mShapeDiameters = this->compute_inscribed_circle_shape_diameter( aAgglom, aAngle, aRelativeChord, aNumCircles );
        return this->integrate_shape_diameter_over_surface();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_QI::compute_global_shape_diameter_shortest_distance(
            const Agglomeration_Parameters& aAgglom,
            real                            aAngle,
            uint                            aNumDistances )
    {
        mShapeDiameters = this->compute_shortest_distance_shape_diameter( aAgglom, aAngle, aNumDistances );
        return this->integrate_shape_diameter_over_surface();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_QI::compute_ddiameter_dvertex( const uint aVertexIndex ) const
    {
        // Get the facets connected to this vertex
        const Vector< moris_index >& tConnectedFacets = mSurfaceMesh.get_vertexs_facet_indices( aVertexIndex );

        Matrix< DDRMat > tIntegratedSensitivity( 1, mSurfaceMesh.get_spatial_dimension(), 0.0 );

        // Loop over connected facets to compute their contribution to the sensitivity as their area is changing
        for ( moris_index iF : tConnectedFacets )
        {
            Matrix< DDRMat > tdMeasuredVertex = mSurfaceMesh.compute_dfacet_measure_dvertex( iF, aVertexIndex, true );
            tIntegratedSensitivity += tdMeasuredVertex * mShapeDiameters( iF );
        }

        // Add shape diameter sensitivity for this vertex (precomputed during shape diameter forward computation, includes summed sensitivities from all facets)
        tIntegratedSensitivity += mdShapeDiameterdVertex.get_row( aVertexIndex );

        return tIntegratedSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_QI::compute_dvolume_dvertex( uint aVertexIndex ) const
    {
        uint tDim = mSurfaceMesh.get_spatial_dimension();

        Matrix< DDRMat > tDVolumeDVertex( 1, tDim, MORIS_REAL_MAX );

        if ( mSurfaceMesh.get_spatial_dimension() == 2 )
        {
            // Get the vertices connected to this vertex
            const Vector< moris_index >& tNeighbors = mSurfaceMesh.get_vertex_neighbors( aVertexIndex );

            moris_index tLowNeighbor;
            moris_index tHighNeighbor;
            if ( aVertexIndex > 1 and std::any_of( tNeighbors.begin(), tNeighbors.end(), []( moris_index i ) { return i == 0; } ) )
            {
                tLowNeighbor  = tNeighbors( 0 ) == 0 ? tNeighbors( 1 ) : tNeighbors( 0 );
                tHighNeighbor = 0;
            }
            else
            {    // Determine the order of the vertices
                tLowNeighbor  = tNeighbors( 0 ) < tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );
                tHighNeighbor = tNeighbors( 0 ) > tNeighbors( 1 ) ? tNeighbors( 0 ) : tNeighbors( 1 );
            }

            // Get vertex coordinates
            const Matrix< DDRMat > tLowCoords  = mSurfaceMesh.get_vertex_coordinates( tLowNeighbor );
            const Matrix< DDRMat > tHighCoords = mSurfaceMesh.get_vertex_coordinates( tHighNeighbor );

            tDVolumeDVertex( 0 ) = 0.5 * ( tHighCoords( 1 ) - tLowCoords( 1 ) );
            tDVolumeDVertex( 1 ) = 0.5 * ( tLowCoords( 0 ) - tHighCoords( 0 ) );
        }
        else
        {
            MORIS_ERROR( false, "Surface Mesh dVolume/dVertex computation only implemented for 2D (lines) meshes" );
        }

        return tDVolumeDVertex;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Surface_Mesh_QI::get_facet_shape_diameter_sensitivities() const
    {
        return mdShapeDiameterdVertex;
    }

    //--------------------------------------------------------------------------------------------------------------

    const real Surface_Mesh_QI::get_facet_shape_diameter( uint aFacetIndex ) const
    {
        return mShapeDiameters( aFacetIndex );
    }

}    // namespace moris::mtk
