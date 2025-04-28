/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh_Geometry.cpp
 *
 */

#include <string>

#include "cl_Tracer.hpp"

#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Floating_Node_Surface_Mesh.hpp"
#include "cl_GEN_Parent_Node.hpp"

#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Stored_Field.hpp"
#include "cl_GEN_Constant_Field.hpp"

#include "fn_cross.hpp"
#include "fn_eye.hpp"
#include "fn_trans.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "fn_MTK_Load_External_Surface_Mesh.hpp"
#include "fn_MTK_mesh_flood_fill.hpp"
#include "fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_SOL_Dist_Map.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Parameters::Surface_Mesh_Parameters( const Parameter_List& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mOffsets( aParameterList.get< Vector< real > >( "offset" ) )
            , mScale( aParameterList.get< Vector< real > >( "scale" ) )
            , mFilePath( aParameterList.get< std::string >( "file_path" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
            , mDiscretizationFactorFunctionName( aParameterList.get< std::string >( "discretization_factor_function_name" ) )
            , mAnalyticADVFunctionName( aParameterList.get< std::string >( "field_function_name" ) )
            , mAnalyticADVSensitivityFunctionName( aParameterList.get< std::string >( "sensitivity_function_name" ) )
            , mOutputFileName( aParameterList.get< std::string >( "output_file_name" ) )
    {
        MORIS_ASSERT( (uint)( mDiscretizationIndex < -1 + not mAnalyticADVFunctionName.empty() ) < 2, "Both a discretization index and an analytical function are provided. Pick only one or neither!" );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry(
            Surface_Mesh_Parameters       aParameters,
            Node_Manager&                 aNodeManager,
            const Vector< ADV >&          aADVs,
            ADV_Manager&                  aADVManager,
            std::shared_ptr< Library_IO > aLibrary )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Surface_Mesh( mtk::load_vertices_from_object_file( aParameters.mFilePath, aParameters.mOffsets, aParameters.mScale ),
                      mtk::load_facets_from_object_file( aParameters.mFilePath ),
                      aParameters.mIntersectionTolerance )
            , mParameters( aParameters )
            , mNodeManager( &aNodeManager )
            , mADVHandler( aADVs )
            , mPerturbationFields( 0 )
            , mOriginalVertexBases( 0, 0 )
            , mOriginalVertexBackgroundElements( 0 )
            , mVertexParametricCoordinates( Surface_Mesh::get_spatial_dimension(), Surface_Mesh::get_number_of_vertices() )
    {
        uint tDim = Surface_Mesh::get_spatial_dimension();

        // parse the file path and extract the file name
        mName = aParameters.mFilePath.substr( aParameters.mFilePath.find_last_of( "/" ) + 1,
                aParameters.mFilePath.find_last_of( "." ) - aParameters.mFilePath.find_last_of( "/" ) - 1 );

        // If this surface mesh is being optimized via B-spline fields, construct seeding fields and check for a scaling function
        if ( this->intended_discretization() )
        {
            // STEP 1: Determine which facet vertices are fixed
            if ( not mParameters.mDiscretizationFactorFunctionName.empty() )
            {
                // Get pointer to function
                get_discretization_scaling_user_defined = aLibrary->load_function< Discretization_Factor_Function >( mParameters.mDiscretizationFactorFunctionName );
            }

            // STEP 2: Construct perturbation fields
            mPerturbationFields.resize( tDim );

            // Build a dummy field for every spatial dimension
            for ( uint iFieldIndex = 0; iFieldIndex < tDim; iFieldIndex++ )
            {
                Vector< uint > tADVIndices;
                Vector< real > tConstants( 1 );
                Vector< uint > tFieldVariableIndices;

                // Build field
                mPerturbationFields( iFieldIndex ) = std::make_shared< Constant_Field >(
                        aADVManager.mADVs,
                        tFieldVariableIndices,
                        tADVIndices,
                        tConstants,
                        mName + "_PERT_" + std::to_string( iFieldIndex ) );
            }
        }
        // Otherwise, check if this surface mesh is being optimized via an analytical function
        else if ( not mParameters.mAnalyticADVFunctionName.empty() )
        {
            // Get pointer to function if so
            get_vertex_adv_dependency_user_defined = aLibrary->load_function< Perturbation_Function >( mParameters.mAnalyticADVFunctionName );
            get_dvertex_dadv_user_defined          = aLibrary->load_function< Sensitivity_Function >( mParameters.mAnalyticADVSensitivityFunctionName );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = &aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Surface_Mesh_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // check if this node has been determined previously
        auto             tNodeIt        = mNodeMeshRegions.find( aNodeIndex );
        bool             tRegionUnknown = tNodeIt == mNodeMeshRegions.end();
        mtk::Mesh_Region tRegion        = tRegionUnknown ? mNodeMeshRegions[ aNodeIndex ] = this->get_region_from_raycast( aNodeCoordinates ) : tNodeIt->second;

        switch ( tRegion )
        {
            case mtk::Mesh_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
            }
            case mtk::Mesh_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
            }
            case mtk::Mesh_Region::INTERFACE:
            {
                return Geometric_Region::INTERFACE;
            }
            default:
            {
                MORIS_ERROR( false, "Unexpected sdf::Object_Region of %d returned from raycast.", tRegion );
                return Geometric_Region::INTERFACE;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region
    Surface_Mesh_Geometry::disambiguate_geometric_region( const Matrix< DDRMat >& aNodeCoordinates )
    {
        // Have to raycast to determine the region
        mtk::Mesh_Region tRegion = this->get_region_from_raycast( aNodeCoordinates );

        switch ( tRegion )
        {
            case mtk::Mesh_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
            }
            case mtk::Mesh_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
            }
            case mtk::Mesh_Region::INTERFACE:
            {
                return Geometric_Region::INTERFACE;
            }
            default:
            {
                MORIS_ERROR( false, "Unexpected sdf::Object_Region of %d returned from raycast.", tRegion );
                return Geometric_Region::INTERFACE;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::has_surface_points( mtk::Cell* aCell )
    {
        // Check if Delaunay triangulation is turned on for this geometry
        if ( not mParameters.mDelaunay )
        {
            return false;
        }

        // Try and see if the pointer is in the stored list
        auto tIt = std::find( mCurrentVertexBackgroundElements.begin(), mCurrentVertexBackgroundElements.end(), aCell );

        while ( tIt != mCurrentVertexBackgroundElements.end() )
        {
            // Get the index of the vertex in this cell
            uint tVertexIndex = std::distance( mCurrentVertexBackgroundElements.begin(), tIt );

            // Get the basis of this vertex
            Matrix< DDRMat > tBasis = mCurrentVertexBases.get_column( tVertexIndex );

            // check if any of the bases for the vertex are zero. If so, the point is on the edge of the cell and we should skip it
            if ( std::any_of( tBasis.cbegin(), tBasis.cend(), [ this ]( real aValue ) { return std::abs( aValue ) < Surface_Mesh::mIntersectionTolerance; } ) )
            {
                // Check for another vertex in this element
                tIt = std::find( std::next( tIt ), mCurrentVertexBackgroundElements.end(), aCell );
            }
            else
            {
                // The surface point is valid
                return true;
            }
        }

        // No points inside the cell that are not on the cell edge
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::get_surface_points( mtk::Cell* aCell )
    {
        // Get the spatial dimension of the surface mesh
        uint tDim = this->get_spatial_dimension();

        // Get the number of vertices for this surface mesh
        uint tNumNodes = this->get_number_of_vertices();

        // Initialize matrix
        Matrix< DDRMat > tSurfacePoints( tNumNodes, tDim );

        // Initialize counter for number of vertices added
        uint iDelaunayPointIndex = 0;

        // Get the coordinates of the background nodes
        Matrix< DDRMat > tBackgroundCellCoords = aCell->get_vertex_coords();

        // Get the length of each dimension of the element
        Vector< real > tElementLengths( tDim );
        for ( uint iDim = 0; iDim < tDim; iDim++ )
        {
            tElementLengths( iDim ) = tBackgroundCellCoords.get_column( iDim ).max() - tBackgroundCellCoords.get_column( iDim ).min();
        }

        // Loop over all surface mesh vertices and check if they are in the background element
        // FIXME: brute force search could be optimized
        for ( uint iSurfaceMeshVertex = 0; iSurfaceMeshVertex < tNumNodes; iSurfaceMeshVertex++ )
        {
            if ( mCurrentVertexBackgroundElements( iSurfaceMeshVertex ) == aCell )
            {
                // Get the basis for this vertex
                Matrix< DDRMat > tBasis = mCurrentVertexBases.get_column( iSurfaceMeshVertex );

                // If the vertex is on the edge of the background element, skip it as this will cause hanging nodes
                if ( std::any_of( tBasis.begin(), tBasis.end(), [ this ]( real aVal ) { return std::abs( aVal ) < Surface_Mesh::mIntersectionTolerance; } ) )
                {
                    continue;
                }

                // convert to parametric coordinates
                Matrix< DDRMat > tParametricCoordinates = trans( this->get_vertex_coordinates( iSurfaceMeshVertex ) ) - tBackgroundCellCoords.get_row( 0 );
                for ( uint iDim = 0; iDim < tDim; iDim++ )
                {
                    tParametricCoordinates( iDim ) = 2.0 / tElementLengths( iDim ) * tParametricCoordinates( iDim ) - 1.0;
                }

                // Add the vertex to the matrix if lies in the cell and isnt on the edge
                tSurfacePoints.set_row( iDelaunayPointIndex++, mVertexParametricCoordinates.get_column( iSurfaceMeshVertex ) );
            }
        }

        // Trim output matrix
        tSurfacePoints.resize( iDelaunayPointIndex, tDim );

        return tSurfacePoints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Floating_Node* Surface_Mesh_Geometry::create_floating_node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Matrix< DDRMat >&           aParametricCoordinates,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
    {
        uint tNumVertices = this->get_number_of_vertices();
        uint tDims        = this->get_spatial_dimension();

        mtk::Interpolation_Function_Factory tInterpolationFactory;
        mtk::Interpolation_Function_Base*   tInterpolation = tInterpolationFactory.create_interpolation_function(
                aBackgroundGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                aBackgroundInterpolationOrder );

        // Perform interpolation using parametric coordinates
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( aParametricCoordinates, tBasis );

        // Size global coordinates based on first locator
        Matrix< DDRMat > tGlobalCoordinate = Matrix< DDRMat >( 1, this->get_spatial_dimension(), 0.0 );
        delete tInterpolation;

        // Add contributions from all locators
        for ( uint iNode = 0; iNode < aBackgroundNodes.size(); iNode++ )
        {
            tGlobalCoordinate += aBackgroundNodes( iNode )->get_global_coordinates() * tBasis( iNode );
        }

        Matrix< DDRMat > tVertexGlobalCoords = Surface_Mesh::get_all_vertex_coordinates();

        // Determine the parent vertex that lies in the background cell
        uint tParentVertex = MORIS_UINT_MAX;
        for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
        {
            bool tVertexFound = true;
            for ( uint iDim = 0; iDim < tDims; iDim++ )
            {
                if ( std::abs( tVertexGlobalCoords( iDim, iVertex ) - tGlobalCoordinate( iDim ) ) > Surface_Mesh::mIntersectionTolerance )
                {
                    tVertexFound = false;
                    break;
                }
            }

            if ( tVertexFound )
            {
                tParentVertex = iVertex;
                break;
            }
        }

        if ( tParentVertex == MORIS_UINT_MAX )
        {
            this->write_to_file( "failed.obj" );
            MORIS_ERROR( false, "Floating node %d does not lie on a vertex of surface mesh \"%s\"", aNodeIndex, this->get_name().c_str() );
        }

        // Create surface mesh floating node
        return new Floating_Node_Surface_Mesh(
                aNodeIndex,
                aBackgroundNodes,
                aParametricCoordinates,
                tParentVertex,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Surface_Mesh_Geometry::create_intersection_node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
    {
        // Determine the local coordinate of the intersection and the facet that intersects the parent edge
        std::pair< uint, real > tIntersection = this->compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode );

        if ( tIntersection.second > 1.0 + Surface_Mesh::mIntersectionTolerance or std::isnan( tIntersection.second ) )
        {
            std::cout << "First parent node index :" << aFirstParentNode.get_index() << std::endl;
            PRINT( aFirstParentNode.get_global_coordinates() );
            std::cout << "Second parent node index :" << aSecondParentNode.get_index() << std::endl;
            PRINT( aSecondParentNode.get_global_coordinates() );
            std::cout << "1st Region by raycasting" << get_region_from_raycast( aFirstParentNode.get_global_coordinates() ) << std::endl;
            std::cout << "2nd Region by raycasting" << get_region_from_raycast( aSecondParentNode.get_global_coordinates() ) << std::endl;
            std::cout << "1st region from stored data " << mNodeMeshRegions.at( aFirstParentNode.get_index() ) << std::endl;
            std::cout << "2nd region from stored data " << mNodeMeshRegions.at( aSecondParentNode.get_index() ) << std::endl;
            tIntersection = this->compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode );
            Surface_Mesh::write_to_file( "failed.obj" );
        }

        MORIS_ERROR( tIntersection.first != MORIS_UINT_MAX and ( tIntersection.second < ( 1.0 + Surface_Mesh::mIntersectionTolerance ) and tIntersection.second > ( -1.0 - Surface_Mesh::mIntersectionTolerance ) ),
                "Intersection node %d on surface mesh %s has local coordinate %f. Should be [-1, 1]",
                aNodeIndex,
                this->get_name().c_str(),
                tIntersection.second );

        // Create surface mesh intersection node
        return new Intersection_Node_Surface_Mesh(
                aNodeIndex,
                aBackgroundNodes,
                aFirstParentNode,
                aSecondParentNode,
                tIntersection,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::pair< uint, real > Surface_Mesh_Geometry::compute_intersection_local_coordinate(
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode )
    {
        // Get the parent node global coordinates. The origin of the ray will be the first parent node
        Matrix< DDRMat > tFirstParentNodeCoordinates  = aFirstParentNode.get_global_coordinates();
        Matrix< DDRMat > tSecondParentNodeCoordinates = aSecondParentNode.get_global_coordinates();

        // Need to cast along the edge to determine where the edge is intersected
        Matrix< DDRMat > tRayDirection = tSecondParentNodeCoordinates - tFirstParentNodeCoordinates;

        // Do a raycast to get the locations of the ray/facet intersections
        bool                     tWarning;
        mtk::Intersection_Vector tLocalCoordinate = this->cast_single_ray( tFirstParentNodeCoordinates, tRayDirection, tWarning );

        // no intersections detected or multiple along parent edge
        if ( tLocalCoordinate.size() == 0 )
        {
            // Check to see if either parent node is on the interface and we can return this
            // FIXME: the raycast should be able to detect this case, but ArborX may not reliably return the candidate facets in this case. Adjust bounding boxes somehow?
            if ( this->get_geometric_region( aFirstParentNode.get_index(), tFirstParentNodeCoordinates ) == Geometric_Region::INTERFACE )
            {
                return std::make_pair( 0, -1.0 );    // FIXME: since we don't know the facet that it intersects, we cant return it. This should get ignored by XTK though
            }
            else if ( this->get_geometric_region( aSecondParentNode.get_index(), tSecondParentNodeCoordinates ) == Geometric_Region::INTERFACE )
            {
                return std::make_pair( 0, 1.0 );    // FIXME: since we don't know the facet that it intersects, we cant return it. This should get ignored by XTK though
            }
        }

        // Process the raycast output and if no valid intersection was found, recast ray with looser tolerance
        return this->process_raycast_for_local_coordinate( Surface_Mesh::mIntersectionTolerance, tFirstParentNodeCoordinates, tRayDirection, tLocalCoordinate );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::pair< uint, real >
    Surface_Mesh_Geometry::process_raycast_for_local_coordinate(
            real                      aOriginalTolerance,
            Matrix< DDRMat >&         aOrigin,
            Matrix< DDRMat >&         aDirection,
            mtk::Intersection_Vector& aRaycastResult )
    {
        // Put the intersections in the local coordinate frame
        for ( uint iIntersection = 0; iIntersection < aRaycastResult.size(); iIntersection++ )
        {
            aRaycastResult( iIntersection ).second = 2.0 * aRaycastResult( iIntersection ).second - 1.0;
        }

        // Iterate through the intersections to find the first value that is within the exclusive range (-1, 1)
        for ( auto& tIntersection : aRaycastResult )
        {
            if ( std::abs( tIntersection.second ) < ( 1.0 + Surface_Mesh::mIntersectionTolerance ) )
            {
                // reset intersection tolerance
                Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;

                // snap if needed
                if ( std::abs( tIntersection.second - 1.0 ) < Surface_Mesh::mIntersectionTolerance )
                {
                    tIntersection.second = 1.0;
                }
                else if ( std::abs( tIntersection.second + 1.0 ) < Surface_Mesh::mIntersectionTolerance )
                {
                    tIntersection.second = -1.0;
                }
                return tIntersection;
            }
        }

        // If no intersection was found, loosen the intersection tolerance and try again
        Surface_Mesh::mIntersectionTolerance *= 10.0;

        // Recast the ray
        bool                     tWarning;
        mtk::Intersection_Vector tNewIntersections = this->cast_single_ray( aOrigin, aDirection, tWarning );

        // Process new intersections
        return this->process_raycast_for_local_coordinate( aOriginalTolerance, aOrigin, aDirection, tNewIntersections );
    }

    //--------------------------------------------------------------------------------------------------------------


#if MORIS_HAVE_ARBORX
    void Surface_Mesh_Geometry::flood_fill_mesh_regions()
    {
        Tracer tTracer( "GEN", "Surface_Mesh_Geometry", "Flood fill mesh nodes" );

        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // ----------------------------------------------------------------------------------------------
        // Initial setup for flood fill
        // ----------------------------------------------------------------------------------------------

        uint tNumNodes = mMesh->get_num_nodes();

        // Matrix< DDRMat >      tNodeCoordinates( Surface_Mesh::get_spatial_dimension(), tNumNodes );
        Vector< moris_index > tCellIndices( tNumNodes );

        // Build the element connectivity from the interpolation mesh
        Matrix< IndexMat > tNodeConnectivity;
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            // Get the connectivity for this node
            Matrix< IndexMat > tSingleNodeConnectivity = mMesh->get_entity_connected_to_entity_loc_inds( iNode, mtk::EntityRank::NODE, mtk::EntityRank::NODE );

            // If this node has more neighbors than any others, resize the connectivity matrix
            if ( tSingleNodeConnectivity.length() > tNodeConnectivity.n_cols() )
            {
                // Get the previous length
                uint tPreviousLength = tNodeConnectivity.n_cols();

                // Resize the connectivity matrix for the new max number of neighbors
                tNodeConnectivity.resize( tNumNodes, tSingleNodeConnectivity.length() );

                // Backfill the previous neighbors with gNoIndex NOTE: This will only work if the Matrix is column major
                std::fill( tNodeConnectivity.begin() + tNumNodes * tPreviousLength, tNodeConnectivity.end(), gNoIndex );
            }

            // Set the connectivity for this node
            for ( uint iConnection = 0; iConnection < tNodeConnectivity.n_cols(); iConnection++ )
            {
                tNodeConnectivity( iNode, iConnection ) = iConnection < tSingleNodeConnectivity.length() ? tSingleNodeConnectivity( iConnection, 0 ) : gNoIndex;
            }
        }

        // Dummy phase index for the entire mesh
        Vector< moris_index > tNodesToInclude( tNumNodes, 1 );

        // Active nodes (all nodes)
        Vector< moris_index > tNodes( tNumNodes );
        std::iota( tNodes.begin(), tNodes.end(), 0 );

        // Return variable for max number of subphases
        moris_index tMaxSubphase;

        // ----------------------------------------------------------------------------------------------
        // Determine which nodes do not lie in a surface mesh bounding box, and add them to the elements to include
        // ----------------------------------------------------------------------------------------------
        // Initialize to include all nodes
        Vector< moris_index > tPhase( tNumNodes, 1 );

        // Initialize the results and offsets
        Kokkos::View< ElementQueryResult*, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int*, MemorySpace >                tOffsets( "offsets", 0 );

        // Build the struct for the bounding boxes for the mesh
        QueryElements< MemorySpace > tQueryElements = this->construct_query_elements< MemorySpace >( tExecutionSpace );

        // Query the bounding volume hierarchy for the surface mesh
        mBVH.query( tExecutionSpace, tQueryElements, ElementIntersectionCallback< MemorySpace >{ tQueryElements }, tResults, tOffsets );

        for ( size_t iIntersection = 0; iIntersection < tResults.extent( 0 ); ++iIntersection )
        {
            // Get the element that was intersected
            moris_index iElement = tResults( iIntersection ).mElementIndex;

            Matrix< IndexMat > tNodesConnectedToElement = mMesh->get_nodes_connected_to_element_loc_inds( iElement );

            // Flag all of the nodes associated with the element
            for ( auto iNode : tNodesConnectedToElement )
            {
                tPhase( iNode ) = 0;
            }
        }

        // ----------------------------------------------------------------------------------------------
        // Run flood fill algorithm and get the subphases
        // ----------------------------------------------------------------------------------------------

        Matrix< IndexMat > tSubphases = mtk::flood_fill(
                tNodeConnectivity,
                tPhase,
                tNodes,
                tNodesToInclude,
                gNoIndex,
                tMaxSubphase,
                true );

        // ----------------------------------------------------------------------------------------------
        // Assign the correct geometric region to each subphase
        // ----------------------------------------------------------------------------------------------

        // For every subphase, determine the geometric region via a raycast (or leave as undefined if the node is in a bounding box)
        Vector< mtk::Mesh_Region > tSubPhaseMeshRegions( tMaxSubphase + 1, mtk::Mesh_Region::UNDEFINED );
        for ( int iSubphase = 0; iSubphase < tMaxSubphase + 1; iSubphase++ )
        {
            // Get the first index of this subphase for seeding
            moris_index iIndex = std::distance( tSubphases.cbegin(), std::find( tSubphases.cbegin(), tSubphases.cend(), iSubphase ) );

            tSubPhaseMeshRegions( iSubphase ) = tPhase( iIndex ) == 0 ? mtk::Mesh_Region::UNDEFINED : Surface_Mesh::get_region_from_raycast( mMesh->get_node_coordinate( iIndex ) );
        }

        // Loop through the nodes and assign their regions
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            if ( tSubPhaseMeshRegions( tSubphases( iNode, 0 ) ) != mtk::Mesh_Region::UNDEFINED )
            {
                mNodeMeshRegions[ iNode ] = tSubPhaseMeshRegions( tSubphases( iNode, 0 ) );
            }
        }
    }


#endif

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::raycast_remaining_unknown_nodes()
    {
        Tracer tTracer( "GEN", "Surface_Mesh_Geometry", "Raycast remaining unknown nodes" );

        // Get the number of nodes in the mesh and the spatial dimension
        uint tDims     = Surface_Mesh::get_spatial_dimension();
        uint tNumNodes = mMesh->get_num_nodes();

        // Initialize a vector to store unknown node indices and their coordinates
        Vector< uint > tUnknownNodes;
        tUnknownNodes.reserve( tNumNodes );    // Reserve space to avoid multiple allocations

        // Initialize a matrix to store the coordinates of unknown nodes
        Matrix< DDRMat > tUnknownNodeCoordinates( tDims, tNumNodes );

        // Fill the vector and matrix with unknown node indices and their coordinates
        uint tNumUnknownNodes = 0;
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            if ( mNodeMeshRegions.find( iNode ) == mNodeMeshRegions.end() )
            {
                tUnknownNodes.emplace_back( iNode );
                tUnknownNodeCoordinates.set_column( tNumUnknownNodes++, trans( mMesh->get_node_coordinate( iNode ) ) );
            }
        }
        tUnknownNodeCoordinates.resize( tDims, tNumUnknownNodes );    // trim the matrix to the correct size

        // Batch cast all of the points
        Vector< mtk::Mesh_Region > tUnknownNodeRegions = Surface_Mesh::batch_get_region_from_raycast( tUnknownNodeCoordinates );

        // Assign the regions to the unknown nodes
        for ( uint iNode = 0; iNode < tNumUnknownNodes; ++iNode )
        {
            mNodeMeshRegions[ tUnknownNodes( iNode ) ] = tUnknownNodeRegions( iNode );
        }
    }


    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        // Store this mesh
        mMesh = aInterpolationMesh;

        // update the perturbation fields with the new mesh
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->reset_nodal_data( aInterpolationMesh );
        }

        this->update_vertex_basis_data();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // since the shape is updating, clear any region information that was stored from the previous shape
        mNodeMeshRegions.clear();

        const uint tDims        = Surface_Mesh::get_spatial_dimension();
        uint       tNumVertices = Surface_Mesh::get_number_of_vertices();

        // Displacements of the owned vertices to communicate to other processors (first rows = coordinates, last row = owned flag)
        Matrix< DDRMat > tOwnedVertexDisplacements( tDims + 1, tNumVertices );

        // Check if this surface mesh is being optimized via B-spline fields
        if ( this->intended_discretization() )
        {
            for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
            {
                // Import advs to field
                mPerturbationFields( iFieldIndex )->import_advs( aOwnedADVs );
            }

            // Add this vertex's movement to the owned vertex coordinates
            for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
            {
                // Get the factor that scales this vertex's movement
                Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 ) : get_discretization_scaling_user_defined( Surface_Mesh::get_original_vertex_coordinates( iVertexIndex ) );

                for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
                {
                    // update the facet vertex if it can move and is owned by this processor
                    if ( this->facet_vertex_depends_on_advs( iVertexIndex ) and mOriginalVertexBackgroundElements( iVertexIndex )->get_owner() == par_rank() )
                    {
                        // Interpolate the bspline field value at the facet vertex location
                        real tInterpolatedPerturbation = this->interpolate_perturbation_from_background_element(
                                mOriginalVertexBackgroundElements( iVertexIndex ),
                                iFieldIndex,
                                iVertexIndex );

                        // build the matrix for new coordinates
                        tOwnedVertexDisplacements( tDims, iVertexIndex )       = 1.0;    // says that this vertex is owned bwwy this proc
                        tOwnedVertexDisplacements( iFieldIndex, iVertexIndex ) = tFactor( iFieldIndex ) * tInterpolatedPerturbation;
                    }
                }
            }
        }
        // The ADVs control some analytic function
        else if ( get_vertex_adv_dependency_user_defined != nullptr )
        {
            mADVHandler.import_advs( aOwnedADVs );

            for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
            {
                // TODO: Implement parallel updating of vertex coordinates
                // Compute the perturbation to this vertex from the user defined function
                Vector< real > tInterpolatedPerturbation = get_vertex_adv_dependency_user_defined( Surface_Mesh::get_original_vertex_coordinates( iVertexIndex ), mADVHandler.get_values() );

                // build the matrix for new coordinates
                for ( uint iFieldIndex = 0; iFieldIndex < tDims; iFieldIndex++ )
                {
                    tOwnedVertexDisplacements( tDims, iVertexIndex )       = 1.0;    // says that this vertex is owned by this proc
                    tOwnedVertexDisplacements( iFieldIndex, iVertexIndex ) = tInterpolatedPerturbation( iFieldIndex );
                }
            }
        }

        // Get the vertex coordinates from all processors and put in a vector of mats on base proc
        Vector< Matrix< DDRMat > > tAllVertexDisplacements;
        gatherv_mats( tOwnedVertexDisplacements, tAllVertexDisplacements );

        // Build matrix with all vertex coordinates on base proc
        Matrix< DDRMat > tCombinedVertexCoordinates( tDims + 1, Surface_Mesh::get_number_of_vertices() );
        if ( par_rank() == 0 )
        {
            for ( uint iProcessor = 0; iProcessor < tAllVertexDisplacements.size(); iProcessor++ )
            {
                for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
                {
                    // TODO: check if vertices are shared and if so that the coordinates are the same

                    // Check to see if the vertex was owned by proc iProcessor
                    if ( (uint)tAllVertexDisplacements( iProcessor )( tDims, iVertexIndex ) == 1 )
                    {
                        // If so, set the node coordinates in the combined matrix
                        tCombinedVertexCoordinates.set_column( iVertexIndex, tAllVertexDisplacements( iProcessor ).get_column( iVertexIndex ) );
                    }
                }
            }
        }

        // Give all the processors the new combined vertex coordinates assembled by base proc
        broadcast_mat( tCombinedVertexCoordinates );

        // Update the surface mesh vertex coordinates
        for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
        {
            // Check if this vertex was updated by any processor
            if ( (uint)tCombinedVertexCoordinates( tDims, iVertexIndex ) == 1 )
            {
                Surface_Mesh::set_vertex_displacement( iVertexIndex, tCombinedVertexCoordinates( { 0, tDims - 1 }, { iVertexIndex, iVertexIndex } ) );
            }
        }

        // Output the updated surface mesh to a file if needed
        if ( not mParameters.mOutputFileName.empty() )
        {
            // Initialize file extension as an obj file as default
            std::string tFileExt = ".obj";

            // check if there is a file extension provided
            if ( mParameters.mOutputFileName.find_last_of( "." ) != std::string::npos )
            {
                // if so, get the file extension
                tFileExt = mParameters.mOutputFileName.substr( mParameters.mOutputFileName.find_last_of( "." ), mParameters.mOutputFileName.length() );
            }

            // write the updated surface mesh to a file
            this->write_to_file( mParameters.mOutputFileName + "_proc_" + std::to_string( par_rank() ) + "_iter_" + std::to_string( gLogger.get_opt_iteration() ) + tFileExt );
        }

        // Update the facet's information based on the new vertex coordinates
        Surface_Mesh::initialize_facet_normals();

        // Determine new region information for the nodes
#if MORIS_HAVE_ARBORX
        Surface_Mesh::construct_bvh();
        this->flood_fill_mesh_regions();
#endif
        this->raycast_remaining_unknown_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::set_advs( sol::Dist_Vector* aADVs )
    {
        mADVHandler.set_advs( aADVs );

        // Have each field import the advs
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->set_advs( aADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    sint Surface_Mesh_Geometry::append_adv_info(
            mtk::Interpolation_Mesh* aMesh,
            Vector< sint >&          aOwnedADVIds,
            Matrix< IdMat >&         aOwnedijklIDs,
            sint                     aOffsetID,
            Vector< real >&          aLowerBounds,
            Vector< real >&          aUpperBounds,
            uint                     aFieldIndex )
    {
        // Get the original offset ID
        sint tOriginalOffsetID = aOffsetID;

        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // Append the ADV info for this field
            aOffsetID = Design::append_adv_info(
                    aMesh,
                    aOwnedADVIds,
                    aOwnedijklIDs,
                    aOffsetID,
                    aLowerBounds,
                    aUpperBounds,
                    iFieldIndex );
        }

        // reset the offset back to the offset for the first perturabtion field (mOffsetID was changed in the above loop)
        mOffsetID = tOriginalOffsetID;

        return aOffsetID;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::depends_on_advs() const
    {
        return this->intended_discretization() or get_vertex_adv_dependency_user_defined != nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::discretize(
            mtk::Mesh_Pair    aMeshPair,
            sol::Dist_Vector* aOwnedADVs )
    {
        MORIS_ASSERT( Design::mSharedADVIDs.size() == Surface_Mesh::get_spatial_dimension() or Design::mSharedADVIDs.size() == 0,
                "mSharedADVIDs should have as many entries as dimensions. Size = %ld",
                mSharedADVIDs.size() );

        if ( mParameters.mDiscretizationIndex >= 0 )
        {
            for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
            {
                // Create a B-spline field
                mPerturbationFields( iFieldIndex ) = std::make_shared< BSpline_Field >(
                        aMeshPair,
                        aOwnedADVs,
                        mSharedADVIDs( iFieldIndex ),
                        mOffsetID + iFieldIndex * aMeshPair.get_interpolation_mesh()->get_max_entity_id( mtk::EntityRank::BSPLINE, mParameters.mDiscretizationIndex ),
                        mParameters.mDiscretizationIndex,
                        mPerturbationFields( iFieldIndex ) );

                // Set analytic field index, for now
                Field::gDiscretizationIndex = mParameters.mDiscretizationIndex;

                mPerturbationFields( iFieldIndex )->mMeshPairForAnalytic = aMeshPair;
            }
        }
        else if ( mParameters.mDiscretizationIndex == -1 )
        {
            for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
            {
                // Just store nodal values
                mPerturbationFields( iFieldIndex ) = std::make_shared< Stored_Field >(
                        aMeshPair.get_interpolation_mesh(),
                        mPerturbationFields( iFieldIndex ) );

                mPerturbationFields( iFieldIndex )->mMeshPairForAnalytic = aMeshPair;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs )
    {
        MORIS_ERROR( false, "Surface mesh bspline fields cannot be remeshed for now" );
        for ( uint iFieldIndex = 0; iFieldIndex < Surface_Mesh::get_spatial_dimension(); iFieldIndex++ )
        {
            if ( mPerturbationFields( iFieldIndex )->get_name() == aMTKField->get_label() )
            {
                if ( mParameters.mDiscretizationIndex >= 0 )
                {
                    // Create a B-spline field
                    mPerturbationFields( iFieldIndex ) = std::make_shared< BSpline_Field >(
                            aOwnedADVs,
                            mSharedADVIDs( iFieldIndex ),
                            mOffsetID + iFieldIndex * aMeshPair.get_interpolation_mesh()->get_max_entity_id( mtk::EntityRank::BSPLINE, mParameters.mDiscretizationIndex ),
                            mParameters.mDiscretizationIndex,
                            aMTKField,
                            aMeshPair );
                }
                else if ( mPerturbationFields( iFieldIndex )->get_name() == aMTKField->get_label() && mParameters.mDiscretizationIndex == -1 )
                {
                    // TODO
                    MORIS_ERROR( false, "Stored field cannot be remeshed for now" );
                }
            }
            mPerturbationFields( iFieldIndex )->mMeshPairForAnalytic = aMeshPair;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_facet_center( const uint aFacetIndex )
    {
        // Get the coordinates of the vertices for this facet
        Matrix< DDRMat > tVertexCoordinates = Surface_Mesh::get_all_vertex_coordinates_of_facet( aFacetIndex );

        // Initialize the output matrix for the center
        Matrix< DDRMat > tFacetCenter( Surface_Mesh::get_spatial_dimension(), 1, 0.0 );

        // Loop over vertices and sum the coordinates
        for ( uint iVertex = 0; iVertex < tVertexCoordinates.n_cols(); iVertex++ )
        {
            tFacetCenter += tVertexCoordinates.get_column( iVertex );
        }

        // Divide by the number of vertices
        tFacetCenter = tFacetCenter / tVertexCoordinates.n_cols();

        return tFacetCenter;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::facet_vertex_depends_on_advs( const uint aFacetVertexIndex )
    {
        // Get the factor that scales this vertex's movement
        Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( Surface_Mesh::get_spatial_dimension(), 1.0 )
                                                                                    : get_discretization_scaling_user_defined( this->get_original_vertex_coordinates( aFacetVertexIndex ) );

        // Return true if this surface mesh can move, its movement was either defined by the user or is discretized
        // and not fixed in all directions by the user, and it lies in the Lagrange mesh domain
        return this->depends_on_advs()
           and ( get_vertex_adv_dependency_user_defined != nullptr
                   or ( mOriginalVertexBackgroundElements( aFacetVertexIndex ) != nullptr
                           and std::any_of( tFactor.cbegin(), tFactor.cend(),    //
                                   []( const real tDirectionFactor ) { return tDirectionFactor != 0.0; } ) ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_dvertex_dadv( const uint aFacetVertexIndex )
    {
        // Initialize sensitivity matrix
        Matrix< DDRMat > tVertexSensitivity;
        if ( this->intended_discretization() )
        {
            // Get the spatial dimension
            const uint tDims = Surface_Mesh::get_spatial_dimension();

            // Determine which directions the vertex can move in by quering the scaling function
            Vector< real > tFactor              = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 )
                                                                                                     : get_discretization_scaling_user_defined( this->get_original_vertex_coordinates( aFacetVertexIndex ) );
            uint           tNumDimsDependOnADVs = std::count_if( tFactor.cbegin(), tFactor.cend(), []( const real aFactor ) { return aFactor != 0.0; } );

            // Get the vertex indices and coordinates of the background element
            Matrix< DDRMat >   tVertexCoordinates = mOriginalVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
            Matrix< IndexMat > tVertexIndices     = mOriginalVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

            // Loop over background nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
            {
                // Get length before adding sensitivities for this node
                uint tNumVertexSensitivities = tVertexSensitivity.n_cols();

                bool tVertexSensitivitySizeDetermined = false;
                uint tDimensionSensitivitiesAdded     = 0;

                // Loop over spatial dimension
                for ( uint iDimensionIndex = 0; iDimensionIndex < tDims; iDimensionIndex++ )
                {
                    // Check that the vertex depends on ADVs in this direction
                    if ( tFactor( iDimensionIndex ) != 0.0 )
                    {
                        Matrix< DDRMat > tNodeSensitivity = tFactor( iDimensionIndex ) * mOriginalVertexBases( iNodeIndex, aFacetVertexIndex ) * mPerturbationFields( iDimensionIndex )->get_dfield_dadvs( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

                        // set size of sensitivity matrix
                        if ( not tVertexSensitivitySizeDetermined )
                        {
                            tVertexSensitivity.resize( tDims, tNumVertexSensitivities + tNumDimsDependOnADVs * tNodeSensitivity.numel() );
                            tVertexSensitivitySizeDetermined = true;
                        }

                        // Each sensitivity is a separate index
                        for ( uint iADVIndex = 0; iADVIndex < tNodeSensitivity.numel(); iADVIndex++ )
                        {
                            tVertexSensitivity( iDimensionIndex, tNumVertexSensitivities + tNodeSensitivity.length() * tDimensionSensitivitiesAdded + iADVIndex ) = tNodeSensitivity( iADVIndex );
                        }

                        tDimensionSensitivitiesAdded++;
                    }
                }
            }
        }
        else if ( get_vertex_adv_dependency_user_defined != nullptr )
        {
            get_dvertex_dadv_user_defined( this->get_original_vertex_coordinates( aFacetVertexIndex ), mADVHandler.get_values(), tVertexSensitivity );

            // Check that the user gave the correct size matrix
            MORIS_ASSERT( tVertexSensitivity.n_cols() == mADVHandler.get_determining_adv_ids().size(), "User defined function for vertex sensitivity needs to have as many columns as ADVs" );
            MORIS_ASSERT( tVertexSensitivity.n_rows() == Surface_Mesh::get_spatial_dimension(), "User defined function for vertex sensitivity needs to have as many rows as spatial dimensions" );
        }

        return tVertexSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Surface_Mesh_Geometry::get_vertex_adv_ids( const uint aFacetVertexIndex )
    {
        // Get spatial dimension of the problem
        const uint tDims = Surface_Mesh::get_spatial_dimension();

        // Initialize vector to be filled
        Vector< sint > tVertexADVIds;

        // The surface mesh is being optimized via B-spline fields
        if ( this->intended_discretization() )
        {
            // Determine which directions the vertex can move in by quering the scaling function
            Vector< real > tFactor              = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 )
                                                                                                     : get_discretization_scaling_user_defined( this->get_original_vertex_coordinates( aFacetVertexIndex ) );
            uint           tNumDimsDependOnADVs = std::count_if( tFactor.cbegin(), tFactor.cend(), []( const real aFactor ) { return aFactor != 0.0; } );

            // Get the vertex indices and coordinates of the background element
            Matrix< DDRMat >   tVertexCoordinates = mOriginalVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
            Matrix< IndexMat > tVertexIndices     = mOriginalVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

            // Loop over background nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
            {
                // Get the ADV IDs for this node
                Vector< sint > tNodeIDs = this->get_determining_adv_ids( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

                // Join the ADV IDs to the output
                // Get the original length
                uint tIDLength = tVertexADVIds.size();

                // Resize to add new ADV IDs
                tVertexADVIds.resize( tVertexADVIds.size() + ( tNodeIDs.size() * tNumDimsDependOnADVs ) / tDims );

                // Join the IDs
                uint tADVsAdded = 0;
                for ( uint iDimension = 0; iDimension < tDims; iDimension++ )
                {
                    if ( tFactor( iDimension ) != 0.0 )
                    {
                        for ( uint iADVIndex = 0; iADVIndex < tNodeIDs.size() / tDims; iADVIndex++ )
                        {
                            tVertexADVIds( tIDLength + tADVsAdded ) = tNodeIDs( ( iDimension * tNodeIDs.size() ) / tDims + iADVIndex );
                            tADVsAdded++;
                        }
                    }
                }
            }
        }
        // The surface mesh is being optimized by the user defined function
        else if ( get_vertex_adv_dependency_user_defined != nullptr )
        {
            // Get all the ADV IDs for this geometry
            tVertexADVIds = mADVHandler.get_determining_adv_ids();
        }

        return tVertexADVIds;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Surface_Mesh_Geometry::get_determining_adv_ids(
            const uint              aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        Vector< sint > tADVIDs;
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // Get the ADV IDs for this field
            Vector< sint > tFieldADVIDs;
            if ( mNodeManager->is_background_node( aNodeIndex ) )
            {
                tFieldADVIDs = mPerturbationFields( iFieldIndex )->get_determining_adv_ids( aNodeIndex, aCoordinates );
            }
            else
            {
                const Node_Manager& tNodeManager = *mNodeManager;
                const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );
                mPerturbationFields( iFieldIndex )->get_determining_adv_ids( tFieldADVIDs, tDerivedNode, *mNodeManager );
            }

            // Append the ADV IDs to the output matrix
            // Resize IDs
            uint tJoinedIDLength = tADVIDs.size();
            tADVIDs.resize( tJoinedIDLength + tFieldADVIDs.size() );

            // Join IDs
            for ( uint tAddedSensitivity = 0; tAddedSensitivity < tFieldADVIDs.size(); tAddedSensitivity++ )
            {
                tADVIDs( tJoinedIDLength + tAddedSensitivity ) = tFieldADVIDs( tAddedSensitivity );
            }
        }


        return tADVIDs;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< real > >
    Surface_Mesh_Geometry::determine_mtk_cell_bounding_box( const mtk::Cell* aElement )
    {
        // Initialize a bounding box
        Vector< Vector< real > > tBoundingBox( 2, Vector< real >( Surface_Mesh::get_spatial_dimension() ) );

        Matrix< DDRMat > tCurrentSearchElementVertexCoordinates = aElement->get_vertex_coords();

        // Build bounding box, set the box as the coordinates for the first vertex
        for ( uint iDimensionIndex = 0; iDimensionIndex < tCurrentSearchElementVertexCoordinates.n_cols(); iDimensionIndex++ )
        {
            tBoundingBox( 0 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( 0, iDimensionIndex );
            tBoundingBox( 1 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( 0, iDimensionIndex );

            // Loop over the rest of the vertices
            for ( uint iVertexIndex = 1; iVertexIndex < tCurrentSearchElementVertexCoordinates.n_rows(); iVertexIndex++ )
            {
                // check if the entry is less than the minimum
                if ( tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex ) < tBoundingBox( 0 )( iDimensionIndex ) )
                {
                    tBoundingBox( 0 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex );
                }
                // check if the entry is greater than the maximum
                if ( tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex ) > tBoundingBox( 1 )( iDimensionIndex ) )
                {
                    tBoundingBox( 1 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex );
                }
            }
        }

        return tBoundingBox;
    }

    //--------------------------------------------------------------------------------------------------------------

    const mtk::Cell* Surface_Mesh_Geometry::find_background_element_from_global_coordinates(
            const Matrix< DDRMat >& aCoordinate )
    {
        uint tDim = Surface_Mesh::get_spatial_dimension();

        // Loop through each mtk::Cell
        for ( uint iCellIndex = 0; iCellIndex < mMesh->get_num_elems(); iCellIndex++ )
        {
            // get this Cell's bounding box
            Vector< Vector< real > > tBoundingBox = this->determine_mtk_cell_bounding_box( &mMesh->get_mtk_cell( iCellIndex ) );

            // assume the point is in this bounding box
            bool tCoordinateInCell = true;

            // Loop over the bounding box and determine if this is true
            for ( uint iDimensionIndex = 0; iDimensionIndex < tDim; iDimensionIndex++ )
            {
                // check if the point is outside the bounding box in this dimension
                if ( aCoordinate( iDimensionIndex ) < tBoundingBox( 0 )( iDimensionIndex )
                        or aCoordinate( iDimensionIndex ) > tBoundingBox( 1 )( iDimensionIndex ) )
                {
                    tCoordinateInCell = false;
                }
            }

            // The coordinate is in the cell if it is within the bounding box in every dimension
            if ( tCoordinateInCell == true )
            {
                // If so, return this element's index
                return &mMesh->get_mtk_cell( iCellIndex );
            }
        }

        // if no element found, return -1
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::compute_vertex_basis(
            const mtk::Cell*        aBackgroundElement,
            const Matrix< DDRMat >& aParametricCoordinates )
    {
        // build interpolator
        mtk::Interpolation_Function_Factory tFactory;
        mtk::Interpolation_Function_Base*   tInterpolation;

        // create interpolation function based on spatial dimension of problem
        tInterpolation = tFactory.create_interpolation_function(
                aBackgroundElement->get_geometry_type(),
                mtk::Interpolation_Type::LAGRANGE,
                aBackgroundElement->get_interpolation_order() );

        // compute basis function at the vertices
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( aParametricCoordinates, tBasis );

        // clean up
        delete tInterpolation;

        return tBasis;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::update_vertex_basis_data()
    {
        // Get the spatial dimension
        uint tSpatialDimension = Surface_Mesh::get_spatial_dimension();

        if ( mCurrentVertexBackgroundElements.size() != Surface_Mesh::get_number_of_vertices() )
        {
            mCurrentVertexBases.resize( mMesh->get_mtk_cell( 0 ).get_number_of_vertices(), Surface_Mesh::get_number_of_vertices() );
            mCurrentVertexBackgroundElements.resize( Surface_Mesh::get_number_of_vertices() );
        }

        // Compute the bases for all facet vertices FIXME boilerplate
        for ( uint iVertexIndex = 0; iVertexIndex < Surface_Mesh::get_number_of_vertices(); iVertexIndex++ )
        {
            // Get this vertex's coordinates
            Matrix< DDRMat > tVertexCoordinates = Surface_Mesh::get_vertex_coordinates( iVertexIndex );

            // Determine which element this vertex lies in, will be the same for every field)
            mCurrentVertexBackgroundElements( iVertexIndex ) = this->find_background_element_from_global_coordinates( tVertexCoordinates );

            // check if the vertex is inside the mesh domain
            if ( mCurrentVertexBackgroundElements( iVertexIndex ) != nullptr )
            {
                // Get the bounding box for this element
                Vector< Vector< real > > tElementBoundingBox = this->determine_mtk_cell_bounding_box( mCurrentVertexBackgroundElements( iVertexIndex ) );

                // determine the local coordinates of the vertex inside the mtk::Cell
                for ( uint iDimensionIndex = 0; iDimensionIndex < tSpatialDimension; iDimensionIndex++ )
                {
                    mVertexParametricCoordinates( iDimensionIndex, iVertexIndex ) = 2.0 * ( tVertexCoordinates( iDimensionIndex ) - tElementBoundingBox( 0 )( iDimensionIndex ) )
                                                                                          / ( tElementBoundingBox( 1 )( iDimensionIndex ) - tElementBoundingBox( 0 )( iDimensionIndex ) )
                                                                                  - 1.0;
                }

                // Get the basis function values at the vertex location
                Matrix< DDRMat > tBasis = this->compute_vertex_basis( mCurrentVertexBackgroundElements( iVertexIndex ), mVertexParametricCoordinates.get_column( iVertexIndex ) );
                mCurrentVertexBases.set_column( iVertexIndex, trans( tBasis ) );
            }
        }

        // Store the first iteration of the bases and the background elements as they are needed for optimization
        if ( not mBasesComputed and this->depends_on_advs() )
        {
            mOriginalVertexBases              = mCurrentVertexBases;
            mOriginalVertexBackgroundElements = mCurrentVertexBackgroundElements;

            mBasesComputed = true;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::interpolate_perturbation_from_background_element(
            const mtk::Cell* aBackgroundElement,
            const uint       aFieldIndex,
            const uint       aFacetVertexIndex )
    {
        // get the indices and coordinates of the background element vertices
        Matrix< IndexMat > tVertexIndices        = aBackgroundElement->get_vertex_inds();
        Matrix< DDRMat >   tAllVertexCoordinates = aBackgroundElement->get_vertex_coords();

        // initialize field value at the node location
        real tPerturbation = 0.0;

        // get perturbation values at the vertices
        for ( uint iBackgroundNodeIndex = 0; iBackgroundNodeIndex < aBackgroundElement->get_number_of_vertices(); iBackgroundNodeIndex++ )
        {
            // add this vertex's field value to the value
            tPerturbation += mOriginalVertexBases.get_column( aFacetVertexIndex )( iBackgroundNodeIndex ) * mPerturbationFields( aFieldIndex )->get_field_value( tVertexIndices( iBackgroundNodeIndex ), { {} } );
        }

        return tPerturbation;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Surface_Mesh_Geometry::get_name()
    {
        return mName;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::string >
    Surface_Mesh_Geometry::get_field_names()
    {
        Vector< std::string > tFieldNames( mPerturbationFields.size() );
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // Should never be empty, as they are created by the geometry with a name
            tFieldNames( iFieldIndex ) = mPerturbationFields( iFieldIndex )->get_name();
        }

        return tFieldNames;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< Field > >
    Surface_Mesh_Geometry::get_fields()
    {
        return mPerturbationFields;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::get_design_info(
            const uint              aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Vector< real >&         aOutputDesignInfo )
    {
        // fit the output to the number of fields the surface mesh has
        aOutputDesignInfo.resize( mPerturbationFields.size() );

        // store the displacement value in every direction in the output
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            aOutputDesignInfo( iFieldIndex ) = mPerturbationFields( iFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::intended_discretization() const
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Surface_Mesh_Geometry::get_discretization_mesh_index() const
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::get_discretization_lower_bound() const
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::get_discretization_upper_bound() const
    {
        return mParameters.mDiscretizationUpperBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::update_dependencies( const Vector< std::shared_ptr< Design > >& aAllUpdatedDesigns )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen