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
#include "fn_join_horiz.hpp"

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
            , mName( aParameterList.get< std::string >( "name" ) )
            , mRegularizationType( aParameterList.get< Regularization_Type >( "regularization_type" ) )
            , mRegularizationFunctionName( aParameterList.get< std::string >( "regularization_function_name" ) )
            , mRegularizationSensitivityFunctionName( aParameterList.get< std::string >( "regularization_sensitivity_function_name" ) )
            , mRegularizationVertexIndsFunctionName( aParameterList.get< std::string >( "regularization_vertex_inds_function_name" ) )
            , mRegularizationFactors( aParameterList.get< Vector< real > >( "regularization_factors" ) )
            , mRegularizationIterations( aParameterList.get< moris_index >( "regularization_iterations" ) )
    {
        MORIS_ASSERT( (uint)( mDiscretizationIndex < -1 + not mAnalyticADVFunctionName.empty() ) < 2, "Both a discretization index and an analytical function are provided. Pick at most 1!" );
        MORIS_ASSERT( mRegularizationType != Regularization_Type::USER_DEFINED or ( mRegularizationType == Regularization_Type::USER_DEFINED and not mRegularizationFunctionName.empty() ),
                "Regularization type likely set to user defined but a regularization function name was not provided. Please set the \"regularization_function_name\" parameter under the surface mesh geometry." );
        MORIS_ASSERT( mRegularizationType == Regularization_Type::NONE or mRegularizationFactors.size() > 0, "If a regularization type is set, you must specify regularization factors corresponding to the regularization type." );
        MORIS_ASSERT( mRegularizationType == Regularization_Type::NONE or mRegularizationIterations > 0, "If a regularization type is set, you must specify the number of regularization iterations." );
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
        mParameters.mName = mParameters.mName.empty() ? mParameters.mFilePath.substr( mParameters.mFilePath.find_last_of( "/" ) + 1,
                                                                mParameters.mFilePath.find_last_of( "." ) - mParameters.mFilePath.find_last_of( "/" ) - 1 )
                                                      : mParameters.mName;

        // Get the regularization functions, build vertex connectivity for regularization if needed
        if ( this->do_regularization() )
        {
            this->load_regularization_function( aLibrary );
            this->build_vertex_connectivity();
        }

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
                        mParameters.mName + "_PERT_" + std::to_string( iFieldIndex ) );
            }
        }
        // Otherwise, check if this surface mesh is being optimized via an analytic function
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

        // Try to see if the pointer is in the stored list
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

        if ( std::abs( tIntersection.second ) > 1.0 + Surface_Mesh::mIntersectionTolerance or std::isnan( tIntersection.second ) )
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

        // Process the raycast output and if no valid intersection was found, recast ray with looser tolerance
        return this->process_raycast_for_local_coordinate( aFirstParentNode, aSecondParentNode, tRayDirection, Surface_Mesh::mIntersectionTolerance, tLocalCoordinate );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::pair< uint, real >
    Surface_Mesh_Geometry::process_raycast_for_local_coordinate(
            const Parent_Node&        aFirstParentNode,
            const Parent_Node&        aSecondParentNode,
            Matrix< DDRMat >&         aDirection,
            real                      aOriginalTolerance,
            mtk::Intersection_Vector& aRaycastResult )
    {
        // Iterate through the intersections to find the first value that is within the exclusive range (-1, 1)
        for ( auto& tIntersection : aRaycastResult )
        {
            // Put the intersection in the local coordinate frame
            tIntersection.second = 2.0 * tIntersection.second - 1.0;

            if ( std::abs( tIntersection.second ) < ( 1.0 - Surface_Mesh::mIntersectionTolerance ) )
            {
                // reset intersection tolerance
                Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;

                // Found a valid intersection that will NOT snap to either parent node
                return tIntersection;
            }
        }

        // If no middle intersection was found, check for intersections that are very close to the parent nodes
        for ( auto& tIntersection : aRaycastResult )
        {
            if ( std::abs( tIntersection.second - 1.0 ) < Surface_Mesh::mIntersectionTolerance )
            {
                // reset intersection tolerance
                Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;

                // Found an intersection that is very close to the second parent node
                // This will be ignored in the intersection node creation
                tIntersection.second = 1.0;
                return tIntersection;
            }
            else if ( std::abs( tIntersection.second + 1.0 ) < Surface_Mesh::mIntersectionTolerance )
            {
                // reset intersection tolerance
                Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;

                // Found an intersection that is very close to the first parent node
                // This will be ignored in the intersection node creation
                tIntersection.second = -1.0;
                return tIntersection;
            }
        }

        // ---------------------------------------------------------
        // No intersection found, check for pathological cases
        // ---------------------------------------------------------

        // Check to see if either parent node is on the interface. Needed if node will snap or raycast tolerancing is off
        if ( this->get_geometric_region( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() ) == Geometric_Region::INTERFACE )
        {
            Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;    // reset intersection tolerance
            return std::make_pair( 0, -1.0 );                             // FIXME: XTK will ignore this node, but ideally this should identify the intersecting facet
        }
        else if ( this->get_geometric_region( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() ) == Geometric_Region::INTERFACE )
        {
            Surface_Mesh::mIntersectionTolerance = aOriginalTolerance;    // reset intersection tolerance
            return std::make_pair( 0, 1.0 );                              // FIXME: XTK will ignore this node, but ideally this should identify the intersecting facet
        }

        // If no intersection was found (but there should be one), loosen the intersection tolerance and try again
        Surface_Mesh::mIntersectionTolerance *= 10.0;

        MORIS_ASSERT( Surface_Mesh::mIntersectionTolerance < 0.1, "Surface mesh intersection tolerance is too large. No valid intersection found." );

        // Recast the ray
        bool                     tWarning;
        mtk::Intersection_Vector tNewIntersections = this->cast_single_ray( aFirstParentNode.get_global_coordinates(), aDirection, tWarning );

        // Recursively call this function to process new intersections
        return this->process_raycast_for_local_coordinate( aFirstParentNode, aSecondParentNode, aDirection, aOriginalTolerance, tNewIntersections );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::load_regularization_function( std::shared_ptr< Library_IO > aLibrary )
    {
        switch ( mParameters.mRegularizationType )
        {
            case Regularization_Type::NONE:
            {
                // Pointers are already assigned to functions that return empty containers
                break;
            }
            case Regularization_Type::ISOTROPIC_LAPLACIAN:
            {
                regularize_mesh            = this->isotropic_laplacian_regularization;
                regularization_sensitivity = this->isotropic_laplacian_regularization_sensitivity;
                regularization_vertex_inds = this->get_determining_vertex_inds_isotropic_laplacian;
                break;
            }
            case Regularization_Type::ANSIOTROPIC_LAPLACIAN:
            {
                regularize_mesh            = this->anisotropic_laplacian_regularization;
                regularization_sensitivity = this->anisotropic_laplacian_regularization_sensitivity;
                regularization_vertex_inds = this->get_determining_vertex_inds_anisotropic_laplacian;
                break;
            }
            case Regularization_Type::TAUBIN:
            {
                regularize_mesh            = this->taubin_regularization;
                regularization_sensitivity = this->taubin_regularization_sensitivity;
                regularization_vertex_inds = this->get_determining_vertex_inds_taubin;
                break;
            }
            case Regularization_Type::USER_DEFINED:
            {
                regularize_mesh            = aLibrary->load_function< Regularization_Function >( mParameters.mRegularizationFunctionName );
                regularization_sensitivity = aLibrary->load_function< Regularization_Sensitivity_Function >( mParameters.mRegularizationSensitivityFunctionName );
                regularization_vertex_inds = aLibrary->load_function< Regularization_Vertex_Inds_Function >( mParameters.mRegularizationVertexIndsFunctionName );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Unknown Regularization type provided to surface mesh geometry %s.", mParameters.mName.c_str() );
                break;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::build_vertex_connectivity()
    {
        // Initialize a set to store unique vertex connectivity for every vertex
        uint                              tNumVertices = this->get_number_of_vertices();
        Vector< std::set< moris_index > > tVertexConnectivity( tNumVertices );

        uint tNumVertsPerFacet = Surface_Mesh::get_facets_vertex_indices( 0 ).size();

        // Loop through all facets of the surface mesh
        for ( auto& tFacet : Surface_Mesh::get_facet_connectivity() )
        {
            // Loop through all vertices of the facet
            for ( uint iVertex = 0; iVertex < tNumVertsPerFacet; iVertex++ )
            {
                // Get the current vertex and the next vertex (wrapping around)
                moris_index tCurrentVertex = tFacet( iVertex );
                moris_index tNextVertex    = tFacet( ( iVertex + 1 ) % tNumVertsPerFacet );

                // Add both to each others sets
                tVertexConnectivity( tCurrentVertex ).insert( tNextVertex );
                tVertexConnectivity( tNextVertex ).insert( tCurrentVertex );
            }
        }

        // Convert the sets to vectors and store them in the vertex connectivity
        mVertexConnectivity.resize( tNumVertices );
        for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
        {
            // Convert the set to a vector
            mVertexConnectivity( iVertex ) = Vector< moris_index >( tVertexConnectivity( iVertex ).begin(), tVertexConnectivity( iVertex ).end() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

#if MORIS_HAVE_ARBORX
    void Surface_Mesh_Geometry::flood_fill_mesh_regions()
    {
        Tracer tTracer( "GEN", "Surface Mesh Geometry", "Flood fill mesh nodes" );

        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // ----------------------------------------------------------------------------------------------
        // Initial setup for flood fill
        // ----------------------------------------------------------------------------------------------

        uint tNumNodes = mMesh->get_num_nodes();

        // Matrix< DDRMat >      tNodeCoordinates( Surface_Mesh::get_spatial_dimension(), tNumNodes );
        Vector< moris_index > tCellIndices( tNumNodes );

        // Build the node connectivity from the interpolation mesh
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
        Tracer tTracer( "GEN", "Surface Mesh Geometry", "Raycast remaining unknown nodes" );

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
                    // update the facet vertex if it moves explicitly via advs and is owned by this processor
                    if ( mOriginalVertexBackgroundElements( iVertexIndex ) != nullptr and mOriginalVertexBackgroundElements( iVertexIndex )->get_owner() == par_rank() )
                    {
                        // Interpolate the bspline field value at the facet vertex location
                        real tInterpolatedPerturbation = this->interpolate_perturbation_from_background_element(
                                mOriginalVertexBackgroundElements( iVertexIndex ),
                                iFieldIndex,
                                iVertexIndex );

                        // build the matrix for new coordinates
                        tOwnedVertexDisplacements( tDims, iVertexIndex )       = 1.0;    // says that this vertex is owned by this proc
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

                // Get the factor that scales this vertex's movement
                Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 ) : get_discretization_scaling_user_defined( Surface_Mesh::get_original_vertex_coordinates( iVertexIndex ) );

                // build the matrix for new coordinates
                for ( uint iFieldIndex = 0; iFieldIndex < tDims; iFieldIndex++ )
                {
                    tOwnedVertexDisplacements( tDims, iVertexIndex )       = 1.0;    // says that this vertex is owned by this proc
                    tOwnedVertexDisplacements( iFieldIndex, iVertexIndex ) = tFactor( iFieldIndex ) * tInterpolatedPerturbation( iFieldIndex );
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

        // Apply any surface mesh regularization if needed
        this->regularize();

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
        // this->flood_fill_mesh_regions();
#endif
        this->raycast_remaining_unknown_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::regularize()
    {
        if ( this->do_regularization() )
        {
            uint tDims = this->get_spatial_dimension();
            uint tNumV = this->get_number_of_vertices();

            // Get the regularization from whatever function was loaded
            Matrix< DDRMat > tRegularizationDisp = regularize_mesh( this->get_all_vertex_coordinates(), this->get_all_original_vertex_coordinates(), this->get_facet_connectivity(), mVertexConnectivity, mParameters.mRegularizationFactors, get_discretization_scaling_user_defined );

            MORIS_ASSERT( tRegularizationDisp.n_rows() == this->get_spatial_dimension(), "Regularization function for surface mesh %s does not return the correct number of rows. Expected %d, got %lu", mParameters.mName.c_str(), this->get_spatial_dimension(), tRegularizationDisp.n_rows() );
            MORIS_ASSERT( tRegularizationDisp.n_cols() == this->get_number_of_vertices(), "Regularization function for surface mesh %s does not return the correct number of cols. Expected %d, got %lu", mParameters.mName.c_str(), this->get_number_of_vertices(), tRegularizationDisp.n_cols() );

            // Scale the regularization displacement for all vertices
            for ( uint iVertexIndex = 0; iVertexIndex < this->get_number_of_vertices(); iVertexIndex++ )
            {
                // Get the factor that scales this vertex's movement
                Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( this->get_spatial_dimension(), 1.0 ) : get_discretization_scaling_user_defined( Surface_Mesh::get_original_vertex_coordinates( iVertexIndex ) );

                for ( uint iDim = 0; iDim < tDims; iDim++ )
                {
                    tRegularizationDisp( iDim, iVertexIndex ) *= tFactor( iDim );
                }
            }

            // Apply the regularization to the vertex displacements
            this->set_all_displacements( this->get_vertex_displacements() + tRegularizationDisp );

            // Compute regularization sensitivities for the first iteration
            if ( this->depends_on_advs() )
            {
                this->initialize_regularization_sensitivities();
            }

            std::map< uint, uint > tADVIdMap;    // Map to store the full sensitivities for this iteration, cleared for each vertex

            // See if we need to apply regularization again
            for ( moris_index iRegIter = 1; iRegIter < mParameters.mRegularizationIterations; iRegIter++ )
            {
                // Get the regularization from whatever function was loaded
                Matrix< DDRMat > tRegularizationDisp = regularize_mesh( this->get_all_vertex_coordinates(), this->get_all_original_vertex_coordinates(), this->get_facet_connectivity(), mVertexConnectivity, mParameters.mRegularizationFactors, get_discretization_scaling_user_defined );

                MORIS_ASSERT( tRegularizationDisp.n_rows() == this->get_spatial_dimension(), "Regularization function for surface mesh %s does not return the correct number of rows. Expected %d, got %lu", mParameters.mName.c_str(), this->get_spatial_dimension(), tRegularizationDisp.n_rows() );
                MORIS_ASSERT( tRegularizationDisp.n_cols() == tNumV, "Regularization function for surface mesh %s does not return the correct number of cols. Expected %d, got %lu", mParameters.mName.c_str(), tNumV, tRegularizationDisp.n_cols() );

                // Compute sensitivities and scale the regularization displacement by the user defined scaling factor
                for ( uint iV = 0; iV < tNumV; iV++ )
                {
                    // Get the factor that scales this vertex's movement
                    Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( this->get_spatial_dimension(), 1.0 ) : get_discretization_scaling_user_defined( Surface_Mesh::get_original_vertex_coordinates( iV ) );

                    // Scale the regularization displacement for this vertex
                    for ( uint iDim = 0; iDim < tDims; iDim++ )
                    {
                        tRegularizationDisp( iDim, iV ) *= tFactor( iDim );
                    }

                    if ( this->depends_on_advs() )
                    {
                        // Initialize new sensitivities
                        tADVIdMap.clear();

                        // Populate the ADV ID map for existing sensitivities for my vertex
                        uint tColumn = 0;
                        for ( const auto iADV : mRegularizationSensitivities( iV ).first )
                        {
                            tADVIdMap[ iADV ] = tColumn++;
                        }

                        // Initialize matrix for new sensitivities
                        Matrix< DDRMat > tVertexSensitivity( tDims, tColumn );

                        // Get surface mesh vertex dependencies
                        Vector< moris_index > tVertexDependencies = regularization_vertex_inds(
                                this->get_all_vertex_coordinates(),
                                this->get_facet_connectivity(),
                                mVertexConnectivity,
                                mParameters.mRegularizationFactors,
                                iV );

                        // Loop over all vertex dependencies and compute dv_i/dADV
                        for ( auto iDependency : tVertexDependencies )
                        {
                            // Check if this vertex depends on any ADVs
                            if ( this->facet_vertex_depends_on_advs( iDependency ) )
                            {
                                // Get dx_i/dx_k for this vertex
                                Matrix< DDRMat > tDviDvk = regularization_sensitivity(
                                        this->get_all_vertex_coordinates(),
                                        this->get_all_original_vertex_coordinates(),
                                        this->get_facet_connectivity(),
                                        mVertexConnectivity,
                                        mParameters.mRegularizationFactors,
                                        iV,
                                        iDependency,
                                        get_discretization_scaling_user_defined );

                                // Get the sensitivity and ADV IDs for this vertex
                                Vector< moris_index > tDependencyADVIds      = this->get_vertex_adv_ids( iDependency );
                                Matrix< DDRMat >      tDependencySensitivity = tDviDvk * this->get_dvertex_dadv( iDependency );

                                MORIS_ASSERT( tDependencyADVIds.size() == tDependencySensitivity.n_cols(),
                                        "Surface mesh %s has %lu ADV IDs for vertex %d, but regularization sensitivities have %lu columns.",
                                        mParameters.mName.c_str(),
                                        tDependencyADVIds.size(),
                                        iDependency,
                                        tDependencySensitivity.n_cols() );

                                // Loop through all the ADV IDs for this new sensitivity and update the sensitivities
                                for ( uint iADV = 0; iADV < tDependencyADVIds.size(); iADV++ )
                                {
                                    // Get the ADV ID for this sensitivity
                                    moris_index tADVID = tDependencyADVIds( iADV );

                                    // Check if this ADV ID is already in the sensitivities map for this iteration
                                    if ( tADVIdMap.find( tADVID ) == tADVIdMap.end() )
                                    {
                                        // If not, add it to the map with a new column index
                                        moris_index tColumnIndex = tADVIdMap.size();
                                        tADVIdMap[ tADVID ]      = tColumnIndex;

                                        // Resize the sensitivity matrix to accommodate the new column
                                        tVertexSensitivity.resize( tDims, tColumnIndex + 1 );
                                        tVertexSensitivity.set_column( tColumnIndex, tDependencySensitivity.get_column( iADV ) );
                                    }
                                    else
                                    {    // Get the column index for this ADV ID
                                        moris_index tColumnIndex = tADVIdMap[ tADVID ];

                                        // Update the sensitivity for this ADV ID
                                        tVertexSensitivity.set_column( tColumnIndex, tVertexSensitivity.get_column( tColumnIndex ) + tDependencySensitivity.get_column( iADV ) );
                                    }
                                }
                            }
                        }

                        // Convert ADV ID map to vector
                        Vector< moris_index > tKeys( tADVIdMap.size() );
                        for ( const auto& tPair : tADVIdMap )
                        {
                            // Store the ADV ID in the vector
                            tKeys( tPair.second ) = tPair.first;
                        }

                        // Set the sensitivities for this vertex
                        mRegularizationSensitivities( iV ) = std::make_pair( tKeys, tVertexSensitivity );

                        MORIS_ASSERT( mRegularizationSensitivities( iV ).first.size() == mRegularizationSensitivities( iV ).second.n_cols(),
                                "Surface mesh %s has %lu ADV IDs for vertex %d, but regularization sensitivities have %lu columns.",
                                mParameters.mName.c_str(),
                                mRegularizationSensitivities( iV ).first.size(),
                                iV,
                                mRegularizationSensitivities( iV ).second.n_cols() );
                    }
                }

                // Apply the regularization to the vertex displacements
                this->set_all_displacements( this->get_vertex_displacements() + tRegularizationDisp );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::initialize_regularization_sensitivities()
    {
        uint tDims = this->get_spatial_dimension();
        uint tNumV = this->get_number_of_vertices();

        mRegularizationSensitivities.resize( tNumV );

        // Initialize map for this vertex to keep track of unique ADV IDs
        std::map< uint, uint > tADVIdMap;

        // Loop through all vertices
        for ( uint iV = 0; iV < tNumV; iV++ )
        {
            tADVIdMap.clear();    // Clear the ADV ID map for this vertex

            // Initialize sensitivity matrix for this vertex
            Matrix< DDRMat > tVertexSensitivity;

            // Get the factor that scales this vertex's movement
            Vector< real > tFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 ) : get_discretization_scaling_user_defined( Surface_Mesh::get_original_vertex_coordinates( iV ) );

            // Get the regularization dependencies for this vertex
            Vector< moris_index > tVertexDependencies = regularization_vertex_inds( this->get_all_vertex_coordinates(), this->get_facet_connectivity(), mVertexConnectivity, mParameters.mRegularizationFactors, iV );

            // Loop through all the vertex dependencies and set the initial sensitivity to be identity
            for ( auto iDependency : tVertexDependencies )
            {
                // Check if this vertex depends on any ADVs explicitly
                if ( this->facet_vertex_explicitly_depends_on_advs( iDependency ) )
                {
                    // Get dx_i/dx_k regularization sensitivity for this dependency vertex
                    Matrix< DDRMat > tDviDvk = regularization_sensitivity(
                            this->get_all_vertex_coordinates(),
                            this->get_all_original_vertex_coordinates(),
                            this->get_facet_connectivity(),
                            mVertexConnectivity,
                            mParameters.mRegularizationFactors,
                            iV,
                            iDependency,
                            get_discretization_scaling_user_defined );

                    // Get the sensitivity and ADV IDs for this vertex
                    Vector< moris_index > tDependencyADVIds      = this->get_vertex_adv_ids_explicit( iDependency );
                    Matrix< DDRMat >      tDependencySensitivity = tDviDvk * this->get_dvertex_dadv_explicit( iDependency );

                    // Loop through all ADV IDs and update the sensitivity matrix
                    for ( uint iADV = 0; iADV < tDependencyADVIds.size(); iADV++ )
                    {
                        // Get the ADV ID for this sensitivity
                        moris_index tADVID = tDependencyADVIds( iADV );

                        // // Scale the sensitivity by the user defined scaling factor
                        // for ( uint iDim = 0; iDim < tDims; iDim++ )
                        // {
                        //     tDependencySensitivity( iDim, iADV ) *= tFactor( iDim );
                        // }

                        // Check if this ADV ID is already in the sensitivities map
                        if ( tADVIdMap.find( tADVID ) == tADVIdMap.end() )
                        {
                            // If not, add it to the map with a new column index
                            moris_index tColumnIndex = tADVIdMap.size();
                            tADVIdMap[ tADVID ]      = tColumnIndex;

                            // Resize the sensitivity matrix to accommodate the new column
                            tVertexSensitivity.resize( tDims, tColumnIndex + 1 );
                            tVertexSensitivity.set_column( tColumnIndex, tDependencySensitivity.get_column( iADV ) );
                        }
                        else
                        {    // Get the column index for this ADV ID
                            moris_index tColumnIndex = tADVIdMap[ tADVID ];

                            // Update the sensitivity for this ADV ID
                            tVertexSensitivity.set_column( tColumnIndex, tVertexSensitivity.get_column( tColumnIndex ) + tDependencySensitivity.get_column( iADV ) );
                        }
                    }
                }
            }

            // Store the sensitivities for this vertex
            Vector< moris_index > tKeys( tADVIdMap.size() );
            for ( const auto& tPair : tADVIdMap )
            {
                // Fill the vector with the keys from the map
                tKeys( tPair.second ) = tPair.first;
            }
            mRegularizationSensitivities( iV ) = std::make_pair( tKeys, tVertexSensitivity );

            MORIS_ASSERT( mRegularizationSensitivities( iV ).first.size() == mRegularizationSensitivities( iV ).second.n_cols(),
                    "Surface mesh %s has %lu ADV IDs for vertex %d, but regularization sensitivities have %lu columns.",
                    mParameters.mName.c_str(),
                    mRegularizationSensitivities( iV ).first.size(),
                    iV,
                    mRegularizationSensitivities( iV ).second.n_cols() );
        }
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

    bool Surface_Mesh_Geometry::facet_vertex_explicitly_depends_on_advs( const uint aFacetVertexIndex ) const
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

    bool Surface_Mesh_Geometry::facet_vertex_depends_on_advs( const uint aFacetVertexIndex ) const
    {
        return this->do_regularization() ? mRegularizationSensitivities( aFacetVertexIndex ).first.size() > 0 : this->facet_vertex_explicitly_depends_on_advs( aFacetVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Vector< moris_index > >&
    Surface_Mesh_Geometry::get_all_vertex_connectivity() const
    {
        return mVertexConnectivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< moris_index >&
    Surface_Mesh_Geometry::get_vertex_connectivity( uint aVertexIndex ) const
    {
        return mVertexConnectivity( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_dvertex_dadv( const uint aFacetVertexIndex ) const
    {
        return this->do_regularization() ? mRegularizationSensitivities( aFacetVertexIndex ).second : this->get_dvertex_dadv_explicit( aFacetVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_dvertex_dadv_explicit( const uint aFacetVertexIndex ) const
    {
        // Get the spatial dimension
        const uint tDims = Surface_Mesh::get_spatial_dimension();

        // Initialize sensitivity matrix
        Matrix< DDRMat > tVertexSensitivity( tDims, 0 );
        if ( this->intended_discretization() and this->facet_vertex_explicitly_depends_on_advs( aFacetVertexIndex ) )
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
    Surface_Mesh_Geometry::get_vertex_adv_ids( const uint aFacetVertexIndex ) const
    {
        return this->do_regularization() ? mRegularizationSensitivities( aFacetVertexIndex ).first : this->get_vertex_adv_ids_explicit( aFacetVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index >
    Surface_Mesh_Geometry::get_vertex_adv_ids_explicit( const uint aFacetVertexIndex ) const
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
            const Matrix< DDRMat >& aCoordinates ) const
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
        if ( mParameters.mDelaunay or this->intended_discretization() )
        {
            Tracer tTracer( "GEN", "Surface Mesh Geometry", "Find background element brute force" );

            // Get the spatial dimension
            uint tSpatialDimension = Surface_Mesh::get_spatial_dimension();

            if ( mCurrentVertexBackgroundElements.size() != Surface_Mesh::get_number_of_vertices() )
            {
                mCurrentVertexBases.resize( mMesh->get_mtk_cell( 0 ).get_number_of_vertices(), Surface_Mesh::get_number_of_vertices() );
                mCurrentVertexBackgroundElements.resize( Surface_Mesh::get_number_of_vertices() );
            }

            // Compute the bases for all facet vertices FIXME brute force is so slow
            for ( uint iVertexIndex = 0; iVertexIndex < Surface_Mesh::get_number_of_vertices(); iVertexIndex++ )
            {
                // Get this vertex's coordinates
                Matrix< DDRMat > tVertexCoordinates = Surface_Mesh::get_vertex_coordinates( iVertexIndex );

                // Determine which element this vertex lies in, will be the same for every field
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
            if ( not mBasesComputed and this->intended_discretization() )
            {
                mOriginalVertexBases              = mCurrentVertexBases;
                mOriginalVertexBackgroundElements = mCurrentVertexBackgroundElements;

                mBasesComputed = true;
            }
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

    Matrix< DDRMat > Surface_Mesh_Geometry::isotropic_laplacian_regularization(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        if ( aFactor.size() > 1 )
        {
            MORIS_LOG_WARNING( "Isotropic laplacian regularization only supports a single factor, using the first value." );
        }

        uint tNumVertices = aVertexCoordinates.n_cols();
        uint tDims        = aVertexCoordinates.n_rows();

        Matrix< DDRMat > tRegularizationDisp( tDims, tNumVertices, 0.0 );

        // Loop over all vertices
        for ( uint iV = 0; iV < tNumVertices; iV++ )
        {
            Matrix< DDRMat > tVertexRegularization( tDims, 1, 0.0 );

            Vector< real > tVertexFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 )
                                                                                              : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( iV ) );

            // Loop over all edges associated with this vertex
            for ( uint iE = 0; iE < aVertexConnectivity( iV ).size(); iE++ )
            {
                // Get the vertex index of the edge
                moris_index tDependencyIndex = aVertexConnectivity( iV )( iE );

                // Compute the edge vector scaled by the factor
                Matrix< DDRMat > tEdge( tDims, 1, 0.0 );
                for ( uint iDimension = 0; iDimension < tDims; iDimension++ )
                {
                    // Scale the edge vector by the factor
                    tEdge( iDimension ) = tVertexFactor( iDimension ) * ( aVertexCoordinates( iDimension, tDependencyIndex ) - aVertexCoordinates( iDimension, iV ) );
                }

                // Add the edge vector to the regularization for this vertex
                tVertexRegularization += tEdge;
            }

            // Divide by the number of edges to get the average, set into output matrix
            tRegularizationDisp.set_column( iV, aFactor( 0 ) / static_cast< real >( aVertexConnectivity( iV ).size() ) * tVertexRegularization );
        }

        return tRegularizationDisp;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::isotropic_laplacian_regularization_sensitivity(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aQueryVertexIndex,
            const uint                             aDependencyVertexIndex,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        uint tDim = aVertexCoordinates.n_rows();

        // Get the factor that scales this vertex's movement
        Vector< real > tVertexFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDim, 1.0 ) : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( aQueryVertexIndex ) );

        // Regularization sensitivity for the vertex itself
        if ( aQueryVertexIndex == aDependencyVertexIndex )
        {
            Matrix< DDRMat > tI = eye( tDim, tDim );

            // Multiply by the factor
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                tI( iDim, iDim ) *= ( 1.0 - aFactor( 0 ) * tVertexFactor( iDim ) );
            }

            return tI;
        }
        // Regularization sensitivity for all other vertices in the dependencies
        else
        {
            Matrix< DDRMat > tI = eye( tDim, tDim );

            // Multiply by the factor
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                tI( iDim, iDim ) *= tVertexFactor( iDim );
            }

            // Get the connectivity for this vertex
            real tNumConnectedVertices = static_cast< real >( aVertexConnectivity( aQueryVertexIndex ).size() );

            return ( aFactor( 0 ) / tNumConnectedVertices ) * tI;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Surface_Mesh_Geometry::get_determining_vertex_inds_isotropic_laplacian(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aVertexIndex )
    {
        Vector< moris_index > tDependencies = aVertexConnectivity( aVertexIndex );

        // Add the requested vertex index to the dependencies
        tDependencies.push_back( aVertexIndex );

        return tDependencies;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::anisotropic_laplacian_regularization(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        MORIS_ERROR( false, "Surface_Mesh_Geometry::anisotropic_laplacian_regularization: implement me" );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::anisotropic_laplacian_regularization_sensitivity(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aQueryVertexIndex,
            const uint                             aDependencyVertexIndex,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        MORIS_ERROR( false, "Surface_Mesh_Geometry::isotropic_laplacian_regularization_sensitivity: implement me" );

        return { { {} } };
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Surface_Mesh_Geometry::get_determining_vertex_inds_anisotropic_laplacian(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aVertexIndex )
    {
        Vector< moris_index > tDependencies = aVertexConnectivity( aVertexIndex );

        // Add the requested vertex index to the dependencies
        tDependencies.push_back( aVertexIndex );

        return tDependencies;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::taubin_regularization(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        if ( aFactor.size() > 2 )
        {
            MORIS_LOG_WARNING( "Taubin regularization only supports using two factors, using the first two." );
        }
        MORIS_ASSERT( aFactor( 0 ) > 0.0 and aFactor( 1 ) < 0.0, "Taubin regularization uses two factors, the first must be positive and the second negative." );
        MORIS_ASSERT( std::abs( aFactor( 1 ) ) > aFactor( 0 ), "Taubin regularization uses a negative factor that is larger than the positive factor." );

        uint tNumVertices = aVertexCoordinates.n_cols();
        uint tDims        = aVertexCoordinates.n_rows();

        Matrix< DDRMat > tRegularizationDisp( tDims, tNumVertices, 0.0 );

        // Loop over all vertices - shrink
        for ( uint iV = 0; iV < tNumVertices; iV++ )
        {
            Matrix< DDRMat > tVertexRegularization( tDims, 1, 0.0 );

            // Get the user defined scaling factor for this vertex
            Vector< real > tVertexFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 )
                                                                                              : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( iV ) );

            // Loop over all edges associated with this vertex
            for ( uint iE = 0; iE < aVertexConnectivity( iV ).size(); iE++ )
            {
                uint tDependencyIndex = aVertexConnectivity( iV )( iE );

                // Compute the edge vector scaled by the factor
                Matrix< DDRMat > tEdge( tDims, 1, 0.0 );
                for ( uint iDimension = 0; iDimension < tDims; iDimension++ )
                {
                    // Scale the edge vector by the factor
                    tEdge( iDimension ) = tVertexFactor( iDimension ) * ( aVertexCoordinates( iDimension, tDependencyIndex ) - aVertexCoordinates( iDimension, iV ) );
                }

                // Add the edge vector to the regularization for this vertex
                tVertexRegularization += tEdge;
            }

            // Divide by the number of edges to get the average, set into output matrix
            tRegularizationDisp.set_column( iV, aFactor( 0 ) / static_cast< real >( aVertexConnectivity( iV ).size() ) * tVertexRegularization );
        }

        // Store the first iteration displacements
        Matrix< DDRMat > tFirstIterationDisplacements = tRegularizationDisp;

        // Loop over all vertices
        for ( uint iV = 0; iV < tNumVertices; iV++ )
        {
            Matrix< DDRMat > tVertexRegularization( tDims, 1, 0.0 );

            // Get the user defined scaling factor for this vertex
            Vector< real > tVertexFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDims, 1.0 )
                                                                                              : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( iV ) );

            // Loop over all edges associated with this vertex
            for ( uint iE = 0; iE < aVertexConnectivity( iV ).size(); iE++ )
            {
                uint tDependencyIndex = aVertexConnectivity( iV )( iE );

                // Compute the edge vector scaled by the factor
                Matrix< DDRMat > tEdge( tDims, 1, 0.0 );
                for ( uint iDimension = 0; iDimension < tDims; iDimension++ )
                {
                    // Scale the edge vector by the factor
                    tEdge( iDimension ) = tVertexFactor( iDimension ) * ( aVertexCoordinates( iDimension, tDependencyIndex ) + tFirstIterationDisplacements( iDimension, tDependencyIndex ) - aVertexCoordinates( iDimension, iV ) - tFirstIterationDisplacements( iDimension, iV ) );
                }

                // Add the edge vector to the regularization for this vertex
                tVertexRegularization += tEdge;
            }

            // Divide by the number of edges to get the average, set into output matrix
            tRegularizationDisp.set_column( iV, aFactor( 1 ) / static_cast< real >( aVertexConnectivity( iV ).size() ) * tVertexRegularization );
        }

        return tRegularizationDisp + tFirstIterationDisplacements;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Surface_Mesh_Geometry::taubin_regularization_sensitivity(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Matrix< DDRMat >&                aOriginalVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aQueryVertexIndex,
            const uint                             aDependencyVertexIndex,
            const Discretization_Factor_Function&  get_discretization_scaling_user_defined )
    {
        uint tDim = aVertexCoordinates.n_rows();

        // Get the user defined scaling factor for this vertex
        Vector< real > tVertexFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDim, 1.0 ) : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( aQueryVertexIndex ) );

        // Get connectivity for this vertex
        const Vector< moris_index >& tVertexConnectivity = aVertexConnectivity( aQueryVertexIndex );

        Matrix< DDRMat > tSens( tDim, tDim );

        // Sensitivity for the vertex itself
        if ( aQueryVertexIndex == aDependencyVertexIndex )
        {
            // Loop over dimensions
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                tSens( iDim, iDim ) = ( 1.0 - tVertexFactor( iDim ) * aFactor( 0 ) ) * ( 1.0 - tVertexFactor( iDim ) * aFactor( 1 ) );
            }
        }
        // Vertex is a direct dependency of the query vertex
        else if ( std::any_of( tVertexConnectivity.cbegin(), tVertexConnectivity.cend(), [ aDependencyVertexIndex ]( moris_index aIndex ) { return (uint)aIndex == aDependencyVertexIndex; } ) )
        {
            // Get the connectivity for this vertex
            real tNumConnectedVertices = static_cast< real >( tVertexConnectivity.size() );

            // Get the user defined scaling factor for this vertex
            Vector< real > tDependencyFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDim, 1.0 )
                                                                                                  : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( aDependencyVertexIndex ) );

            // Loop over dimensions
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                tSens( iDim, iDim ) = ( 1.0 - tVertexFactor( iDim ) * aFactor( 1 ) ) * tVertexFactor( iDim ) * aFactor( 0 ) / tNumConnectedVertices
                                    + ( 1.0 - tDependencyFactor( iDim ) * aFactor( 0 ) ) * tVertexFactor( iDim ) * aFactor( 1 ) / tNumConnectedVertices;
            }
        }

        // Loop over the connected vertices and check if the dependency vertex is in the connectivity of the connected vertex (check for second order dependencies)
        for ( auto iDependency : tVertexConnectivity )
        {
            if ( std::find( aVertexConnectivity( iDependency ).cbegin(), aVertexConnectivity( iDependency ).cend(), aDependencyVertexIndex ) != aVertexConnectivity( iDependency ).cend() )
            {
                // Get the user defined scaling factor for this vertex
                Vector< real > tDependencyFactor = get_discretization_scaling_user_defined == nullptr ? Vector< real >( tDim, 1.0 )
                                                                                                      : get_discretization_scaling_user_defined( aOriginalVertexCoordinates.get_column( iDependency ) );

                // Get the connectivity for this vertex
                real tNumDependencyConnectedVertices = static_cast< real >( aVertexConnectivity( iDependency ).size() );
                real tNumConnectedVertices           = static_cast< real >( tVertexConnectivity.size() );

                // Loop over dimensions
                for ( uint iDim = 0; iDim < tDim; iDim++ )
                {
                    tSens( iDim, iDim ) += tVertexFactor( iDim ) * tDependencyFactor( iDim ) * aFactor( 0 ) * aFactor( 1 ) / tNumConnectedVertices / tNumDependencyConnectedVertices;
                }
            }
        }

        return tSens;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Surface_Mesh_Geometry::get_determining_vertex_inds_taubin(
            const Matrix< DDRMat >&                aVertexCoordinates,
            const Vector< Vector< moris_index > >& aFacetConnectivity,
            const Vector< Vector< moris_index > >& aVertexConnectivity,
            const Vector< real >&                  aFactor,
            const uint                             aVertexIndex )
    {
        // Output size: number of connected vertices + number of connections for each connected vertex + 1 for the vertex itself
        uint tNumFirstConnections = aVertexConnectivity( aVertexIndex ).size();

        // Loop through each connection and determine the number of connected vertices to determine total output size
        uint tNumSecondConnections = 0;
        for ( auto iV : aVertexConnectivity( aVertexIndex ) )
        {
            tNumSecondConnections += aVertexConnectivity( iV ).size();
        }
        uint tNumSensitivities = tNumFirstConnections + tNumSecondConnections + 1;    // + 1 for the vertex itself

        // Initialize output vector
        Vector< moris_index > tDependencies( tNumSensitivities );

        // Loop through each connection and get that vertex's connection indices
        uint iDependency = 0;
        for ( uint iV = 0; iV < aVertexConnectivity( aVertexIndex ).size(); iV++ )
        {
            // Set this vertex index as a dependency
            tDependencies( iDependency++ ) = aVertexConnectivity( aVertexIndex )( iV );

            // Loop through connections to this other vertex
            for ( uint iSecondV = 0; iSecondV < aVertexConnectivity( aVertexConnectivity( aVertexIndex )( iV ) ).size(); iSecondV++ )
            {
                // Get the vertex index of this connection
                tDependencies( iDependency++ ) = aVertexConnectivity( aVertexConnectivity( aVertexIndex )( iV ) )( iSecondV );
            }
        }

        // Add the vertex index itself to the dependencies
        tDependencies( iDependency++ ) = aVertexIndex;

        // Make unique
        unique( tDependencies );

        MORIS_ASSERT( std::find( tDependencies.cbegin(), tDependencies.cend(), aVertexIndex ) != tDependencies.cend(),
                "Surface_Mesh_Geometry::get_determining_vertex_inds_taubin: The vertex index should be in the dependencies." );

        return tDependencies;
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
        return mParameters.mName;
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

    bool Surface_Mesh_Geometry::do_regularization() const
    {
        return mParameters.mRegularizationType != Regularization_Type::NONE;
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

    real Surface_Mesh_Geometry::compute_GQI_curvature() const
    {
        // TO IMPLEMENT
        return 0.0;
    }

}    // namespace moris::gen