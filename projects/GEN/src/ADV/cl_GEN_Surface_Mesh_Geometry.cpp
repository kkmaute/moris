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

#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
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
#include "cl_MTK_Enums.hpp"

#include "cl_SOL_Dist_Map.hpp"

// BRENDAN FORMAT
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
            , mDiscretizationFactorFuntionName( aParameterList.get< std::string >( "discretization_factor_function_name" ) )
            , mAnalyticPerturbationFunctionName( aParameterList.get< std::string >( "field_function_name" ) )
            , mPerturbationSensitivityFunctionName( aParameterList.get< std::string >( "sensitivity_function_name" ) )
    {
        MORIS_ASSERT( ( mDiscretizationIndex > -2 ) + ( not mAnalyticPerturbationFunctionName.empty() ) < 2, "Surface mesh should be optimized by either a B-spline field or an analytic function - not both!" );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry(
            mtk::Mesh*                     aMesh,
            ADV_Manager&                   aADVManager,
            const Surface_Mesh_Parameters& aParameters,
            const Vector< ADV >&           aADVs,
            Node_Manager&                  aNodeManager,
            std::shared_ptr< Library_IO >  aLibrary )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Surface_Mesh( mtk::load_vertices_from_object_file( aParameters.mFilePath, aParameters.mOffsets, aParameters.mScale ),
                      mtk::load_facets_from_object_file( aParameters.mFilePath ),
                      aParameters.mIntersectionTolerance )
            , mParameters( aParameters )
            , mADVHandler( aADVs )
            , mNodeManager( &aNodeManager )
            , mMesh( aMesh )
            , mPerturbationFields( 0 )
            , mVertexBases( 0, 0 )
            , mVertexBackgroundElements( 0 )
    {
        uint tDim = Surface_Mesh::get_spatial_dimension();

        // parse the file path and extract the file name
        mName = aParameters.mFilePath.substr( aParameters.mFilePath.find_last_of( "/" ) + 1,
                aParameters.mFilePath.find_last_of( "." ) - aParameters.mFilePath.find_last_of( "/" ) - 1 );

        // If this surface mesh is being optimized, construct fields,
        // store original vertex coordinates, determine which facet vertices are fixed, and compute the bases for all vertices
        if ( this->intended_discretization() )
        {
            // STEP 1: Determine which facet vertices are fixed
            if ( not mParameters.mDiscretizationFactorFuntionName.empty() )
            {
                // Get pointer to function
                get_discretization_factor_user_defined = aLibrary->load_function< DISCRETIZATION_FACTOR_FUNCTION >( mParameters.mDiscretizationFactorFuntionName );
            }

            // STEP 2: Construct perturbation fields
            mPerturbationFields.resize( tDim );

            // build perturbation fields
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
        else if ( not mParameters.mAnalyticPerturbationFunctionName.empty() )
        {
            // Get pointer to perturbations and sensitivities
            get_analytic_perturbation_user_defined = aLibrary->load_function< ANALYTIC_PERTURBATION_FUNCTION >( mParameters.mAnalyticPerturbationFunctionName );

            get_perturbation_sensitivity_user_defined = aLibrary->load_function< PERTURBATION_SENSITIVITY_FUNCTION >( mParameters.mPerturbationSensitivityFunctionName );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::~Surface_Mesh_Geometry()
    {
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
        // Raycast from the point
        mtk::Mesh_Region tRegion = this->raycast_point( aNodeCoordinates );

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
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Unexpected sdf::Object_Region of %d returned from raycast.", tRegion );
                return Geometric_Region::INTERFACE;
            }
        }
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
        uint tParentFacet     = MORIS_UINT_MAX;
        real tLocalCoordinate = this->compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode, tParentFacet );

        MORIS_ERROR( tParentFacet != MORIS_UINT_MAX or ( tLocalCoordinate < 1.0 + this->get_intersection_tolerance() and tLocalCoordinate > -1.0 - this->get_intersection_tolerance() ),
                "Intersection node local coordinate is not between -1 and 1 or parent facet is null. Local coordinate = %f",
                tLocalCoordinate );


        // Create surface mesh intersection node
        return new Intersection_Node_Surface_Mesh(
                aNodeIndex,
                aBackgroundNodes,
                aFirstParentNode,
                aSecondParentNode,
                tLocalCoordinate,
                tParentFacet,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::compute_intersection_local_coordinate(
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            uint&                             aParentFacet )
    {
        // Get the direction for the raycast
        Matrix< DDRMat > tRayDirection = aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates();

        //  Compute the distance from the first parent to all the facets
        Vector< uint > tIntersectionFacets;
        Vector< real > tLocalCoordinate = this->compute_ray_facet_intersections( aFirstParentNode.get_global_coordinates(), tRayDirection, tIntersectionFacets );

        // Put the intersections in the local coordinate frame
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            tLocalCoordinate( iIntersection ) = 2.0 / norm( aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates() ) * ( tLocalCoordinate( iIntersection ) ) - 1.0;
        }

        // -------------------------------------------------------------------------------------
        // STEP 3: Process the intersection information and determine if the surface mesh intersects the parent edge
        // -------------------------------------------------------------------------------------

        // check number of intersections along parent edge
        uint tNumberOfParentEdgeIntersections = 0;
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            if ( tLocalCoordinate( iIntersection ) < 1.0 + this->get_intersection_tolerance() and tLocalCoordinate( iIntersection ) > -1.0 - this->get_intersection_tolerance() )
            {
                tNumberOfParentEdgeIntersections++;
            }
        }

        // no intersections detected or multiple along parent edge
        if ( tLocalCoordinate.size() == 0 )
        {
            aParentFacet = MORIS_UINT_MAX;
            return MORIS_REAL_MAX;
        }

        // Set return values for intersection location and associated facet
        MORIS_ASSERT( tIntersectionFacets.size() == tLocalCoordinate.size(), "Inconsistent size of facet vector (size %lu) and local coordinate vector (size %lu)", tIntersectionFacets.size(), tLocalCoordinate.size() );
        aParentFacet = tIntersectionFacets( 0 );

        return tLocalCoordinate( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        uint tDim         = Surface_Mesh::get_spatial_dimension();
        uint tNumVertices = Surface_Mesh::get_number_of_vertices();

        // Get the coordinates of the owned vertices to communicate to other processors (first rows = coordinates, last row = owned flag)
        Matrix< DDRMat > tOwnedVertexDisplacements( tDim + 1, tNumVertices );

        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // Import advs to field
            mPerturbationFields( iFieldIndex )->import_advs( aOwnedADVs );
        }

        // Add this vertex's movement to the owned vertex coordinates
        for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
        {
            if ( this->intended_discretization() )
            {
                // Get the factor that scales this vertex's movement
                Vector< real > tFactor = get_discretization_factor_user_defined == nullptr ? Vector< real >( tDim, 1.0 ) : get_discretization_factor_user_defined( Surface_Mesh::get_vertex_coordinates( iVertexIndex ) );

                for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
                {
                    // update the facet vertex if it can move and is owned by this processor
                    if ( this->facet_vertex_depends_on_advs( iVertexIndex ) and mVertexBackgroundElements( iVertexIndex )->get_owner() == par_rank() )
                    {
                        // Interpolate the bspline field value at the facet vertex location
                        real tInterpolatedPerturbation = this->interpolate_perturbation_from_background_element(
                                mVertexBackgroundElements( iVertexIndex ),
                                iFieldIndex,
                                iVertexIndex );

                        // build the matrix for new coordinates
                        tOwnedVertexDisplacements( tDim, iVertexIndex )        = 1.0;    // says that this vertex is owned by this proc
                        tOwnedVertexDisplacements( iFieldIndex, iVertexIndex ) = tFactor( iFieldIndex ) * tInterpolatedPerturbation;
                    }
                }
            }
            else if ( get_analytic_perturbation_user_defined != nullptr )
            {
                // Get the perturbation value at the vertex
                Vector< real > tPerturbation = get_analytic_perturbation_user_defined( Surface_Mesh::get_vertex_coordinates( iVertexIndex ), mADVHandler.get_values() );

                // build the matrix for new coordinates
                tOwnedVertexDisplacements( tDim, iVertexIndex ) = 1.0;    // says that this vertex is owned by this proc
                for ( uint iFieldIndex = 0; iFieldIndex < tDim; iFieldIndex++ )
                {
                    tOwnedVertexDisplacements( iFieldIndex, iVertexIndex ) = tPerturbation( iFieldIndex );
                }
            }
        }

        // Get the vertex coordinates from all processors and put in a vector of mats on base proc
        Vector< Matrix< DDRMat > > tAllVertexDisplacements;
        gatherv_mats( tOwnedVertexDisplacements, tAllVertexDisplacements );

        // Build matrix with all vertex coordinates on base proc
        Matrix< DDRMat > tCombinedVertexCoordinates( tDim + 1, tNumVertices );
        if ( par_rank() == 0 )
        {
            for ( uint iProcessor = 0; iProcessor < tAllVertexDisplacements.size(); iProcessor++ )
            {
                for ( uint iVertexIndex = 0; iVertexIndex < tNumVertices; iVertexIndex++ )
                {
                    // TODO: check if vertices are shared and if so that the coordinates are the same

                    // Check to see if the vertex was owned by proc iProcessor
                    if ( (uint)tAllVertexDisplacements( iProcessor )( tDim, iVertexIndex ) == 1 )
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
            if ( (uint)tCombinedVertexCoordinates( tDim, iVertexIndex ) == 1 )
            {
                Surface_Mesh::append_vertex_displacement( iVertexIndex, tCombinedVertexCoordinates.get_column( iVertexIndex ) );
            }
        }

        // Update the facet's information based on the new vertex coordinates
        Surface_Mesh::initialize_facet_normals();
        Surface_Mesh::construct_bvh();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::set_advs( sol::Dist_Vector* aADVs )
    {
        // Have each field import the advs
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->set_advs( aADVs );
        }

        // Let the ADV handler know about the ADVs as well
        mADVHandler.set_advs( aADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    sint Surface_Mesh_Geometry::append_adv_info(
            mtk::Interpolation_Mesh* aMesh,
            Vector< sint >&          aOwnedADVIds,
            Matrix< IdMat >&         aOwnedijklIDs,
            sint                     aOffsetID,
            Vector< real >&          aLowerBounds,
            Vector< real >&          aUpperBounds )
    {
        // Get the original offset ID
        sint tOriginalOffsetID = aOffsetID;

        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            aOffsetID = Design::append_adv_info(
                    aMesh,
                    aOwnedADVIds,
                    aOwnedijklIDs,
                    aOffsetID,
                    aLowerBounds,
                    aUpperBounds );
        }

        // reset the offset back to the offset for the first perturabtion field (mOffsetID was changed in the above loop)
        mOffsetID = tOriginalOffsetID;

        return aOffsetID;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::depends_on_advs() const
    {
        return mParameters.mDiscretizationIndex > -1 or get_analytic_perturbation_user_defined != nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        // update the perturbation fields with the new mesh
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->reset_nodal_data( aInterpolationMesh );
        }

        // update the stored mtk interpolation mesh with the new mesh
        mMesh = aInterpolationMesh;

        if ( !mBasesComputed and this->depends_on_advs() )
        {
            this->update_vertex_basis_data();
            mBasesComputed = true;
        }
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

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_facet_center( uint aFacetIndex )
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

    bool Surface_Mesh_Geometry::facet_vertex_depends_on_advs( uint aFacetVertexIndex )
    {
        if ( get_analytic_perturbation_user_defined != nullptr )
        {
            return true;
        }
        else
        {
            // Get the factor that scales this vertex's movement
            Vector< real > tFactor = get_discretization_factor_user_defined == nullptr ? Vector< real >( Surface_Mesh::get_spatial_dimension(), 1.0 ) : get_discretization_factor_user_defined( Surface_Mesh::get_vertex_coordinates( aFacetVertexIndex ) );

            // Return true if this surface mesh can move, its movement was not fixed in all directions by the user, and it lies in the Lagrange mesh domain
            return this->depends_on_advs()
               and mVertexBackgroundElements( aFacetVertexIndex ) != nullptr    //
               and std::any_of( tFactor.cbegin(), tFactor.cend(), []( real tFac ) { return tFac != 0.0; } );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_dvertex_dadv( uint aFacetVertexIndex )
    {
        uint tDim = Surface_Mesh::get_spatial_dimension();

        // Initialize sensitivity matrix
        Matrix< DDRMat > tVertexSensitivity;

        if ( this->intended_discretization() )
        {
            // Determine which directions the vertex can move in
            Vector< bool > tVertexDependsOnADVs( tDim );

            Vector< real > tFactor              = get_discretization_factor_user_defined == nullptr ? Vector< real >( tDim, 1.0 ) : get_discretization_factor_user_defined( Surface_Mesh::get_vertex_coordinates( aFacetVertexIndex ) );
            uint           tNumDimsDependOnADVs = std::count_if( tFactor.cbegin(), tFactor.cend(), []( real tFac ) { return tFac != 0.0; } );

            // Get the vertex indices and coordinates of the background element
            Matrix< DDRMat >   tVertexCoordinates = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
            Matrix< IndexMat > tVertexIndices     = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

            // Loop over background nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
            {
                // Get length before adding sensitivities for this node
                uint tNumVertexSensitivities = tVertexSensitivity.n_cols();

                bool tVertexSensitivitySizeDetermined = false;
                uint tDimensionSensitivitiesAdded     = 0;

                // Get the sensitivity factor of the node in this direction
                Vector< real > tFactor = get_discretization_factor_user_defined == nullptr ? Vector< real >( tDim ) : get_discretization_factor_user_defined( this->get_vertex_coordinates( aFacetVertexIndex ) );
                // Loop over spatial dimension
                for ( uint iDimensionIndex = 0; iDimensionIndex < tDim; iDimensionIndex++ )
                {
                    // Check that the vertex depends on ADVs in this direction
                    if ( tFactor( iDimensionIndex ) != 0.0 )
                    {
                        Matrix< DDRMat > tNodeSensitivity = tFactor( iDimensionIndex ) * mVertexBases( iNodeIndex, aFacetVertexIndex ) * mPerturbationFields( iDimensionIndex )->get_dfield_dadvs( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

                        // set size of sensitivity matrix
                        if ( not tVertexSensitivitySizeDetermined )
                        {
                            tVertexSensitivity.resize( tDim, tNumVertexSensitivities + tNumDimsDependOnADVs * tNodeSensitivity.numel() );
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
        else if ( get_analytic_perturbation_user_defined != nullptr )
        {
            get_perturbation_sensitivity_user_defined( this->get_vertex_coordinates( aFacetVertexIndex ), mADVHandler.get_values(), tVertexSensitivity );    // BRENDAN NEED TO PASS ADVS TO FUNCTION
        }

        return tVertexSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Surface_Mesh_Geometry::get_vertex_adv_ids( uint aFacetVertexIndex )
    {
        if ( this->intended_discretization() )
        {
            uint tDim = Surface_Mesh::get_spatial_dimension();

            // Determine which directions the vertex can move in
            Vector< bool > tVertexDependsOnADVs( tDim );
            Vector< real > tFactor              = get_discretization_factor_user_defined == nullptr ? Vector< real >( tDim, 1.0 ) : get_discretization_factor_user_defined( Surface_Mesh::get_vertex_coordinates( aFacetVertexIndex ) );
            uint           tNumDimsDependOnADVs = std::count_if( tFactor.cbegin(), tFactor.cend(), []( real tFac ) { return tFac != 0.0; } );

            // Initialize matrix to be filled
            Vector< sint > tVertexADVIds;

            // Get the vertex indices and coordinates of the background element
            Matrix< DDRMat >   tVertexCoordinates = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
            Matrix< IndexMat > tVertexIndices     = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

            // Loop over background nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
            {
                // Get the ADV IDs for this node
                Vector< sint > tNodeIDs = this->get_determining_adv_ids( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

                // Join the ADV IDs to the output
                // Get the original length
                uint tIDLength = tVertexADVIds.size();

                // Resize to add new ADV IDs
                tVertexADVIds.resize( tVertexADVIds.size() + ( tNodeIDs.size() * tNumDimsDependOnADVs ) / tDim );

                // Join the IDs
                uint tADVsAdded = 0;
                for ( uint iDimension = 0; iDimension < tDim; iDimension++ )
                {
                    if ( tVertexDependsOnADVs( iDimension ) )
                    {
                        for ( uint iADVIndex = 0; iADVIndex < tNodeIDs.size() / tDim; iADVIndex++ )
                        {
                            tVertexADVIds( tIDLength + tADVsAdded ) = tNodeIDs( ( iDimension * tNodeIDs.size() ) / tDim + iADVIndex );
                            tADVsAdded++;
                        }
                    }
                }
            }

            return tVertexADVIds;
        }
        else
        {
            return mADVHandler.get_determining_adv_ids();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Surface_Mesh_Geometry::get_determining_adv_ids(
            uint                    aNodeIndex,
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

    void Surface_Mesh_Geometry::get_design_info(
            uint                    aNodeIndex,
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

    bool Surface_Mesh_Geometry::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Surface_Mesh_Geometry::get_discretization_mesh_index()
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::update_dependencies( const Vector< std::shared_ptr< Design > >& aAllUpdatedDesigns )
    {
    }


    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< real > >
    Surface_Mesh_Geometry::determine_mtk_cell_bounding_box( mtk::Cell* aElement )
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

    mtk::Cell*
    Surface_Mesh_Geometry::find_background_element_from_global_coordinates( const Matrix< DDRMat >& aCoordinate )
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
            mtk::Cell*              aBackgroundElement,
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
        uint tDim = Surface_Mesh::get_spatial_dimension();

        // Set size if it has not been set already
        if ( mVertexBases.n_cols() != Surface_Mesh::get_number_of_vertices() )
        {
            mVertexBases.resize( mMesh->get_mtk_cell( 0 ).get_number_of_vertices(), Surface_Mesh::get_number_of_vertices() );
            mVertexBackgroundElements.resize( Surface_Mesh::get_number_of_vertices() );
        }

        // Compute the bases for all facet vertices
        for ( uint iVertexIndex = 0; iVertexIndex < Surface_Mesh::get_number_of_vertices(); iVertexIndex++ )
        {
            // Get this vertex's coordinates
            Matrix< DDRMat > tVertexCoordinates = Surface_Mesh::get_vertex_coordinates( iVertexIndex );

            Matrix< DDRMat > tVertexParametricCoordinates( tDim, 1 );

            // Determine which element this vertex lies in, will be the same for every field
            mVertexBackgroundElements( iVertexIndex ) = this->find_background_element_from_global_coordinates( Surface_Mesh::get_vertex_coordinates( iVertexIndex ) );

            // check if the vertex is inside the mesh domain
            if ( mVertexBackgroundElements( iVertexIndex ) != nullptr )
            {
                // Get the bounding box for this element
                Vector< Vector< real > > tElementBoundingBox = this->determine_mtk_cell_bounding_box( mVertexBackgroundElements( iVertexIndex ) );

                // determine the local coordinates of the vertex inside the mtk::Cell
                for ( uint iDimensionIndex = 0; iDimensionIndex < tDim; iDimensionIndex++ )
                {
                    tVertexParametricCoordinates( iDimensionIndex, 0 ) = 2.0 * ( tVertexCoordinates( iDimensionIndex ) - tElementBoundingBox( 0 )( iDimensionIndex ) )
                                                                               / ( tElementBoundingBox( 1 )( iDimensionIndex ) - tElementBoundingBox( 0 )( iDimensionIndex ) )
                                                                       - 1.0;
                }

                // Get the basis function values at the vertex location
                Matrix< DDRMat > tBasis = this->compute_vertex_basis( mVertexBackgroundElements( iVertexIndex ), tVertexParametricCoordinates );
                mVertexBases.set_column( iVertexIndex, trans( tBasis ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::interpolate_perturbation_from_background_element(
            mtk::Cell* aBackgroundElement,
            uint       aFieldIndex,
            uint       aFacetVertexIndex )
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
            tPerturbation += mVertexBases.get_column( aFacetVertexIndex )( iBackgroundNodeIndex ) * mPerturbationFields( aFieldIndex )->get_field_value( tVertexIndices( iBackgroundNodeIndex ), { {} } );
        }

        return tPerturbation;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen