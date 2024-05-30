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
#include "cl_MTK_Enums.hpp"

#include "cl_SOL_Dist_Map.hpp"
namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Parameters::Surface_Mesh_Parameters( const Parameter_List& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mFilePath( aParameterList.get< std::string >( "file_path" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
            , mVertexFactorFunctionName( aParameterList.get< std::string >( "vertex_factor_function_name" ) )
    {
        string_to_cell( aParameterList.get< std::string >( "offset" ), mOffsets );
        string_to_cell( aParameterList.get< std::string >( "scale" ), mScale );
        string_to_cell( aParameterList.get< std::string >( "adv_indices" ), mADVIndices );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry(
            mtk::Mesh*                    aMesh,
            Matrix< DDRMat >              aADVs,
            Surface_Mesh_Parameters       aParameters,
            Node_Manager&                 aNodeManager,
            std::shared_ptr< Library_IO > aLibrary )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Object( aParameters.mFilePath, aParameters.mIntersectionTolerance, aParameters.mOffsets, aParameters.mScale )
            , mParameters( aParameters )
            , mNodeManager( &aNodeManager )
            , mMesh( aMesh )
            , mPerturbationFields( 0 )
            , mVertexBases( 0, 0 )
            , mVertexBackgroundElements( 0 )
            , mOriginalVertexCoordinates( Object::mVertices.size(), Object::mDimension )
    {
        // Check the correct number of ADVs are provided (either as many as dimensions or zero)
        MORIS_ERROR( aParameters.mADVIndices.size() == Object::mDimension
                             or aParameters.mADVIndices.size() == 0,
                "GEN - %ld ADV Indices provided to surface mesh. Should be %d or 0 (for now).",
                aParameters.mADVIndices.size(),
                Object::mDimension );

        // parse the file path and extract the file name
        mName = aParameters.mFilePath.substr( aParameters.mFilePath.find_last_of( "/" ) + 1,
                aParameters.mFilePath.find_last_of( "." ) - aParameters.mFilePath.find_last_of( "/" ) - 1 );

        // copy original vertex coords from Object base class for perturbation
        for ( uint iVertexIndex = 0; iVertexIndex < Object::mVertices.size(); iVertexIndex++ )
        {
            for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
            {
                mOriginalVertexCoordinates( iVertexIndex )( iDimension ) = Object::mVertices( iVertexIndex )->get_coord( iDimension );
            }
        }

        // If this surface mesh is being optimized, construct fields,
        // store original vertex coordinates, determine which facet vertices are fixed, and compute the bases for all vertices
        if ( this->depends_on_advs() )
        {
            // STEP 1: Determine which facet vertices are fixed
            if ( not mParameters.mVertexFactorFunctionName.empty() )
            {
                // Get pointer to function
                mVertexFactorFunction = aLibrary->load_function< VERTEX_FACTOR_FUNCTION >( mParameters.mVertexFactorFunctionName );
            }

            // STEP 2: Construct perturbation fields
            mPerturbationFields.resize( Object::mDimension );

            // build perturbation fields
            for ( uint iFieldIndex = 0; iFieldIndex < Object::get_dimension(); iFieldIndex++ )
            {
                Matrix< DDUMat > tADVIndices;
                Matrix< DDRMat > tConstants;
                Matrix< DDUMat > tFieldVariableIndices;
                // construct field to be discretized into a bspline field eventually
                if ( aParameters.mDiscretizationIndex > -1 )
                {
                    // allocate space for contants information
                    tConstants.resize( 1, 1 );

                    // set constants
                    tConstants( 0, 0 ) = 0.0;
                }
                // construct constant field with an ADV to rigidly displace the surface mesh
                else
                {
                    // allocate space for ADV information
                    tADVIndices.resize( 1, 1 );
                    tFieldVariableIndices.resize( 1, 1 );

                    // set ADV information
                    tADVIndices( 0, 0 )           = mParameters.mADVIndices( iFieldIndex );
                    tFieldVariableIndices( 0, 0 ) = 0;
                }
                // Build field
                mPerturbationFields( iFieldIndex ) = std::make_shared< Constant_Field >(
                        aADVs,
                        tFieldVariableIndices,
                        tADVIndices,
                        tConstants,
                        mName + "_PERT_" + std::to_string( iFieldIndex ) );
            }
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
        sdf::Object_Region tRegion = raycast_point( *this, aNodeCoordinates );

        switch ( tRegion )
        {
            case sdf::Object_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
            }
            case sdf::Object_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
            }
            case sdf::Object_Region::INTERFACE:
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
        sdf::Facet* tParentFacet     = nullptr;
        real        tLocalCoordinate = this->compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode, tParentFacet );

        MORIS_ERROR( tParentFacet == nullptr or ( tLocalCoordinate < 1 + this->get_intersection_tolerance() and tLocalCoordinate > -1 - this->get_intersection_tolerance() ),
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
            sdf::Facet*&                      aParentFacet )
    {
        // transform the interface geometry to local coordinates
        this->transform_surface_mesh_to_local_coordinate( aFirstParentNode, aSecondParentNode );

        // Compute the distance to the facets
        Matrix< DDRMat > tCastPoint( Object::mDimension, 1 );
        tCastPoint.fill( 0.0 );
        Vector< sdf::Facet* > tIntersectionFacets;
        Vector< real >        tLocalCoordinate = sdf::compute_distance_to_facets( *this, tCastPoint, 0, tIntersectionFacets );

        // shift local coordinate to be between -1 and 1
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            tLocalCoordinate( iIntersection ) += -1.0;
        }

        // reset the object to the vertex coordinates at the current design iteration
        this->reset_coordinates();

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
            aParentFacet = nullptr;
            return MORIS_REAL_MAX;
        }
        else if ( tNumberOfParentEdgeIntersections > 1 )
        {
            MORIS_LOG_WARNING( "Multiple facet intersections detected along parent edge. Using first intersection." );
        }

        // Set return values for intersection location and associated facet
        MORIS_ASSERT( tIntersectionFacets.size() == tLocalCoordinate.size(), "Inconsistent size of facet vector (size %lu) and local coordinate vector (size %lu)", tIntersectionFacets.size(), tLocalCoordinate.size() );
        aParentFacet = tIntersectionFacets( 0 );
        return tLocalCoordinate( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::transform_surface_mesh_to_local_coordinate(
            const Parent_Node& aFirstParentNode,
            const Parent_Node& aSecondParentNode )
    {
        // STEP 1: shift the object so the first parent is at the origin
        Matrix< DDRMat > tFirstParentNodeGlobalCoordinates = aFirstParentNode.get_global_coordinates();
        Vector< real >   tShift( Object::mDimension );
        MORIS_ASSERT( tFirstParentNodeGlobalCoordinates.numel() == tShift.size(),
                "Intersection Node Surface Mesh::transform_mesh_to_local_coordinates() inconsistent parent node and interface geometry dimensions." );
        for ( uint iCoord = 0; iCoord < tShift.size(); iCoord++ )
        {
            tShift( iCoord ) = -1.0 * tFirstParentNodeGlobalCoordinates( iCoord );
        }
        this->shift( tShift );

        // STEP 2: rotate the object
        // get unit axis to rotate to
        Matrix< DDRMat > tParentVector = aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates();

        // augment with zero if 2D
        if ( tParentVector.numel() == 2 )
        {
            tParentVector.reshape( 3, 1 );
            tParentVector( 2, 0 ) = 0.0;
        }

        real tParentVectorNorm = norm( tParentVector );

        tParentVector = tParentVector / tParentVectorNorm;

        // create vector orthogonal to parent vector and cast axis
        // in 2D, this vector is the z axis
        Matrix< DDRMat > tRotationMatrix( 3, 1 );
        Matrix< DDRMat > tCastAxis = { { 1.0 }, { 0.0 }, { 0.0 } };

        if ( norm( tParentVector + tCastAxis ) < this->get_intersection_tolerance() )
        {
            tRotationMatrix = { { -1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
        }
        else
        {
            Matrix< DDRMat > tAntisymmetricCrossProduct = { { 0, tParentVector( 1 ), tParentVector( 2 ) },
                { -tParentVector( 1 ), 0.0, 0.0 },
                { -tParentVector( 2 ), 0.0, 0.0 } };

            Matrix< DDRMat > tAntisymmetricCrossProductSquared = { { -std::pow( tParentVector( 1 ), 2 ) - std::pow( tParentVector( 2 ), 2 ), 0.0, 0.0 },
                { 0.0, -std::pow( tParentVector( 1 ), 2 ), -tParentVector( 1 ) * tParentVector( 2 ) },
                { 0.0, -tParentVector( 1 ) * tParentVector( 2 ), -std::pow( tParentVector( 2 ), 2 ) } };

            tRotationMatrix = eye( 3, 3 ) + tAntisymmetricCrossProduct + ( 1 / ( 1 + tParentVector( 0 ) ) ) * tAntisymmetricCrossProductSquared;
        }

        // check that the rotation matrix is correct by ensuring the parent vector was rotated to the x axis
        MORIS_ASSERT( norm( tRotationMatrix * tParentVector - tCastAxis ) < this->get_intersection_tolerance(),
                "Rotation matrix should rotate the parent vector to the x axis." );

        // trim the transformation matrix if 2D
        if ( Object::mDimension == 2 )
        {
            tRotationMatrix.resize( 2, 2 );
        }

        // rotate the object
        this->rotate( tRotationMatrix );

        // STEP 3: scale the object
        Vector< real > tScaling( Object::mDimension, 2.0 / tParentVectorNorm );
        this->scale( tScaling );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // STEP 1: Import advs to field
            mPerturbationFields( iFieldIndex )->import_advs( aOwnedADVs );

            // STEP 2: Apply field value to surface mesh nodes
            for ( uint iVertexIndex = 0; iVertexIndex < Object::mVertices.size(); iVertexIndex++ )
            {
                // Get the factor that scales this vertex's movement
                real tFactor = mVertexFactorFunction == nullptr ? 1.0 : mVertexFactorFunction( iVertexIndex, Object::mVertices( iVertexIndex )->get_coords(), iFieldIndex );

                // Move vertex if needed
                if ( tFactor != 0.0 and mVertexBackgroundElements( iVertexIndex ) != nullptr )
                {
                    // Interpolate the bspline field value at the facet vertex location
                    real tInterpolatedPerturbation = this->interpolate_perturbation_from_background_element(
                            mVertexBackgroundElements( iVertexIndex ),
                            iFieldIndex,
                            iVertexIndex );

                    // Displace the vertex by the total perturbation
                    Object::mVertices( iVertexIndex )->set_node_coord( mOriginalVertexCoordinates( iVertexIndex )( iFieldIndex ) + tFactor * tInterpolatedPerturbation, iFieldIndex );
                }
            }
        }

        // STEP 3: Update all facet data
        this->update_all_facets();

        // this->write_to_obj_file( mName + "_" + std::to_string( mIteration ) + ".txt" );
        // mIteration++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::set_advs( sol::Dist_Vector* aADVs )
    {
        // Have each field import the advs
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->set_advs( aADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Surface_Mesh_Geometry::append_adv_info(
            mtk::Interpolation_Mesh* aMesh,
            Matrix< DDSMat >&        aOwnedADVIds,
            Matrix< IdMat >&         aOwnedijklIDs,
            sint                     aOffsetID,
            Matrix< DDRMat >&        aLowerBounds,
            Matrix< DDRMat >&        aUpperBounds )
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

    bool
    Surface_Mesh_Geometry::depends_on_advs() const
    {
        return mParameters.mADVIndices.size() > 0 or mParameters.mDiscretizationIndex > -1;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        // update the perturbation fields with the new mesh
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            mPerturbationFields( iFieldIndex )->reset_nodal_data( aInterpolationMesh );
        }

        // update the stored mtk interpolation mesh with the new mesh
        mMesh = aInterpolationMesh;

        if ( !mBasesComputed )
        {
            this->update_vertex_basis_data();
            mBasesComputed = true;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::discretize(
            mtk::Mesh_Pair    aMeshPair,
            sol::Dist_Vector* aOwnedADVs )
    {
        MORIS_ASSERT( Design::mSharedADVIDs.size() == Object::mDimension or Design::mSharedADVIDs.size() == 0,
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
        for ( uint iFieldIndex = 0; iFieldIndex < Object::mDimension; iFieldIndex++ )
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

    bool
    Surface_Mesh_Geometry::facet_vertex_depends_on_advs( uint aFacetVertexIndex )
    {
        // check if the whole surface mesh depends on advs or the vertex is outside the domain
        if ( !this->depends_on_advs() or mVertexBackgroundElements( aFacetVertexIndex ) == nullptr )
        {
            return false;
        }

        // check if this node was specified to be fixed
        if ( mVertexFactorFunction == nullptr )
        {
            return true;
        }
        else
        {
            bool tVertexDependsOnADVs = false;
            for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
            {
                if ( mVertexFactorFunction( aFacetVertexIndex, Object::mVertices( aFacetVertexIndex )->get_coords(), iDimension ) != 0.0 )
                {
                    tVertexDependsOnADVs = true;
                    break;
                }
            }

            return tVertexDependsOnADVs;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Surface_Mesh_Geometry::get_dvertex_dadv( uint aFacetVertexIndex )
    {
        // Determine which directions the vertex can move in
        Vector< bool > tVertexDependsOnADVs( Object::mDimension );
        uint           tNumDimsDependOnADVs = 0;
        for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
        {
            tVertexDependsOnADVs( iDimension ) = mVertexFactorFunction == nullptr ? true : mVertexFactorFunction( aFacetVertexIndex, Object::mVertices( aFacetVertexIndex )->get_coords(), iDimension ) != 0.0;
            if ( tVertexDependsOnADVs( iDimension ) )
            {
                tNumDimsDependOnADVs++;
            }
        }

        // Get the vertex indices and coordinates of the background element
        Matrix< DDRMat >   tVertexCoordinates = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
        Matrix< IndexMat > tVertexIndices     = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

        // Initialize sensitivity matrix
        Matrix< DDRMat > tVertexSensitivity;

        // Loop over background nodes
        for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
        {
            // Get length before adding sensitivities for this node
            uint tNumVertexSensitivities = tVertexSensitivity.n_cols();

            bool tVertexSensitivitySizeDetermined = false;
            uint tDimensionSensitivitiesAdded = 0;

            // Loop over spatial dimension
            for ( uint iDimensionIndex = 0; iDimensionIndex < Object::mDimension; iDimensionIndex++ )
            {
                // Get the sensitivity factor of the node in this direction
                real tFactor = mVertexFactorFunction == nullptr ? 1.0 : mVertexFactorFunction( aFacetVertexIndex, Object::mVertices( aFacetVertexIndex )->get_coords(), iDimensionIndex );

                // Check that the vertex depends on ADVs in this direction
                if ( tFactor != 0.0 )
                {
                    Matrix< DDRMat > tNodeSensitivity = tFactor * mVertexBases( iNodeIndex, aFacetVertexIndex ) * mPerturbationFields( iDimensionIndex )->get_dfield_dadvs( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

                    // set size of sensitivity matrix
                    if ( not tVertexSensitivitySizeDetermined )
                    {
                        tVertexSensitivity.resize( Object::mDimension, tNumVertexSensitivities + tNumDimsDependOnADVs * tNodeSensitivity.numel()  );
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

        return tVertexSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Surface_Mesh_Geometry::get_vertex_adv_ids( uint aFacetVertexIndex )
    {
        // Determine which directions the vertex can move in
        Vector< bool > tVertexDependsOnADVs( Object::mDimension );
        uint           tNumDimsDependOnADVs = 0;
        for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
        {
            tVertexDependsOnADVs( iDimension ) = mVertexFactorFunction == nullptr ? true : mVertexFactorFunction( aFacetVertexIndex, Object::mVertices( aFacetVertexIndex )->get_coords(), iDimension ) != 0.0;
            if ( tVertexDependsOnADVs( iDimension ) )
            {
                tNumDimsDependOnADVs++;
            }
        }

        // Initialize matrix to be filled
        Matrix< DDSMat > tVertexADVIds;

        // Get the vertex indices and coordinates of the background element
        Matrix< DDRMat >   tVertexCoordinates = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_coords();
        Matrix< IndexMat > tVertexIndices     = mVertexBackgroundElements( aFacetVertexIndex )->get_vertex_inds();

        // Loop over background nodes
        for ( uint iNodeIndex = 0; iNodeIndex < tVertexCoordinates.n_rows(); iNodeIndex++ )
        {
            // Get the ADV IDs for this node
            Matrix< DDSMat > tNodeIDs = this->get_determining_adv_ids( tVertexIndices( iNodeIndex ), tVertexCoordinates.get_row( iNodeIndex ) );

            // Join the ADV IDs to the output
            // Get the original length
            uint tIDLength = tVertexADVIds.length();

            // Resize to add new ADV IDs
            tVertexADVIds.resize( 1, tVertexADVIds.length() + ( tNodeIDs.length() * tNumDimsDependOnADVs ) / Object::mDimension );

            // Join the IDs
            uint tADVsAdded = 0;
            for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
            {
                if ( tVertexDependsOnADVs( iDimension ) )
                {
                    for ( uint iADVIndex = 0; iADVIndex < tNodeIDs.length() / Object::mDimension; iADVIndex++ )
                    {
                        tVertexADVIds( tIDLength + tADVsAdded ) = tNodeIDs( ( iDimension * tNodeIDs.length() ) / Object::mDimension + iADVIndex );
                        tADVsAdded++;
                    }
                }
            }
        }

        return tVertexADVIds;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Surface_Mesh_Geometry::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        Matrix< DDSMat > tADVIDs;
        for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
        {
            // Get the ADV IDs for this field
            Matrix< DDSMat > tFieldADVIDs;
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
            uint tJoinedIDLength = tADVIDs.n_cols();
            tADVIDs.resize( 1, tJoinedIDLength + tFieldADVIDs.length() );

            // Join IDs
            for ( uint tAddedSensitivity = 0; tAddedSensitivity < tFieldADVIDs.length(); tAddedSensitivity++ )
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

    Vector< Vector< real > >
    Surface_Mesh_Geometry::determine_mtk_cell_bounding_box( mtk::Cell* aElement )
    {
        // Initialize a bounding box
        Vector< Vector< real > > tBoundingBox( 2, Vector< real >( Object::mDimension ) );

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
        // Loop through each mtk::Cell
        for ( uint iCellIndex = 0; iCellIndex < mMesh->get_num_elems(); iCellIndex++ )
        {
            // get this Cell's bounding box BRENDAN CAN MAKE THIS EASIER ASSUMING COORDINATES IN THE CORRECT ORDER
            Vector< Vector< real > > tBoundingBox = this->determine_mtk_cell_bounding_box( &mMesh->get_mtk_cell( iCellIndex ) );

            // assume the point is in this bounding box
            bool tCoordinateInCell = true;

            // Loop over the bounding box and determine if this is true
            for ( uint iDimensionIndex = 0; iDimensionIndex < Object::mDimension; iDimensionIndex++ )
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

    void
    Surface_Mesh_Geometry::update_vertex_basis_data()
    {
        // Set size if it has not been set already
        if ( mVertexBases.n_cols() != Object::mVertices.size() )
        {
            mVertexBases.resize( mMesh->get_mtk_cell( 0 ).get_number_of_vertices(), Object::mVertices.size() );
            mVertexBackgroundElements.resize( Object::mVertices.size() );
        }

        // Compute the bases for all facet vertices
        for ( uint iVertexIndex = 0; iVertexIndex < Object::mVertices.size(); iVertexIndex++ )
        {
            Matrix< DDRMat > tVertexParametricCoordinates( Object::mDimension, 1 );

            // Determine which element this vertex lies in, will be the same for every field
            mVertexBackgroundElements( iVertexIndex ) = this->find_background_element_from_global_coordinates( Object::mVertices( iVertexIndex )->get_coords() );

            // check if the vertex is inside the mesh domain
            if ( mVertexBackgroundElements( iVertexIndex ) != nullptr )
            {
                // Get the bounding box for this element
                Vector< Vector< real > > tElementBoundingBox = this->determine_mtk_cell_bounding_box( mVertexBackgroundElements( iVertexIndex ) );

                // determine the local coordinates of the vertex inside the mtk::Cell
                for ( uint iDimensionIndex = 0; iDimensionIndex < Object::mDimension; iDimensionIndex++ )
                {
                    tVertexParametricCoordinates( iDimensionIndex, 0 ) = 2.0 * ( Object::mVertices( iVertexIndex )->get_coord( iDimensionIndex ) - tElementBoundingBox( 0 )( iDimensionIndex ) )
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

        // Vector< mtk::Vertex_Interpolation* > tBsplineToLagrange = aBackgroundElement->get_vertex_interpolations( mParameters.mDiscretizationIndex ); BRENDAN this is how to get T matrices for a cell, idk if needed

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