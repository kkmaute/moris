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

#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Enums.hpp"
namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Parameters::Surface_Mesh_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mFilePath( aParameterList.get< std::string >( "file_path" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
    {
        string_to_cell( aParameterList.get< std::string >( "offset" ), mOffsets );
        string_to_cell( aParameterList.get< std::string >( "scale" ), mScale );
        string_to_cell( aParameterList.get< std::string >( "adv_indices" ), mADVIndices );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry( mtk::Mesh* aMesh, Matrix< DDRMat > aADVs, Surface_Mesh_Parameters aParameters )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Object( aParameters.mFilePath, aParameters.mIntersectionTolerance, aParameters.mOffsets, aParameters.mScale )
            , mParameters( aParameters )
            , mMesh( aMesh )
            , mOriginalVertexCoordinates( 0,0 )
            , mPerturbationFields( 0 )
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


        // If this surface mesh is being optimized, construct fields and store original vertex coordinates
        if ( aParameters.mADVIndices.size() > 0 or aParameters.mDiscretizationIndex > -1 )
        {
            mOriginalVertexCoordinates.resize( Object::mVertices.size(), Object::mDimension );
            // copy original vertex coords from Object base class for perturbation
            for ( uint iVertexIndex = 0; iVertexIndex < Object::mVertices.size(); iVertexIndex++ )
            {
                for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
                {
                    mOriginalVertexCoordinates( iVertexIndex )( iDimension ) = Object::mVertices( iVertexIndex )->get_coord( iDimension );
                }
            }

            // Allocate memory for perturbation fields
            mPerturbationFields.resize( Object::mDimension );

            // build perturbation fields
            Matrix< DDUMat > tFieldVariableIndices = { { 0 } };
            for ( uint iFieldIndex = 0; iFieldIndex < Object::get_dimension(); iFieldIndex++ )
            {
                Matrix< DDUMat > tADVIndices;
                Matrix< DDRMat > tADVs;
                Matrix< DDRMat > tConstants;
                // construct field to be discretized into a bspline field eventually
                if ( aParameters.mDiscretizationIndex > -1 )
                {
                    tConstants.resize( 1, 1 );
                    tConstants( 0, 0 ) = 1.0;
                }
                // construct constant field with an ADV to rigidly displace the surface mesh
                else
                {
                    tADVIndices.resize( 1, 1 );
                    tADVs.resize( 1, 1 );
                    tADVIndices = { { iFieldIndex } };
                    tADVs( 0, 0 ) = mParameters.mADVIndices( iFieldIndex );
                }
                // Build field
                mPerturbationFields( iFieldIndex ) = std::make_shared< Constant_Field >(
                        tADVs,
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
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        // Create linear intersection node
        return new Intersection_Node_Surface_Mesh(
                aNodeIndex,
                aBackgroundNodes,
                aFirstParentNode,
                aSecondParentNode,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Surface_Mesh_Geometry::compute_intersection_local_coordinate(
            const Cell< Node* >& aBackgroundNodes,
            const Parent_Node&   aFirstParentNode,
            const Parent_Node&   aSecondParentNode )
    {
        // transform the interface geometry to local coordinates
        this->transform_surface_mesh_to_local_coordinate( aFirstParentNode, aSecondParentNode );

        // Compute the distance to the facets
        Matrix< DDRMat > tCastPoint( Object::mDimension, 1 );
        tCastPoint.fill( 0.0 );
        Cell< real > tLocalCoordinate = sdf::compute_distance_to_facets( *this, tCastPoint, 0 );

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
        if ( tLocalCoordinate.size() == 0 or tNumberOfParentEdgeIntersections > 1 )
        {
            return MORIS_REAL_MAX;
        }

        return tLocalCoordinate( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::transform_surface_mesh_to_local_coordinate(
            const Parent_Node& aFirstParentNode,
            const Parent_Node& aSecondParentNode )
    {
        // step 1: shift the object so the first parent is at the origin
        Matrix< DDRMat > tFirstParentNodeGlobalCoordinates = aFirstParentNode.get_global_coordinates();
        Cell< real >     tShift( Object::mDimension );
        MORIS_ASSERT( tFirstParentNodeGlobalCoordinates.numel() == tShift.size(),
                "Intersection Node Surface Mesh::transform_mesh_to_local_coordinates() inconsistent parent node and interface geometry dimensions." );
        for ( uint iCoord = 0; iCoord < tShift.size(); iCoord++ )
        {
            tShift( iCoord ) = -1.0 * tFirstParentNodeGlobalCoordinates( iCoord );
        }
        this->shift( tShift );

        // step 2: rotate the object
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

        // step 3: scale the object
        Cell< real > tScaling( Object::mDimension, 2.0 / tParentVectorNorm );
        this->scale( tScaling );
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        // TODO BRENDAN: maybe?
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // Have each field import the advs
        for ( auto tPerturbationField : mPerturbationFields )
        {
            tPerturbationField->import_advs( aOwnedADVs );
        }

        // update the surface mesh with the new field data
        // FIXME: Put this in the above loop so it doesn't get called unless necessary
        Cell< Cell< real > > tElementBoundingBox( 2, Object::mDimension );
        Matrix< DDRMat >     tVertexParametricCoordinates( Object::mDimension, 1 );
        for ( uint iVertexIndex = 0; iVertexIndex < Object::mVertices.size(); iVertexIndex++ )
        {
            // Initialize perturbation vector
            Cell< real > tNewVertexCoordinates = mOriginalVertexCoordinates( iVertexIndex );

            // Determine which element this vertex lies in, will be the same for every field
            int tElementIndex = this->find_background_element_from_global_coordinates(
                    Object::mVertices( iVertexIndex )->get_coords(),
                    tElementBoundingBox );

            if ( tElementIndex != -1 )
            {    // determine the local coordinates of the vertex inside the mtk::Cell
                for ( uint iDimensionIndex = 0; iDimensionIndex < Object::mDimension; iDimensionIndex++ )
                {
                    tVertexParametricCoordinates( iDimensionIndex, 0 ) = ( Object::mVertices( iVertexIndex )->get_coord( iDimensionIndex )
                                                                                 - tElementBoundingBox( 0 )( iDimensionIndex ) )
                                                                       / ( tElementBoundingBox( 1 )( iDimensionIndex )
                                                                               - tElementBoundingBox( 0 )( iDimensionIndex ) );
                }

                // Loop through each field and interpolate its displacement value at the vertex's location
                for ( uint iFieldIndex = 0; iFieldIndex < mPerturbationFields.size(); iFieldIndex++ )
                {
                    // Interpolate the bspline field value at the facet vertex location
                    real tInterpolatedPerturbation = this->interpolate_perturbation_from_background_element( &mMesh->get_mtk_cell( tElementIndex ), iFieldIndex, tVertexParametricCoordinates );

                    // Add the perturbation to the list
                    tNewVertexCoordinates( iFieldIndex ) += tInterpolatedPerturbation;
                }

                // Displace the vertex by the total perturbation
                Object::mVertices( iVertexIndex )->set_node_coords( tNewVertexCoordinates );
            }
        }

        // Update all facet data
        this->Object::update_all_facets();

        // BRENDAN
        std::cout << "surface mesh imported ADVs\n";
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
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::discretize(
            mtk::Mesh_Pair          aMeshPair,
            sol::Dist_Vector*       aOwnedADVs,
            const Matrix< DDSMat >& aSharedADVIds,
            uint                    aADVOffsetID )
    {
        if ( mParameters.mDiscretizationIndex >= 0 )
        {
            // Allocate memory for perturbation fields
            mPerturbationFields.resize( Object::mDimension );
            for ( uint iFieldIndex = 0; iFieldIndex < Object::mDimension; iFieldIndex++ )
            {
                // Create a B-spline field
                mPerturbationFields( iFieldIndex ) = std::make_shared< BSpline_Field >(
                        aMeshPair,
                        aOwnedADVs,
                        aSharedADVIds,
                        aADVOffsetID,
                        mParameters.mDiscretizationIndex,
                        mPerturbationFields( iFieldIndex ) );

                // Set analytic field index, for now
                Field::gDiscretizationIndex = mParameters.mDiscretizationIndex;
            }
        }
        else if ( mParameters.mDiscretizationIndex == -1 )
        {
            // Allocate memory for perturbation fields
            mPerturbationFields.resize( Object::mDimension );
            for ( uint iFieldIndex = 0; iFieldIndex < Object::mDimension; iFieldIndex++ )
            {
                // Just store nodal values
                mPerturbationFields( iFieldIndex ) = std::make_shared< Stored_Field >(
                        aMeshPair.get_interpolation_mesh(),
                        mPerturbationFields( iFieldIndex ) );
            }
        }
        for ( auto iPerturbationField : mPerturbationFields )
        {
            iPerturbationField->mMeshPairForAnalytic = aMeshPair;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs,
            const Matrix< DDSMat >&       aSharedADVIds,
            uint                          aADVOffsetID )
    {
        for ( uint iFieldIndex = 0; iFieldIndex < Object::mDimension; iFieldIndex++ )
        {
            // Allocate memory for perturbation fields
            mPerturbationFields.resize( Object::mDimension );
            if ( mParameters.mDiscretizationIndex >= 0 )
            {
                // Create a B-spline field
                mPerturbationFields( iFieldIndex ) = std::make_shared< BSpline_Field >(
                        aOwnedADVs,
                        aSharedADVIds,
                        aADVOffsetID,
                        mParameters.mDiscretizationIndex,
                        aMTKField,
                        aMeshPair );
            }
            else if ( mParameters.mDiscretizationIndex == -1 )
            {
                // TODO
                MORIS_ERROR( false, "Stored field cannot be remeshed for now" );
            }
            mPerturbationFields( iFieldIndex )->mMeshPairForAnalytic = aMeshPair;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Surface_Mesh_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Cell< real >&           aOutputDesignInfo )
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

    moris_index
    Surface_Mesh_Geometry::find_background_element_from_global_coordinates(
            const Matrix< DDRMat >& aCoordinate,
            Cell< Cell< real > >&   aBoundingBox )
    {
        // Loop through each mtk::Cell
        for ( uint iElementIndex = 0; iElementIndex < mMesh->get_num_elems(); iElementIndex++ )
        {
            // Get the vertices of the mtk::Cell at this index
            Matrix< DDRMat > tCurrentSearchElementVertexCoordinates = mMesh->get_mtk_cell( iElementIndex ).get_vertex_coords();

            // Build bounding box for this mtk::Cell
            // Loop over dimensions
            for ( uint iDimensionIndex = 0; iDimensionIndex < tCurrentSearchElementVertexCoordinates.n_cols(); iDimensionIndex++ )
            {
                // Loop over vertices
                for ( uint iVertexIndex = 0; iVertexIndex < tCurrentSearchElementVertexCoordinates.n_rows(); iVertexIndex++ )
                {
                    // check if the entry is less than the minimum
                    if ( tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex ) < aBoundingBox( 0 )( iDimensionIndex ) )
                    {
                        aBoundingBox( 0 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex );
                    }
                    // check if the entry is greater than the maximum
                    if ( tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex ) > aBoundingBox( 1 )( iDimensionIndex ) )
                    {
                        aBoundingBox( 1 )( iDimensionIndex ) = tCurrentSearchElementVertexCoordinates( iVertexIndex, iDimensionIndex );
                    }
                }
            }

            // Check if the query point is inside the bounding box
            bool tInsideBoundingBox = true;

            // Loop over dimensions
            for ( uint iDimensionIndex = 0; iDimensionIndex < tCurrentSearchElementVertexCoordinates.n_cols(); iDimensionIndex++ )
            {
                // Check if the point is outside the bounding box
                if ( aCoordinate( iDimensionIndex ) < aBoundingBox( 0 )( iDimensionIndex )
                        or aCoordinate( iDimensionIndex ) > aBoundingBox( 1 )( iDimensionIndex ) )
                {
                    // if outside, set flag to false
                    tInsideBoundingBox = false;
                }
            }

            // if the element has been found, return its index
            if ( tInsideBoundingBox )
            {
                return iElementIndex;
            }
        }

        // if no element found, return -1
        return -1;
    }

    real Surface_Mesh_Geometry::interpolate_perturbation_from_background_element(
            mtk::Cell*              aBackgroundElement,
            uint                    aFieldIndex,
            const Matrix< DDRMat >& aParametricCoordinates )
    {
        // number of nodes to be used for interpolation
        uint tNumBases = aBackgroundElement->get_number_of_vertices();

        // build interpolator
        mtk::Interpolation_Function_Factory tFactory;
        mtk::Interpolation_Function_Base*   tInterpolation;

        // create interpolation function based on spatial dimension  of problem
        tInterpolation = tFactory.create_interpolation_function(
                aBackgroundElement->get_geometry_type(),
                mtk::Interpolation_Type::LAGRANGE,
                aBackgroundElement->get_interpolation_order() );

        // compute basis function at the vertices
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( aParametricCoordinates, tBasis );

        // get the indices and coordinates of the background element vertices
        Matrix< IndexMat > tVertexIndices        = aBackgroundElement->get_vertex_inds();
        Matrix< DDRMat >   tAllVertexCoordinates = aBackgroundElement->get_vertex_coords();

        // initialize field value at the node location
        real tPerturbation = 0.0;

        // get perturbation values at the vertices
        for ( uint iBackgroundNodeIndex = 0; iBackgroundNodeIndex < tNumBases; ++iBackgroundNodeIndex )
        {
            Matrix< DDRMat > tVertexCoordinates( &tVertexCoordinates( iBackgroundNodeIndex, 0 ), 0, Object::mDimension );

            // // FIXME: get vertex coordinates for this vertex. This copy is slow to have to do for every surface mesh vertex
            // Matrix< DDRMat > tVertexCoordinates( 1, Object::mDimension );
            // for ( uint iDimension = 0; iDimension < Object::mDimension; iDimension++ )
            // {
            //     tVertexCoordinates( 0, iDimension ) = tAllVertexCoordinates( iBackgroundNodeIndex, iDimension );
            // }

            // add this vertex's field value to the value
            tPerturbation += tBasis( iBackgroundNodeIndex ) * mPerturbationFields( aFieldIndex )->get_field_value( tVertexIndices( iBackgroundNodeIndex ), tVertexCoordinates );
        }

        // clean up
        delete tInterpolation;

        return tPerturbation;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::ge
