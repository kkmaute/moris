/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Delaunay_Subdivision_Interface.cpp
 *
 */

// BRENDAN check for const correctness where possible

#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Delaunay_Subdivision_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Interpolation_Function.hpp"
#include "cl_XTK_Octree_Interface.hpp"

// Fortran function used during the triangulation of 3D physical element. Will return an organized array of the triangulation (number of tetrahedra-by-4 vertices)
extern "C" {
void tetlst_( uint&, uint*, uint*, uint&, uint* );
}

// Fortran function that finds the closest prime number to argument provided. Only used during the 3D triangulation routine
extern "C" {
int prime_( uint& );
}

// Fortran function to perform triangulation of 2D physical element
extern "C" {
void dtris2_( uint&, uint&, real*, uint*, uint&, uint*, uint*, uint*, uint& );
}

// Fortran function to perform triangulation of 3D physical element
extern "C" {
void dtris3_( uint&, uint&, uint&, uint&, real*, uint*, uint&, uint&, uint&, uint&, uint*, uint*, uint*, uint& );
}

namespace moris::xtk
{
    Delaunay_Subdivision_Interface::Delaunay_Subdivision_Interface( Parameter_List& aParameterList )
    {
    }

    //--------------------------------------------------------------------------------------------------

    Vector< moris_index > Delaunay_Subdivision_Interface::get_decomposed_cell_indices()
    {
        return mMeshGenerationData->mDelaunayBgCellInds;
    }

    //--------------------------------------------------------------------------------------------------

    Vector< Vector< moris_index > >
    Delaunay_Subdivision_Interface::triangulation()
    {
        // Reset the number of cells created by this algorithm
        mNumNewCells += mNumTotalCells;
        mNumTotalCells = 0;

        // Get the number of child mesh points - these are all the child meshes, not just the ones with surface points
        uint tNumCMPoints = mMeshGenerationData->mDelaunayBgCellInds.size();

        // Initialize return vector for connectivities
        Vector< Vector< moris_index > > tConnectivities( tNumCMPoints );

        // Allocate variables used for geompack3d call
        uint  tMaxNumTriangles = 24;                                                    // Unused but required by the function
        uint* tStack           = (uint*)alloca( sizeof( uint ) * tMaxNumTriangles );    // Unused but required by the function
        uint  tNumTriangles    = 0;                                                     // Total number of cells after triangulation
        uint  tError           = 0;                                                     // Error flag

        // Loop through all the surface points
        for ( uint iCM = 0; iCM < tNumCMPoints; iCM++ )
        {
            // Get the number of surface points for this child mesh to estimate the number of potential triangles
            uint tNumAllSurfacePoints = mAllSurfacePoints( iCM ).size() / mBackgroundMesh->get_spatial_dim();

            // Check if delaunay triangulation is needed for this cell
            if ( tNumAllSurfacePoints > 0 )
            {
                // Allocate outputs for delaunay triangulation
                Vector< uint > tDelaunayTriangulation( tNumAllSurfacePoints * tMaxNumTriangles );       // Vector containing the connectivity of the triangles
                Vector< uint > tDelaunayTriangulationNBR( tNumAllSurfacePoints * tMaxNumTriangles );    // I don't know what this is for and I don't use it but the function requires it
                Vector< uint > tLocalIndices( tNumAllSurfacePoints );                                   // Indices of the points in the triangulation
                std::iota( tLocalIndices.begin(), tLocalIndices.end(), 1 );

                Vector< real > tSurfacePoints = mAllSurfacePoints( iCM );

                // Perform delaunay triangulation to get the connectivity for the new cells BRENDAN add compiler directive
                dtris2_(
                        tNumAllSurfacePoints,
                        tMaxNumTriangles,
                        &tSurfacePoints( 0 ),
                        &tLocalIndices( 0 ),
                        tNumTriangles,
                        &tDelaunayTriangulation( 0 ),
                        &tDelaunayTriangulationNBR( 0 ),
                        tStack,
                        tError );

                MORIS_ERROR( not tError, "geompack3d reported an error during delaunay triangulation for child mesh %d.", mMeshGenerationData->mDelaunayBgCellInds( iCM ) );

                // Add to the number of new cells created by this algorithm
                mNumTotalCells += tNumTriangles;

                // Annoying conversion to moris_index
                size_t                tSize = tNumTriangles * 3;
                Vector< moris_index > tDelaunayTriangulationIndex( tSize );
                std::transform(
                        tDelaunayTriangulation.begin(),
                        tDelaunayTriangulation.begin() + tSize,
                        tDelaunayTriangulationIndex.begin(),
                        []( uint aValue ) { return static_cast< moris_index >( aValue ); } );

                // Store the connectivity for this child mesh
                tConnectivities( iCM ) = tDelaunayTriangulationIndex;
            }
        }

        // return success if finished
        return tConnectivities;
    }

    //--------------------------------------------------------------------------------------------------

    bool Delaunay_Subdivision_Interface::has_geometric_dependent_vertices() const
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------

    void Delaunay_Subdivision_Interface::perform_impl_vertex_requests(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        // keep track of some useful classes (avoid passing to every function)
        mGeometryEngine     = aMeshGenerator->get_geom_engine();
        mMeshGenerationData = aMeshGenerationData;
        mDecompositionData  = aDecompositionData;
        mCutIntegrationMesh = aCutIntegrationMesh;
        mBackgroundMesh     = aBackgroundMesh;
        mGenerator          = aMeshGenerator;

        // Get the total number of cell in the background mesh as well as the number of delaunay cells
        uint tNumBackgroundCells = mCutIntegrationMesh->get_num_ig_cell_groups();
        uint tNumDelaunayCells   = aMeshGenerationData->mDelaunayBgCellInds.size();

        // Setup decomposition data for this algorithm
        mDecompositionData->mHasSecondaryIdentifier = true;
        mDecompositionData->tCMNewNodeLoc           = Vector< Vector< moris_index > >( tNumBackgroundCells );
        mDecompositionData->tCMNewNodeParamCoord    = Vector< Vector< Matrix< DDRMat > > >( tNumBackgroundCells );

        // Stores the surface points for each child mesh as a flattened matrix
        mAllSurfacePoints.resize( tNumDelaunayCells );

        // Loop through the list of cells with surface points
        for ( uint iCell = 0; iCell < tNumDelaunayCells; iCell++ )
        {
            // Get the index of the cell
            uint tBgCellIndex = aMeshGenerationData->mDelaunayBgCellInds( iCell );

            // Get the intersected cell from the background mesh
            mtk::Cell& tCell = aBackgroundMesh->get_mtk_cell( tBgCellIndex );

            // Get the type of cell (tri, quad, tet, hex)
            mtk::Geometry_Type tCellType = tCell.get_geometry_type();

            // Get the global coordinates of all the nodes of aCell
            Matrix< DDRMat > tBackgroundCellCoords = tCell.get_vertex_coords();
            uint             tNumCornerNodes       = this->get_num_geometric_nodes( tCellType );
            uint             tDim                  = tBackgroundCellCoords.n_cols();

            // Get the surface points for this cell - parametric coordinates
            Matrix< DDRMat > tSurfacePoints = aMeshGenerationData->mDelaunayPoints( tBgCellIndex );

            // Get the associated global coordinates for the surface points
            Matrix< DDRMat > tGlobalSurfacePoints = this->get_surface_point_global_coordinates( tSurfacePoints, tCell );

            // Make requests for these vertices
            for ( uint iPoint = 0; iPoint < tSurfacePoints.n_rows(); iPoint++ )
            {
                // get the parent rank of the new node
                mtk::EntityRank tParentRank = mtk::EntityRank::ELEMENT;

                // get a unique ID for this vertex? BRENDAN: need to figure out how to do this
                moris_index tSecondaryID = std::stoul( std::to_string( tBgCellIndex ) + std::to_string( iPoint ) );

                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                // check if new node for current edge has already been requested ...
                bool tRequestExist = mDecompositionData->request_exists(
                        tBgCellIndex,
                        tSecondaryID,
                        tParentRank,
                        tNewNodeIndexInSubdivision );

                // ... if not, request it
                if ( not tRequestExist )
                {
                    // find out which processor owns parent entity of currently treated element
                    moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tBgCellIndex, tParentRank );

                    // Register new node request
                    tNewNodeIndexInSubdivision = mDecompositionData->register_new_request(
                            tBgCellIndex,
                            tSecondaryID,
                            tOwningProc,
                            (mtk::EntityRank)tParentRank,
                            tGlobalSurfacePoints.get_row( iPoint ) );

                    std::shared_ptr< Matrix< DDRMat > > tNewNodeXiPtr = std::make_shared< Matrix< DDRMat > >( tSurfacePoints.get_row( iPoint ) );

                    // Associate this node with the child mesh
                    mDecompositionData->tCMNewNodeLoc( tBgCellIndex ).push_back( tNewNodeIndexInSubdivision );
                    mDecompositionData->tCMNewNodeParamCoord( tBgCellIndex ).push_back( tSurfacePoints.get_row( iPoint ) );

                    // Register new floating node in GEN
                    mGeometryEngine->create_floating_node(
                            mMeshGenerationData->mDelaunayGeometryIndices( tBgCellIndex )( iPoint ),
                            tCell,
                            tSurfacePoints.get_row( iPoint ),
                            tCellType,
                            tCell.get_interpolation_order() );
                }
            }    // end for: iterate through surface points

            // Compute the total number of points for this child mesh
            uint tNumPointsForDelaunay = tSurfacePoints.n_rows() + tNumCornerNodes;

            // Get all of the points for this child mesh into a single vector - background points
            Vector< real > tSurfacePointsVector( tNumPointsForDelaunay * tDim );
            for ( uint iPoint = 0; iPoint < tNumCornerNodes; iPoint++ )
            {
                for ( uint iDim = 0; iDim < tDim; iDim++ )
                {
                    tSurfacePointsVector( iDim + iPoint * tDim ) = tBackgroundCellCoords( iPoint, iDim );
                }
            }

            // Get all of the points for this child mesh into a single vector - surface points
            for ( uint iPoint = 0; iPoint < tSurfacePoints.n_rows(); iPoint++ )
            {
                for ( uint iDim = 0; iDim < tDim; iDim++ )
                {
                    tSurfacePointsVector( iDim + ( iPoint + tNumCornerNodes ) * tDim ) = tGlobalSurfacePoints( iPoint, iDim );
                }
            }

            // Store the surface points for this child mesh so we can access it to make ig cells later
            mAllSurfacePoints( iCell ) = tSurfacePointsVector;

        }    // end for: iterate through background cells
    }

    //--------------------------------------------------------------------------------------------------

    void Delaunay_Subdivision_Interface::perform_impl_generate_mesh(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        // First, update the newly created nodes in the geometry engine
        for ( uint iNodeRequest = 0; iNodeRequest < aDecompositionData->tNewNodeIndex.size(); iNodeRequest++ )
        {
            moris_index tNodeIndex = aDecompositionData->tNewNodeIndex( iNodeRequest );
            moris_id    tNodeId    = aDecompositionData->tNewNodeId( iNodeRequest );
            moris_index tNodeOwner = aDecompositionData->tNewNodeOwner( iNodeRequest );

            mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
        }

        // Number of background cells
        uint tNumDelaunayCells = aMeshGenerationData->mDelaunayBgCellInds.size();

        // new cell info
        moris::mtk::Cell_Info_Factory            tFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tIgCellInfo = tFactory.create_cell_info_sp( this->get_ig_cell_topology() );

        // Perform a delaunay triangulation on the surface points to get new cell connectivity in the cells local indexing
        Vector< Vector< moris_index > > tCellToVertexConnectivity = this->triangulation();

        // This is the global node indices and is used in the integration mesh generator, so we must set it here
        mNewCellToVertexConnectivity = Vector< Vector< moris::moris_index > >( mNumTotalCells, Vector< moris::moris_index >( 3 ) );

        // Since each child mesh should just be a background cell, there are no cell indices we need to replace
        mNewCellCellIndexToReplace = Vector< moris::moris_index >( mNumTotalCells, MORIS_INDEX_MAX );

        // All the new cells will be the same type
        mNewCellCellInfo = Vector< std::shared_ptr< moris::mtk::Cell_Info > >( mNumTotalCells, tIgCellInfo );

        // Set the correct size for the number of new cells
        mNewCellChildMeshIndex.resize( mNumTotalCells );

        // Counter for the new cells
        uint tNewCellIndex = 0;

        // Loop through the child meshes
        for ( uint iCell = 0; iCell < tNumDelaunayCells; iCell++ )
        {
            // Get the child mesh index
            uint iCM = aMeshGenerationData->mDelaunayBgCellInds( iCell );

            // Get the child mesh for this background cell
            std::shared_ptr< Child_Mesh_Experimental > tChildMesh = aCutIntegrationMesh->get_child_mesh( iCM );

            // Get the number of new cells to be created for this child mesh
            uint tNumNewCells = tCellToVertexConnectivity( iCell ).size() / 3;

            // Loop through the new cells to be created from the connectivity table
            for ( uint iNewCell = 0; iNewCell < tNumNewCells; iNewCell++ )
            {
                mNewCellChildMeshIndex( tNewCellIndex ) = iCM;

                // Add the global indices of the vertices to the new cell connectivity
                for ( uint iV = 0; iV < 3; iV++ )    // brendan hard coded
                {
                    // Get the local index of the vertex to be added to the cell from the local connectivity
                    moris_index tNewVertexCMOrdinal = tCellToVertexConnectivity( iCell )( iNewCell * 3 + iV ) - 1;
                    MORIS_ERROR( tNewVertexCMOrdinal < (moris::moris_index)tChildMesh->mIgVerts->size(), "Template ordinal out of bounds" );

                    // Get the global index of the vertex and add it to the new cell connectivity
                    mNewCellToVertexConnectivity( tNewCellIndex )( iV ) = tChildMesh->mIgVerts->get_vertex( tNewVertexCMOrdinal )->get_index();
                }

                tNewCellIndex++;
            }
        }

        mNumNewCells += mNumTotalCells;
    }    // end function: Delaunay_Subdivision_Interface::perform_impl_generate_mesh()

    //--------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Delaunay_Subdivision_Interface::get_surface_point_global_coordinates( const Matrix< DDRMat >& aParametricCoordinates, const mtk::Cell& aCell )
    {
        // Create interpolator
        mtk::Interpolation_Function_Factory tInterpolationFactory;
        mtk::Interpolation_Function_Base*   tInterpolation = tInterpolationFactory.create_interpolation_function(
                aCell.get_geometry_type(),
                mtk::Interpolation_Type::LAGRANGE,
                aCell.get_interpolation_order() );

        // Get number of bases
        uint tNumberOfBases = tInterpolation->get_number_of_bases();

        // Size output matrix
        Matrix< DDRMat > tGlobalCoordinates( aParametricCoordinates.n_rows(), aParametricCoordinates.n_cols() );

        // Get the global coordinates of the background nodes
        Matrix< DDRMat > tBackgroundCellCoords = aCell.get_vertex_coords();

        for ( uint iPoint = 0; iPoint < aParametricCoordinates.n_rows(); iPoint++ )
        {

            // Perform interpolation using parametric coordinates
            Matrix< DDRMat > tBasis;
            tInterpolation->eval_N( aParametricCoordinates.get_row( iPoint ), tBasis );

            // Size global coordinates based on parametric coordinates
            Matrix< DDRMat > tPointGlobalCoordinates = Matrix< DDRMat >( 1, aParametricCoordinates.n_cols(), 0.0 );

            // Add contributions from all locators
            for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBases; iBasisIndex++ )
            {
                tPointGlobalCoordinates += tBackgroundCellCoords.get_row( iBasisIndex ) * tBasis( iBasisIndex );
            }

            tGlobalCoordinates.set_row( iPoint, tPointGlobalCoordinates );
        }

        // Clean up
        delete tInterpolation;

        return tGlobalCoordinates;
    }

    //--------------------------------------------------------------------------------------------------

    enum Decomposition_Algorithm_Type
    Delaunay_Subdivision_Interface::get_algorithm_type() const
    {
        return Decomposition_Algorithm_Type::DELAUNAY;
    }

    //--------------------------------------------------------------------------------------------------

    moris_index
    Delaunay_Subdivision_Interface::get_signature() const
    {
        return 333;
    }

    //--------------------------------------------------------------------------------------------------

    mtk::CellTopology
    Delaunay_Subdivision_Interface::get_ig_cell_topology() const
    {
        // BRENDAN UPDATE
        return mtk::CellTopology::TRI3;
    }

    //--------------------------------------------------------------------------------------------------

    uint
    Delaunay_Subdivision_Interface::get_num_geometric_nodes( mtk::Geometry_Type aCellType ) const
    {
        switch ( aCellType )
        {
            case mtk::Geometry_Type::POINT:
                return 1;
            case mtk::Geometry_Type::LINE:
                return 2;
            case mtk::Geometry_Type::TRI:
                return 3;
            case mtk::Geometry_Type::QUAD:
                return 4;
            case mtk::Geometry_Type::TET:
                return 4;
            case mtk::Geometry_Type::HEX:
                return 8;
            case mtk::Geometry_Type::PENTA:
                return 6;
            default:
                MORIS_ASSERT( false, "Unknown geometry type for cell" );
                return 0;
        }
    }

}    // namespace moris::xtk
