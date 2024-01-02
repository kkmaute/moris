/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh_Cleanup.cpp
 *
 */

#include "cl_XTK_Integration_Mesh_Cleanup.hpp"

#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "moris_typedefs.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_MTK_Cell.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_det.hpp"
#include <typeinfo>
#include <numeric> 
#include "cl_Communication_Tools.hpp"
#include "cl_Stopwatch.hpp"

using namespace moris;

namespace xtk
{

    // separate line
    Integration_Mesh_Cleanup::Integration_Mesh_Cleanup(
            Cut_Integration_Mesh*     aMTKMesh,
            Facet_Based_Connectivity* aFacetConnectivity )
            : mInputMesh( aMTKMesh )
            , mInputFacetConnectivity( aFacetConnectivity )
    {
        mNumVerts      = mInputMesh->mIntegrationVertices.size();
        mNumVertGroups = mInputMesh->mIntegrationVertexGroups.size();
        mNumCellGroups = mInputMesh->mIntegrationCellGroups.size();

        // mChildMeshBoundary.resize(mNumVerts, false);
        // mBlkPhaseBoundary.resize(mNumVerts, false);
        // mVertIndToCells.resize(mNumVerts);

        mVertexGroups = mInputMesh->mIntegrationVertexGroups;
        mCellGroups   = mInputMesh->mIntegrationCellGroups;
    }

    // ----------------------------------------------------------------------------

    Integration_Mesh_Cleanup::~Integration_Mesh_Cleanup()
    {
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::set_vertex_cell_connectivity()
    {
        // Each vertex index corresponds to list of attached cells
        mVertIndToCells.resize( mNumVerts );

        // Go through controlled cell indices
        for ( uint iC = mInputMesh->mFirstControlledCellIndex; iC < mInputMesh->mIntegrationCells.size(); iC++ )
        {
            // matrix of vertex indices
            moris::Matrix< IndexMat > VertInds = mInputMesh->mIntegrationCells( iC )->get_vertex_inds();

            // go through cell vertex indices, and add cell index to mVertIndToCells, if it is not already in the list
            for ( uint iV = 0; iV < VertInds.numel(); iV++ )
            {
                if ( std::find( mVertIndToCells( VertInds( iV ) ).begin(), mVertIndToCells( VertInds( iV ) ).end(), iC ) == mVertIndToCells( VertInds( iV ) ).end() )
                {
                    mVertIndToCells( VertInds( iV ) ).push_back( iC );
                }
            }
        }

        // sort cell indices
        for ( uint iV = 0; iV < mVertIndToCells.size(); iV++ )
        {
            std::sort( mVertIndToCells( iV ).begin(), mVertIndToCells( iV ).end() );
        }
    }

    // ------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::set_vertex_vertex_connectivity()
    {
        // Each vertex index corresponds to list of attached vertices
        mVertIndToVerts.resize( mNumVerts );

        // Go through controlled vertex indices
        for ( moris_index iV = 0; (uint)iV < mInputMesh->mIntegrationVertices.size(); iV++ )
        {
            // go through attached cells
            for ( moris_index iC = 0; (uint)iC < mVertIndToCells( iV ).size(); iC++ )
            {
                // vertices of attached cell
                moris::Matrix< IndexMat > tVertInds = mInputMesh->mIntegrationCells( mVertIndToCells( iV )( iC ) )->get_vertex_inds();

                for ( moris_index iVv = 0; (uint)iVv < tVertInds.numel(); iVv++ )    // go through vertices of attached cell
                {
                    // if not adding itself or a duplicate, add vertex
                    if ( tVertInds( iVv ) != iV && std::find( mVertIndToVerts( iV ).begin(), mVertIndToVerts( iV ).end(), tVertInds( iVv ) ) == mVertIndToVerts( iV ).end() )
                    {
                        mVertIndToVerts( iV ).push_back( tVertInds( iVv ) );
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::set_blk_phase_vertex_flags()
    {

        // vertex index corresponds to the bulk phase of one of its attached cells
        Vector< moris_index > VertBlkPhase;
        VertBlkPhase.resize( mNumVerts, -1 );

        mBlkPhaseBoundary.resize( mNumVerts, false );

        // Go through controlled cell indices
        for ( uint iC = mInputMesh->mFirstControlledCellIndex; iC < mInputMesh->mIntegrationCells.size(); iC++ )
        {
            // vertex indices of cell
            moris::Matrix< IndexMat > VertInds = mInputMesh->mIntegrationCells( iC )->get_vertex_inds();

            // save the bulk phase of the current cell
            moris_index tBlkPhase = mInputMesh->mIntegrationCellBulkPhase( iC );

            // for each vertex of the current cell, assign bulk phase,
            // or set flag true if a different bulk phase is present than one already assigned
            for ( uint iV = 0; iV < VertInds.numel(); iV++ )
            {
                if ( VertBlkPhase( VertInds( iV ) ) == -1 )
                {
                    VertBlkPhase( VertInds( iV ) ) = tBlkPhase;
                }
                else if ( VertBlkPhase( VertInds( iV ) ) != tBlkPhase )
                {
                    mBlkPhaseBoundary( VertInds( iV ) ) = true;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::set_owned_verts_flag()
    {
        mOwnedVertices.resize( mNumVerts, true );


        Vector< moris_index > tNotOwnedIgCells;
        tNotOwnedIgCells.reserve( mInputMesh->mIntegrationCells.size() );

        mNotOwnedIgCellIds.reserve( mInputMesh->mIntegrationCells.size() );

        for ( moris_index iG = 0; (uint)iG < mInputMesh->mNotOwnedIntegrationCellGroups.size(); iG++ )
        {
            moris_index tGroupInd = mInputMesh->mNotOwnedIntegrationCellGroups( iG );

            for ( moris_index iC = 0; (uint)iC < mInputMesh->mIntegrationCellGroups( tGroupInd )->mIgCellGroup.size(); iC++ )
            {
                moris_index tCellInd = mInputMesh->mIntegrationCellGroups( tGroupInd )->mIgCellGroup( iC )->get_index();

                if ( std::find( tNotOwnedIgCells.begin(), tNotOwnedIgCells.end(), tCellInd ) == tNotOwnedIgCells.end() )
                {
                    tNotOwnedIgCells.push_back( mInputMesh->mIntegrationCellGroups( tGroupInd )->mIgCellGroup( iC )->get_index() );
                    mNotOwnedIgCellIds.push_back( mInputMesh->mIntegrationCellGroups( tGroupInd )->mIgCellGroup( iC )->get_id() );
                }
            }
        }

        for ( moris_index iV = 0; (uint)iV < mNumVerts; iV++ )
        {
            for ( moris_index iC = 0; (uint)iC < mVertIndToCells( iV ).size(); iC++ )
            {
                if ( std::find( tNotOwnedIgCells.begin(), tNotOwnedIgCells.end(), mVertIndToCells( iV )( iC ) ) != tNotOwnedIgCells.end() )
                {
                    mOwnedVertices( iV ) = false;
                }
            }
        }


        tNotOwnedIgCells.shrink_to_fit();
        mNotOwnedIgCellIds.shrink_to_fit();

        std::sort( mNotOwnedIgCellIds.begin(), mNotOwnedIgCellIds.end() );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::set_child_mesh_boundary_flag()
    {
        mChildMeshBoundary.resize( mNumVerts, false );


        for ( moris_index iV = 0; (uint)iV < mNumVerts; iV++ )
        {
            // go through groups, to use IG_Vertex_Group::get_vertex_local_coords()
            // can be made more efficient if get_vertex_local_coords is implemented for mtk::Vertex
            for ( moris::uint iGroup = 0; iGroup < mVertexGroups.size(); iGroup++ )
            {
                if ( mVertexGroups( iGroup )->vertex_is_in_group( iV ) )
                {
                    moris::Matrix< DDRMat >* tVertLocCoord = mVertexGroups( iGroup )->get_vertex_local_coords( iV ).get();

                    // If any of the local coordinates are +/- 1, then flag is set true
                    if ( tVertLocCoord->numel() == 2 )    // x,y
                    {
                        if ( std::abs( tVertLocCoord->operator()( 0 ) ) == 1 || std::abs( tVertLocCoord->operator()( 1 ) ) == 1 )
                        {
                            mChildMeshBoundary( iV ) = true;
                            break;    // next Vertex
                        }
                    }
                    else if ( tVertLocCoord->numel() == 3 )    // x,y,z
                    {
                        if ( std::abs( tVertLocCoord->operator()( 0 ) ) == 1 || std::abs( tVertLocCoord->operator()( 1 ) ) == 1 || std::abs( tVertLocCoord->operator()( 2 ) ) == 1 )
                        {
                            mChildMeshBoundary( iV ) = true;
                            break;    // next Vertex
                        }
                    }
                    else
                    {
                        MORIS_ERROR( false, "Dimension out of bounds" );
                    }
                }
            }
        }
    }


    // ----------------------------------------------------------------------------

    bool_t
    Integration_Mesh_Cleanup::calcInvert( moris_index aVert1, moris_index aVert2 )
    {
        // cell indices attached to aVert1
        Vector< moris_index > tCells = this->mVertIndToCells( aVert1 );
        tCells.remove( MORIS_INDEX_MAX );

        // flag, if merging aVert1 to aVert2 results in the inversion of a cell
        Vector< bool_t > tinverts;
        tinverts.reserve( tCells.size() );

        // sign of determinant of tri/tet vertices
        // describes if the cell is convex/concave
        bool_t tDetSign;

        for ( moris_index iC = 0; (uint)iC < tCells.size(); iC++ )    // go through Cells attached to aVert1
        {
            // vertex indices of attached cell
            moris::Matrix< IndexMat > tVerts = mInputMesh->mIntegrationCells( tCells( iC ) )->get_vertex_inds();

            bool_t duplicate = true;

            // if cell contains aVert1 and aVert2 => duplicate = true
            // if cell only contains aVert1 => duplicate = false
            if ( mInputMesh->get_spatial_dim() == 2 )
            {
                if ( tVerts( 0 ) != aVert2 && tVerts( 1 ) != aVert2 && tVerts( 2 ) != aVert2 ) { duplicate = false; }
            }
            else if ( mInputMesh->get_spatial_dim() == 3 )
            {
                if ( tVerts( 0 ) != aVert2 && tVerts( 1 ) != aVert2 && tVerts( 2 ) != aVert2 && tVerts( 3 ) != aVert2 ) { duplicate = false; }
            }

            // if cell only contains aVert1, and not aVert2
            //  reason is, transformation of cells is analyzed, but cells containing both aVert1 and aVert2 are deleted in the merge,
            //  so they are omitted in the analysis
            if ( !duplicate )
            {

                Vector< moris::Matrix< DDRMat > > tCoords;
                tCoords.resize( mInputMesh->get_spatial_dim() + 1 );

                // get vertex coordinates of cell being analyzed
                for ( moris_index iP = 0; (uint)iP < tCoords.size(); iP++ )
                {
                    tCoords( iP ) = mInputMesh->mIntegrationVertices( tVerts( iP ) )->get_coords();
                }

                // allocate matrix of vertex coordinates
                moris::Matrix< DDRMat > mat1( mInputMesh->get_spatial_dim(), mInputMesh->get_spatial_dim() );

                // fill matrix
                if ( mInputMesh->get_spatial_dim() == 3 )
                {
                    mat1( { 0, 0 }, { 0, 2 } ) = tCoords( 1 ) - tCoords( 0 );
                    mat1( { 1, 1 }, { 0, 2 } ) = tCoords( 2 ) - tCoords( 0 );
                    mat1( { 2, 2 }, { 0, 2 } ) = tCoords( 3 ) - tCoords( 0 );
                }
                else if ( mInputMesh->get_spatial_dim() == 2 )
                {
                    mat1( { 0, 0 }, { 0, 1 } ) = tCoords( 1 ) - tCoords( 0 );
                    mat1( { 1, 1 }, { 0, 1 } ) = tCoords( 2 ) - tCoords( 0 );
                }

                // calculate determinant of matrix
                if ( det( mat1 ) >= 0 )
                {
                    tDetSign = true;
                }
                else
                {
                    tDetSign = false;
                }

                // replace aVert1 with aVert2 in tCoords, to analyze if the new determinant has changed signs
                for ( moris_index iP = 0; (uint)iP < tCoords.size(); iP++ )
                {
                    if ( mInputMesh->mIntegrationVertices( tVerts( iP ) )->get_index() == aVert1 )
                    {
                        tCoords( iP ) = mInputMesh->mIntegrationVertices( aVert2 )->get_coords();
                    }
                    else
                    {
                        tCoords( iP ) = mInputMesh->mIntegrationVertices( tVerts( iP ) )->get_coords();
                    }
                }

                // refill matrix
                if ( mInputMesh->get_spatial_dim() == 2 )
                {
                    mat1( { 0, 0 }, { 0, 1 } ) = tCoords( 1 ) - tCoords( 0 );
                    mat1( { 1, 1 }, { 0, 1 } ) = tCoords( 2 ) - tCoords( 0 );
                }
                else if ( mInputMesh->get_spatial_dim() == 3 )
                {
                    // moris::Matrix<DDRMat> mat1(3,3);
                    mat1( { 0, 0 }, { 0, 2 } ) = tCoords( 1 ) - tCoords( 0 );
                    mat1( { 1, 1 }, { 0, 2 } ) = tCoords( 2 ) - tCoords( 0 );
                    mat1( { 2, 2 }, { 0, 2 } ) = tCoords( 3 ) - tCoords( 0 );
                }

                // if determinant sign has not changed, tinverts = true, otherwise false
                if ( det( mat1 ) >= 0 )
                {
                    if ( tDetSign == true ) { tinverts.push_back( true ); }
                    else { tinverts.push_back( false ); }
                }
                else
                {
                    if ( tDetSign == true ) { tinverts.push_back( false ); }
                    else { tinverts.push_back( true ); }
                }
            }
        }

        // if no cells were inverted (change in det sign) as a result of merging aVert1 and aVert2, return true
        // otherwise return false
        if ( find( std::begin( tinverts ), std::end( tinverts ), false ) == std::end( tinverts ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }


    // ----------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::delauny( moris_index aVert1, moris_index aVert2 )
    {
        // indices of cells attached to vertex of index aVert1
        Vector< moris_index > tCells = this->mVertIndToCells( aVert1 );

        tCells.remove( MORIS_INDEX_MAX );

        double_t tMinAngle = 99;
        double_t tAngle;

        for ( moris_index iC = 0; (uint)iC < tCells.size(); iC++ )    // go through Cells attached to aVert1
        {
            // vertex indices of cell
            moris::Matrix< IndexMat > tVerts = mInputMesh->mIntegrationCells( tCells( iC ) )->get_vertex_inds();

            bool_t duplicate = true;

            // if cell contains aVert1 and aVert2 => duplicate = true
            // if cell only contains aVert1 => duplicate = false
            if ( mInputMesh->get_spatial_dim() == 2 )
            {
                if ( tVerts( 0 ) != aVert2 && tVerts( 1 ) != aVert2 && tVerts( 2 ) != aVert2 ) { duplicate = false; }
            }
            else if ( mInputMesh->get_spatial_dim() == 3 )
            {
                if ( tVerts( 0 ) != aVert2 && tVerts( 1 ) != aVert2 && tVerts( 2 ) != aVert2 && tVerts( 3 ) != aVert2 ) { duplicate = false; }
            }

            // if cell only contains aVert1, and not aVert2
            //  reason is, transformation of cells is analyzed, but cells containing both aVert1 and aVert2 are deleted in the merge,
            //  so they are omitted in the analysis
            if ( !duplicate )
            {

                Vector< moris::Matrix< DDRMat > > tCoords;
                tCoords.resize( mInputMesh->get_spatial_dim() + 1 );

                // replace aVert1 with aVert2 in tCoords, to analyze the resulting triangle/tet
                for ( moris_index iP = 0; (uint)iP < tCoords.size(); iP++ )
                {
                    if ( mInputMesh->mIntegrationVertices( tVerts( iP ) )->get_index() == aVert1 )
                    {
                        tCoords( iP ) = mInputMesh->mIntegrationVertices( aVert2 )->get_coords();
                    }
                    else
                    {
                        tCoords( iP ) = mInputMesh->mIntegrationVertices( tVerts( iP ) )->get_coords();
                    }
                }

                // calculate the minimum angle/solid angle of the resulting triangle/tet
                if ( mInputMesh->get_spatial_dim() == 2 )
                {
                    tAngle = minAngle2D( tCoords( 0 ), tCoords( 1 ), tCoords( 2 ) );
                }
                else if ( mInputMesh->get_spatial_dim() == 3 )
                {
                    tAngle = minAngle3D( tCoords( 0 ), tCoords( 1 ), tCoords( 2 ), tCoords( 3 ) );
                }
                else
                {
                    MORIS_ERROR( false, "Integration_Mesh_Cleanup::delauny() - unknown number of spatial dimensions." );
                    tAngle = 0.0;
                }

                // save minimum angle of all the transformed cells
                if ( tAngle < tMinAngle ) { tMinAngle = tAngle; }
            }
        }

        return tMinAngle;
    }

    // ----------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::minAngle2D( moris::Matrix< DDRMat > C1, moris::Matrix< DDRMat > C2, moris::Matrix< DDRMat > C3 )
    {
        double_t A1 = angle( C2 - C1, C3 - C1 );
        double_t A2 = angle( C3 - C2, C1 - C2 );
        double_t A3 = angle( C1 - C3, C2 - C3 );

        return std::min( std::min( A1, A2 ), A3 );
    }
    // ----------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::minAngle3D( moris::Matrix< DDRMat > C1,
            moris::Matrix< DDRMat >                               C2,
            moris::Matrix< DDRMat >                               C3,
            moris::Matrix< DDRMat >                               C4 )
    {

        double_t sA1 = solidAngle( C1, C2, C3, C4 );
        double_t sA2 = solidAngle( C2, C3, C4, C1 );
        double_t sA3 = solidAngle( C3, C4, C1, C2 );
        double_t sA4 = solidAngle( C4, C1, C2, C3 );

        return std::min( std::min( sA1, sA2 ), std::min( sA3, sA4 ) );
    }

    // --------------------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::solidAngle( moris::Matrix< DDRMat > nodeO,
            moris::Matrix< DDRMat >                               nodeA,
            moris::Matrix< DDRMat >                               nodeB,
            moris::Matrix< DDRMat >                               nodeC )
    {
        MORIS_ASSERT( nodeO.numel() == 3, "solid angle only for 3d points" );

        double_t TAOB = angle( nodeA - nodeO, nodeB - nodeO );
        double_t TAOC = angle( nodeA - nodeO, nodeC - nodeO );
        double_t TBOC = angle( nodeB - nodeO, nodeC - nodeO );

        double_t Ts = ( TAOB + TAOC + TBOC ) / 2;

        if ( Ts - TAOB < 0 && Ts - TAOB > -1e-6 )
        {
            return 0;
        }
        else if ( Ts - TAOC < 0 && Ts - TAOC > -1e-6 )
        {
            return 0;
        }
        else if ( Ts - TBOC < 0 && Ts - TBOC > -1e-6 )
        {
            return 0;
        }
        if ( Ts - TAOB < 0 || Ts - TAOC < 0 || Ts - TBOC < 0 ) { MORIS_ASSERT( false, "Negative solid angle outside tolerance." ); }

        return 4 * atan( sqrt( tan( Ts / 2 ) * tan( ( Ts - TAOB ) / 2 ) * tan( ( Ts - TAOC ) / 2 ) * tan( ( Ts - TBOC ) / 2 ) ) );
    }

    // ---------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::angle( moris::Matrix< DDRMat > nodeA,
            moris::Matrix< DDRMat >                          nodeB )
    {

        double_t ab    = dot( nodeA, nodeB );
        double_t mag_a = norm( nodeA );
        double_t mag_b = norm( nodeB );

        return safeAcos( ab / ( mag_a * mag_b ) );
    }

    // ---------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::safeAcos( double_t x )
    {
        if ( x < -1.0 )
        {
            if ( x < -1 - 1e-6 ) { MORIS_ASSERT( false, "Acos function out of angle tolerance range" ); }
            x = -1.0;
        }
        else if ( x > 1.0 )
        {
            if ( x > 1 + 1e-6 ) { MORIS_ASSERT( false, "Acos function out of angle tolerance range" ); }
            x = 1.0;
        }
        return acos( x );
    }

    // ----------------------------------------------------------------------------

    Vector< moris_index >
    Integration_Mesh_Cleanup::check_flats()
    {

        Vector< moris_index > tFlats;
        tFlats.reserve( mInputMesh->mIntegrationCells.size() );

        // go through controlled ig cells
        for ( uint iC = mInputMesh->mFirstControlledCellIndex; iC < mInputMesh->mIntegrationCells.size(); iC++ )
        {
            // save vertex indices
            moris::Matrix< DDRMat > tCoords = mInputMesh->mIntegrationCells( iC )->get_vertex_coords();
            MORIS_ASSERT( tCoords.n_rows() == mInputMesh->get_spatial_dim() + 1, "Error: only checking subdivided integration cells" );

            // area/volume of triangle/tet
            double_t tn_area = n_area( tCoords );

            // if smaller than tolerance, add index
            if ( tn_area < 1e-6 && tn_area > -1e-6 )
            {
                tFlats.push_back( mInputMesh->mIntegrationCells( iC )->get_index() );
            }
            else if ( tn_area < -1e-6 )
            {
                MORIS_ASSERT( false, "negative area/volume" );
            }
        }

        return tFlats;
    }

    // ---------------------------------------------------------------------------

    // merges all coincident vertices
    void
    Integration_Mesh_Cleanup::check_coinc_verts( Vector< moris_index > aFlats, Vector< moris::mtk::Cell* >& aActiveIgCells )
    {
        // Vector<Vector<moris_index>> MergeVerts;
        // MergeVerts.reserve(mInputMesh->mIntegrationCells.size());
        int iC = 0;
        while ( (uint)iC < mInputMesh->mIntegrationCells.size() )
        // for (moris::uint iC = 0; iC < mInputMesh->mIntegrationCells.size(); iC++)
        {
            moris::Matrix< DDRMat > tVertexCoords = mInputMesh->mIntegrationCells( iC )->get_vertex_coords();
            for ( uint iV1 = 0; iV1 < 4; iV1++ )
            {
                for ( uint iV2 = iV1 + 1; iV2 < 4; iV2++ )
                {
                    double_t tDist = dist( tVertexCoords( { iV1, iV1 }, { 0, 2 } ), tVertexCoords( { iV2, iV2 }, { 0, 2 } ) );
                    if ( tDist < 1e-6 )
                    {
                        moris_index Ind1 = mInputMesh->mIntegrationCells( iC )->get_vertex_inds()( iV1 );
                        moris_index Ind2 = mInputMesh->mIntegrationCells( iC )->get_vertex_inds()( iV2 );
                        this->merge( Ind1, Ind2, aActiveIgCells );
                        iV1 = 4;
                        iV2 = 4;
                        iC  = -1;
                    }
                }
            }

            iC++;
        }
    }

    // ---------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::dist( moris::Matrix< DDRMat > aV1, moris::Matrix< DDRMat > aV2 )
    {
        return norm( aV1 - aV2 );
    }

    // ---------------------------------------------------------------------------

    double_t
    Integration_Mesh_Cleanup::n_area( moris::Matrix< DDRMat > aCoords )
    {
        if ( aCoords.numel() == 6 )    // 2D
        {
            // vectors spanning triangle
            moris::Matrix< DDRMat > tVec_a = aCoords( { 1, 1 }, { 0, 1 } ) - aCoords( { 0, 0 }, { 0, 1 } );
            moris::Matrix< DDRMat > tVec_b = aCoords( { 2, 2 }, { 0, 1 } ) - aCoords( { 0, 0 }, { 0, 1 } );

            double_t mag_a = norm( tVec_a );
            double_t mag_b = norm( tVec_b );

            return 0.5 * angle( tVec_a, tVec_b ) * mag_a * mag_b;
        }
        else if ( aCoords.numel() == 12 )    // 3D
        {
            // vectors spanning tet
            moris::Matrix< DDRMat > tVec_a = aCoords( { 1, 1 }, { 0, 2 } ) - aCoords( { 0, 0 }, { 0, 2 } );
            moris::Matrix< DDRMat > tVec_b = aCoords( { 2, 2 }, { 0, 2 } ) - aCoords( { 0, 0 }, { 0, 2 } );
            moris::Matrix< DDRMat > tVec_c = aCoords( { 3, 3 }, { 0, 2 } ) - aCoords( { 0, 0 }, { 0, 2 } );

            moris::Matrix< DDRMat > tCrossProd = cross( tVec_a, tVec_b );

            return ( 1.0 / 6 ) * dot( tCrossProd, tVec_c );
        }
        else
        {
            MORIS_ERROR( false, "incorrect number of elements" );
            return -1;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::num_merges()
    {

        // first count how many vertices are eligible for merge
        mMergeNum = 0;
        for ( moris_index iV = 0; (uint)iV < mNumVerts; iV++ )
        {
            if ( !( mBlkPhaseBoundary( iV ) || mChildMeshBoundary( iV ) ) )
            {
                mMergeNum++;
            }
        }

        // reserve list of vertex indices
        mMergeInds.reserve( mMergeNum );

        // assign vertex indices to merge
        for ( moris_index iV = 0; (uint)iV < mNumVerts; iV++ )
        {
            if ( !( mBlkPhaseBoundary( iV ) || mChildMeshBoundary( iV ) ) )
            {
                mMergeInds.push_back( iV );
            }
        }


        // sort list of indices to merge by id
        // sorting by id ensures that different processors will merge vertices in the same order
        std::sort( mMergeInds.begin(), mMergeInds.end(), [ & ]( const auto& i, const auto& j ) { return mInputMesh->mIntegrationVertices( i )->get_id() < mInputMesh->mIntegrationVertices( j )->get_id(); } );
    }


    void
    Integration_Mesh_Cleanup::merge( moris_index Vert1, moris_index Vert2, Vector< moris::mtk::Cell* >& aActiveIgCells )
    {


        // Info on vertices to be merged ------------------------------------------------------------
        moris_index V1Ind = Vert1;
        moris_id    V1Id  = mInputMesh->mIntegrationVertices( Vert1 )->get_id();

        moris_index V2Ind = Vert2;
        // moris_id V2Id = mInputMesh->mIntegrationVertices(Vert2)->get_id();

        // Cell indices in mControlledIgVerts
        moris_index ControlledVertsInd1 = -1;
        moris_index ControlledVertsInd2 = -1;
        for ( moris_index iV = 0; (uint)iV < mInputMesh->mControlledIgVerts.size(); iV++ )
        {
            if ( V1Ind == mInputMesh->mControlledIgVerts( iV )->get_index() )
            {
                ControlledVertsInd1 = iV;
                if ( ControlledVertsInd2 > -1 ) { break; }
            }
            else if ( V2Ind == mInputMesh->mControlledIgVerts( iV )->get_index() )
            {
                ControlledVertsInd2 = iV;
                if ( ControlledVertsInd1 > -1 ) { break; }
            }
        }


        // replace vertex pointers of cells attached to Vert1 ---------------------------------------

        mtk::Vertex* mVertex2 = mInputMesh->mIntegrationVertices( V2Ind );

        for ( moris_index iC = 0; (uint)iC < mVertIndToCells( V1Ind ).size(); iC++ )    // go through cells attached to Vert1
        {
            if ( mVertIndToCells( V1Ind )( iC ) != MORIS_INDEX_MAX )
            {
                // index of cell in mControlledIgCells list
                moris_index ControlledInd = mInputMesh->get_integration_cell_controlled_index( mVertIndToCells( V1Ind )( iC ) );

                // vertex indices of cell
                moris::Matrix< IndexMat > CellVertInds = mInputMesh->mControlledIgCells( ControlledInd )->get_vertex_inds();

                for ( moris_index iV = 0; (uint)iV < CellVertInds.numel(); iV++ )    // go through cell vertices
                {
                    if ( V1Ind == CellVertInds( iV ) )
                    {
                        mInputMesh->mControlledIgCells( ControlledInd )->replace_vertex_pointer( mVertex2, iV );
                        break;    // next cell
                    }
                }
            }
        }


        // merge vertices in connectivity
        Facet_Based_Connectivity* tOldFacetConnectivity = new Facet_Based_Connectivity;
        *tOldFacetConnectivity                          = *mInputFacetConnectivity;
        mInputFacetConnectivity->merge_vertices( V1Ind, V2Ind, mFacetMergeInds, tOldFacetConnectivity, mVertex2 );


        // Merge cells erased during merge -------------------------------------------
        for ( moris_index iC = 0; (uint)iC < mVertIndToCells( V1Ind ).size(); iC++ )
        {
            if ( mVertIndToCells( V1Ind )( iC ) != MORIS_INDEX_MAX )
            {
                // index of cell
                moris_index Ind = mVertIndToCells( V1Ind )( iC );

                // id of cell
                // moris_id Id = mInputMesh->mIntegrationCells(Ind)->get_id();

                // index of cell in mControlledIgCells list
                moris_index ControlledInd = mInputMesh->get_integration_cell_controlled_index( Ind );


                // if cell contains two identical vertex indices
                if ( !mInputMesh->mControlledIgCells( ControlledInd )->check_unique_vertex_inds() )
                {
                    mCellMergeInds.push_back( Ind );

                    // merge coinciding facets
                    mInputFacetConnectivity->merge_facets( Ind, V1Ind, mFacetMergeInds, tOldFacetConnectivity );

                    // set removed cell index in mControlledIgCells to max
                    mInputMesh->mControlledIgCells( ControlledInd )->set_index( MORIS_INDEX_MAX );


                    // set cell index in mVertIndToCells to max
                    moris::Matrix< IndexMat > tVertInds = mInputMesh->mIntegrationCells( Ind )->get_vertex_inds();
                    for ( int i = 0; (uint)i < tVertInds.numel(); i++ )
                    {
                        if ( tVertInds( i ) == V2Ind )
                        {
                            tVertInds( i ) = V1Ind;
                            break;
                        }
                    }
                    for ( moris_index iV = 0; (uint)iV < tVertInds.numel(); iV++ )
                    {
                        for ( moris_index icC = 0; (uint)icC < mVertIndToCells( tVertInds( iV ) ).size(); icC++ )
                        {
                            if ( mVertIndToCells( tVertInds( iV ) )( icC ) == Ind )
                            {
                                mVertIndToCells( tVertInds( iV ) )( icC ) = MORIS_INDEX_MAX;
                                break;
                            }
                        }
                    }


                    // set cell index in mIntegrationCellBulkPhase to max
                    mInputMesh->mIntegrationCellBulkPhase( Ind ) = MORIS_INDEX_MAX;


                    // set index in mIntegrationCellToCellGroupIndex to max
                    mInputMesh->mIntegrationCellToCellGroupIndex( Ind ) = { MORIS_INDEX_MAX };
                }
            }
        }


        // set vert1 index in mControlledIgVerts to max
        mInputMesh->mControlledIgVerts( ControlledVertsInd1 )->set_index( MORIS_INDEX_MAX );


        // transfer vertex-cell connectivity and set vertex index in mVertIndToCells to max
        for ( moris_index iC = 0; (uint)iC < mVertIndToCells( V1Ind ).size(); iC++ )
        {    // if(mVertIndToCells(V1Ind)(iC) != MORIS_INDEX_MAX)
            mVertIndToCells( V2Ind ).push_back( mVertIndToCells( V1Ind )( iC ) );
            // check for duplicates?
        }
        mVertIndToCells( V1Ind ) = { MORIS_INDEX_MAX };

        // correct vertex-vertex connectivity to account for merge
        this->adjust_vertex_vertex_connectivity( V1Ind, V2Ind );


        // set vertex index from mIgVertexParentEntityIndex to max
        mInputMesh->mIgVertexParentEntityIndex( V1Ind ) = MORIS_INDEX_MAX;

        // set vertex index from mIgVertexParentEntityRank to max
        mInputMesh->mIgVertexParentEntityRank( V1Ind ) = MORIS_INDEX_MAX;

        // set vertex id in mIntegrationVertexIdToIndexMap to max
        mInputMesh->mIntegrationVertexIdToIndexMap[ V1Id ] = MORIS_INDEX_MAX;    // use .at() instead? returns error if key V1Id not found
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::shift_indices( Vector< moris::mtk::Cell* >& aActiveIgCells )
    {


        // Cells ----------------------------------------------------------------------------------

        // get merged vertex ids and nullify vertex pointers in mIntegrationVertices
        Vector< moris_id > tMergeIds;
        tMergeIds.resize( mMergeInds.size(), -1 );
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            tMergeIds( iV )                                      = mInputMesh->mIntegrationVertices( mMergeInds( iV ) )->get_id();
            mInputMesh->mIntegrationVertices( mMergeInds( iV ) ) = NULL;
        }

        // get merged cell ids
        Vector< moris_id > tCellMergeIds;
        tCellMergeIds.resize( mCellMergeInds.size(), -1 );
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            tCellMergeIds( iC ) = mInputMesh->mIntegrationCells( mCellMergeInds( iC ) )->get_id();
        }


        // remove cell pointers from mActiveIgCells
        for ( moris_index iC = 0; (uint)iC < aActiveIgCells.size(); iC++ )
        {
            if ( aActiveIgCells( iC )->get_index() == MORIS_INDEX_MAX )
            {
                aActiveIgCells.erase( iC );
                iC--;
                // break off once all cells are deleted? or not worth it?
            }
        }

        // erase deleted cells from mVertIndToCells
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            moris::Matrix< IndexMat > tVertInds = mInputMesh->mIntegrationCells( mCellMergeInds( iC ) )->get_vertex_inds();
            for ( moris_index iV = 0; (uint)iV < tVertInds.numel(); iV++ )
            {
                if ( tVertInds( iV ) != MORIS_INDEX_MAX )
                {
                    for ( moris_index icC = 0; (uint)icC < mVertIndToCells( tVertInds( iV ) ).size(); icC++ )
                    {
                        if ( mVertIndToCells( tVertInds( iV ) )( icC ) == MORIS_INDEX_MAX )
                        {
                            mVertIndToCells( tVertInds( iV ) ).erase( icC );
                            icC--;
                        }
                    }
                }
            }
        }
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            for ( moris_index icC = 0; (uint)icC < mVertIndToCells( mMergeInds( iV ) ).size(); icC++ )
            {
                if ( mVertIndToCells( mMergeInds( iV ) )( icC ) == MORIS_INDEX_MAX )    // assert that cell has only 1 member?
                {
                    mVertIndToCells( mMergeInds( iV ) ).erase( icC );
                    icC--;
                }
            }
        }


        // nullify cell pointer in mIntegrationCells
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            mInputMesh->mIntegrationCells( mCellMergeInds( iC ) ) = NULL;
        }


        // adjust indices in mVertIndToCells
        moris_index tCellMergeMax = mCellMergeInds( mCellMergeInds.size() - 1 );
        for ( moris_index iV = 0; (uint)iV < mVertIndToCells.size(); iV++ )    // start iV at first integration cell instead of 0?
        {
            for ( moris_index iC = 0; (uint)iC < mVertIndToCells( iV ).size(); iC++ )
            {
                if ( mVertIndToCells( iV )( iC ) > tCellMergeMax )
                {
                    // indices are shifted by the total amount of merged cells
                    mVertIndToCells( iV )( iC ) -= mCellMergeInds.size();
                }
                else if ( mVertIndToCells( iV )( iC ) > mCellMergeInds( 0 ) )
                {
                    // indices are shifted based on their position in the merged cells
                    Vector< moris_index > tvector = mCellMergeInds;    // maybe implement a cell of index to new index instead of repeating this
                    tvector.push_back( mVertIndToCells( iV )( iC ) );
                    std::sort( tvector.begin(), tvector.end() );

                    std::vector< moris_index >::iterator itr = std::find( tvector.begin(), tvector.end(), mVertIndToCells( iV )( iC ) );

                    mVertIndToCells( iV )( iC ) -= std::distance( tvector.begin(), itr );
                }
            }
        }

        // remove deleted cells from mIntegrationCellBulkPhase
        mInputMesh->mIntegrationCellBulkPhase.remove( MORIS_INDEX_MAX );

        // remove deleted cells from mIntegrationCellToCellGroupIndex
        for ( moris_index iC = 0; (uint)iC < mInputMesh->mIntegrationCellToCellGroupIndex.size(); iC++ )
        {
            if ( mInputMesh->mIntegrationCellToCellGroupIndex( iC ).size() > 0 )
            {
                if ( mInputMesh->mIntegrationCellToCellGroupIndex( iC )( 0 ) == MORIS_INDEX_MAX )
                {
                    MORIS_ASSERT( mInputMesh->mIntegrationCellToCellGroupIndex( iC ).size() == 1, "Error" );
                    mInputMesh->mIntegrationCellToCellGroupIndex.erase( iC );
                    iC--;
                }
            }
        }


        // remove deallocated cell pointers from mIntegrationCellGroups
        moris_index tIndOffset = 0;
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            for ( moris_index iG = 0; (uint)iG < mInputMesh->mIntegrationCellGroups.size(); iG++ )
            {
                if ( mInputMesh->mIntegrationCellGroups( iG )->cell_is_in_group( mCellMergeInds( iC ) - tIndOffset ) )
                {
                    mInputMesh->mIntegrationCellGroups( iG )->remove_cell( mCellMergeInds( iC ) - tIndOffset );    // also adjusts indices in group
                }
                else
                {
                    mInputMesh->mIntegrationCellGroups( iG )->shift_indices( mCellMergeInds( iC ) - tIndOffset );
                }
            }
            tIndOffset++;
        }


        // remove cell from mNotOwnedIgCellIds
        for ( moris_index iC = 0; (uint)iC < mNotOwnedIgCellIds.size(); iC++ )
        {
            if ( std::find( tMergeIds.begin(), tMergeIds.end(), mNotOwnedIgCellIds( iC ) ) != tMergeIds.end() )
            {
                mNotOwnedIgCellIds.erase( iC );
                iC--;
            }
        }

        // erase deallocated cells from mIntegrationCells
        tIndOffset = 0;
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            mInputMesh->mIntegrationCells.erase( mCellMergeInds( iC ) - tIndOffset );
            tIndOffset++;
        }

        // fix mIntegrationCellIdToIndexMap
        for ( moris_index iC = 0; (uint)iC < tCellMergeIds.size(); iC++ )
        {
            mInputMesh->mIntegrationCellIdToIndexMap.erase( tCellMergeIds( iC ) );
        }
        for ( moris_index icC = 0; (uint)icC < mInputMesh->mIntegrationCells.size(); icC++ )
        {
            mInputMesh->mIntegrationCellIdToIndexMap[ mInputMesh->mIntegrationCells( icC )->get_id() ] = (moris_index)icC;
        }
        MORIS_ASSERT( mInputMesh->mIntegrationCellIdToIndexMap.size() == mInputMesh->mIntegrationCells.size(), "Size issue" );


        // shift indices in mControlledIgCells
        moris_index offset = 0;
        for ( moris_index iC = 1; (uint)iC < mInputMesh->mControlledIgCells.size(); iC++ )
        {
            // assumed mControlledIgCells is sorted by index
            // MORIS_ASSERT(mInputMesh->mControlledIgCells(iC)->get_index() > mInputMesh->mControlledIgCells(iC-1)->get_index(),"Not sorted");

            if ( mInputMesh->mControlledIgCells( iC - 1 )->get_index() == MORIS_INDEX_MAX )
            {
                offset++;
            }

            if ( mInputMesh->mControlledIgCells( iC )->get_index() != MORIS_INDEX_MAX )
            {
                mInputMesh->mControlledIgCells( iC )->set_index( mInputMesh->mControlledIgCells( iC )->get_index() - offset );
            }
        }


        // deallocate and erase cells from mControlledIgCells
        for ( moris_index iC = 0; (uint)iC < mInputMesh->mControlledIgCells.size(); iC++ )
        {
            if ( mInputMesh->mControlledIgCells( iC )->get_index() == MORIS_INDEX_MAX )
            {
                mInputMesh->mControlledIgCells( iC ).reset();
                mInputMesh->mControlledIgCells.erase( iC );
                iC--;
            }
        }


        // Vertices ----------------------------------------------------------------------------------

        // remove vertex pointer from mIntegrationVertices
        tIndOffset = 0;
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            mInputMesh->mIntegrationVertices.erase( mMergeInds( iV ) - tIndOffset );
            tIndOffset++;
        }


        // remove vertex pointer from mIntegrationVertexGroups
        tIndOffset = 0;
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            for ( moris_index iG = 0; (uint)iG < mInputMesh->mIntegrationVertexGroups.size(); iG++ )
            {
                if ( mInputMesh->mIntegrationVertexGroups( iG )->vertex_is_in_group( mMergeInds( iV ) - tIndOffset ) )
                {
                    mInputMesh->mIntegrationVertexGroups( iG )->remove_vertex( mMergeInds( iV ) - tIndOffset );    // also adjusts indices in group
                }
                else
                {
                    mInputMesh->mIntegrationVertexGroups( iG )->shift_indices( mMergeInds( iV ) - tIndOffset );
                }
            }
            tIndOffset++;
        }


        // remove deleted vertices from mVertIndToCells
        tIndOffset = 0;
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            if ( mVertIndToCells( mMergeInds( iV ) - tIndOffset ).size() == 0 )
            {
                mVertIndToCells.erase( mMergeInds( iV ) - tIndOffset );
                tIndOffset++;
            }
            else
            {
                MORIS_ASSERT( false, "Error, list should be empty" );
            }
        }


        // remove deleted vertices from mVertIndToVerts
        for ( moris_index iV = 0; (uint)iV < mMergeIndsV2.size(); iV++ )
        {
            for ( moris_index ivV = 0; (uint)ivV < mVertIndToVerts( mMergeIndsV2( iV ) ).size(); ivV++ )
            {
                if ( mVertIndToVerts( mMergeIndsV2( iV ) )( ivV ) != MORIS_INDEX_MAX )
                {
                    for ( moris_index ivvV = 0; (uint)ivvV < mVertIndToVerts( mVertIndToVerts( mMergeIndsV2( iV ) )( ivV ) ).size(); ivvV++ )
                    {
                        if ( mVertIndToVerts( mVertIndToVerts( mMergeIndsV2( iV ) )( ivV ) )( ivvV ) == MORIS_INDEX_MAX )
                        {
                            mVertIndToVerts( mVertIndToVerts( mMergeIndsV2( iV ) )( ivV ) ).erase( ivvV );
                            ivvV--;
                        }
                    }
                }
                else
                {
                    mVertIndToVerts( mMergeIndsV2( iV ) ).erase( ivV );
                    ivV--;
                }
            }
        }
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            // MORIS_ASSERT(mVertIndToVerts(mMergeInds(iV)).size() == 1, "Error");
            for ( moris_index ivV = 0; (uint)ivV < mVertIndToVerts( mMergeInds( iV ) ).size(); ivV++ )    // remove this loop? all mMergeInds should just be deleted?
            {
                if ( mVertIndToVerts( mMergeInds( iV ) )( ivV ) == MORIS_INDEX_MAX )
                {
                    mVertIndToVerts( mMergeInds( iV ) ).erase( ivV );
                    ivV--;
                }
            }
        }

        // shift indices of mVertIndToVerts
        moris_index tVertMergeMax = mMergeInds( mMergeInds.size() - 1 );       // check if sorted, otherwise just max?
        for ( moris_index iV = 0; (uint)iV < mVertIndToVerts.size(); iV++ )    // start iV at first controlled index instead of 0?
        {
            for ( moris_index ivV = 0; (uint)ivV < mVertIndToVerts( iV ).size(); ivV++ )
            {
                if ( mVertIndToVerts( iV )( ivV ) > tVertMergeMax )
                {
                    mVertIndToVerts( iV )( ivV ) -= mMergeInds.size();
                }
                else if ( mVertIndToVerts( iV )( ivV ) > mMergeInds( 0 ) )    // check if mMergeInds is sorted? or use minimum here
                {
                    Vector< moris_index > tvector = mMergeInds;    // maybe a cell of index to new index instead of repeating this
                    tvector.push_back( mVertIndToVerts( iV )( ivV ) );
                    std::sort( tvector.begin(), tvector.end() );

                    std::vector< moris_index >::iterator itr = std::find( tvector.begin(), tvector.end(), mVertIndToVerts( iV )( ivV ) );

                    mVertIndToVerts( iV )( ivV ) -= std::distance( tvector.begin(), itr );
                }
            }
        }

        tIndOffset = 0;
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            if ( mVertIndToVerts( mMergeInds( iV ) - tIndOffset ).size() == 0 )
            {
                mVertIndToVerts.erase( mMergeInds( iV ) - tIndOffset );
                tIndOffset++;
            }
            else
            {
                MORIS_ASSERT( false, "Error, list should be empty" );
            }
        }


        // remove vertex index from mIgVertexParentEntitiyIndex
        mInputMesh->mIgVertexParentEntityIndex.remove( MORIS_INDEX_MAX );    // delete elems from mMergeInds instead?

        // remove vertex index from mIgVertexParentEntityRank
        mInputMesh->mIgVertexParentEntityRank.remove( MORIS_INDEX_MAX );    // delete elems from mMergeInds instead?

        // deallocate and remove vertex from mVertexCoordinates
        tIndOffset = 0;
        for ( moris_index iV = 0; (uint)iV < mMergeInds.size(); iV++ )
        {
            mInputMesh->mVertexCoordinates( mMergeInds( iV ) - tIndOffset ).reset();
            mInputMesh->mVertexCoordinates.erase( mMergeInds( iV ) - tIndOffset );
            tIndOffset++;
        }


        // adjust vertex id to index map
        for ( auto it = mInputMesh->mIntegrationVertexIdToIndexMap.begin(); it != mInputMesh->mIntegrationVertexIdToIndexMap.end(); )
        {
            if ( it->second == MORIS_INDEX_MAX )
            {
                mInputMesh->mIntegrationVertexIdToIndexMap.erase( it++ );
            }
            else { ++it; }
        }


        MORIS_ASSERT( mInputMesh->mIntegrationVertexIdToIndexMap.size() == mInputMesh->mIntegrationVertices.size(), "Error size" );
        for ( moris_index iV = 0; (uint)iV < mInputMesh->mIntegrationVertices.size(); iV++ )
        {
            mInputMesh->mIntegrationVertexIdToIndexMap[ mInputMesh->mIntegrationVertices( iV )->get_id() ] = (moris_index)iV;
        }
        MORIS_ASSERT( mInputMesh->mIntegrationVertexIdToIndexMap.size() == mInputMesh->mIntegrationVertices.size(), "Error size" );

        // shift indices in mControlledIgVerts
        offset = 0;
        for ( moris_index iV = 1; (uint)iV < mInputMesh->mControlledIgVerts.size(); iV++ )
        {
            // assumed mControlledIgVerts is sorted by index
            // MORIS_ASSERT(mInputMesh->mControlledIgVerts(iV)->get_index() > mInputMesh->mControlledIgVerts(iV-1)->get_index(),"Not sorted");

            if ( mInputMesh->mControlledIgVerts( iV - 1 )->get_index() == MORIS_INDEX_MAX )
            {
                offset++;
            }

            if ( mInputMesh->mControlledIgVerts( iV )->get_index() != MORIS_INDEX_MAX )
            {
                mInputMesh->mControlledIgVerts( iV )->set_index( mInputMesh->mControlledIgVerts( iV )->get_index() - offset );
            }
        }

        // deallocate and erase verts from mControlledIgVerts
        for ( moris_index iV = 0; (uint)iV < mInputMesh->mControlledIgVerts.size(); iV++ )
        {
            if ( mInputMesh->mControlledIgVerts( iV )->get_index() == MORIS_INDEX_MAX )
            {
                mInputMesh->mControlledIgVerts( iV ).reset();
                mInputMesh->mControlledIgVerts.erase( iV );    // break after all verts have been deleted? or just go through controlled indices of verts to be deleted
                iV--;
            }
        }

        // adjust mFirstControlledVertexIndex
        mInputMesh->mFirstControlledVertexIndex = mInputMesh->mControlledIgVerts( 0 )->get_index();


        // Facets ----------------------------------------------------------------------------------


        // new map of mInputFacetConnectivity->mCellIndexToCellOrdinal
        for ( moris_index iC = 0; (uint)iC < mCellMergeInds.size(); iC++ )
        {
            MORIS_ASSERT( mInputFacetConnectivity->mCellIndexToCellOrdinal.at( mCellMergeInds( iC ) ) == MORIS_INDEX_MAX, "Error" );
            mInputFacetConnectivity->mCellIndexToCellOrdinal.erase( mCellMergeInds( iC ) );
        }
        Vector< moris_index > oldOrds;
        oldOrds.reserve( mInputFacetConnectivity->mCellIndexToCellOrdinal.size() );
        for ( std::pair< moris_id, moris_index > el : mInputFacetConnectivity->mCellIndexToCellOrdinal )
        {
            oldOrds.push_back( el.second );
        }

        Vector< moris_index > transOrds;
        transOrds.resize( oldOrds.size(), -1 );
        std::iota( transOrds.begin(), transOrds.end(), 0 );
        std::sort( transOrds.begin(), transOrds.end(), [ & ]( moris_index i, moris_index j ) { return oldOrds( i ) < oldOrds( j ); } );

        Vector< moris_index > newOrds;
        newOrds.resize( oldOrds.size(), -1 );
        for ( moris_index iO = 0; (uint)iO < oldOrds.size(); iO++ )
        {
            newOrds( transOrds( iO ) ) = iO;
        }

        std::unordered_map< moris_id, moris_index > tNewMap;

        Vector< moris_id > newKeys;
        newKeys.resize( newOrds.size(), -1 );


        moris_index i = 0;
        for ( std::pair< moris_id, moris_index > el : mInputFacetConnectivity->mCellIndexToCellOrdinal )    // todo: this is slow
        {

            if ( el.first > tCellMergeMax )    //
            {
                newKeys( i ) = el.first - mCellMergeInds.size();
                i++;
            }
            else if ( el.first > mCellMergeInds( 0 ) )    // check if mCellMergeInds is sorted? or use minimum here
            {
                Vector< moris_index > tvector = mCellMergeInds;    // maybe a cell of index to new index instead of repeating this
                tvector.push_back( el.first );
                std::sort( tvector.begin(), tvector.end() );

                std::vector< moris_index >::iterator itr = std::find( tvector.begin(), tvector.end(), el.first );

                newKeys( i ) = el.first - std::distance( tvector.begin(), itr );
                i++;
            }
            else
            {
                newKeys( i ) = el.first;
                i++;
            }
        }

        for ( moris_index iC = (uint)newKeys.size() - 1; iC >= 0; iC-- )
        {
            tNewMap[ newKeys( iC ) ] = newOrds( iC );
        }

        mInputFacetConnectivity->mCellIndexToCellOrdinal = tNewMap;


        // erase cells from mCellToFacet
        // todo: could do this directly with mFacetMergeInds? mFacetMergeInds might also simplify the whole mCellIndexToCellOrdinal rewrite
        for ( moris_index iC = 0; (uint)iC < mInputFacetConnectivity->mCellToFacet.size(); iC++ )
        {
            if ( mInputFacetConnectivity->mCellToFacet( iC ).size() == 1 )
            {
                MORIS_ASSERT( mInputFacetConnectivity->mCellToFacet( iC )( 0 ) == MORIS_INDEX_MAX, "Error in deletion" );
                mInputFacetConnectivity->mCellToFacet.erase( iC );
                iC--;
            }
        }

        // shift indices in mCellToFacet
        moris_index tFacetMergeMax = *max_element( mFacetMergeInds.begin(), mFacetMergeInds.end() );
        moris_index tFacetMergeMin = *min_element( mFacetMergeInds.begin(), mFacetMergeInds.end() );
        for ( moris_index iC = 0; (uint)iC < mInputFacetConnectivity->mCellToFacet.size(); iC++ )    // todo: also slow, figure out index shift once, then apply here
        {
            for ( moris_index iF = 0; (uint)iF < mInputFacetConnectivity->mCellToFacet( iC ).size(); iF++ )
            {
                if ( mInputFacetConnectivity->mCellToFacet( iC )( iF ) > tFacetMergeMax )
                {
                    mInputFacetConnectivity->mCellToFacet( iC )( iF ) -= mFacetMergeInds.size();
                }
                else if ( mInputFacetConnectivity->mCellToFacet( iC )( iF ) > tFacetMergeMin )
                {
                    Vector< moris_index > tvector = mFacetMergeInds;    // maybe a cell of index to new index instead of repeating this
                    tvector.push_back( mInputFacetConnectivity->mCellToFacet( iC )( iF ) );
                    std::sort( tvector.begin(), tvector.end() );

                    std::vector< moris_index >::iterator itr = std::find( tvector.begin(), tvector.end(), mInputFacetConnectivity->mCellToFacet( iC )( iF ) );

                    mInputFacetConnectivity->mCellToFacet( iC )( iF ) -= std::distance( tvector.begin(), itr );
                }
            }
        }


        // mFacetToCellEdgeOrdinal
        for ( moris_index iF = 0; (uint)iF < mInputFacetConnectivity->mFacetToCellEdgeOrdinal.size(); iF++ )
        {
            if ( mInputFacetConnectivity->mFacetToCellEdgeOrdinal( iF ).size() == 1 )
            {
                if ( mInputFacetConnectivity->mFacetToCellEdgeOrdinal( iF )( 0 ) == MORIS_INDEX_MAX )
                {
                    mInputFacetConnectivity->mFacetToCellEdgeOrdinal.erase( iF );
                    iF--;
                }
            }
            else
            {
                for ( moris_index iE = 0; (uint)iE < mInputFacetConnectivity->mFacetToCellEdgeOrdinal( iF ).size(); iE++ )
                {
                    if ( mInputFacetConnectivity->mFacetToCellEdgeOrdinal( iF )( iE ) == MORIS_INDEX_MAX )
                    {
                        mInputFacetConnectivity->mFacetToCellEdgeOrdinal( iF ).erase( iE );
                        iE--;
                    }
                }
            }
        }

        // mFacetToCell
        std::sort( mFacetMergeInds.begin(), mFacetMergeInds.end() );    // todo: put this elsewhere
        tIndOffset = 0;
        for ( moris_index iF = 0; (uint)iF < mFacetMergeInds.size(); iF++ )    // todo: these types of loops are slow, instead copy list, leaving out the indices to be removed.
        {
            mInputFacetConnectivity->mFacetToCell.erase( mFacetMergeInds( iF ) - tIndOffset );
            tIndOffset++;
        }

        for ( moris_index iF = 0; (uint)iF < mInputFacetConnectivity->mFacetToCell.size(); iF++ )
        {
            for ( moris_index iC = 0; (uint)iC < mInputFacetConnectivity->mFacetToCell( iF ).size(); iC++ )
            {
                if ( mInputFacetConnectivity->mFacetToCell( iF )( iC ) == NULL )
                {
                    mInputFacetConnectivity->mFacetToCell( iF ).erase( iC );
                    iC--;
                }
            }
        }

        // mFacetVertices
        for ( moris_index iF = 0; (uint)iF < mInputFacetConnectivity->mFacetVertices.size(); iF++ )
        {
            if ( mInputFacetConnectivity->mFacetVertices( iF )( 0 ) == NULL )    // todo: replace with mFacetVertices(iF) == {NULL, NULL}
            {
                MORIS_ASSERT( mInputFacetConnectivity->mFacetVertices( iF )( 1 ) == NULL, "Error" );
                mInputFacetConnectivity->mFacetVertices.erase( iF );
                iF--;
            }
            else
            {
                for ( moris_index iV = 0; (uint)iV < mInputFacetConnectivity->mFacetVertices( iF ).size(); iV++ )
                {
                    if ( mInputFacetConnectivity->mFacetVertices( iF )( iV ) == NULL )
                    {
                        MORIS_ASSERT( false, "Error" );
                        mInputFacetConnectivity->mFacetVertices( iF ).erase( iV );
                        iV--;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh_Cleanup::get_attached_vertex( moris_index aVert1 )
    {
        // list of attached vertices (indices)
        Vector< moris_index > attachedNodes = mVertIndToVerts( aVert1 );

        // list of minimum angles created by merging to attached vertex
        Vector< double > tAttachedVertsAngle;
        tAttachedVertsAngle.reserve( attachedNodes.size() );

        // flag if merge inverts a cell, true if no cell is inverted
        Vector< bool_t > tAttachedVertsFlag;
        tAttachedVertsFlag.reserve( attachedNodes.size() );

        // fill angle list and flags for all attached vertices
        for ( moris_index iN = 0; (uint)iN < attachedNodes.size(); iN++ )
        {
            if ( attachedNodes( iN ) != MORIS_INDEX_MAX )
            {
                tAttachedVertsAngle.push_back( delauny( aVert1, attachedNodes( iN ) ) );
                tAttachedVertsFlag.push_back( calcInvert( aVert1, attachedNodes( iN ) ) );
            }
            else
            {
                tAttachedVertsAngle.push_back( -1 );
                tAttachedVertsFlag.push_back( false );
            }
        }

        // set merge in attachedVertsAngle to -1 if the merge creates an inverted tri/tet or double merge
        for ( moris_index iN = 0; (uint)iN < attachedNodes.size(); iN++ )
        {
            if ( !tAttachedVertsFlag( iN ) )
            {
                tAttachedVertsAngle( iN ) = -1;
            }
        }


        // index of max value in tAttachedVertsAngle
        moris_index maxAngleInd = std::max_element( tAttachedVertsAngle.begin(), tAttachedVertsAngle.end() ) - tAttachedVertsAngle.begin();


        if ( tAttachedVertsAngle( maxAngleInd ) == 0 )
        {
            // std::cout<<"Merge creates flat tet"<<std::endl;
        }
        if ( tAttachedVertsAngle( maxAngleInd ) == -1 )
        {
            // std::cout<<"Merge creates inverted tet"<<std::endl;
            return -1;
        }


        return attachedNodes( maxAngleInd );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::adjust_vertex_vertex_connectivity( moris_index aVert1, moris_index aVert2 )
    {

        if ( aVert2 == 784 )
        {
            std::cout << "stp" << std::endl;
        }


        // add verts from aVert1 to aVert2 without duplicates
        for ( moris_index iV = 0; (uint)iV < mVertIndToVerts( aVert1 ).size(); iV++ )
        {
            // check for itself
            if ( mVertIndToVerts( aVert1 )( iV ) != aVert2 )
            {
                // check for duplicates
                if ( std::find( mVertIndToVerts( aVert2 ).begin(), mVertIndToVerts( aVert2 ).end(), mVertIndToVerts( aVert1 )( iV ) ) == mVertIndToVerts( aVert2 ).end() )
                {
                    moris_index tAddInd = mVertIndToVerts( aVert1 )( iV );
                    mVertIndToVerts( aVert2 ).push_back( tAddInd );    // todo: fix?
                }
            }
        }

        // set aVert1 in all relevant lists to max
        for ( moris_index iV = 0; (uint)iV < mVertIndToVerts( aVert1 ).size(); iV++ )
        {
            if ( mVertIndToVerts( aVert1 )( iV ) != MORIS_INDEX_MAX )
            {
                for ( moris_index ivV = 0; (uint)ivV < mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) ).size(); ivV++ )
                {
                    if ( mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) )( ivV ) == aVert1 )
                    {
                        mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) )( ivV ) = MORIS_INDEX_MAX;
                    }
                }
            }
        }

        // add aVert2 to verts connected to aVert1
        for ( moris_index iV = 0; (uint)iV < mVertIndToVerts( aVert1 ).size(); iV++ )
        {
            if ( mVertIndToVerts( aVert1 )( iV ) != MORIS_INDEX_MAX )
            {
                for ( moris_index ivV = 0; (uint)ivV < mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) ).size(); ivV++ )
                {
                    Vector< moris_index > tattachedCells = mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) );
                    if ( std::find( tattachedCells.begin(), tattachedCells.end(), aVert2 ) == tattachedCells.end() && mVertIndToVerts( aVert1 )( iV ) != aVert2 )
                    {
                        mVertIndToVerts( mVertIndToVerts( aVert1 )( iV ) ).push_back( aVert2 );
                    }
                }
            }
        }


        // set vertex itself to single element with max
        mVertIndToVerts( aVert1 ) = { MORIS_INDEX_MAX };
    }

    // ----------------------------------------------------------------------------


    void
    Integration_Mesh_Cleanup::communicate_merged_cells( moris::Matrix< IdMat > aNotOwnedCellVerts )
    {

        // get current proc's rank
        moris_id tMyRank   = par_rank();
        moris_id tCommSize = par_size();


        if ( tCommSize == 1 )    // serial
        {
            // check that all cells are owned in serial
            MORIS_ASSERT( mInputMesh->mOwnedChildMeshes.size() == mInputMesh->mChildMeshes.size(),
                    "Integration_Mesh_Cleanup::communicate_merged_cells() - Code running in serial; not all cells are owned by proc 0." );
        }
        else    // parallel
        {

            // prepare communication list
            Matrix< IdMat > tCommunicationList( tCommSize, 1 );
            for ( moris_id iProc = 0; iProc < tCommSize; iProc++ )
            {
                tCommunicationList( iProc ) = iProc;
            }

            Vector< Matrix< IdMat > > tRequestNotOwnedCellIds( tCommSize );
            Vector< Matrix< IdMat > > tRequestNotOwnedCellVerts( tCommSize );


            int dim = mInputMesh->get_spatial_dim();

            // resize and fill communication data
            for ( moris_index iProc = 0; iProc < tCommSize; iProc++ )
            {
                tRequestNotOwnedCellIds( iProc ).set_size( 1, mNotOwnedIgCellIds.size() );
                tRequestNotOwnedCellVerts( iProc ).set_size( mNotOwnedIgCellIds.size(), dim + 1 );

                for ( moris_index iC = 0; (uint)iC < mNotOwnedIgCellIds.size(); iC++ )
                {
                    tRequestNotOwnedCellIds( iProc )( iC )                       = (moris_id)mNotOwnedIgCellIds( iC );
                    tRequestNotOwnedCellVerts( iProc )( { iC, iC }, { 0, dim } ) = aNotOwnedCellVerts( { iC, iC }, { 0, dim } );
                }
            }

            // communication data being received
            Vector< Matrix< IdMat > > tReceivedNotOwnedCellIds;
            Vector< Matrix< IdMat > > tReceivedNotOwnedCellVerts;

            // communicate between processors
            communicate_mats( tCommunicationList, tRequestNotOwnedCellIds, tReceivedNotOwnedCellIds );
            communicate_mats( tCommunicationList, tRequestNotOwnedCellVerts, tReceivedNotOwnedCellVerts );


            // compare
            for ( moris_index iProc = 0; iProc < tCommSize; iProc++ )
            {
                // go through all processors, skip the current active processor
                if ( iProc != tMyRank )
                {
                    // go through received not owned cell ids from other processor (ids are in rising order)
                    for ( moris_index iC = 0; (uint)iC < tReceivedNotOwnedCellIds( iProc ).numel(); iC++ )
                    {
                        // if received cell id is in current processor, compare vertex ids
                        if ( mInputMesh->mIntegrationCellIdToIndexMap.find( tReceivedNotOwnedCellIds( iProc )( iC ) ) != mInputMesh->mIntegrationCellIdToIndexMap.end() )
                        {
                            moris_index CellInd = mInputMesh->mIntegrationCellIdToIndexMap[ tReceivedNotOwnedCellIds( iProc )( iC ) ];

                            for ( int iDim = 0; iDim <= dim; iDim++ )
                            {
                                if ( mInputMesh->mIntegrationCells( CellInd )->get_vertex_ids()( iDim ) != tReceivedNotOwnedCellVerts( iProc )( iC, iDim ) )
                                {
                                    MORIS_ASSERT( false, "Error, aura cell mismatch" );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::make_GenMeshMap( moris::ge::Geometry_Engine* aGeometryEngine )
    {
        Vector< moris_index > GenMeshMap;
        GenMeshMap.reserve( mNumVerts );

        for ( moris_index iV = 0; (uint)iV < mNumVerts; iV++ )
        {
            // this may assume mMergeInds is sorted by index? May be issues if sorting by id does not result in a list sorted by index as well
            if ( std::find( mMergeInds.begin(), mMergeInds.end(), iV ) == mMergeInds.end() )
            {
                GenMeshMap.push_back( iV );
            }
        }

        aGeometryEngine->get_pdv_host_manager()->set_GenMeshMap( GenMeshMap );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::merge_list( Vector< moris::mtk::Cell* >& aActiveIgCells )
    {

        // allocate lists of indices (cells, facets, vertices) being merged
        mMergeIndsV2.reserve( mMergeNum );
        mCellMergeInds.reserve( mMergeNum * 3 );
        mFacetMergeInds.reserve( mMergeNum * 3 );


        for ( uint iM = 0; iM < mMergeInds.size(); iM++ )
        {
            // set vertex being merged (Vert1), and vertex it is being "merged to" (aVert2)
            moris_index Vert1 = mMergeInds( iM );
            moris_index Vert2 = this->get_attached_vertex( mMergeInds( iM ) );

            if ( Vert2 != -1 )    // merge only if a suitable vertex to merge to was found
            {
                mMergeIndsV2.push_back( Vert2 );

                // merges Vert1->Vert2, keeping index values
                this->merge( Vert1, Vert2, aActiveIgCells );
            }
            else
            {
                mMergeInds.erase( iM );
                iM--;
            }
        }

        mCellMergeInds.shrink_to_fit();
        mFacetMergeInds.shrink_to_fit();

        std::sort( mCellMergeInds.begin(), mCellMergeInds.end() );

        // correct cell, facet, vertex indices
        this->shift_indices( aActiveIgCells );
    }


    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Cleanup::perform( Vector< moris::mtk::Cell* >& aActiveIgCells, moris::ge::Geometry_Engine* aGeometryEngine )
    {

        // mFlats = this->check_flats(); //checks for degenerate cells
        // std::cout<<"There are "<<mFlats.size()<<" flat tets before merge."<<std::endl;


        // sets mVertIndToCells, in: index of vertex || out: list of cells (indices) attached to vertex
        this->set_vertex_cell_connectivity();

        // sets mVertIndToVerts, in: index of vertex || out: list of vertices (indices) attached to vertex
        this->set_vertex_vertex_connectivity();

        // sets mBlkPhaseBoundary, in: index of vertex || out: bool that is true if vertex is on the boundary of a blk phase
        this->set_blk_phase_vertex_flags();

        // sets mChildMeshBoundary, in: index of vertex || out: bool that is true if vertex is on the boundary of a child mesh
        this->set_child_mesh_boundary_flag();

        // sets mOwnedVertices, in: index of vertex || out: bool that is true if vertex is owned by processor
        // sets mNotOwnedIgCellIds, list of not owned cell ids
        this->set_owned_verts_flag();


        // sets mMergeInds, indices of vertices not on bulk phase boundary and not on child mesh boundary (controls which verts are merged)
        this->num_merges();

        // translate vertex indices from before merge to after merge (used in FEM)
        this->make_GenMeshMap( aGeometryEngine );

        // merge the list of vertices set in num_merges()
        this->merge_list( aActiveIgCells );


        // store not owned cell vertex ids for parallel cross check
        int dim = mInputMesh->get_spatial_dim();

        moris::Matrix< IdMat > tNotOwnedCellVerts;
        tNotOwnedCellVerts.set_size( mNotOwnedIgCellIds.size(), dim + 1 );

        for ( moris_index iC = 0; (uint)iC < mNotOwnedIgCellIds.size(); iC++ )
        {
            tNotOwnedCellVerts( { iC, iC }, { 0, dim } ) = mInputMesh->mIntegrationCells( mInputMesh->mIntegrationCellIdToIndexMap[ mNotOwnedIgCellIds( iC ) ] )->get_vertex_ids().matrix_data();
        }

        // parallel cross check
        this->communicate_merged_cells( tNotOwnedCellVerts );


        // check for degenerate cells
        // mFlats = this->check_flats();
        // std::cout<<"There are "<<mFlats.size()<<" flat tets after merge."<<std::endl;
    }

}    // namespace xtk