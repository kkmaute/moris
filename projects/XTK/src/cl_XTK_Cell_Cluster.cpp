/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_Cluster.cpp
 *
 */

#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_MTK_Cluster_Group.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "fn_stringify_matrix.hpp"
#include "HDF5_Tools.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "fn_linsolve.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "fn_linsolve.hpp"
#include "fn_det.hpp"
#include "fn_cross.hpp"


// namespace moris
// {
namespace moris::xtk
{
    //----------------------------------------------------------------

    Cell_Cluster::Cell_Cluster()
            : mTrivial( true )
            , mVoid( false )
            , mInvalid( false )
            , mInterpolationCell( nullptr )
            , mChildMesh( nullptr )
            , mPrimaryIntegrationCells( 0, nullptr )
            , mVoidIntegrationCells( 0, nullptr )
            , mVerticesInCluster( 0, nullptr )
            , mClusterGroups( 0 )

    {
    }

    Cell_Cluster::Cell_Cluster( bool aOnlyForVisualization )
            : mTrivial( true )
            , mVoid( false )
            , mInvalid( false )
            , mInterpolationCell( nullptr )
            , mChildMesh( nullptr )
            , mPrimaryIntegrationCells( 0, nullptr )
            , mVoidIntegrationCells( 0, nullptr )
            , mVerticesInCluster( 0, nullptr )
            , mClusterGroups( 0 )
    {
        mOnlyForVis = aOnlyForVisualization;
    }

    //----------------------------------------------------------------

    Cell_Cluster::~Cell_Cluster() {}

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_trivial( const mtk::Leader_Follower aIsLeader ) const
    {
        return mTrivial;
    }

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_full() const
    {
        // Cluster is full if it is trivial and not void
        return ( mTrivial && !mVoid );
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Cell const * > const &
    Cell_Cluster::get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mPrimaryIntegrationCells;
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Cell const * > const &
    Cell_Cluster::get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }

    //----------------------------------------------------------------

    moris::mtk::Cell const &
    Cell_Cluster::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        return *mInterpolationCell;
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Cell_Cluster::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mVerticesInCluster;
    }

    //----------------------------------------------------------------

    Matrix< DDRMat >
    Cell_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        if ( !mTrivial )
        {
            return mLocalCoords;
        }
        else
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const *tCellInfo = mInterpolationCell->get_cell_info();

            // local coordinate matrix
            Matrix< DDRMat > tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coords_of_cell( tXi );

            return tXi;
        }
    }

    //----------------------------------------------------------------

    Matrix< DDRMat >
    Cell_Cluster::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const *aVertex,
            const mtk::Leader_Follower                                                   aIsLeader ) const
    {
        MORIS_ERROR( !mTrivial, "Accessing local coordinates on a trivial cell cluster is not allowed" );
        return *mVertexGroup->get_vertex_local_coords( aVertex->get_index() );
    }

    //----------------------------------------------------------------

    moris_index
    Cell_Cluster::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
    {
        return this->get_vertices_local_coordinates_wrt_interp_cell( aIsLeader ).n_cols();
    }

    Matrix< DDRMat >
    Cell_Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const
    {
        if ( mTrivial )
        {
            return mLocalCoords;
        }
        else
        {
            // MORIS_ERROR(!this->is_trivial(),"get_primary_cell_local_coords_on_side_wrt_interp_cell on trivial cluster is not allowed");
            MORIS_ASSERT( aPrimaryCellClusterIndex < (moris_index)this->get_num_primary_cells(), "Integration Cell Cluster index out of bounds" );

            // get the integration cell of interest
            moris::mtk::Cell const *tIntegrationCell = this->get_primary_cells_in_cluster()( aPrimaryCellClusterIndex );

            // get the vertex pointers on the side - for the bulk this is all vertices on the integration cell
            Vector< moris::mtk::Vertex * > tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (n_node x dim_xsi)
            Matrix< DDRMat > tVertexParamCoords( tVerticesOnCell.size(), this->get_dim_of_param_coord() );

            // iterate through vertices and collect local coordinates
            for ( moris::uint i = 0; i < tVerticesOnCell.size(); i++ )
            {
                tVertexParamCoords.get_row( i ) = this->get_vertex_local_coordinate_wrt_interp_cell( tVerticesOnCell( i ) ).get_row( 0 );
            }

            return tVertexParamCoords;
        }
    }

    //----------------------------------------------------------------
    Interpolation_Cell_Unzipped const *
    Cell_Cluster::get_xtk_interpolation_cell() const
    {
        return mInterpolationCell;
    }

    //----------------------------------------------------------------

    Matrix< IndexMat > Cell_Cluster::get_hanging_nodes() const
    {
        MORIS_ERROR( 0, "FIXME" );
        return mChildMesh->get_hanging_nodes();
    }

    //----------------------------------------------------------------

    size_t
    Cell_Cluster::capacity()
    {
        size_t tTotalSize = 0;
        tTotalSize += sizeof( mTrivial );
        tTotalSize += sizeof( mInterpolationCell );
        tTotalSize += sizeof( mChildMesh );
        tTotalSize += mPrimaryIntegrationCells.capacity();
        tTotalSize += mVoidIntegrationCells.capacity();
        tTotalSize += mVerticesInCluster.capacity();
        return tTotalSize;
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_primary_integration_cell_group( const std::shared_ptr< IG_Cell_Group > &aPrimaryIgCells )
    {
        mPrimaryIgCellGroup = { aPrimaryIgCells };

        mPrimaryIntegrationCells.resize( aPrimaryIgCells->mIgCellGroup.size() );

        for ( moris::uint i = 0; i < aPrimaryIgCells->mIgCellGroup.size(); i++ )
        {
            mPrimaryIntegrationCells( i ) = aPrimaryIgCells->mIgCellGroup( i );
        }
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_primary_integration_cell_groups( Vector< std::shared_ptr< IG_Cell_Group > > aPrimaryIgCells )
    {
        // store IG cell groups with cell cluster
        mPrimaryIgCellGroup = aPrimaryIgCells;

        // count total number of IG cells in all groups passed into function
        moris::uint tCount = 0;
        for ( moris::uint iCellGroup = 0; iCellGroup < aPrimaryIgCells.size(); iCellGroup++ )
        {
            tCount = tCount + aPrimaryIgCells( iCellGroup )->mIgCellGroup.size();
        }

        // initialize list of IG cells
        mPrimaryIntegrationCells.resize( tCount );

        // reset counter to track position in list
        tCount = 0;

        // store IG cells in list
        for ( moris::uint iCellGroup = 0; iCellGroup < aPrimaryIgCells.size(); iCellGroup++ )
        {
            for ( moris::uint jCellInGroup = 0; jCellInGroup < aPrimaryIgCells( iCellGroup )->mIgCellGroup.size(); jCellInGroup++ )
            {
                mVoidIntegrationCells( tCount ) = aPrimaryIgCells( iCellGroup )->mIgCellGroup( jCellInGroup );
                tCount++;
            }
        }
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_void_integration_cell_groups( Vector< std::shared_ptr< IG_Cell_Group > > &aVoidIgCells )
    {
        mVoidIgCellGroup = aVoidIgCells;

        moris::uint tCount = 0;
        for ( moris::uint i = 0; i < aVoidIgCells.size(); i++ )
        {
            tCount = tCount + aVoidIgCells( i )->mIgCellGroup.size();
        }

        mVoidIntegrationCells.resize( tCount );

        tCount = 0;

        for ( moris::uint i = 0; i < aVoidIgCells.size(); i++ )
        {
            for ( moris::uint j = 0; j < aVoidIgCells( i )->mIgCellGroup.size(); j++ )
            {
                mVoidIntegrationCells( tCount++ ) = aVoidIgCells( i )->mIgCellGroup( j );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster::set_ig_vertex_group( const std::shared_ptr< IG_Vertex_Group > &aVertexGroup )
    {
        mVertexGroup = aVertexGroup;

        mVerticesInCluster.resize( mVertexGroup->size() );
        mLocalCoords.resize( mVertexGroup->size(), mVertexGroup->get_vertex_local_coords( aVertexGroup->get_vertex( 0 )->get_index() )->n_cols() );

        for ( moris::uint i = 0; i < mVertexGroup->size(); i++ )
        {
            mVerticesInCluster( i ) = mVertexGroup->get_vertex( i );
            mLocalCoords.set_row( i, *mVertexGroup->get_vertex_local_coords( mVerticesInCluster( i )->get_index() ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster::set_quadrature_points( const uint aOrder, const uint aDim )
    {
    
        //MORIS_ASSERT( aDim < 3, "Currently moment fitting only works for 2D problems" );
        if ( aDim == 2 )
        {
            if ( aOrder == 1 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else if ( aOrder == 2 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_3x3 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else if ( aOrder == 3 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_4x4 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else
            {
                MORIS_ASSERT( aOrder < 4, "Only 3rd order supported currently");
            }

            


        }

        else if ( aDim == 3 )
        {
            if ( aOrder == 1 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_2x2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else if ( aOrder == 2 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_3x3x3 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else if ( aOrder == 3 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_4x4x4 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKPoints;
        
                tIntData.get_points( tMTKPoints );

                mQuadraturePoints = tMTKPoints;

            }
            else
            {
                MORIS_ASSERT( aOrder < 4, "Only 3rd order supported currently");
            }

        
        
        }
        //MORIS_ASSERT( aOrder < 2, "Currently moment fitting only works for linear problems" );

        

        // Get weights from mtk::integrator


            
        /*Vector< real > mOneDPoints;

        if (mOrder == 1)
        {
            mOneDPoints = {1.0/std::sqrt(3.0) , -1.0/std::sqrt(3.0)};
        }
        else
        {
            mOneDPoints =  {0.0, std::sqrt(3.0/5.0) , -std::sqrt(3.0/5.0)};
        }
        
        MORIS_ASSERT(mOrder < 2, "Only 2nd order supported currently");
        
        mQuadraturePoints.resize( mOneDPoints.size()*mOneDPoints.size() , 2 );

        for (uint iRowIndex = 0; iRowIndex < mOneDPoints.size() ; iRowIndex ++ )
        {
            for (uint iColIndex = 0; iColIndex < mOneDPoints.size() ; iColIndex ++ )
            {
                uint iMatIndex = 2*(iRowIndex) + iColIndex;
                
                mQuadraturePoints( iMatIndex , 0 ) = mOneDPoints( iRowIndex );

                mQuadraturePoints( iMatIndex , 1 ) = mOneDPoints( iColIndex );

            }
        }*/

        
        
    

    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cell_Cluster::get_quadrature_points() const
    {
        return mQuadraturePoints;
    }

    //------------------------------------------------------------------------------


    // This function for the trivial case of the volume fraction being 1.

    void
    Cell_Cluster::set_quadrature_weights( const uint aOrder, const uint aDim )
    {
        // if 3D problem throw exception

        //MORIS_ASSERT(aDim < 3, "Currently moment fitting only works for 2D problems");

        if ( aDim == 2 )
        {

            if ( aOrder == 1 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else if ( aOrder == 2 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_3x3 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else if ( aOrder == 3 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_4x4 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat > tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else
            {
                MORIS_ASSERT( aOrder < 4, "Only 3rd order supported currently");
            }


        }
        else if ( aDim == 3 )
        {

            if ( aOrder == 1 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_2x2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else if ( aOrder == 2 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_3x3x3 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat >  tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else if ( aOrder == 3 )
            {
                mtk::Integration_Rule tIntObj( mtk::Geometry_Type::HEX , mtk::Integration_Type::GAUSS , mtk::Integration_Order::HEX_4x4x4 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
                const mtk::Integrator tIntData( tIntObj );

                Matrix< DDRMat > tMTKWeights;
        
                tIntData.get_weights( tMTKWeights );

                mQuadratureWeights = tMTKWeights;

            }
            else
            {
                MORIS_ASSERT( aOrder < 4, "Only 3rd order supported currently");
            }


        }
        
        
        //mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
        //mtk::Integrator tIntData( tIntObj );

        // Get weights from mtk::integrator

        //Matrix< DDRMat >  tMTKWeights;
        
        //tIntData.get_weights( tMTKWeights );
        
        //Vector< real > mOneDWeights;

        // if (mOrder == 1)
        // {
        //     mOneDWeights = { 1.0 , 1.0 };
        // }
        // else
        // {
        //     mOneDWeights =  { 8.0/9.0, 5.0/9.0 , 5.0/9.0 };
        // }
        
        //MORIS_ASSERT( aOrder < 2, "Only 1st order supported currently");
        
        /*mQuadratureWeights.resize( mOneDWeights.size()*mOneDWeights.size() , 2 );

        for (uint iRowIndex = 0; iRowIndex < mOneDWeights.size() ; iRowIndex ++ )
        {
            for (uint iColIndex = 0; iColIndex < mOneDWeights.size() ; iColIndex ++ )
            {
                uint iMatIndex = 2*(iRowIndex) + iColIndex;
                
                mQuadratureWeights( iMatIndex ) = mOneDWeights( iRowIndex ) * mOneDWeights( iColIndex );

            }
        }*/
        

        

    }

    //------------------------------------------------------------------------------

    void 
    Cell_Cluster::find_subphase_boundary_vertices(
                const std::shared_ptr< IG_Cell_Group >  aSubphaseIGCells,
                const std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
                const uint aDim
                )
    {
        // Get the primary subphase IG cells first
        Vector< mtk::Cell* > tSubphaseIgCellsPtr = aSubphaseIGCells->mIgCellGroup;

        // Initialize the vector containing the facets
        Vector< moris_index > tSubphaseFacets;

        // Now run a for loop to extract the facets corresponding to the subphase
        for (uint iCellIndex = 0; iCellIndex < tSubphaseIgCellsPtr.size(); iCellIndex++)
        {
            // Get the global cell index
            moris_index tCellIndexGlobal = tSubphaseIgCellsPtr( iCellIndex )->get_index();

            // Get facet connectivity specific cell index
            moris_index tCellIndexLocal = aFacetConnectivity->get_cell_ordinal( tCellIndexGlobal );

            // Get facets attached to this cell
            tSubphaseFacets.append(aFacetConnectivity->mCellToFacet( tCellIndexLocal ));

        }

        // Now filter out the facets only belonging to one cell
        for (uint iFacetIndex =0; iFacetIndex < tSubphaseFacets.size(); iFacetIndex++ )
        {
            // Get facet
            moris_index tSingleFacet = tSubphaseFacets( iFacetIndex );

            // Now check how many elements does the facet belong to
            Vector< mtk::Cell* > tOwningElements = aFacetConnectivity->mFacetToCell( tSingleFacet );

            // Only retain cells in the subphase from the facet to cell map
            Vector< mtk::Cell* > tOwningElementsInSubphase;
            
            // For facet to cell edge ordinal map
            uint tCellIndexForFacet = 0;

            // Now check if the facet belongs to only one cell in the subphase
            for (uint iSubphaseCellIndex = 0; iSubphaseCellIndex < tSubphaseIgCellsPtr.size(); iSubphaseCellIndex++)
            {
                // Get cell
                mtk::Cell* tSubphaseCellToCheck = tSubphaseIgCellsPtr( iSubphaseCellIndex );

                //fprintf( stdout , "Subphase cell index %d\n", tSubphaseCellToCheck->get_index() );

                // Get local cell index
                //moris_index tSubphaseCellIndexLocal = aFacetConnectivity->get_cell_ordinal( tSubphaseCellToCheck->get_index() );

                for (uint iOwningCells = 0; iOwningCells < tOwningElements.size() ; iOwningCells++) 
                {                    
                    // Check if subphase cell part of facet owning cells
                    if ( tOwningElements( iOwningCells )->get_index() == tSubphaseCellToCheck->get_index() )
                    {
                        // Add to vector containing only subphase cells sharing facet
                        tOwningElementsInSubphase.push_back( tOwningElements( iOwningCells ) );

                        // Get cell index value for facet to cell ordinal map
                        tCellIndexForFacet = iOwningCells;
                               
                    } 

                }

            } 

            // If belonging to one cell then add it into the vector containing pointer to facet indices.
            if ( tOwningElementsInSubphase.size()  == 1 )
            {
                mFacetVerticesOnSubphaseBoundary.push_back( aFacetConnectivity->mFacetVertices( tSingleFacet ) );

                // Get all possible cell ordinals for this facet
                Vector< moris_index > tFacetOrdinals = aFacetConnectivity->mFacetToCellEdgeOrdinal( tSingleFacet ) ;

                // Get ordinal corresponding to given subphase cell _and_ facet
                moris_index tCellOrdinal = tFacetOrdinals( tCellIndexForFacet );

                // Compute normal for this facet based on the side ordinal data
                Matrix< DDRMat > tFacetNormal = tOwningElementsInSubphase( 0 )->compute_outward_side_normal( tCellOrdinal );

                // Add to vector containing facet normals
                mFacetNormals.push_back( tFacetNormal );

                // Get the (ordered) vertices for this facet
                Vector< moris::mtk::Vertex const * > tVerticesOnFacet = tOwningElementsInSubphase( 0 )->get_vertices_on_side_ordinal( tCellOrdinal );

                Matrix< DDRMat > tVertexCoordinatesFacetOrdered;

                if ( aDim == 2 )
                {
                    tVertexCoordinatesFacetOrdered.reshape( tVerticesOnFacet.size() , 2 );
                }
                else if ( aDim == 3 )
                {
                    tVertexCoordinatesFacetOrdered.reshape( tVerticesOnFacet.size() , 3 );
                }
                

                // Get the vertex local coordinates and add to mFacetCoordinates
                for ( uint iVertex = 0; iVertex < tVerticesOnFacet.size(); iVertex++ )
                {
                    // Get vertex
                    moris::mtk::Vertex const * tVertexInFacet = tVerticesOnFacet( iVertex );

                    // Get local coordinates
                    Matrix< DDRMat > tVertexCoords = get_vertex_local_coordinate_wrt_interp_cell( tVertexInFacet );

                    // Add to Matrix
                    tVertexCoordinatesFacetOrdered.set_row( iVertex , tVertexCoords.get_row( 0 ) );

                }

                mFacetVertexCoordinates.push_back( tVertexCoordinatesFacetOrdered );

                // get the vertex coordinates
                //moris::Matrix< moris::DDRMat > tVertexCoords = tOwningElementsInSubphase( 0 )->get_vertex_coords();

                // Get the nodes which need to be used
                //moris::Matrix< moris::IndexMat > tEdgeNodesForNormal = tOwningElementsInSubphase( 0 )->get_cell_info()->get_node_map_outward_normal( tCellOrdinal );

                // Declare node matrix
                //Matrix < DDRMat > tFacetVertexCoordinates;
                //tFacetVertexCoordinates.reshape( 2 , 2 );
                
                // Order the nodes
                //moris_index tFirstNode  = tEdgeNodesForNormal( 0 );
                //moris_index tSecondNode = tEdgeNodesForNormal( 1 );
                //tFacetVertexCoordinates( 0 , 0 ) = tVertexCoords( tFirstNode, 0 );
                //tFacetVertexCoordinates( 1 , 0 ) = tVertexCoords( tSecondNode, 0 );
                //tFacetVertexCoordinates( 0 , 1 ) = tVertexCoords( tFirstNode, 1 );
                //tFacetVertexCoordinates( 1 , 1 ) = tVertexCoords( tSecondNode, 1 );

                // Add the vertex coordinate matrix to the vector
                //mFacetVertexCoordinates.push_back( tFacetVertexCoordinates );
                
            }

            

        }  

        // Code below for testing the facet normal computation function. Not needed otherwise, but still kept.
        
        Matrix < DDRMat > tFacetNormals;

        if ( aDim == 2 )
        {
            tFacetNormals.reshape( mFacetNormals.size() , 2 );
        }
        else if ( aDim == 3 )
        {
            tFacetNormals.reshape( mFacetNormals.size() , 3 );
        }


        


        for (uint iFacetNormalIndex = 0; iFacetNormalIndex < mFacetNormals.size(); iFacetNormalIndex++ )
        {
            Matrix < DDRMat > tFacetNormal = mFacetNormals( iFacetNormalIndex );

            Matrix < DDRMat > tFacetNormalReshaped  = trans( tFacetNormal );
            
            tFacetNormals.set_row( iFacetNormalIndex , tFacetNormalReshaped );
                
        }

        Matrix < DDRMat > tFacetCoords;

        if ( aDim == 2 )
        {
            tFacetCoords.reshape( 2*mFacetVertexCoordinates.size() , 2 );
        }
        else if ( aDim == 3 )
        {
            tFacetCoords.reshape( 3*mFacetVertexCoordinates.size() , 3 );
        }

        

        for (uint iFacetNormalIndex = 0; iFacetNormalIndex < mFacetVertexCoordinates.size(); iFacetNormalIndex++ )
        {
            Matrix < DDRMat > tFacetCoord = mFacetVertexCoordinates( iFacetNormalIndex );
            
            if ( aDim == 2 )
            {
                tFacetCoords.set_row( 2*iFacetNormalIndex , tFacetCoord.get_row( 0 ) );

                tFacetCoords.set_row( 2*iFacetNormalIndex + 1 , tFacetCoord.get_row( 1 ) );

            }
            else if ( aDim == 3 )
            {
                tFacetCoords.set_row( 3*iFacetNormalIndex , tFacetCoord.get_row( 0 ) );

                tFacetCoords.set_row( 3*iFacetNormalIndex + 1 , tFacetCoord.get_row( 1 ) );

                tFacetCoords.set_row( 3*iFacetNormalIndex + 2 , tFacetCoord.get_row( 2 ) );

            }
            
                
        }



        // The code below was all for testing the facet filtering function. Not needed otherwise, but still kept.

        std::vector< double > tFacetCoordinatesFileVectorX;

        std::vector< double > tFacetCoordinatesFileVectorY;
        
        std::vector< double > tFacetCoordinatesFileVectorZ;
        


        for(uint iVertInd = 0; iVertInd < mFacetVerticesOnSubphaseBoundary.size(); iVertInd++)
        {
             // Get vertex
             
             moris::Vector< moris::mtk::Vertex* > tFacetIndexBdry = mFacetVerticesOnSubphaseBoundary( iVertInd ) ;

             for (uint iFacetVertInd = 0; iFacetVertInd < tFacetIndexBdry.size() ; iFacetVertInd++ )
             
             {
                // Get one facet vertex
                moris::mtk::Vertex* tBdryFacet = tFacetIndexBdry( iFacetVertInd );
                
                // Get coords
                Matrix<DDRMat> tBdryCoords = tBdryFacet->get_coords();

                tFacetCoordinatesFileVectorX.push_back(tBdryCoords( 0 , 0 ));

                tFacetCoordinatesFileVectorY.push_back(tBdryCoords( 0 , 1 ));

                if ( aDim == 3 )
                {
                    tFacetCoordinatesFileVectorZ.push_back(tBdryCoords( 0 , 2 ));;
                }

             }

                   
                   
        
        }  

        hid_t  tFileID = create_hdf5_file( "FacetVertices_X.hdf5" );
        herr_t tStatus = 0;
        save_vector_to_hdf5_file( tFileID, std::string("Coords"), tFacetCoordinatesFileVectorX, tStatus );

        close_hdf5_file( tFileID );


        hid_t  tFileID1 = create_hdf5_file( "FacetVertices_Y.hdf5" );
        herr_t tStatus1 = 0;
        save_vector_to_hdf5_file( tFileID1, std::string("Coords"), tFacetCoordinatesFileVectorY, tStatus1 );

        close_hdf5_file( tFileID1 );

        hid_t  tFileID2 = create_hdf5_file( "FacetNormals_new.hdf5" );
        herr_t tStatus2 = 0;
        save_matrix_to_hdf5_file( tFileID2, std::string("Coords"), tFacetNormals, tStatus2 );

        close_hdf5_file( tFileID2 );

        hid_t  tFileID3 = create_hdf5_file( "FacetVertices_New.hdf5" );
        herr_t tStatus3 = 0;
        save_matrix_to_hdf5_file( tFileID3, std::string("Coords"), tFacetCoords, tStatus3 );

        close_hdf5_file( tFileID3 );

        hid_t  tFileID4 = create_hdf5_file( "FacetVertices_Z.hdf5" );
        herr_t tStatus4 = 0;
        save_vector_to_hdf5_file( tFileID4, std::string("Coords"), tFacetCoordinatesFileVectorZ, tStatus4 );

        close_hdf5_file( tFileID4 );

        

    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::compute_quadrature_weights( const uint aOrder, const uint aDim )
    {
        //MORIS_ASSERT( aDim < 3 , "Only 2D supported at present");

        // Get midpoint of facet coordinates -

        /*Matrix< DDRMat > tFacetMidpoints;

        tFacetMidpoints.reshape( mFacetVerticesOnSubphaseBoundary.size() , aDim );

        for (uint iFacetIndex = 0; iFacetIndex < mFacetVerticesOnSubphaseBoundary.size(); iFacetIndex++)
        {
            // Get coordinates corresponding to one facet -
            Vector< mtk::Vertex* > tFacetCoords = mFacetVerticesOnSubphaseBoundary( iFacetIndex );

            // Get coordinates and sum them to get midpoint
            mtk::Vertex* tFirstFacetCoord = tFacetCoords( 0 );

            mtk::Vertex* tSecondFacetCoord = tFacetCoords( 1 );

            mtk::Vertex* tThirdFacetCoord;
            
            if ( aDim == 3 )
            {
                tThirdFacetCoord = tFacetCoords( 2 );                
            }

            Matrix< DDRMat > tFacetMidpointCoords;
            
            if ( aDim == 2 )
            {
                tFacetMidpointCoords = 0.5*(tFirstFacetCoord->get_coords() + tSecondFacetCoord->get_coords()) ;

            }
            else if ( aDim == 3 )
            {
                tFacetMidpointCoords = ( 1.0 / 3.0 )*(tFirstFacetCoord->get_coords() + tSecondFacetCoord->get_coords() + tThirdFacetCoord->get_coords() ) ;
            }

            tFacetMidpoints.set_row( iFacetIndex, tFacetMidpointCoords ) ;
            
        }*/

        // Now, create interpolation objects for specifying the basis function for each facet. One for the IP and one for the geometry.
        //mtk::Interpolation_Rule tInterpolationRule( mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

        if ( aDim == 2 )
        {
            if ( aOrder == 1 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::QUAD , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else if ( aOrder == 2 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::QUAD , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::QUADRATIC , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else if ( aOrder == 3 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::QUAD , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::CUBIC , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else
            {

              MORIS_ASSERT( aOrder < 4, "Only upto 3rd order supported currently");

            }

        }
        else if ( aDim == 3 )
        {
            if ( aOrder == 1 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::HEX , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else if ( aOrder == 2 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::HEX , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::QUADRATIC , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else if ( aOrder == 3 )
            {
                mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::HEX , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::CUBIC , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

                mIPInterp = tIPInterpolationRule.create_space_interpolation_function();

            }
            else
            {

              MORIS_ASSERT( aOrder < 4, "Only upto 3rd order supported currently");

            }

        }
        

        

        // Declare cellshape
        //mtk::CellShape tCellShape = mtk::CellShape::GENERAL;

        // Set space sideset to true for interpolation to start
        //bool tSpaceSideset = true;

        // Create space interpolator object
        //mtk::Space_Interpolator tSpaceInterpolationObject( tInterpolationRule , tIPInterpolationRule , tCellShape , tSpaceSideset );

        //Matrix< DDRMat > tFacetNormals;

        /*tFacetNormals.reshape( mFacetVerticesOnSubphaseBoundary.size() , 2 );
        
        // compute the normals of all the boundary facets.
        for (uint iFacetMidpointIndex = 0; iFacetMidpointIndex < mFacetVerticesOnSubphaseBoundary.size(); iFacetMidpointIndex++ )
        {
            // Location at which to compute the normal in parent coordinates (midpoint of facet)
            Matrix< DDRMat > aXiLocal = {{0}};

            // Define parent element coordinate
            tSpaceInterpolationObject.set_space( aXiLocal );

            // Get facet end coordinates (physical coordinates of line element)
            Vector< mtk::Vertex* > tFacetCoords = mFacetVerticesOnSubphaseBoundary( iFacetMidpointIndex );

            // Get coordinates 
            Matrix< DDRMat > tFirstFacetCoord = tFacetCoords( 0 )->get_coords();

            Matrix< DDRMat > tSecondFacetCoord = tFacetCoords( 1 )->get_coords();

            Matrix< DDRMat > tXiHat;

            tXiHat.reshape( 2 , 2 );

            tXiHat.set_row( 0 , tFirstFacetCoord );

            tXiHat.set_row( 1 , tSecondFacetCoord );
            
            // Specify physical coordinate
            tSpaceInterpolationObject.set_space_coeff( tXiHat );

            // Declare facet normal vector
            Matrix< DDRMat > tFacetNormal;
            
            // Compute normal
            tSpaceInterpolationObject.get_normal( tFacetNormal );

            // Store the computed normals
            tFacetNormals.set_row( iFacetMidpointIndex ,  trans( tFacetNormal ) );


        }

        // For inspection purposes only - check the value of the computed facet normals
        hid_t  tFileID2 = create_hdf5_file( "FacetNormals.hdf5" );
        herr_t tStatus2 = 0;
        save_matrix_to_hdf5_file( tFileID2, std::string("Coords"), tFacetNormals, tStatus2 );

        close_hdf5_file( tFileID2 );*/

        //Declare LHS matrix
        Matrix < DDRMat > tMomentFittingLHS;

        // Determine the number of moments
        uint tNmoments = std::pow( aOrder + 1 , aDim ); 

        tMomentFittingLHS.reshape( tNmoments , tNmoments );
            
        // Generate LHS
        for (uint iQuadPointIndex = 0; iQuadPointIndex < tMomentFittingLHS.n_cols() ; iQuadPointIndex++)
        {
            // Declare matrix for basis function values
            Matrix< DDRMat > tN;

            // Get quad point
            Matrix< DDRMat > tXi = mQuadraturePoints.get_column( iQuadPointIndex );

            // Get value of basis functions at quad point
            mIPInterp->eval_N( tXi , tN );

            // Place it in LHS 
            tMomentFittingLHS.set_column( iQuadPointIndex , trans( tN ) );

        }

        // create RHS vector
        Matrix < DDRMat > tMomentFittingRHS;
        
        // Allocate memory for RHS
        tMomentFittingRHS.reshape( tNmoments , 1 );
        
        // Compute the RHS
        for (uint iFacetIndex = 0; iFacetIndex < mFacetVerticesOnSubphaseBoundary.size() ; iFacetIndex++)
        {
            // Get subphase facet coordinates
            Matrix< DDRMat > tFacetCoords = mFacetVertexCoordinates( iFacetIndex );

            // Get facet normal
            Matrix< DDRMat > tFacetNormal = mFacetNormals( iFacetIndex );

            // Get individual facet normal coordinates
            real tNx  = tFacetNormal( 0 , 0 );
            real tNy  = tFacetNormal( 1 , 0 );
            real tNz;
            if ( aDim == 3 )
            {
                tNz = tFacetNormal( 2 , 0 );
            }

            // Get quadrature object
            mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

            mtk::Integration_Order tIntegrationOrder = mtk::Integration_Order::BAR_1;

            if ( aDim == 2 )
            {
                tGeometryType = mtk::Geometry_Type::LINE ;

                tIntegrationOrder = mtk::Integration_Order::BAR_32 ;

            }
            else if ( aDim == 3 )
            {
                tGeometryType = mtk::Geometry_Type::TRI ;

                tIntegrationOrder = mtk::Integration_Order::TRI_12 ;

            }

            mtk::Integration_Rule tIntObjLine( tGeometryType , mtk::Integration_Type::GAUSS , tIntegrationOrder , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
            
            mtk::Integrator tIntDataLine( tIntObjLine );
            
            // Get quadrature points (line or surface)
            Matrix< DDRMat > tIntPointsLine;

            tIntDataLine.get_points( tIntPointsLine );

            // Get line quadrature weights
            Matrix< DDRMat > tIntWeightsLine;

            tIntDataLine.get_weights( tIntWeightsLine );

            // Get geometric jacobian
            real tD ;
            
            if ( aDim == 2 )
            {
                tD  = std::sqrt(std::pow( tFacetCoords( 1, 0 ) - tFacetCoords( 0, 0 ) , 2) + std::pow( tFacetCoords( 1, 1 ) - tFacetCoords( 0, 1 ) , 2 ));
            }

            
            // Loop over line quadrature points to compute the integral
            for( uint iQuadPtIndex = 0 ; iQuadPtIndex < tIntWeightsLine.numel() ; iQuadPtIndex++ )
            {
                // Get quad point
                Matrix< DDRMat > tQuadPoint = tIntPointsLine.get_column( iQuadPtIndex );

                // Get quad weight
                real tQuadWeight = tIntWeightsLine( iQuadPtIndex );

                real tXm;
                real tYm;
                real tZm;

                Matrix< DDRMat > tXvector ; 
                Matrix< DDRMat > tYvector ;
                Matrix< DDRMat > tZvector ; 
                if ( aDim == 2 )
                {
                    tXvector  = {{ tFacetCoords( 0 , 0 ) },{ tFacetCoords( 1 , 0 ) }};
                    tYvector  = {{ tFacetCoords( 0 , 1 ) },{ tFacetCoords( 1 , 1 ) }};

                }
                if ( aDim == 3 )
                {
                    tXvector = {{ tFacetCoords( 0 , 0 ) },{ tFacetCoords( 1 , 0 ) },{ tFacetCoords( 2 , 0 ) }};
                    tYvector = {{ tFacetCoords( 0 , 1 ) },{ tFacetCoords( 1 , 1 ) },{ tFacetCoords( 2 , 1 ) }};
                    tZvector = {{ tFacetCoords( 0 , 2 ) },{ tFacetCoords( 1 , 2 ) },{ tFacetCoords( 2 , 2 ) }};

                } 
                


                if ( aDim == 2 )
                {
                    // Get the mapped quad point value x
                    tXm = 0.5 * ( 1.0 - tQuadPoint( 0 ) ) * tFacetCoords( 0 , 0 ) + 0.5 * ( 1.0 + tQuadPoint( 0 ) ) * tFacetCoords( 1 , 0 );

                    // Get the mapped quad point value x
                    tYm = 0.5 * ( 1.0 - tQuadPoint( 0 ) ) * tFacetCoords( 0 , 1 ) + 0.5 * ( 1.0 + tQuadPoint( 0 ) ) * tFacetCoords( 1 , 1 );
                }
                else if ( aDim == 3 )
                {
                    // Define interpolation rule
                    mtk::Interpolation_Rule tIGInterpolationRule( mtk::Geometry_Type::TRI , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );
                    mtk::Interpolation_Function_Base* tIGInterp = tIGInterpolationRule.create_space_interpolation_function();

                    // Declare matrix for basis functions
                    Matrix< DDRMat > tNmapping;

                    // Get the mapped quad point value x
                    tIGInterp->eval_N( tQuadPoint , tNmapping );                    

                    // get mapped quad point value x
                    Matrix< DDRMat > tXmat = tNmapping * tXvector;
                    tXm = tXmat( 0 );

                    // get mapped quad point value y
                    Matrix< DDRMat > tYmat = tNmapping * tYvector;
                    tYm = tYmat( 0 );

                    // get mapped quad point value z
                    Matrix< DDRMat > tZmat = tNmapping * tZvector;
                    tZm = tZmat( 0 );

                    // Allocate matrix for storing mapping
                    Matrix< DDRMat > tNximapping;

                    // Compute geometric Jacobian
                    tIGInterp->eval_dNdXi( tQuadPoint , tNximapping );

                    Matrix< DDRMat > tRxi_x = tNximapping.get_row(0) * tXvector;
                    Matrix< DDRMat > tRxi_y = tNximapping.get_row(0) * tYvector;
                    Matrix< DDRMat > tRxi_z = tNximapping.get_row(0) * tZvector;

                    Matrix< DDRMat > tReta_x = tNximapping.get_row(1) * tXvector;
                    Matrix< DDRMat > tReta_y = tNximapping.get_row(1) * tYvector;
                    Matrix< DDRMat > tReta_z = tNximapping.get_row(1) * tZvector;
                    

                    Matrix< DDRMat > tRxi = {{ tRxi_x( 0 ) }, { tRxi_y( 0 ) }, { tRxi_z( 0 ) }};

                    Matrix< DDRMat > tReta = {{ tReta_x( 0 ) }, { tReta_y( 0 ) }, { tReta_z( 0 ) }};

                    Matrix < DDRMat > tCrossProdRxiReta = cross( tRxi , tReta );

                    tD = std::sqrt( std::pow(tCrossProdRxiReta( 0 ), 2) + std::pow(tCrossProdRxiReta( 1 ), 2) + std::pow(tCrossProdRxiReta( 2 ),2) ); 

                }

                

                // Evaluate the moments
                if ( aOrder == 1 && aDim == 2)
                {
                    tMomentFittingRHS( 0 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm - 0.5 * ( tXm * tXm ) ) * ( 1.0 - tYm ) + 0.25 * tNy * ( tYm - 0.5 * ( tYm * tYm ) ) * ( 1.0 - tXm ) ) ;

                    tMomentFittingRHS( 1 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm + 0.5 * ( tXm * tXm ) ) * ( 1.0 - tYm ) + 0.25 * tNy * ( tYm - 0.5 * ( tYm * tYm ) ) * ( 1.0 + tXm ) ) ;

                    tMomentFittingRHS( 2 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm + 0.5 * ( tXm * tXm ) ) * ( 1.0 + tYm ) + 0.25 * tNy * ( tYm + 0.5 * ( tYm * tYm ) ) * ( 1.0 + tXm ) ) ;

                    tMomentFittingRHS( 3 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( tXm - 0.5 * ( tXm * tXm ) ) * ( 1.0 + tYm ) + 0.25 * tNy * ( tYm + 0.5 * ( tYm * tYm ) ) * ( 1.0 - tXm ) ) ;

                }
                else if ( aOrder == 2 && aDim == 2 )
                {
                    tMomentFittingRHS( 0 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( -0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm ) * ( -tYm ) + 0.25 * tNy * ( -0.5 * std::pow( tYm, 2 ) + (1.0/3.0) * std::pow( tYm, 3 ) ) * ( 1.0 - tXm ) * ( -tXm ) ) ;

                    tMomentFittingRHS( 1 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * (  0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm ) * ( -tYm ) + 0.25 * tNy * ( -0.5 * std::pow( tYm, 2 ) + (1.0/3.0) * std::pow( tYm, 3 ) ) * ( 1.0 + tXm ) * (  tXm ) ) ;

                    tMomentFittingRHS( 2 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * (  0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) * ( 1.0 + tYm ) * (  tYm ) + 0.25 * tNy * (  0.5 * std::pow( tYm, 2 ) + (1.0/3.0) * std::pow( tYm, 3 ) ) * ( 1.0 + tXm ) * (  tXm ) ) ;

                    tMomentFittingRHS( 3 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.25 * tNx * ( -0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) * ( 1.0 + tYm ) * (  tYm ) + 0.25 * tNy * (  0.5 * std::pow( tYm, 2 ) + (1.0/3.0) * std::pow( tYm, 3 ) ) * ( 1.0 - tXm ) * ( -tXm ) ) ;

                    tMomentFittingRHS( 4 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( 1.0 - tYm ) * ( -tYm ) * ( tXm - ( 1.0/3.0 ) * std::pow( tXm, 3 ) ) + 0.5 * tNy * ( 1.0 - tXm * tXm )*( -0.5 * std::pow( tYm, 2 ) + (1.0/3.0) * std::pow( tYm, 3 ) ) ) ;

                    tMomentFittingRHS( 5 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * (1.0 - tYm * tYm) * (  0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) + 0.5 * tNy * (tXm + 1.0 ) * ( tXm ) * ( tYm - ( 1.0/3.0 )*( std::pow(tYm, 3) ) ) ) ;

                    tMomentFittingRHS( 6 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( tXm - ( 1.0 / 3.0 ) * ( std::pow( tXm , 3 ) ) ) * ( 1.0 + tYm ) * ( tYm ) + 0.5 * tNy * ( 1.0 - tXm * tXm ) * ( 0.5 * std::pow( tYm , 2 ) + ( 1.0 / 3.0 ) * ( std::pow( tYm , 3 ) ) ) ) ;

                    tMomentFittingRHS( 7 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 0.5 * tNx * ( -0.5 * std::pow( tXm, 2 ) + (1.0/3.0) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm * tYm ) + 0.5 * tNy * ( 1.0 - tXm ) * ( -tXm ) * ( tYm - ( 1.0/3.0 )*( std::pow(tYm, 3) ) ) ) ; 

                    tMomentFittingRHS( 8 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * ( 1.0 * tNx * ( tXm - ( 1.0/3.0 ) * std::pow( tXm, 3 ) ) * ( 1.0 - tYm * tYm ) + 1.0 * tNy * ( tYm - ( 1.0/3.0 )*( std::pow(tYm, 3) ) ) * ( 1.0 - tXm * tXm ) ) ;

                } 
                else if ( aOrder == 3 && aDim == 2 )
                {
                    real tt2 = std::pow(tXm,2);
                    real tt3 = std::pow(tXm,3);
                    real tt4 = std::pow(tYm,2);
                    real tt5 = std::pow(tYm,3);

                    tMomentFittingRHS( 0 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tXm,2)*(tt4*(9.0/5.12e+2)-tt5*(9.0/5.12e+2)+tYm/5.12e+2-1.0/5.12e+2)-std::pow(tXm,4)*(tt4*7.91015625e-2-tt5*7.91015625e-2+tYm*8.7890625e-3-8.7890625e-3)-tXm*(tt4*(9.0/2.56e+2)-tt5*(9.0/2.56e+2)+tYm/2.56e+2-1.0/2.56e+2) ) ;
                    tMomentFittingRHS( 0 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tYm,2)*(tt2*(9.0/5.12e+2)-tt3*(9.0/5.12e+2)+tXm/5.12e+2-1.0/5.12e+2)-std::pow(tYm,4)*(tt2*7.91015625e-2-tt3*7.91015625e-2+tXm*8.7890625e-3-8.7890625e-3)-tYm*(tt2*(9.0/2.56e+2)-tt3*(9.0/2.56e+2)+tXm/2.56e+2-1.0/2.56e+2)) ;

                    tMomentFittingRHS( 1 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tXm,2)*(tt4*(9.0/5.12e+2)-tt5*(9.0/5.12e+2)+tYm/5.12e+2-1.0/5.12e+2)+std::pow(tXm,4)*(tt4*7.91015625e-2-tt5*7.91015625e-2+tYm*8.7890625e-3-8.7890625e-3)-tXm*(tt4*(9.0/2.56e+2)-tt5*(9.0/2.56e+2)+tYm/2.56e+2-1.0/2.56e+2) ) ;
                    tMomentFittingRHS( 1 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tYm,2)*(tt2*(9.0/5.12e+2)+tt3*(9.0/5.12e+2)-tXm/5.12e+2-1.0/5.12e+2)-std::pow(tYm,4)*(tt2*7.91015625e-2+tt3*7.91015625e-2-tXm*8.7890625e-3-8.7890625e-3)-tYm*(tt2*(9.0/2.56e+2)+tt3*(9.0/2.56e+2)-tXm/2.56e+2-1.0/2.56e+2));

                    tMomentFittingRHS( 2 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tXm,2)*(tt4*(9.0/5.12e+2)+tt5*(9.0/5.12e+2)-tYm/5.12e+2-1.0/5.12e+2)+std::pow(tXm,4)*(tt4*7.91015625e-2+tt5*7.91015625e-2-tYm*8.7890625e-3-8.7890625e-3)-tXm*(tt4*(9.0/2.56e+2)+tt5*(9.0/2.56e+2)-tYm/2.56e+2-1.0/2.56e+2) ) ;
                    tMomentFittingRHS( 2 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tYm,2)*(tt2*(9.0/5.12e+2)+tt3*(9.0/5.12e+2)-tXm/5.12e+2-1.0/5.12e+2)+std::pow(tYm,4)*(tt2*7.91015625e-2+tt3*7.91015625e-2-tXm*8.7890625e-3-8.7890625e-3)-tYm*(tt2*(9.0/2.56e+2)+tt3*(9.0/2.56e+2)-tXm/2.56e+2-1.0/2.56e+2));

                    tMomentFittingRHS( 3 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tXm,2)*(tt4*(9.0/5.12e+2)+tt5*(9.0/5.12e+2)-tYm/5.12e+2-1.0/5.12e+2)-std::pow(tXm,4)*(tt4*7.91015625e-2+tt5*7.91015625e-2-tYm*8.7890625e-3-8.7890625e-3)-tXm*(tt4*(9.0/2.56e+2)+tt5*(9.0/2.56e+2)-tYm/2.56e+2-1.0/2.56e+2) );
                    tMomentFittingRHS( 3 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tYm,2)*(tt2*(9.0/5.12e+2)-tt3*(9.0/5.12e+2)+tXm/5.12e+2-1.0/5.12e+2)+std::pow(tYm,4)*(tt2*7.91015625e-2-tt3*7.91015625e-2+tXm*8.7890625e-3-8.7890625e-3)-tYm*(tt2*(9.0/2.56e+2)-tt3*(9.0/2.56e+2)+tXm/2.56e+2-1.0/2.56e+2) ) ;

                    tMomentFittingRHS( 4 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)-tt5*(2.43e+2/5.12e+2)+tYm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)+std::pow(tXm,4)*(tt4*2.373046875e-1-tt5*2.373046875e-1+tYm*2.63671875e-2-2.63671875e-2)+tXm*(tt4*(8.1e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(9.0/2.56e+2)-9.0/2.56e+2) );
                    tMomentFittingRHS( 4 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tYm,2)*(tt2*(9.0/5.12e+2)-tt3*(2.7e+1/5.12e+2)+tXm*(2.7e+1/5.12e+2)-9.0/5.12e+2)+std::pow(tYm,4)*(tt2*7.91015625e-2-tt3*2.373046875e-1+tXm*2.373046875e-1-7.91015625e-2)+tYm*(tt2*(9.0/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(2.7e+1/2.56e+2)-9.0/2.56e+2) );
                    
                    tMomentFittingRHS( 5 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)-tt5*(2.43e+2/5.12e+2)+tYm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)-std::pow(tXm,4)*(tt4*2.373046875e-1-tt5*2.373046875e-1+tYm*2.63671875e-2-2.63671875e-2)+tXm*(tt4*(8.1e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(9.0/2.56e+2)-9.0/2.56e+2) )  ;
                    tMomentFittingRHS( 5 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tYm,2)*(tt2*(9.0/5.12e+2)+tt3*(2.7e+1/5.12e+2)-tXm*(2.7e+1/5.12e+2)-9.0/5.12e+2)+std::pow(tYm,4)*(tt2*7.91015625e-2+tt3*2.373046875e-1-tXm*2.373046875e-1-7.91015625e-2)+tYm*(tt2*(9.0/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(2.7e+1/2.56e+2)-9.0/2.56e+2) ) ;
                    
                    tMomentFittingRHS( 6 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tXm,2)*(tt4*(9.0/5.12e+2)-tt5*(2.7e+1/5.12e+2)+tYm*(2.7e+1/5.12e+2)-9.0/5.12e+2)-std::pow(tXm,4)*(tt4*7.91015625e-2-tt5*2.373046875e-1+tYm*2.373046875e-1-7.91015625e-2)+tXm*(tt4*(9.0/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(2.7e+1/2.56e+2)-9.0/2.56e+2) )  ;
                    tMomentFittingRHS( 6 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)+tt3*(2.43e+2/5.12e+2)-tXm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)+std::pow(tYm,4)*(tt2*2.373046875e-1+tt3*2.373046875e-1-tXm*2.63671875e-2-2.63671875e-2)+tYm*(tt2*(8.1e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(9.0/2.56e+2)-9.0/2.56e+2) );

                    tMomentFittingRHS( 7 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tXm,2)*(tt4*(9.0/5.12e+2)+tt5*(2.7e+1/5.12e+2)-tYm*(2.7e+1/5.12e+2)-9.0/5.12e+2)-std::pow(tXm,4)*(tt4*7.91015625e-2+tt5*2.373046875e-1-tYm*2.373046875e-1-7.91015625e-2)+tXm*(tt4*(9.0/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(2.7e+1/2.56e+2)-9.0/2.56e+2) ) ;
                    tMomentFittingRHS( 7 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)+tt3*(2.43e+2/5.12e+2)-tXm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)-std::pow(tYm,4)*(tt2*2.373046875e-1+tt3*2.373046875e-1-tXm*2.63671875e-2-2.63671875e-2)+tYm*(tt2*(8.1e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(9.0/2.56e+2)-9.0/2.56e+2) ); 
                
                    tMomentFittingRHS( 8 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)+tt5*(2.43e+2/5.12e+2)-tYm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)-std::pow(tXm,4)*(tt4*2.373046875e-1+tt5*2.373046875e-1-tYm*2.63671875e-2-2.63671875e-2)+tXm*(tt4*(8.1e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(9.0/2.56e+2)-9.0/2.56e+2) )  ;
                    tMomentFittingRHS( 8 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tYm,2)*(tt2*(9.0/5.12e+2)+tt3*(2.7e+1/5.12e+2)-tXm*(2.7e+1/5.12e+2)-9.0/5.12e+2)-std::pow(tYm,4)*(tt2*7.91015625e-2+tt3*2.373046875e-1-tXm*2.373046875e-1-7.91015625e-2)+tYm*(tt2*(9.0/2.56e+2)+tt3*(2.7e+1/2.56e+2)-tXm*(2.7e+1/2.56e+2)-9.0/2.56e+2) )  ;

                    tMomentFittingRHS( 9 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)+tt5*(2.43e+2/5.12e+2)-tYm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)+std::pow(tXm,4)*(tt4*2.373046875e-1+tt5*2.373046875e-1-tYm*2.63671875e-2-2.63671875e-2)+tXm*(tt4*(8.1e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(9.0/2.56e+2)-9.0/2.56e+2) ) ;
                    tMomentFittingRHS( 9 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tYm,2)*(tt2*(9.0/5.12e+2)-tt3*(2.7e+1/5.12e+2)+tXm*(2.7e+1/5.12e+2)-9.0/5.12e+2)-std::pow(tYm,4)*(tt2*7.91015625e-2-tt3*2.373046875e-1+tXm*2.373046875e-1-7.91015625e-2)+tYm*(tt2*(9.0/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(2.7e+1/2.56e+2)-9.0/2.56e+2) );

                    tMomentFittingRHS( 10 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tXm,2)*(tt4*(9.0/5.12e+2)+tt5*(2.7e+1/5.12e+2)-tYm*(2.7e+1/5.12e+2)-9.0/5.12e+2)+std::pow(tXm,4)*(tt4*7.91015625e-2+tt5*2.373046875e-1-tYm*2.373046875e-1-7.91015625e-2)+tXm*(tt4*(9.0/2.56e+2)+tt5*(2.7e+1/2.56e+2)-tYm*(2.7e+1/2.56e+2)-9.0/2.56e+2) ) ;
                    tMomentFittingRHS( 10 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(3.0/2.56e+2)-3.0/2.56e+2)+std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)-tt3*(2.43e+2/5.12e+2)+tXm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)-std::pow(tYm,4)*(tt2*2.373046875e-1-tt3*2.373046875e-1+tXm*2.63671875e-2-2.63671875e-2)+tYm*(tt2*(8.1e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(9.0/2.56e+2)-9.0/2.56e+2) );

                    tMomentFittingRHS( 11 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( -std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tXm,2)*(tt4*(9.0/5.12e+2)-tt5*(2.7e+1/5.12e+2)+tYm*(2.7e+1/5.12e+2)-9.0/5.12e+2)+std::pow(tXm,4)*(tt4*7.91015625e-2-tt5*2.373046875e-1+tYm*2.373046875e-1-7.91015625e-2)+tXm*(tt4*(9.0/2.56e+2)-tt5*(2.7e+1/2.56e+2)+tYm*(2.7e+1/2.56e+2)-9.0/2.56e+2) ) ;
                    tMomentFittingRHS( 11 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( -std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(2.7e+1/2.56e+2)+tXm*(3.0/2.56e+2)-3.0/2.56e+2)-std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)-tt3*(2.43e+2/5.12e+2)+tXm*(2.7e+1/5.12e+2)-2.7e+1/5.12e+2)+std::pow(tYm,4)*(tt2*2.373046875e-1-tt3*2.373046875e-1+tXm*2.63671875e-2-2.63671875e-2)+tYm*(tt2*(8.1e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(9.0/2.56e+2)-9.0/2.56e+2) ) ;

                    tMomentFittingRHS( 12 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)-tt5*(7.29e+2/5.12e+2)+tYm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)-std::pow(tXm,4)*(tt4*2.373046875e-1-tt5*7.119140625e-1+tYm*7.119140625e-1-2.373046875e-1)-tXm*(tt4*(8.1e+1/2.56e+2)-tt5*(2.43e+2/2.56e+2)+tYm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2)) ;
                    tMomentFittingRHS( 12 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)-tt3*(7.29e+2/5.12e+2)+tXm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)-std::pow(tYm,4)*(tt2*2.373046875e-1-tt3*7.119140625e-1+tXm*7.119140625e-1-2.373046875e-1)-tYm*(tt2*(8.1e+1/2.56e+2)-tt3*(2.43e+2/2.56e+2)+tXm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) );

                    tMomentFittingRHS( 13 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)-tt5*(8.1e+1/2.56e+2)+tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)-tt5*(7.29e+2/5.12e+2)+tYm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)+std::pow(tXm,4)*(tt4*2.373046875e-1-tt5*7.119140625e-1+tYm*7.119140625e-1-2.373046875e-1)-tXm*(tt4*(8.1e+1/2.56e+2)-tt5*(2.43e+2/2.56e+2)+tYm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) ) ;
                    tMomentFittingRHS( 13 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)+tt3*(7.29e+2/5.12e+2)-tXm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)-std::pow(tYm,4)*(tt2*2.373046875e-1+tt3*7.119140625e-1-tXm*7.119140625e-1-2.373046875e-1)-tYm*(tt2*(8.1e+1/2.56e+2)+tt3*(2.43e+2/2.56e+2)-tXm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) ) ;

                    tMomentFittingRHS( 14 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)+tt5*(7.29e+2/5.12e+2)-tYm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)+std::pow(tXm,4)*(tt4*2.373046875e-1+tt5*7.119140625e-1-tYm*7.119140625e-1-2.373046875e-1)-tXm*(tt4*(8.1e+1/2.56e+2)+tt5*(2.43e+2/2.56e+2)-tYm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) );
                    tMomentFittingRHS( 14 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)+tt3*(8.1e+1/2.56e+2)-tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)+tt3*(7.29e+2/5.12e+2)-tXm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)+std::pow(tYm,4)*(tt2*2.373046875e-1+tt3*7.119140625e-1-tXm*7.119140625e-1-2.373046875e-1)-tYm*(tt2*(8.1e+1/2.56e+2)+tt3*(2.43e+2/2.56e+2)-tXm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) );

                    tMomentFittingRHS( 15 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNx *( std::pow(tXm,3)*(tt4*(2.7e+1/2.56e+2)+tt5*(8.1e+1/2.56e+2)-tYm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)+std::pow(tXm,2)*(tt4*(2.43e+2/5.12e+2)+tt5*(7.29e+2/5.12e+2)-tYm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)-std::pow(tXm,4)*(tt4*2.373046875e-1+tt5*7.119140625e-1-tYm*7.119140625e-1-2.373046875e-1)-tXm*(tt4*(8.1e+1/2.56e+2)+tt5*(2.43e+2/2.56e+2)-tYm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) ) ;
                    tMomentFittingRHS( 15 , 0 ) += 0.5 * 0.5 * tQuadWeight * tD * tNy *( std::pow(tYm,3)*(tt2*(2.7e+1/2.56e+2)-tt3*(8.1e+1/2.56e+2)+tXm*(8.1e+1/2.56e+2)-2.7e+1/2.56e+2)-std::pow(tYm,2)*(tt2*(2.43e+2/5.12e+2)-tt3*(7.29e+2/5.12e+2)+tXm*(7.29e+2/5.12e+2)-2.43e+2/5.12e+2)+std::pow(tYm,4)*(tt2*2.373046875e-1-tt3*7.119140625e-1+tXm*7.119140625e-1-2.373046875e-1)-tYm*(tt2*(8.1e+1/2.56e+2)-tt3*(2.43e+2/2.56e+2)+tXm*(2.43e+2/2.56e+2)-8.1e+1/2.56e+2) );

                }
                else if ( aOrder == 1 && aDim == 3 )
                {
                    tMomentFittingRHS( 0 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm - 1),2)*(tYm - 1.0)*(tZm - 1.0)) / 16.0 + tNy * -((tXm - 1.0)*std::pow((tYm - 1.0),2)*(tZm - 1.0))/16.0 + tNz * -((tXm - 1.0)*(tYm - 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 1 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm + 1.0),2)*(tYm - 1.0)*(tZm - 1.0))/16.0 + tNy * ((tXm + 1.0)*std::pow((tYm - 1.0),2)*(tZm - 1.0))/16.0 + tNz * ((tXm + 1.0)*(tYm - 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 2 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm + 1.0),2)*(tYm + 1.0)*(tZm - 1.0))/16.0 + tNy * -((tXm + 1.0)*std::pow((tYm + 1.0),2)*(tZm - 1.0))/16.0 + tNz * -((tXm + 1.0)*(tYm + 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 3 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm - 1.0),2)*(tYm + 1.0)*(tZm - 1.0))/16.0 + tNy * ((tXm - 1.0)*std::pow((tYm + 1.0),2)*(tZm - 1.0))/16.0 + tNz * ((tXm - 1.0)*(tYm + 1.0)*std::pow((tZm - 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 4 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm - 1.0),2)*(tYm - 1.0)*(tZm + 1.0))/16.0 + tNy * ((tXm - 1.0)*std::pow((tYm - 1.0),2)*(tZm + 1.0))/16.0 + tNz * ((tXm - 1.0)*(tYm - 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 5 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * -(std::pow((tXm + 1.0),2)*(tYm - 1.0)*(tZm + 1.0))/16.0 + tNy * -((tXm + 1.0)*std::pow((tYm - 1.0),2)*(tZm + 1.0))/16.0 + tNz * -((tXm + 1.0)*(tYm - 1.0)*std::pow((tZm + 1.0),2))/16.0 );

                    tMomentFittingRHS( 6 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx * (std::pow((tXm + 1.0),2)*(tYm + 1.0)*(tZm + 1.0))/16.0 + tNy * ((tXm + 1.0)*std::pow((tYm + 1.0),2)*(tZm + 1.0))/16.0 + tNz * ((tXm + 1.0)*(tYm + 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;

                    tMomentFittingRHS( 7 , 0 ) += 0.5 * ( 1.0 / 3.0 ) * tQuadWeight * tD * ( tNx* -(std::pow((tXm - 1),2)*(tYm + 1.0)*(tZm + 1.0))/16.0 + tNy * -((tXm - 1.0)*std::pow((tYm + 1.0),2)*(tZm + 1.0))/16.0 + tNz * -((tXm - 1.0)*(tYm + 1.0)*std::pow((tZm + 1.0),2))/16.0 ) ;


                    


                }
                
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

            //tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5 * ( (std::pow( (tx0 - tx1) , 2)*( ty0 - 1.0 ))/24.0 - (std::pow( (tx0 - tx1) , 2)*( ty0 - ty1 ))/32.0 + ((( 1 / 8.0 ) * std::pow( tx0 , 2 ) -  tx0 / 4.0 )*( ty0 - ty1 ))/2.0 )*/
            
            //tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5*tNx*(( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + 0.5 * (( 1.0 / 8.0 )* (- std::pow( tx0 , 2 )) + ( 1.0 / 4.0 )*( tx0 ) )*( ty0 - ty1 ) + (1.0 / 3.0)*( ty0 - ty1 )*( tx1/4.0 - tx0/4.0 + (tx0/4.0) * (tx0 - tx1) ) - ( -std::pow( tx0 , 2 )/8.0 + tx0/4.0 )*(ty0 - 1) - 0.5*(ty0 - 1.0)*(tx1/4.0 - tx0/4.0 + (tx0/4.0)*(tx0 - tx1) ));  

            //tMomentFittingRHS( 0 , 0 ) = tMomentFittingRHS( 0 , 0 ) + 0.5*tNy*( 0.25*0.5*(tx0 - tx1)*(0.5 * std::pow( ty0 , 2 ) + ty0) - (1.0/32.0)*(tx0-tx1)*(std::pow( (ty0-ty1), 2)) + (1.0/12.0)*(tx0 - tx1)*(ty1 - ty0 + ty0*(ty0 - ty1)) + (1.0/24.0)*(tx0 - 1.0)*(std::pow( (ty0-ty1), 2)) - ( -0.5*std::pow( ty0, 2 ) + ty0 )*0.25*( tx0 - 1.0 ) - 0.5*0.25*( tx0 - 1 )*( ty1- ty0 + ty0*(ty0 - ty1) ) );   

            //tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) + 0.5*tNx*( ( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + (( 1.0/8.0 )*std::pow( tx0 , 2) - tx0/4.0 )*0.5*( ty0 -ty1 ) - (1.0/3.0)*(ty0 - ty1)*( tx0/4.0 - tx1/4.0 + (tx0/4.0)*( tx0 - tx1 )) - (std::pow(tx0,2)/8.0 + tx0/4.0)*(ty0 - 1.0) + 0.5*(ty0 - 1.0)*(tx0/4.0 - tx1/4.0 + (tx0/4.0)*(tx0 - tx1) ) );

            //tMomentFittingRHS( 1 , 0 ) = tMomentFittingRHS( 1 , 0 ) + 0.5*tNy*((1.0/32.0)*(tx0 - tx1)*(std::pow( (ty0 - ty1), 2)) - (1.0/8.0)*(tx0 - tx1)*( 0.5*std::pow( ty0 , 2 ) + ty0 ) - (1.0/12.0)*(tx0 - tx1)*(ty1 - ty0 + ty0*(ty0 - ty1)) -  (1.0/24.0)*(std::pow(ty0 - ty1, 2))*( tx0 + 1.0) + 0.25*(-std::pow(ty0,2) + ty0 ) + 0.25*0.5*(tx0 + 1)*(ty1 - ty0 + ty0*(ty0 - ty1)));

            //tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + 0.5*tNx*( ( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 + 1.0 ) - ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) + 0.5*((1.0/8.0)*std::pow(tx0,2)+ tx0/4.0 )*(ty0 - ty1) + (1.0/3.0)*(ty0 -ty1)*0.25*(tx0 - tx1 + tx0*(tx0 - tx1)) + 0.25*(0.5*std::pow(tx0,2) + tx0)*(ty0 + 1.0) - 0.5*0.25*(ty0 + 1.0)*(tx0 - tx1 + tx0*(tx0 - tx1)));

            //tMomentFittingRHS( 2 , 0 ) = tMomentFittingRHS( 2 , 0 ) + 0.5*tNy*( (1.0/12.0)*(tx0 - tx1)*(ty0 - ty1 + ty0*(ty0 - ty1)) - 0.5*0.25*(tx0 - tx1)*(0.5*std::pow(ty0,2)+ty0) - (1.0/32.0)*(tx0 - tx1)*(std::pow(ty0 - ty1,2)) + (1.0/24.0)*(tx0 + 1.0)*(std::pow(ty0 - ty1, 2)) + 0.25*(0.5*std::pow(ty0,2) + ty0)*(tx0 + 1.0) - 0.5*0.25*(tx0 + 1.0)*(ty0 - ty1 +ty0*(ty0 - ty1)));

            //tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) + 0.5*tNx*( -( 1.0 / 24.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 + 1.0 ) + ( 1.0 / 32.0 ) * std::pow(( tx0 - tx1 ), 2 )* ( ty0 - ty1 ) - 0.25*(-0.5*std::pow(tx0 , 2) + tx0 )*( ty0 - ty1 ) - (1.0/12.0)*(ty0 - 1)*( tx1 - tx0 + tx0*(tx0 - tx1)) + 0.25*(ty0 + 1.0)*( -0.5*std::pow(tx0 , 2) + tx0 ) + (1.0/8.0)*(ty0 + 1)*( tx1 - tx0 + tx0*(tx0 - tx1)) );

            //tMomentFittingRHS( 3 , 0 ) = tMomentFittingRHS( 3 , 0 ) + 0.5*tNy*( (1.0/32.0)*(tx0 -tx1)*(std::pow(ty0 - ty1, 2)) + 0.25*0.5*(tx0 - tx1)*( 0.5*std::pow(ty0, 2) + ty0 ) - (1.0/12.0)*(tx0 - tx1)*(ty0 - ty1 + ty0*(ty0 - ty1)) - (1.0/24.0)*(tx0 - 1.0)*(std::pow(ty0 - ty1, 2)) + 0.25*(tx0 - 1)*(0.25*std::pow(ty0,2)) + 0.5*0.25*(tx0 - 1)*( ty0 - ty1 + ty0*(ty0 - ty1)));

            //fprintf(stdout, "Facet Index %d \n", iFacetIndex );

        }

        // Solve the system
        //mQuadratureWeights = {{0.0}, {0.0}, {0.0}, {0.0}};
        
        if ( ( aOrder == 1 || aOrder == 2 )  && aDim == 2 )
        {
            mQuadratureWeights = ( inv( tMomentFittingLHS ) * tMomentFittingRHS);
        }
        else if ( aOrder == 1 && aDim == 3 )
        {
            mQuadratureWeights = ( inv( tMomentFittingLHS ) * tMomentFittingRHS);
        }
        else if ( aOrder == 2 && aDim == 3 )
        {
            mQuadratureWeights =  solve( tMomentFittingLHS  , tMomentFittingRHS);
        }
        else if ( aOrder == 3 )
        {
            mQuadratureWeights =  solve( tMomentFittingLHS  , tMomentFittingRHS);
        }
        
        // Reshape so as to be compatible with format required downstream
        mQuadratureWeights.reshape( 1 , tNmoments );
        
        
    }
    

    //----------------------------------------------------------------

    Matrix< DDRMat >
    Cell_Cluster::get_quadrature_weights() const
    {
        return mQuadratureWeights;
    }


    //----------------------------------------------------------------

    bool
    Cell_Cluster::has_cluster_group( const moris_index aDiscretizationMeshIndex ) const
    {
        // if the list of B-spline meshes for which the cluster groups have been set isn't large enough ...
        if ( mClusterGroups.size() < (uint)aDiscretizationMeshIndex + 1 )
        {
            // ... the cluster group has not been set yet
            return false;
        }

        // if the entry exists but is empty
        else if ( mClusterGroups( aDiscretizationMeshIndex ).expired() )
        {
            return false;
        }

        // if the entry both exists and is not empty, the cluster has a group
        else
        {
            return true;
        }
    }

    //----------------------------------------------------------------

    std::shared_ptr< mtk::Cluster_Group >
    Cell_Cluster::get_cluster_group( const moris_index aDiscretizationMeshIndex ) const
    {
        // check that the cluster group exists and is set
        MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
                "xtk::Cell_Cluster::get_cluster_group() - Cluster group is not set or does not exist." );

        // return the pointer to the cluster group
        return mClusterGroups( aDiscretizationMeshIndex ).lock();
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster::set_cluster_group(
            const moris_index                            aDiscretizationMeshIndex,
            const std::shared_ptr< mtk::Cluster_Group > &aClusterGroupPtr )
    {
        // check that the cluster group is set to the correct B-spline list index
        MORIS_ASSERT( aClusterGroupPtr->get_discretization_mesh_index_for_cluster_group() == aDiscretizationMeshIndex,
                "xtk::Cell_Cluster::set_cluster_group() - Index which the cluster group lives on is not the list index it gets set to on the cluster." );

        // check if the list of cluster groups is big enough to accommodate the cluster groups for each B-spline mesh
        if ( mClusterGroups.size() < (uint)aDiscretizationMeshIndex + 1 )
        {
            // ... if not increase the size
            mClusterGroups.resize( aDiscretizationMeshIndex + 1 );
        }

        // store pointer to the cluster group associated with this B-spline mesh
        mClusterGroups( aDiscretizationMeshIndex ) = aClusterGroupPtr;
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< IG_Vertex_Group >
    Cell_Cluster::get_ig_vertex_group()
    {
        return mVertexGroup;
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell_Cluster::compute_cluster_group_cell_measure(
            const moris_index          aDiscretizationMeshIndex,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // only compute // FIXME: ghost clusters for visualization should also have their cluster volumes assigned
        if ( !mOnlyForVis )
        {
            // check that the cluster group exists and is set
            MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
                    "xtk::Cell_Cluster::compute_cluster_group_cell_measure() - Cluster group is not set or does not exist." );

            // compute the group measure and return it
            return mClusterGroups( aDiscretizationMeshIndex ).lock()->compute_cluster_group_cell_measure( aPrimaryOrVoid, aIsLeader );
        }
        else    // cluster groups are not defined on clusters that are only for visualization purposes (e.g. for ghost visualization)
        {
            return -1.0;
        }
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell_Cluster::compute_cluster_group_cell_measure_derivative(
            const moris_index          aDiscretizationMeshIndex,
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the cluster group exists and is set
        MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
                "xtk::Cell_Cluster::compute_cluster_group_cell_measure_derivative() - Cluster group is not set or does not exist." );

        // compute the group measure derivative and return it
        return mClusterGroups( aDiscretizationMeshIndex ).lock()->compute_cluster_group_cell_measure_derivative( aPerturbedVertexCoords, aDirection, aPrimaryOrVoid, aIsLeader );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::xtk
// } // namespace moris
