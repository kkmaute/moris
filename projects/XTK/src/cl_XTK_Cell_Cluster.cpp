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
    
        MORIS_ASSERT( aDim < 3, "Currently moment fitting only works for 2D problems" );

        MORIS_ASSERT( aOrder < 2, "Currently moment fitting only works for linear problems" );

        mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
        const mtk::Integrator tIntData( tIntObj );

        // Get weights from mtk::integrator

        Matrix< DDRMat >  tMTKPoints;
        
        tIntData.get_points( tMTKPoints );
            
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

        mQuadraturePoints = tMTKPoints;
        
    

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

        MORIS_ASSERT(aDim < 3, "Currently moment fitting only works for 2D problems");

        mtk::Integration_Rule tIntObj( mtk::Geometry_Type::QUAD , mtk::Integration_Type::GAUSS , mtk::Integration_Order::QUAD_2x2 , mtk::Geometry_Type::LINE , mtk::Integration_Type::GAUSS , mtk::Integration_Order::BAR_1  );
        
        mtk::Integrator tIntData( tIntObj );

        // Get weights from mtk::integrator

        Matrix< DDRMat >  tMTKWeights;
        
        tIntData.get_weights( tMTKWeights );
        
        //Vector< real > mOneDWeights;

        // if (mOrder == 1)
        // {
        //     mOneDWeights = { 1.0 , 1.0 };
        // }
        // else
        // {
        //     mOneDWeights =  { 8.0/9.0, 5.0/9.0 , 5.0/9.0 };
        // }
        
        MORIS_ASSERT( aOrder < 2, "Only 1st order supported currently");
        
        /*mQuadratureWeights.resize( mOneDWeights.size()*mOneDWeights.size() , 2 );

        for (uint iRowIndex = 0; iRowIndex < mOneDWeights.size() ; iRowIndex ++ )
        {
            for (uint iColIndex = 0; iColIndex < mOneDWeights.size() ; iColIndex ++ )
            {
                uint iMatIndex = 2*(iRowIndex) + iColIndex;
                
                mQuadratureWeights( iMatIndex ) = mOneDWeights( iRowIndex ) * mOneDWeights( iColIndex );

            }
        }*/
        mQuadratureWeights = tMTKWeights;

        

    }

    //------------------------------------------------------------------------------

    void 
    Cell_Cluster::find_subphase_boundary_vertices(
                const std::shared_ptr< IG_Cell_Group >  aSubphaseIGCells,
                const std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity
                )
    {
        // Get the primary subphase IG cells first
        Vector< mtk::Cell* > tSubphaseIgCellsPtr = aSubphaseIGCells->mIgCellGroup;

        for (uint iCellIndex = 0; iCellIndex < tSubphaseIgCellsPtr.size(); iCellIndex++)
            {
                // get cell
                mtk::Cell* tCellObj = tSubphaseIgCellsPtr( iCellIndex );

                fprintf( stdout, " Cell ID in subphase %d\n", (moris_id)tCellObj->get_id() );
            }
        
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
            Vector< moris::mtk::Cell* > tOwningElementsInSubphase;

            for (uint iSubphaseCellIndex = 0; iSubphaseCellIndex < tSubphaseIgCellsPtr.size(); iSubphaseCellIndex++)
            {
                // Get cell
                mtk::Cell* tSubphaseCellToCheck = tSubphaseIgCellsPtr( iSubphaseCellIndex );

                for (uint iOwningCells = 0; iOwningCells < tOwningElements.size() ; iOwningCells++) 
                {
                    // Check if subphase cell part of facet owning cells
                    if (tOwningElements( iOwningCells ) == tSubphaseCellToCheck)
                    {
                        tOwningElementsInSubphase.push_back( tOwningElements( iOwningCells ) );
                    } 


                }

            } 

            // If belonging to one cell then add it into the vector containing pointer to facet indices.
            if ( tOwningElementsInSubphase.size()  == 1 )
            {
                mFacetVerticesOnSubphaseBoundary.push_back( aFacetConnectivity->mFacetVertices( tSingleFacet ) );
            }

        }  

        // The code below was all for testing the facet filtering function. Not needed otherwise, but still kept.

        std::vector< double > tFacetCoordinatesFileVectorX;

        std::vector< double > tFacetCoordinatesFileVectorY;


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
        

    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::compute_quadrature_weights( const uint aOrder, const uint aDim )
    {

        // Get midpoint of facet coordinates -

        Matrix< DDRMat > tFacetMidpoints;

        tFacetMidpoints.reshape( mFacetVerticesOnSubphaseBoundary.size() , 2 );

        for (uint iFacetIndex = 0; iFacetIndex < mFacetVerticesOnSubphaseBoundary.size(); iFacetIndex++)
        {
            // Get coordinates corresponding to one facet -
            Vector< mtk::Vertex* > tFacetCoords = mFacetVerticesOnSubphaseBoundary( iFacetIndex );

            // Get coordinates and sum them to get midpoint
            mtk::Vertex* tFirstFacetCoord = tFacetCoords( 0 );

            mtk::Vertex* tSecondFacetCoord = tFacetCoords( 1 );

            Matrix< DDRMat > tFacetMidpointCoords = 0.5*(tFirstFacetCoord->get_coords() + tSecondFacetCoord->get_coords()) ;

            tFacetMidpoints.set_row( iFacetIndex, tFacetMidpointCoords ) ;
            
        }

        // Now, create interpolation objects for specifying the basis function for each facet. One for the IP and one for the geometry.
        mtk::Interpolation_Rule tInterpolationRule( mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

        mtk::Interpolation_Rule tIPInterpolationRule( mtk::Geometry_Type::QUAD , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

        // Declare cellshape
        mtk::CellShape tCellShape = mtk::CellShape::GENERAL;

        // Set space sideset to true for interpolation to start
        bool tSpaceSideset = true;

        // Create space interpolator object
        mtk::Space_Interpolator tSpaceInterpolationObject( tInterpolationRule , tIPInterpolationRule , tCellShape , tSpaceSideset );

        Matrix< DDRMat > tFacetNormals;

        tFacetNormals.reshape( mFacetVerticesOnSubphaseBoundary.size() , 2 );
        
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

        close_hdf5_file( tFileID2 );

        //Declare LHS matrix
        Matrix < DDRMat > tMomentFittingLHS;
        tMomentFittingLHS.reshape( 4 , 4 );
        
        mtk::Interpolation_Function_Base* tIPInterp = tIPInterpolationRule.create_space_interpolation_function();

        // Generate LHS
        for (uint iQuadPointIndex = 0; iQuadPointIndex < mQuadraturePoints.n_rows() ; iQuadPointIndex++)
        {
            // Declare matrix for basis function values
            Matrix< DDRMat > tN;

            // Get quad point
            Matrix< DDRMat > tXi = mQuadraturePoints.get_column( iQuadPointIndex );

            // Get value of basis functions at quad point
            tIPInterp->eval_N( tXi , tN );

            // Place it in LHS 
            tMomentFittingLHS.set_column( iQuadPointIndex , trans( tN ) );

        }
        
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
