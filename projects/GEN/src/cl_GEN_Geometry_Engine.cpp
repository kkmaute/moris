/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Engine.cpp
 *
 */

// MRS
#include "cl_Param_List.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Tracer.hpp"
#include "cl_Library_IO.hpp"

// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"

// MTK
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

// SOL TODO if we move this out of SOL, GEN doesn't have to depend on it
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------
    // PUBLIC
    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine::Geometry_Engine(
            Vector< Vector< ParameterList > >    aParameterLists,
            const std::shared_ptr< Library_IO >& aLibrary,
            mtk::Mesh*                           aMesh )
            : mNodeManager( aMesh )
            , mPhaseTable( create_phase_table( aParameterLists, aLibrary ) )
            , mPDVHostManager( mNodeManager )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create geometry engine" );

        // Requested IQIs
        mRequestedIQIs = string_to_cell< std::string >( aParameterLists( 0 )( 0 ).get< std::string >( "IQI_types" ) );

        // Geometries
        mGeometryFieldFile = aParameterLists( 0 )( 0 ).get< std::string >( "geometry_field_file" );
        mOutputMeshFile    = aParameterLists( 0 )( 0 ).get< std::string >( "output_mesh_file" );
        mTimeOffset        = aParameterLists( 0 )( 0 ).get< real >( "time_offset" );

        // Read ADVs
        if ( aParameterLists( 0 )( 0 ).get< sint >( "advs_size" ) )
        {
            mInitialPrimitiveADVs = Matrix< DDRMat >( aParameterLists( 0 )( 0 ).get< sint >( "advs_size" ), 1, aParameterLists( 0 )( 0 ).get< real >( "initial_advs_fill" ) );
            mLowerBounds          = Matrix< DDRMat >( aParameterLists( 0 )( 0 ).get< sint >( "advs_size" ), 1, aParameterLists( 0 )( 0 ).get< real >( "lower_bounds_fill" ) );
            mUpperBounds          = Matrix< DDRMat >( aParameterLists( 0 )( 0 ).get< sint >( "advs_size" ), 1, aParameterLists( 0 )( 0 ).get< real >( "upper_bounds_fill" ) );
        }
        else
        {
            mInitialPrimitiveADVs = string_to_mat< DDRMat >( aParameterLists( 0 )( 0 ).get< std::string >( "initial_advs" ) );
            mLowerBounds          = string_to_mat< DDRMat >( aParameterLists( 0 )( 0 ).get< std::string >( "lower_bounds" ) );
            mUpperBounds          = string_to_mat< DDRMat >( aParameterLists( 0 )( 0 ).get< std::string >( "upper_bounds" ) );

            // check that advs and bounds are vectors
            MORIS_ERROR( isvector( mInitialPrimitiveADVs ), "ADVs need to be of type vector.\n" );
            MORIS_ERROR( isvector( mLowerBounds ), "ADV lower bounds need to be of type vector.\n" );
            MORIS_ERROR( isvector( mUpperBounds ), "ADV upper bounds need to be of type vector.\n" );

            // ensure that advs and bounds are column vectors
            mInitialPrimitiveADVs = mInitialPrimitiveADVs.n_rows() == 1 ? trans( mInitialPrimitiveADVs ) : mInitialPrimitiveADVs;
            mLowerBounds          = mLowerBounds.n_rows() == 1 ? trans( mLowerBounds ) : mLowerBounds;
            mUpperBounds          = mUpperBounds.n_rows() == 1 ? trans( mUpperBounds ) : mUpperBounds;
        }

        // Create designs with the factory
        for ( uint iParameterIndex = 2; iParameterIndex < aParameterLists.size(); iParameterIndex++ )
        {
            aParameterLists( 1 ).append( aParameterLists( iParameterIndex ) );
        }
        Design_Factory tDesignFactory( aParameterLists( 1 ), mInitialPrimitiveADVs, aLibrary, aMesh, mNodeManager );

        // Get geometries and properties from the factory
        mGeometries = tDesignFactory.get_geometries();
        mProperties = tDesignFactory.get_properties();

        MORIS_ERROR( mGeometries.size() <= MAX_GEOMETRIES,
                "Number of geometries exceeds MAX_GEOMETRIES, please change this in GEN_Data_Types.hpp" );

        // Set requested PDVs
        Vector< std::string > tRequestedPDVNames = string_to_cell< std::string >( aParameterLists( 0 )( 0 ).get< std::string >( "PDV_types" ) );
        Vector< PDV_Type >    tRequestedPDVTypes( tRequestedPDVNames.size() );

        map< std::string, PDV_Type > tPDVTypeMap = get_pdv_type_map();

        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tRequestedPDVTypes.size(); tPDVTypeIndex++ )
        {
            tRequestedPDVTypes( tPDVTypeIndex ) = tPDVTypeMap[ tRequestedPDVNames( tPDVTypeIndex ) ];
        }
        mPDVHostManager.set_requested_interpolation_pdv_types( tRequestedPDVTypes );

        // Initialize PDV type list
        this->initialize_pdv_type_list();

        // Print the phase table if requested
        if ( aParameterLists( 0 )( 0 ).get< bool >( "print_phase_table" ) and par_rank() == 0 )
        {
            mPhaseTable.print();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine::Geometry_Engine(
            mtk::Interpolation_Mesh*   aMesh,
            Geometry_Engine_Parameters aParameters )
            : mNodeManager( aMesh )
            , mGeometries( aParameters.mGeometries )
            , mProperties( aParameters.mProperties )
            , mPhaseTable( create_phase_table( aParameters.mGeometries.size(), aParameters.mBulkPhases ) )
            , mInitialPrimitiveADVs( aParameters.mADVs )
            , mTimeOffset( aParameters.mTimeOffset )
            , mPDVHostManager( mNodeManager )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create geometry engine" );

        // Create integration mesh
        mtk::Integration_Mesh* tIntegrationMesh = create_integration_mesh_from_interpolation_mesh(
                aMesh->get_mesh_type(),
                aMesh );

        // Register mesh pair
        mtk::Mesh_Pair tMeshPair( aMesh, tIntegrationMesh );

        // Distribute ADVs
        this->distribute_advs( tMeshPair );
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine::~Geometry_Engine()
    {
        // Delete stored distributed vectors
        delete mOwnedADVs;
        delete mPrimitiveADVs;

        // Delete queued intersection node, in case it wasn't admitted
        delete mQueuedIntersectionNode;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::set_advs( const Matrix< DDRMat >& aNewADVs )
    {
        // Set new ADVs
        mOwnedADVs->vec_put_scalar( 0 );
        mOwnedADVs->replace_global_values( mFullADVIds, aNewADVs );
        mOwnedADVs->vector_global_assembly();
        mPrimitiveADVs->import_local_to_global( *mOwnedADVs );

        // Import ADVs into fields that need it
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            mGeometries( tGeometryIndex )->import_advs( mOwnedADVs );
        }
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            mProperties( tPropertyIndex )->import_advs( mOwnedADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >&
    Geometry_Engine::get_advs()
    {
        // Create full ADVs
        sol::Matrix_Vector_Factory tDistributedFactory;

        sol::Dist_Map*    tFullMap    = tDistributedFactory.create_map( mFullADVIds );
        sol::Dist_Vector* tFullVector = tDistributedFactory.create_vector( tFullMap, 1, false, true );

        // Import ADVs
        tFullVector->import_local_to_global( *mOwnedADVs );

        // Extract copy
        tFullVector->extract_copy( mADVs );

        // Delete full ADVs/map
        delete tFullVector;

        return mADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::set_phase_function(
            PHASE_FUNCTION               aPhaseFunction,
            uint                         aNumPhases,
            const Vector< std::string >& aPhaseNames )
    {
        mPhaseTable.set_phase_function( aPhaseFunction, aNumPhases, aPhaseNames );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::set_dQIdp(
            const Vector< moris::Matrix< DDRMat >* >& adQIdp,
            moris::Matrix< moris::DDSMat >*           aMap )
    {
        mPDVHostManager.set_dQIdp( adQIdp, aMap );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< IdMat >&
    Geometry_Engine::get_IjklIDs()
    {
        return mFullijklIDs;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >&
    Geometry_Engine::get_lower_bounds()
    {
        return mLowerBounds;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >&
    Geometry_Engine::get_upper_bounds()
    {
        return mUpperBounds;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::communicate_requested_IQIs()
    {
        mPDVHostManager.set_requested_IQIs( mRequestedIQIs );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Geometry_Engine::get_dcriteria_dadv()
    {
        return mPDVHostManager.compute_diqi_dadv( mFullADVIds );
    }

    //--------------------------------------------------------------------------------------------------------------

    MSI::Design_Variable_Interface*
    Geometry_Engine::get_design_variable_interface()
    {
        return &mPDVHostManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::is_intersected_by_active_geometry( const Matrix< IndexMat >& aNodeIndices )
    {
        // Get first geometric region
        Geometric_Region tFirstNodeGeometricRegion = mGeometries( mActiveGeometryIndex )->get_geometric_region( aNodeIndices( 0 ), mNodeManager.get_node( aNodeIndices( 0 ) ).get_global_coordinates() );

        // Test nodes for other geometric regions
        for ( uint iNodeInEntityIndex = 0; iNodeInEntityIndex < aNodeIndices.length(); iNodeInEntityIndex++ )
        {
            // Get test geometric region
            Geometric_Region tTestGeometricRegion = mGeometries( mActiveGeometryIndex )->get_geometric_region( aNodeIndices( iNodeInEntityIndex ), mNodeManager.get_node( aNodeIndices( iNodeInEntityIndex ) ).get_global_coordinates() );

            // Test if it is different from the first region. If so, the entity is intersected
            if ( tTestGeometricRegion != tFirstNodeGeometricRegion )
            {
                return true;
            }
        }

        // If no differences were found, this entity is not intersected
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::is_intersected(
            uint                      aGeometryIndex,
            const Matrix< IndexMat >& aNodeIndices )
    {
        // Get first geometric region
        Geometric_Region tFirstNodeGeometricRegion = mGeometries( aGeometryIndex )->get_geometric_region( aNodeIndices( 0 ), mNodeManager.get_node( aNodeIndices( 0 ) ).get_global_coordinates() );

        // Test nodes for other geometric regions
        for ( uint iNodeInEntityIndex = 0; iNodeInEntityIndex < aNodeIndices.length(); iNodeInEntityIndex++ )
        {
            // Get test geometric region
            Geometric_Region tTestGeometricRegion = mGeometries( aGeometryIndex )->get_geometric_region( aNodeIndices( iNodeInEntityIndex ), mNodeManager.get_node( aNodeIndices( iNodeInEntityIndex ) ).get_global_coordinates() );

            // Test if it is different from the first region. If so, the entity is intersected
            if ( tTestGeometricRegion != tFirstNodeGeometricRegion )
            {
                return true;
            }
        }

        // If no differences were found, this entity is not intersected
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queue_intersection(
            uint                     aEdgeFirstNodeIndex,
            uint                     aEdgeSecondNodeIndex,
            const Matrix< DDRMat >&  aEdgeFirstNodeParametricCoordinates,
            const Matrix< DDRMat >&  aEdgeSecondNodeParametricCoordinates,
            const Matrix< DDUMat >&  aBackgroundElementNodeIndices,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        // If previous intersection node was not admitted, this will delete it
        delete mQueuedIntersectionNode;

        // Get background nodes
        Vector< Background_Node* > tBackgroundNodes( aBackgroundElementNodeIndices.length() );
        for ( uint iNode = 0; iNode < tBackgroundNodes.size(); iNode++ )
        {
            tBackgroundNodes( iNode ) = &( mNodeManager.get_background_node( aBackgroundElementNodeIndices( iNode ) ) );
        }

        // Create parent nodes
        Parent_Node tFirstParentNode( mNodeManager.get_node( aEdgeFirstNodeIndex ), aEdgeFirstNodeParametricCoordinates );
        Parent_Node tSecondParentNode( mNodeManager.get_node( aEdgeSecondNodeIndex ), aEdgeSecondNodeParametricCoordinates );

        // Have the active geometry create a new intersection node
        mQueuedIntersectionNode = mGeometries( mActiveGeometryIndex )->create_intersection_node( mNodeManager.get_total_number_of_nodes(), tBackgroundNodes, tFirstParentNode, tSecondParentNode, aBackgroundGeometryType, aBackgroundInterpolationOrder );

        // Return if queued intersected node is on the parent edge
        return mQueuedIntersectionNode->parent_edge_is_intersected();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queued_intersection_first_parent_on_interface()
    {
        return mQueuedIntersectionNode->is_first_parent_on_interface();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queued_intersection_second_parent_on_interface()
    {
        return mQueuedIntersectionNode->is_second_parent_on_interface();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Geometry_Engine::get_queued_intersection_local_coordinate()
    {
        return mQueuedIntersectionNode->get_local_coordinate();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Geometry_Engine::get_queued_intersection_global_coordinates()
    {
        return mQueuedIntersectionNode->get_global_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Geometry_Engine::admit_queued_intersection()
    {
        // Add new derived node to the node manager
        mNodeManager.add_derived_node( mQueuedIntersectionNode );

        // Set to nullptr to let the geometry engine know it was admitted
        mQueuedIntersectionNode = nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::update_intersection_node(
            uint        aNodeIndex,
            moris_id    aNodeId,
            moris_index aNodeOwner )
    {
        mNodeManager.update_derived_node( aNodeIndex, aNodeId, aNodeOwner );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::create_new_derived_nodes(
            Vector< mtk::Cell* >&                                aNewNodeParentCell,
            const Vector< std::shared_ptr< Matrix< DDRMat > > >& aParametricCoordinates )
    {
        // This function can't be traced; Right now XTK does not always call it from all processors.

        // Get number of new nodes
        uint tNumberOfNewDerivedNodes = aNewNodeParentCell.size();

        // Get vertex indices from parent cell and change parametric coordinate type
        Vector< Matrix< IndexMat > > tVertexIndices( tNumberOfNewDerivedNodes );
        Vector< Matrix< DDRMat > >   tParametricCoordinates( tNumberOfNewDerivedNodes );
        for ( uint iCellIndex = 0; iCellIndex < tNumberOfNewDerivedNodes; iCellIndex++ )
        {
            tVertexIndices( iCellIndex )         = aNewNodeParentCell( iCellIndex )->get_vertex_inds();
            tParametricCoordinates( iCellIndex ) = *aParametricCoordinates( iCellIndex );
        }

        // Call overloaded function
        if ( tNumberOfNewDerivedNodes > 0 )
        {
            this->create_new_derived_nodes(
                    tVertexIndices,
                    tParametricCoordinates,
                    aNewNodeParentCell( 0 )->get_geometry_type(),
                    aNewNodeParentCell( 0 )->get_interpolation_order() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::create_new_derived_nodes(
            const Vector< Matrix< IndexMat > >& aVertexIndices,
            const Vector< Matrix< DDRMat > >&   aParametricCoordinates,
            mtk::Geometry_Type                  aBackgroundGeometryType,
            mtk::Interpolation_Order            aBackgroundInterpolationOrder )
    {
        // This function can't be traced; Right now XTK does not always call it from all processors.

        // Loop over nodes
        for ( uint iNode = 0; iNode < aVertexIndices.size(); iNode++ )
        {
            // Create basis nodes
            Vector< Background_Node* > tBackgroundNodes( aVertexIndices( iNode ).length() );
            for ( uint iBaseNode = 0; iBaseNode < aVertexIndices( iNode ).length(); iBaseNode++ )
            {
                tBackgroundNodes( iBaseNode ) = &( mNodeManager.get_background_node( aVertexIndices( iNode )( iBaseNode ) ) );
            }

            // Create new derived node
            mNodeManager.create_derived_node(
                    tBackgroundNodes,
                    aParametricCoordinates( iNode ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_num_phases()
    {
        return mPhaseTable.get_num_phases();
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_phase_index(
            moris_index             aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // Initialize bitset of geometry signs
        Geometry_Bitset tGeometrySigns( 0 );

        // Flip bits as needed
        for ( uint iGeometryIndex = 0; iGeometryIndex < mGeometries.size(); iGeometryIndex++ )
        {
            Geometric_Region tGeometricRegion = mGeometries( iGeometryIndex )->get_geometric_region( aNodeIndex, aNodeCoordinates );
            tGeometrySigns.set( iGeometryIndex, tGeometricRegion == Geometric_Region::POSITIVE );
        }

        return mPhaseTable.get_phase_index( tGeometrySigns );
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Geometry_Engine::get_geometric_region(
            uint                    aGeometryIndex,
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        return mGeometries( aGeometryIndex )->get_geometric_region( aNodeIndex, aNodeCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Geometry_Engine::get_element_phase_index( const mtk::Cell& aCell )
    {
        // Get the vertices that are a part of this cell
        Vector< mtk::Vertex* > tVertices = aCell.get_vertex_pointers();

        // Start geometry bitset
        Geometry_Bitset tGeometrySigns( 0 );

        // Loop over all geometries
        for ( uint iGeometryIndex = 0; iGeometryIndex < mGeometries.size(); iGeometryIndex++ )
        {
            // Loop over vertices on the cell
            for ( auto iVertex : tVertices )
            {
                // Get geometric region
                Geometric_Region tGeometricRegion = mGeometries( iGeometryIndex )->get_geometric_region( iVertex->get_index(), iVertex->get_coords() );

                // If we can determine the region already, do so
                if ( tGeometricRegion == Geometric_Region::NEGATIVE )
                {
                    tGeometrySigns.set( iGeometryIndex, false );
                    goto region_determined;
                }
                else if ( tGeometricRegion == Geometric_Region::POSITIVE )
                {
                    tGeometrySigns.set( iGeometryIndex, true );
                    goto region_determined;
                }
            }

            // All vertices are (somehow) on the interface; return that this is a problem
            return MORIS_INDEX_MAX;

        region_determined:;
        }

        return mPhaseTable.get_phase_index( tGeometrySigns );
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_number_of_geometries()
    {
        return mGeometries.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_num_bulk_phase()
    {
        return mPhaseTable.get_num_phases();
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_active_geometry_index()
    {
        return mActiveGeometryIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::advance_geometry_index()
    {
        MORIS_ASSERT( mActiveGeometryIndex < mGeometries.size(),
                "Trying to advance past the number of geometries in the geometry engine" );
        mActiveGeometryIndex += 1;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Geometry_Engine::get_num_refinement_fields()
    {
        return mGeometries.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< uint >&
    Geometry_Engine::get_num_refinements( uint aFieldIndex )
    {
        return mGeometries( aFieldIndex )->get_num_refinements();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< uint >&
    Geometry_Engine::get_refinement_mesh_indices( uint aFieldIndex )
    {
        return mGeometries( aFieldIndex )->get_refinement_mesh_indices();
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > >
    Geometry_Engine::get_mtk_fields()
    {
        // Initialize vector of mtk fields
        Vector< std::shared_ptr< mtk::Field > > tMTKFields;

        // Loop over geometries
        for ( const auto& iGeometry : mGeometries )
        {
            // Add MTK fields
            tMTKFields.append( iGeometry->get_mtk_fields() );
        }

        // Loop over properties
        for ( const auto& iProperty : mProperties )
        {
            // Add MTK field, if it exists
            auto tMTKField = iProperty->get_mtk_field();
            if ( tMTKField )
            {
                tMTKFields.push_back( tMTKField );
            }
        }

        // Return final list
        return tMTKFields;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Geometry_Engine::get_refinement_function_index(
            uint aFieldIndex,
            uint aRefinementIndex )
    {
        return mGeometries( aFieldIndex )->get_refinement_function_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::create_pdvs( const mtk::Mesh_Pair& aMeshPair )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create PDVs" );

        // Get meshes using first mesh on mesh manager: Lagrange mesh with numbered aura (default)
        mtk::Integration_Mesh*   tIntegrationMesh   = aMeshPair.get_integration_mesh();
        mtk::Interpolation_Mesh* tInterpolationMesh = aMeshPair.get_interpolation_mesh();

        // Initialize PDV type groups and mesh set info from integration mesh
        Vector< Vector< Vector< PDV_Type > > > tPDVTypes( tIntegrationMesh->get_num_sets() );
        Vector< PDV_Type >                     tPDVTypeGroup( 1 );

        // Loop over properties to create PDVs
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            // Get PDV type from property
            tPDVTypeGroup( 0 ) = mProperties( tPropertyIndex )->get_pdv_type();

            // Get mesh set indices and names
            Vector< uint > tMeshSetIndices = mProperties( tPropertyIndex )->get_pdv_mesh_set_indices( tIntegrationMesh );

            // Assign PDV types to each mesh set
            for ( uint tSetIndexPosition = 0; tSetIndexPosition < tMeshSetIndices.size(); tSetIndexPosition++ )
            {
                uint tMeshSetIndex = tMeshSetIndices( tSetIndexPosition );
                tPDVTypes( tMeshSetIndex ).push_back( tPDVTypeGroup );
            }
        }

        // Set interpolation PDV types in host manager
        mPDVHostManager.set_interpolation_pdv_types( tPDVTypes );

        // Get and save communication map from IP mesh
        Matrix< IdMat > tCommTable = tInterpolationMesh->get_communication_table();
        mPDVHostManager.set_communication_table( tCommTable );

        // Get and save global to local vertex maps from IP and IG meshes
        std::unordered_map< moris_id, moris_index > tIPVertexGlobaToLocalMap =
                tInterpolationMesh->get_vertex_glb_id_to_loc_vertex_ind_map();
        std::unordered_map< moris_id, moris_index > tIGVertexGlobaToLocalMap =
                tIntegrationMesh->get_vertex_glb_id_to_loc_vertex_ind_map();
        mPDVHostManager.set_vertex_global_to_local_maps(
                tIPVertexGlobaToLocalMap,
                tIGVertexGlobaToLocalMap );

        // Create PDV hosts
        this->create_interpolation_pdvs(
                tInterpolationMesh,
                tIntegrationMesh,
                tPDVTypes );

        // Set integration PDV types
        if ( mShapeSensitivities )
        {
            // Set integration PDV types
            this->set_integration_pdv_types( tIntegrationMesh );
        }

        // Create PDV IDs
        mPDVHostManager.create_pdv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::print_gen_vertices(
            const std::string& aFile,
            mtk::Mesh*         aMesh )
    {
        std::ostringstream tStringStream;
        tStringStream.clear();
        tStringStream.str( "" );

        tStringStream << "Vert_Id,";
        tStringStream << "Vert Ind,";
        tStringStream << "Owner,";
        tStringStream << "Prank,";
        for ( moris::uint iHeader = 0; iHeader < aMesh->get_spatial_dim(); iHeader++ )
        {
            tStringStream << "Coords_" + std::to_string( iHeader ) << ",";
        }

        for ( moris::uint iHeader = 0; iHeader < this->get_number_of_geometries(); iHeader++ )
        {
            tStringStream << "gval_" + std::to_string( iHeader ) << ",";
        }
        for ( moris::uint iHeader = 0; iHeader < this->get_number_of_geometries(); iHeader++ )
        {
            tStringStream << "gprox_" + std::to_string( iHeader );
            if ( iHeader != this->get_number_of_geometries() - 1 )
            {
                tStringStream << ",";
            }
        }
        tStringStream << std::endl;
        // iterate through vertices
        for ( moris::uint iV = 0; iV < aMesh->get_num_nodes(); iV++ )
        {
            mtk::Vertex& tVertex = aMesh->get_mtk_vertex( (moris_index)iV );
            tStringStream << tVertex.get_id() << ",";
            tStringStream << tVertex.get_index() << ",";
            tStringStream << tVertex.get_owner() << ",";
            tStringStream << par_rank() << ",";
            moris::Matrix< moris::DDRMat > tCoords = tVertex.get_coords();

            for ( moris::uint iSp = 0; iSp < aMesh->get_spatial_dim(); iSp++ )
            {
                tStringStream << std::scientific << tCoords( iSp ) << ",";
            }
            for ( moris::uint iGeom = 0; iGeom < this->get_number_of_geometries(); iGeom++ )
            {
                Vector< real > tGeometryInfo;
                mGeometries( iGeom )->get_design_info( tVertex.get_index(), tCoords, tGeometryInfo );

                for ( uint iGeomFields = 0; iGeomFields < mGeometries( iGeom )->get_num_fields(); iGeomFields++ )
                {
                    tStringStream << tGeometryInfo( iGeomFields ) << ",";
                }
            }
            for ( moris::uint iGeom = 0; iGeom < this->get_number_of_geometries(); iGeom++ )
            {
                // Get geometric region
                Geometric_Region tGeometricRegion = mGeometries( iGeom )->get_geometric_region( tVertex.get_index(), tVertex.get_coords() );

                // Add to string stream based on region
                if ( tGeometricRegion == Geometric_Region::NEGATIVE )
                {
                    tStringStream << "-";
                }
                else if ( tGeometricRegion == Geometric_Region::INTERFACE )
                {
                    tStringStream << "=";
                }
                else
                {
                    tStringStream << "+";
                }
                if ( iGeom != this->get_number_of_geometries() - 1 )
                {
                    tStringStream << ",";
                }
            }
            tStringStream << std::endl;
        }
        if ( not aFile.empty() )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::setup_diagnostics(
            bool                aDiagnostics,
            std::string const & aDiagnosticPath,
            std::string const & aDiagnosticLabel )
    {
        mDiagnostics = aDiagnostics;

        if ( mDiagnostics )
        {
            mDiagnosticPath = aDiagnosticPath;
            mDiagnosticId   = aDiagnosticLabel;

            MORIS_ERROR( !mDiagnosticPath.empty(), "If diagnostics are turned on, a diagnostics path must be specified" );
            if ( mDiagnosticId.empty() )
            {
                mDiagnosticId = "no_spec";
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Geometry_Engine::get_diagnostic_file_name( std::string const & aLabel ) const
    {
        MORIS_ASSERT( mDiagnostics, "Only callable with diagnostics on" );
        return mDiagnosticPath + "/id_" + mDiagnosticId + "_p_" + std::to_string( moris::par_rank() ) + "_" + aLabel + ".diag";
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::initialize_pdv_type_list()
    {
        // Reserve of temporary pdv type list
        Vector< enum PDV_Type > tTemporaryPDVTypeList;

        tTemporaryPDVTypeList.reserve( static_cast< int >( PDV_Type::UNDEFINED ) + 1 );

        Matrix< DDUMat > tListToCheckIfEnumExist( ( static_cast< int >( PDV_Type::UNDEFINED ) + 1 ), 1, 0 );

        // PDV type map
        map< std::string, PDV_Type > tPDVTypeMap = get_pdv_type_map();

        // Loop over properties to build parallel consistent pdv list
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            // PDV type and mesh set names/indices from parameter list
            PDV_Type tPDVType = mProperties( tPropertyIndex )->get_pdv_type();

            if ( tListToCheckIfEnumExist( static_cast< int >( tPDVType ), 0 ) == 0 )
            {
                // Set 1 at position of the enum value
                tListToCheckIfEnumExist( static_cast< int >( tPDVType ), 0 ) = 1;

                tTemporaryPDVTypeList.push_back( tPDVType );
            }
        }

        // Shrink pdv type list to fit
        tTemporaryPDVTypeList.shrink_to_fit();

        // Communicate dof types so that all processors have the same unique list
        mPDVHostManager.communicate_dof_types( tTemporaryPDVTypeList );

        // Create a map
        mPDVHostManager.create_dv_type_map();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::distribute_advs(
            mtk::Mesh_Pair                          aMeshPair,
            Vector< std::shared_ptr< mtk::Field > > aFields,
            mtk::EntityRank                         aADVEntityRank )
    {
        // Tracer
        Tracer tTracer( "GEN", "Distribute ADVs" );

        // Gather all designs
        Vector< std::shared_ptr< Design > > tDesigns( mGeometries.size() + mProperties.size() );
        std::copy( mGeometries.begin(), mGeometries.end(), tDesigns.begin() );
        std::copy( mProperties.begin(), mProperties.end(), tDesigns.begin() + mGeometries.size() );

        // Get interpolation mesh
        mtk::Interpolation_Mesh* tMesh = aMeshPair.get_interpolation_mesh();

        //------------------------------------//
        // Determine owned and shared ADV IDs //
        //------------------------------------//
        clock_t tStart_Owned_Shared_ADVs = clock();

        // Set primitive IDs
        Matrix< DDSMat > tPrimitiveADVIds( mInitialPrimitiveADVs.length(), 1 );
        for ( uint iADVIndex = 0; iADVIndex < mInitialPrimitiveADVs.length(); iADVIndex++ )
        {
            tPrimitiveADVIds( iADVIndex ) = iADVIndex;
        }

        // Start with primitive IDs for owned IDs on processor 0
        Matrix< DDSMat > tOwnedADVIds( 0, 0 );
        if ( par_rank() == 0 )
        {
            tOwnedADVIds = tPrimitiveADVIds;
        }

        // this is done to initialize primitive adv positions with gNoID
        Matrix< IdMat > tOwnedijklIDs( tPrimitiveADVIds.numel(), 1, gNoID );

        // Owned and shared ADVs per field
        // Vector< Matrix< DDSMat > > tSharedADVIds( tDesigns.size() ); BRENDAN
        Matrix< DDUMat >           tAllOffsetIDs( tDesigns.size(), 1 );

        // Loop over all geometries to get number of new ADVs
        sint tOffsetID = tPrimitiveADVIds.length();
        for ( uint iDesignIndex = 0; iDesignIndex < tDesigns.size(); iDesignIndex++ )
        {
            tOffsetID = tDesigns( iDesignIndex )->append_adv_info(    //
                    aMeshPair.get_interpolation_mesh(),
                    tOwnedADVIds,
                    tOwnedijklIDs,
                    tOffsetID,
                    mLowerBounds,
                    mUpperBounds );
        }

        // Set owned ADV IDs
        mPDVHostManager.set_owned_adv_ids( tOwnedADVIds );

        MORIS_LOG_INFO( "Time to collect owned and shared ADVs: %f sec", ( moris::real )( clock() - tStart_Owned_Shared_ADVs ) / CLOCKS_PER_SEC );

        //----------------------------------------//
        // Create owned ADV vector                //
        //----------------------------------------//
        clock_t tStart_Create_Owned_ADVs = clock();

        // Create factory for distributed ADV vector
        sol::Matrix_Vector_Factory tDistributedFactory;

        // Create owned ADV vector
        sol::Dist_Map*    tOwnedADVMap  = tDistributedFactory.create_map( tOwnedADVIds );
        sol::Dist_Vector* tNewOwnedADVs = tDistributedFactory.create_vector( tOwnedADVMap, 1, false, true );

        // Determine if primitive ADVs have been created already
        if ( mPrimitiveADVs )
        {
            // Assign old primitive ADVs to initial vector
            for ( uint iADVIndex = 0; iADVIndex < mInitialPrimitiveADVs.length(); iADVIndex++ )
            {
                mInitialPrimitiveADVs( iADVIndex ) = mPrimitiveADVs->operator()( iADVIndex, 0 );
            }
        }
        else
        {
            // Create primitive ADV vector
            sol::Dist_Map* tPrimitiveADVMap = tDistributedFactory.create_map( tPrimitiveADVIds );
            mPrimitiveADVs                  = tDistributedFactory.create_vector( tPrimitiveADVMap, 1, false, true );
        }

        // Assign primitive ADVs to the owned vector
        if ( par_rank() == 0 )
        {
            tNewOwnedADVs->replace_global_values( tPrimitiveADVIds, mInitialPrimitiveADVs );
        }

        // Global assembly
        tNewOwnedADVs->vector_global_assembly();

        // Get primitive ADVs from owned vector
        mPrimitiveADVs->import_local_to_global( *tNewOwnedADVs );

        // Set field ADVs using distributed vector
        if ( mInitialPrimitiveADVs.length() > 0 )
        {
            for ( uint tDesignIndex = 0; tDesignIndex < tDesigns.size(); tDesignIndex++ )
            {
                tDesigns( tDesignIndex )->set_advs( mPrimitiveADVs );
            }
        }

        MORIS_LOG_INFO( "Time to create owned ADVs: %f sec", ( moris::real )( clock() - tStart_Create_Owned_ADVs ) / CLOCKS_PER_SEC );

        //----------------------------------------//
        // Convert to B-spline fields             //
        //----------------------------------------//
        clock_t tStart_Convert_to_Bspline_Fields = clock();

        // Loop to discretize geometries when requested
        // FIXME: make the check for the field name inside the geometry, remove the field index parameter,
        for ( uint iGeometryIndex = 0; iGeometryIndex < mGeometries.size(); iGeometryIndex++ )
        {
            for ( uint iGeometryFieldIndex = 0; iGeometryFieldIndex < mGeometries( iGeometryIndex )->get_num_fields(); iGeometryFieldIndex++ )
            {
                // Loop over MTK fields to find a match
                bool tUseMTKField = false;
                for ( const auto& iMTKField : aFields )
                {
                    mGeometries( iGeometryIndex )->discretize( iMTKField, aMeshPair, tNewOwnedADVs );
                    tUseMTKField = true;
                    break;
                }

                // Otherwise discretize with original field
                if ( not tUseMTKField )
                {
                    mGeometries( iGeometryIndex )->discretize( aMeshPair, tNewOwnedADVs );
                }

                // Shape sensitivities logic}
                mShapeSensitivities = ( mShapeSensitivities or mGeometries( iGeometryIndex )->depends_on_advs() );
            }
        }

        // Loop to discretize properties when requested
        for ( uint iPropertyIndex = 0; iPropertyIndex < mProperties.size(); iPropertyIndex++ )
        {
            // Loop over MTK fields to find a match
            bool tUseMTKField = false;
            for ( const auto& iMTKField : aFields )
            {
                if ( mProperties( iPropertyIndex )->get_name() == iMTKField->get_label() )
                {
                    mProperties( iPropertyIndex )->discretize( iMTKField, aMeshPair, tNewOwnedADVs );
                    tUseMTKField = true;
                    break;
                }
            }

            // Otherwise discretize with original field
            if ( not tUseMTKField )
            {
                mProperties( iPropertyIndex )->discretize( aMeshPair, tNewOwnedADVs );
            }
        }

        // Register node manager with each geometry/property TODO figure out a better way to do this; have node manager automatically copy over when discretizing?
        for ( const auto& iGeometry : mGeometries )
        {
            iGeometry->set_node_manager( mNodeManager );
        }
        for ( const auto& iProperty : mProperties )
        {
            iProperty->set_node_manager( mNodeManager );
        }

        // Update dependencies
        if ( tDesigns.size() > 0 )
        {
            std::copy( mGeometries.begin(), mGeometries.end(), tDesigns.begin() );
            std::copy( mProperties.begin(), mProperties.end(), tDesigns.begin() + mGeometries.size() );
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                mProperties( tPropertyIndex )->update_dependencies( tDesigns );
            }
        }

        // Delete old owned ADVs and save ADVs
        delete mOwnedADVs;
        mOwnedADVs = tNewOwnedADVs;

        MORIS_LOG_INFO( "Time to convert to Bspline fields: %f sec", ( moris::real )( clock() - tStart_Convert_to_Bspline_Fields ) / CLOCKS_PER_SEC );

        //----------------------------------------//
        // Communicate all ADV IDs to processor 0 //
        //----------------------------------------//

        clock_t tStart_Communicate_ADV_IDs = clock();

        // Sending mats
        Vector< Matrix< DDSMat > > tSendingIDs( 0 );
        Vector< Matrix< DDRMat > > tSendingLowerBounds( 0 );
        Vector< Matrix< DDRMat > > tSendingUpperBounds( 0 );
        Vector< Matrix< DDSMat > > tSendingijklIDs( 0 );

        // Receiving mats
        Vector< Matrix< DDSMat > > tReceivingIDs( 0 );
        Vector< Matrix< DDRMat > > tReceivingLowerBounds( 0 );
        Vector< Matrix< DDRMat > > tReceivingUpperBounds( 0 );
        Vector< Matrix< DDSMat > > tReceivingjklIDs( 0 );

        // Set up communication list for communicating ADV IDs
        Matrix< IdMat > tCommunicationList( 1, 1, 0 );
        if ( par_rank() == 0 )
        {
            // Resize communication list and sending mats
            tCommunicationList.resize( par_size() - 1, 1 );
            tSendingIDs.resize( par_size() - 1 );
            tSendingLowerBounds.resize( par_size() - 1 );
            tSendingUpperBounds.resize( par_size() - 1 );
            tSendingijklIDs.resize( par_size() - 1 );

            // Assign communication list
            for ( uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++ )
            {
                tCommunicationList( tProcessorIndex - 1 ) = tProcessorIndex;
            }
        }
        else
        {
            tSendingIDs         = { tOwnedADVIds };
            tSendingLowerBounds = { mLowerBounds };
            tSendingUpperBounds = { mUpperBounds };
            tSendingijklIDs     = { tOwnedijklIDs };
        }

        // Communicate mats
        communicate_mats( tCommunicationList, tSendingIDs, tReceivingIDs );
        communicate_mats( tCommunicationList, tSendingLowerBounds, tReceivingLowerBounds );
        communicate_mats( tCommunicationList, tSendingUpperBounds, tReceivingUpperBounds );
        if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
        {
            communicate_mats( tCommunicationList, tSendingijklIDs, tReceivingjklIDs );
        }

        MORIS_LOG_INFO( "Time to communicate ADV IDs: %f sec", ( moris::real )( clock() - tStart_Communicate_ADV_IDs ) / CLOCKS_PER_SEC );

        // Assemble full ADVs/bounds
        clock_t tStart_ADV_Bounds = clock();

        if ( par_rank() == 0 )
        {
            // Start full IDs with owned IDs on processor 0
            mFullADVIds = tOwnedADVIds;
            if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
            {
                mFullijklIDs = tOwnedijklIDs;
            }
            else
            {
                mFullijklIDs.set_size( 0, 0 );
            }

            // Assemble additional IDs/bounds from other processors
            for ( uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++ )
            {
                // Get number of received ADVs
                uint tFullADVsFilled    = mFullADVIds.length();
                uint tFullijklIDsFilled = mFullijklIDs.length();
                uint tNumReceivedADVs   = tReceivingIDs( tProcessorIndex - 1 ).length();

                // Resize full ADV IDs and bounds
                mFullADVIds.resize( tFullADVsFilled + tNumReceivedADVs, 1 );
                mLowerBounds.resize( tFullADVsFilled + tNumReceivedADVs, 1 );
                mUpperBounds.resize( tFullADVsFilled + tNumReceivedADVs, 1 );
                if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
                {
                    mFullijklIDs.resize( tFullijklIDsFilled + tNumReceivedADVs, 1 );
                }

                // Assign received ADV IDs
                for ( uint tADVIndex = 0; tADVIndex < tNumReceivedADVs; tADVIndex++ )
                {
                    mFullADVIds( tFullADVsFilled + tADVIndex ) = tReceivingIDs( tProcessorIndex - 1 )( tADVIndex );
                    mLowerBounds( tFullADVsFilled + tADVIndex ) =
                            tReceivingLowerBounds( tProcessorIndex - 1 )( tADVIndex );
                    mUpperBounds( tFullADVsFilled + tADVIndex ) =
                            tReceivingUpperBounds( tProcessorIndex - 1 )( tADVIndex );

                    if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
                    {
                        mFullijklIDs( tFullijklIDsFilled + tADVIndex ) =
                                tReceivingjklIDs( tProcessorIndex - 1 )( tADVIndex );
                    }
                }
            }
        }
        else
        {
            mLowerBounds.set_size( 0, 0 );
            mUpperBounds.set_size( 0, 0 );
            mFullijklIDs.set_size( 0, 0 );
        }

        MORIS_LOG_INFO( "Time to assemble ADVs and bounds on Proc 0: %f sec", ( moris::real )( clock() - tStart_ADV_Bounds ) / CLOCKS_PER_SEC );

        // Reset mesh information
        clock_t tStart_Reset_Mesh_Info = clock();
        this->reset_mesh_information( tMesh );

        MORIS_LOG_INFO( "Time to reset mesh information: %f sec", ( moris::real )( clock() - tStart_Reset_Mesh_Info ) / CLOCKS_PER_SEC );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::reset_mesh_information( mtk::Interpolation_Mesh* aMesh )
    {
        // Register spatial dimension
        mNumSpatialDimensions = aMesh->get_spatial_dim();

        // Set GEN nodes
        mNodeManager.reset_background_nodes( aMesh );

        // Reset PDV host manager
        mPDVHostManager.reset();

        // Reset info related to the mesh
        mActiveGeometryIndex = 0;
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            mGeometries( tGeometryIndex )->reset_nodal_data( aMesh );
        }
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            mProperties( tPropertyIndex )->reset_nodal_data( aMesh );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::output_fields( mtk::Mesh* aMesh )
    {
        // Tracer
        Tracer tTracer( "GEN", "Output fields" );

        this->output_fields_on_mesh( aMesh, mOutputMeshFile );
        this->write_geometry_fields( aMesh, mGeometryFieldFile );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::output_fields_on_mesh(
            mtk::Mesh*  aMesh,
            std::string aExodusFileName )
    {
        // time shift
        real tTimeShift = 0.0;

        if ( not aExodusFileName.empty() )
        {
            if ( mTimeOffset > 0 )
            {
                // get optimization iteration
                uint tOptIter = gLogger.get_opt_iteration();

                // set name
                std::string tOptIterStrg = std::to_string( tOptIter );
                aExodusFileName += ".e-s." + std::string( 4 - tOptIterStrg.length(), '0' ) + tOptIterStrg;

                // determine time shift
                tTimeShift = tOptIter * mTimeOffset;
            }

            // Write mesh
            mtk::Writer_Exodus tWriter( aMesh );
            tWriter.write_mesh( "./", aExodusFileName, "./", "gen_temp.exo" );

            // Setup field names
            uint tNumGeometryFields = 0;
            for ( uint iGeom = 0; iGeom < mGeometries.size(); iGeom++ )
            {
                tNumGeometryFields += mGeometries( iGeom )->get_num_fields();
            }
            Vector< std::string > tFieldNames(0 );

            // Geometry field names
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                // Get the geometries field names
                Vector< std::string > tGeometryFieldNames = mGeometries( tGeometryIndex )->get_field_names();

                // Append the field names to the list
                tFieldNames.append( tGeometryFieldNames );
            }

            // Property field names
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                // Get the geometries field names
                Vector< std::string > tPropertyFieldNames = mProperties( tPropertyIndex )->get_field_names();

                // Append the field names to the list
                tFieldNames.append( tPropertyFieldNames );
            }

            // write time to file
            tWriter.set_time( tTimeShift );

            // Set nodal fields based on field names
            tWriter.set_nodal_fields( tFieldNames );

            // Get all node coordinates
            Vector< Matrix< DDRMat > > tNodeCoordinates( aMesh->get_num_nodes() );
            for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
            {
                tNodeCoordinates( tNodeIndex ) = aMesh->get_node_coordinate( tNodeIndex );
            }

            // Loop over geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                for ( uint iGeometryFieldIndex = 0; iGeometryFieldIndex < mGeometries( tGeometryIndex )->get_num_fields(); iGeometryFieldIndex++ )
                {
                    // Create field vector
                    Matrix< DDRMat > tFieldData( aMesh->get_num_nodes(), 1 );

                    for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                    {
                        // Get design info from the geometry
                        Vector< real > tGeometryInfo;
                        mGeometries( tGeometryIndex )->get_design_info( tNodeIndex, tNodeCoordinates( tNodeIndex ), tGeometryInfo );

                        // Assign field to vector
                        tFieldData( tNodeIndex ) = tGeometryInfo( iGeometryFieldIndex );
                    }

                    // Create field on mesh
                    tWriter.write_nodal_field( tFieldNames( iGeometryFieldIndex ), tFieldData );
                }
            }

            // Loop over properties
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                for ( uint iPropertyFieldIndex = 0; iPropertyFieldIndex < mProperties( tPropertyIndex )->get_num_fields(); iPropertyFieldIndex++ )
                {
                    // Create field vector
                    Matrix< DDRMat > tFieldData( aMesh->get_num_nodes(), 1 );

                    // Loop over all nodes on the mesh
                    for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                    {
                        // Get design info from the property
                        Vector< real > tPropertyInfo;
                        mProperties( tPropertyIndex )->get_design_info( tNodeIndex, tNodeCoordinates( tNodeIndex ), tPropertyInfo );

                        // Assign field to vector
                        tFieldData( tNodeIndex ) = tPropertyInfo( iPropertyFieldIndex );
                    }

                    // Create field on mesh
                    tWriter.write_nodal_field( tFieldNames( tNumGeometryFields + tPropertyIndex ), tFieldData );
                }
            }

            // Finalize
            tWriter.close_file( true );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Geometry_Engine::write_geometry_fields(
            mtk::Mesh*         aMesh,
            const std::string& aBaseFileName )
    {
        if ( not aBaseFileName.empty() )
        {
            // Get all node coordinates
            Vector< Matrix< DDRMat > > tNodeCoordinates( aMesh->get_num_nodes() );
            for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
            {
                tNodeCoordinates( tNodeIndex ) = aMesh->get_node_coordinate( tNodeIndex );
            }

            // Loop over geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                // Loop over fields for this geometry
                Vector< std::ofstream > tOutFiles( mGeometries( tGeometryIndex )->get_num_fields() );
                for ( uint iGeometryFieldIndex = 0; iGeometryFieldIndex < mGeometries( tGeometryIndex )->get_num_fields(); iGeometryFieldIndex++ )
                {
                    // Create file
                    tOutFiles( iGeometryFieldIndex ).open( aBaseFileName + "_" + std::to_string( tGeometryIndex ) + "field_" + std::to_string( iGeometryFieldIndex ) + ".txt" );
                }

                // Write to file
                for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                {
                    Vector< real > tGeometryInfo;
                    mGeometries( tGeometryIndex )->get_design_info( tNodeIndex, tNodeCoordinates( tNodeIndex ), tGeometryInfo );

                    for ( uint iGeometryFieldIndex = 0; iGeometryFieldIndex < mGeometries( tGeometryIndex )->get_num_fields(); iGeometryFieldIndex++ )
                    {
                        // Coordinates
                        for ( uint tDimension = 0; tDimension < mNumSpatialDimensions; tDimension++ )
                        {
                            tOutFiles( iGeometryFieldIndex ) << tNodeCoordinates( tNodeIndex )( tDimension ) << ", ";
                        }

                        // Fill unused dimensions with zeros
                        for ( uint tDimension = mNumSpatialDimensions; tDimension < 3; tDimension++ )
                        {
                            tOutFiles( iGeometryFieldIndex ) << 0.0 << ", ";
                        }

                        // Level-set field
                        tOutFiles( iGeometryFieldIndex ) << tGeometryInfo( iGeometryFieldIndex ) << std::endl;
                    }
                }

                // Close files
                for ( uint iOutFile = 0; iOutFile < tOutFiles.size(); iOutFile++ )
                {
                    tOutFiles( iOutFile ).close();
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    // PRIVATE
    //--------------------------------------------------------------------------------------------------------------
    void Geometry_Engine::create_interpolation_pdvs(
            mtk::Interpolation_Mesh*               aInterpolationMesh,
            mtk::Integration_Mesh*                 aIntegrationMesh,
            Vector< Vector< Vector< PDV_Type > > > aPDVTypes )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create interpolation PDV hosts" );

        // Get information from integration mesh
        uint tNumSets = aPDVTypes.size();

        // Size node information cells
        Vector< Vector< uint > >   tNodeIndicesPerSet( tNumSets );
        Vector< Vector< sint > >   tNodeIdsPerSet( tNumSets );
        Vector< Vector< uint > >   tNodeOwnersPerSet( tNumSets );
        Vector< Matrix< DDRMat > > tNodeCoordinatesPerSet( tNumSets );

        // Get communication table and map
        Matrix< IdMat >    tCommTable             = aInterpolationMesh->get_communication_table();
        Vector< moris_id > tCommunicationTableMap = build_communication_table_map( tCommTable );

        // TODO change over to just use a cell to begin with
        Vector< moris_index > tCommunicationTable( tCommTable.length() );
        for ( uint iCommTableIndex = 0; iCommTableIndex < tCommunicationTable.size(); iCommTableIndex++ )
        {
            tCommunicationTable( iCommTableIndex ) = tCommTable( iCommTableIndex );
        }

        // Loop through sets in integration mesh
        for ( uint iMeshSetIndex = 0; iMeshSetIndex < tNumSets; iMeshSetIndex++ )
        {
            // Determine number of nodes if there are PDVs on this set
            if ( aPDVTypes( iMeshSetIndex ).size() > 0 )
            {
                // Get set pointer
                mtk::Set* tSet = aIntegrationMesh->get_set_by_index( iMeshSetIndex );

                // Select sides of interpolation cells to get info from
                Vector< mtk::Leader_Follower > tSetSides = mtk::get_leader_follower( tSet->get_set_type() );

                // Get number of clusters on set
                uint tNumberOfClusters = tSet->get_num_clusters_on_set();

                // Number of nodes on this set
                uint tNumberOfNodesInSet = 0;

                // Number of shared nodes on this set per proc
                Vector< uint > tNumSharedNodesPerProc( tCommunicationTable.size() );

                // Loop over clusters on this set to count nodes
                for ( uint tClusterIndex = 0; tClusterIndex < tNumberOfClusters; tClusterIndex++ )
                {
                    // get pointer to cluster
                    const mtk::Cluster* tCluster = tSet->get_clusters_by_index( tClusterIndex );

                    // Loop over leader/follower
                    for ( mtk::Leader_Follower iLeaderFollower : tSetSides )
                    {
                        // Get node indices in cluster
                        mtk::Cell const *  tBaseCell            = tCluster->get_interpolation_cell( iLeaderFollower ).get_base_cell();
                        Matrix< IndexMat > tNodeOwnersInCluster = tBaseCell->get_vertex_owners();

                        // Add to the number of base nodes on this set
                        tNumberOfNodesInSet += tNodeOwnersInCluster.length();

                        // Determine if we need to check for shared nodes
                        if ( par_size() > 1 )
                        {
                            // Get number of shared nodes with each proc
                            for ( uint iNodeInCluster = 0; iNodeInCluster < tNodeOwnersInCluster.length(); iNodeInCluster++ )
                            {
                                // Determine if this proc is not node owner
                                moris_index tNodeOwner = tNodeOwnersInCluster( iNodeInCluster );
                                if ( tNodeOwner not_eq par_rank() )
                                {
                                    // Count up number of shared nodes with that proc
                                    tNumSharedNodesPerProc( tCommunicationTableMap( tNodeOwner ) )++;
                                }
                            }
                        }
                    }
                }

                // Communicate to owning proc about shared nodes
                Vector< uint > tNumOwnedNodesPerProc( tCommunicationTable.size() );
                communicate_scalars( tCommunicationTable, tNumSharedNodesPerProc, tNumOwnedNodesPerProc );

                // Add number of nodes this proc owns that it may not know is on the set
                tNumberOfNodesInSet += std::accumulate( tNumOwnedNodesPerProc.begin(), tNumOwnedNodesPerProc.end(), 0 );

                // Resize node indices, IDs, owners, and coordinates for this set
                tNodeIndicesPerSet( iMeshSetIndex ).resize( tNumberOfNodesInSet, 0 );
                tNodeIdsPerSet( iMeshSetIndex ).resize( tNumberOfNodesInSet, -1 );
                tNodeOwnersPerSet( iMeshSetIndex ).resize( tNumberOfNodesInSet, 0 );
                tNodeCoordinatesPerSet( iMeshSetIndex ).resize( tNumberOfNodesInSet, mNumSpatialDimensions );

                // Loop over clusters on this set
                uint tCurrentNode = 0;
                for ( uint iClusterIndex = 0; iClusterIndex < tNumberOfClusters; iClusterIndex++ )
                {
                    // get pointer to cluster
                    const mtk::Cluster* tCluster = tSet->get_clusters_by_index( iClusterIndex );

                    // Loop over leader/follower
                    for ( mtk::Leader_Follower iLeaderFollower : tSetSides )
                    {
                        // Get base interpolation cell
                        mtk::Cell const * tBaseCell = tCluster->get_interpolation_cell( iLeaderFollower ).get_base_cell();

                        // Indices, IDs, and ownership of base cell nodes in cluster
                        Matrix< IndexMat > tNodeIndicesInCluster     = tBaseCell->get_vertex_inds();
                        Matrix< IdMat >    tNodeIdsInCluster         = tBaseCell->get_vertex_ids();
                        Matrix< IndexMat > tNodeOwnersInCluster      = tBaseCell->get_vertex_owners();
                        Matrix< DDRMat >   tNodeCoordinatesInCluster = tBaseCell->get_vertex_coords();

                        // number of base nodes in cluster
                        uint tNumberOfBaseNodes = tNodeIndicesInCluster.length();

                        // check for consistency
                        MORIS_ASSERT( tNodeIdsInCluster.length() == tNumberOfBaseNodes and tNodeOwnersInCluster.length() == tNumberOfBaseNodes,
                                "Geometry_Engine::create_interpolation_pdvs - inconsistent cluster information.\n" );

                        // FIXME This list has duplicate entries. Functionality is working as is, but is slightly slower. Mesh needs to provide unique nodes per set to fix this.
                        for ( uint iNodeInCluster = 0; iNodeInCluster < tNumberOfBaseNodes; iNodeInCluster++ )
                        {
                            // Set node index/ID/owner info
                            tNodeIndicesPerSet( iMeshSetIndex )( tCurrentNode ) = tNodeIndicesInCluster( iNodeInCluster );
                            tNodeIdsPerSet( iMeshSetIndex )( tCurrentNode )     = tNodeIdsInCluster( iNodeInCluster );
                            tNodeOwnersPerSet( iMeshSetIndex )( tCurrentNode )  = tNodeOwnersInCluster( iNodeInCluster );
                            tNodeCoordinatesPerSet( iMeshSetIndex ).get_row( tCurrentNode ) =
                                    tNodeCoordinatesInCluster.get_row( iNodeInCluster );

                            // Increment index in overall lists
                            tCurrentNode++;
                        }
                    }
                }

                // Determine if we need to check for shared nodes
                if ( par_size() > 1 )
                {
                    // Create lists of shared nodes to communicate
                    Vector< Vector< sint > > tSharedNodeIdsOnSet( tCommunicationTable.size() );
                    for ( uint iProcIndex = 0; iProcIndex < tCommunicationTable.size(); iProcIndex++ )
                    {
                        tSharedNodeIdsOnSet( iProcIndex ).resize( tNumSharedNodesPerProc( iProcIndex ) );
                        tNumSharedNodesPerProc( iProcIndex ) = 0;
                    }

                    // Copy over shared node IDs from currently known nodes
                    for ( uint iNodeInSet = 0; iNodeInSet < tCurrentNode; iNodeInSet++ )
                    {
                        // Determine if this proc is not node owner
                        moris_index tNodeOwner = tNodeOwnersPerSet( iMeshSetIndex )( iNodeInSet );
                        if ( tNodeOwner not_eq par_rank() )
                        {
                            // Get mapped proc index
                            uint tCommunicationProcIndex = tCommunicationTableMap( tNodeOwner );

                            // Add node ID to communication list
                            tSharedNodeIdsOnSet( tCommunicationProcIndex )( tNumSharedNodesPerProc( tCommunicationProcIndex )++ ) = tNodeIdsPerSet( iMeshSetIndex )( iNodeInSet );
                        }
                    }

                    // Create owned ID cell
                    Vector< Vector< sint > > tOwnedNodeIdsOnSet( tCommunicationTable.size() );

                    // Communicate IDs of shared nodes to the owning processor
                    communicate_cells(
                            tCommunicationTable,
                            tSharedNodeIdsOnSet,
                            tOwnedNodeIdsOnSet );

                    // Add new node information to this set
                    for ( uint iProcIndex = 0; iProcIndex < tCommunicationTable.size(); iProcIndex++ )
                    {

                        for ( uint iSharedNodeListIndex = 0; iSharedNodeListIndex < tOwnedNodeIdsOnSet( iProcIndex ).size(); iSharedNodeListIndex++ )
                        {
                            // Get owned node ID
                            moris_id tNodeId = tOwnedNodeIdsOnSet( iProcIndex )( iSharedNodeListIndex );

                            // Get index on this proc
                            moris_index tNodeIndex = aInterpolationMesh->get_background_mesh().get_loc_entity_ind_from_entity_glb_id( tNodeId, mtk::EntityRank::NODE );

                            // Get coordinates
                            Matrix< DDRMat > tNodeCoordinates = aInterpolationMesh->get_node_coordinate( tNodeIndex );

                            // Set node index/ID/owner/coordinates
                            tNodeIndicesPerSet( iMeshSetIndex )( tCurrentNode ) = tNodeIndex;
                            tNodeIdsPerSet( iMeshSetIndex )( tCurrentNode )     = tNodeId;
                            tNodeOwnersPerSet( iMeshSetIndex )( tCurrentNode )  = par_rank();
                            tNodeCoordinatesPerSet( iMeshSetIndex ).set_row( tCurrentNode, tNodeCoordinates );

                            // Increment index in overall lists
                            tCurrentNode++;
                        }
                    }
                }
            }
        }

        // Create hosts of nodes of non-unzipped interpolation nodes
        mPDVHostManager.create_interpolation_pdv_hosts(
                tNodeIndicesPerSet,
                tNodeIdsPerSet,
                tNodeOwnersPerSet,
                tNodeCoordinatesPerSet );

        // Loop over properties to assign PDVs on this set
        for ( const auto& iProperty : mProperties )
        {
            // Check if this is an interpolation PDV
            MORIS_ERROR( iProperty->is_interpolation_pdv(),
                    "Assignment of PDVs is only supported with an interpolation mesh right now." );

            // Get PDV type and all mesh set indices for this property
            PDV_Type       tPDVType        = iProperty->get_pdv_type();
            Vector< uint > tMeshSetIndices = iProperty->get_pdv_mesh_set_indices( aIntegrationMesh );

            // Loop through nodes in these sets
            for ( uint iMeshSetIndex : tMeshSetIndices )
            {
                for ( uint iNodeInSet = 0; iNodeInSet < tNodeIndicesPerSet( iMeshSetIndex ).size(); iNodeInSet++ )
                {
                    // Create interpolation PDV
                    mPDVHostManager.create_interpolation_pdv( tNodeIndicesPerSet( iMeshSetIndex )( iNodeInSet ), tPDVType, iProperty );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Geometry_Engine::set_integration_pdv_types( mtk::Integration_Mesh* aIntegrationMesh )
    {
        // Tracer
        Tracer tTracer( "GEN", "Set integration PDV types" );

        // Get information from integration mesh
        uint tNumSets = aIntegrationMesh->get_num_sets();

        // Cell of IG PDV_Type types
        Vector< PDV_Type > tCoordinatePDVs( mNumSpatialDimensions );

        switch ( mNumSpatialDimensions )
        {
            case 2:
            {
                tCoordinatePDVs( 0 ) = PDV_Type::X_COORDINATE;
                tCoordinatePDVs( 1 ) = PDV_Type::Y_COORDINATE;
                break;
            }
            case 3:
            {
                tCoordinatePDVs( 0 ) = PDV_Type::X_COORDINATE;
                tCoordinatePDVs( 1 ) = PDV_Type::Y_COORDINATE;
                tCoordinatePDVs( 2 ) = PDV_Type::Z_COORDINATE;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Geometry Engine only works for 2D and 3D models." );
            }
        }

        // Loop through sets
        Vector< Vector< Vector< PDV_Type > > > tPDVTypes( tNumSets );
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            // PDV_Type types per set
            tPDVTypes( tMeshSetIndex ).resize( 1 );
            tPDVTypes( tMeshSetIndex )( 0 ) = tCoordinatePDVs;
        }

        // Set PDV types
        mPDVHostManager.set_integration_pdv_types( tPDVTypes );
        mPDVHostManager.set_requested_integration_pdv_types( tCoordinatePDVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table
    Geometry_Engine::create_phase_table(
            const Vector< Vector< ParameterList > >& aParameterLists,
            const std::shared_ptr< Library_IO >&     aLibrary )
    {
        // Get number of geometries
        uint tNumGeometries = aParameterLists( 1 ).size();

        // Recreate phase table via different methods if needed
        std::string tPhaseFunctionName = aParameterLists( 0 )( 0 ).get< std::string >( "phase_function_name" );
        if ( not tPhaseFunctionName.empty() )
        {
            // User-defined phase function
            return { aLibrary->load_function< PHASE_FUNCTION >( tPhaseFunctionName ),
                static_cast< uint >( aParameterLists( 0 )( 0 ).get< sint >( "number_of_phases" ) ) };
        }
        else if ( not aParameterLists( 0 )( 0 ).get< std::string >( "phase_table" ).empty() )
        {
            // User-defined bulk phases
            return { tNumGeometries, string_to_mat< DDUMat >( aParameterLists( 0 )( 0 ).get< std::string >( "phase_table" ) ) };
        }
        else
        {
            // Unique phase per geometry combination
            return { tNumGeometries };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table
    Geometry_Engine::create_phase_table(
            uint                    aNumGeometries,
            const Matrix< DDUMat >& aBulkPhases,
            PHASE_FUNCTION          aPhaseFunction,
            uint                    aNumPhases )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create phase table" );

        if ( aPhaseFunction )
        {
            return { aPhaseFunction, aNumPhases };
        }
        else if ( aBulkPhases.length() > 0 )
        {
            return { aNumGeometries, aBulkPhases };
        }
        else
        {
            return { aNumGeometries };
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::gen
