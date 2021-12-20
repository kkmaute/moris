// MRS
#include "fn_Parsing_Tools.hpp"
#include "cl_Tracer.hpp"

// GEN
#include "cl_GEN_Geometry_Engine.hpp"
#include "GEN_typedefs.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "cl_GEN_BSpline_Geometry.hpp"
#include "cl_GEN_BSpline_Property.hpp"
#include "cl_GEN_Stored_Geometry.hpp"
#include "fn_GEN_create_properties.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Intersection_Node_Bilinear.hpp"

// MTK
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

// XTK FIXME
#include "cl_XTK_Topology.hpp"

// SOL FIXME
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

#include "cl_Stopwatch.hpp"//CHR/src

namespace moris
{
namespace ge
{

    //--------------------------------------------------------------------------------------------------------------
    // PUBLIC
    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine::Geometry_Engine(
        Cell< Cell< ParameterList > > aParameterLists,
        std::shared_ptr< Library_IO > aLibrary,
        mtk::Mesh*                    aMesh ) :
        mPhaseTable( create_phase_table( aParameterLists, aLibrary ) )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create geometry engine" );

        // Level set options
        mIsocontourThreshold   = aParameterLists( 0 )( 0 ).get< real >( "isocontour_threshold" );
        mIsocontourTolerance   = aParameterLists( 0 )( 0 ).get< real >( "isocontour_tolerance" );
        mIntersectionTolerance = aParameterLists( 0 )( 0 ).get< real >( "intersection_tolerance" );

        MORIS_ERROR( mIsocontourTolerance > 1e-14,
            "Geometry_Engine::Geometry_Engine - Isocontour tolerance should be larger than 1e-14" );

        MORIS_ERROR( mIntersectionTolerance > 1e-14,
            "Geometry_Engine::Geometry_Engine - Intersection tolerance should be larger than 1e-14" );

        // Requested IQIs
        mRequestedIQIs = string_to_cell< std::string >( aParameterLists( 0 )( 0 ).get< std::string >( "IQI_types" ) );

        // Set library
        mLibrary = aLibrary;

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

        // Geometries
        mGeometries = create_geometries(
            aParameterLists( 1 ),
            mInitialPrimitiveADVs,
            mLibrary,
            aMesh );

        MORIS_ERROR( mGeometries.size() <= MAX_GEOMETRIES,
            "Number of geometries exceeds MAX_GEOMETRIES, please change this in GEN_typedefs.hpp" );

        // Properties
        mProperties = create_properties(
            aParameterLists( 2 ),
            mInitialPrimitiveADVs,
            mGeometries,
            mLibrary );

        // Get intersection mode
        std::string                           tIntersectionModeString = aParameterLists( 0 )( 0 ).get< std::string >( "intersection_mode" );
        map< std::string, Intersection_Mode > tIntersectionModeMap    = get_intersection_mode_map();
        mIntersectionMode                                             = tIntersectionModeMap[tIntersectionModeString];

        // Set requested PDVs
        Cell< std::string >          tRequestedPdvNames = string_to_cell< std::string >( aParameterLists( 0 )( 0 ).get< std::string >( "PDV_types" ) );
        Cell< PDV_Type >             tRequestedPdvTypes( tRequestedPdvNames.size() );
        map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();
        for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tRequestedPdvTypes.size(); tPdvTypeIndex++ )
        {
            tRequestedPdvTypes( tPdvTypeIndex ) = tPdvTypeMap[tRequestedPdvNames( tPdvTypeIndex )];
        }
        mPDVHostManager.set_requested_interpolation_pdv_types( tRequestedPdvTypes );

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
        Geometry_Engine_Parameters aParameters ) :
        mGeometries( aParameters.mGeometries ),
        mProperties( aParameters.mProperties ), mPhaseTable( create_phase_table( aParameters.mGeometries.size(), aParameters.mBulkPhases ) ), mIntersectionMode( aParameters.mIntersectionMode ), mIsocontourThreshold( aParameters.mIsocontourThreshold ), mIsocontourTolerance( aParameters.mIsocontourTolerance ), mIntersectionTolerance( aParameters.mIntersectionTolerance ), mInitialPrimitiveADVs( aParameters.mADVs ), mTimeOffset( aParameters.mTimeOffset )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create geometry engine" );

        MORIS_ERROR( mIsocontourTolerance > 1e-14,
            "Geometry_Engine::Geometry_Engine - Isocontour tolerance should be larger than 1e-14" );

        MORIS_ERROR( mIntersectionTolerance > 1e-14,
            "Geometry_Engine::Geometry_Engine - Intersection tolerance should be larger than 1e-14" );

        mtk::Integration_Mesh* tIntegrationMesh = create_integration_mesh_from_interpolation_mesh(
            aMesh->get_mesh_type(),
            aMesh );

        mtk::Mesh_Pair tMeshPair( aMesh, tIntegrationMesh );
        this->distribute_advs( tMeshPair, {} );
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine::~Geometry_Engine()
    {
        delete mOwnedADVs;
        delete mPrimitiveADVs;
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

    void
    Geometry_Engine::communicate_requested_IQIs( Cell< std::string > aIQINames )
    {
        mPDVHostManager.set_requested_IQIs( aIQINames );
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

    bool
    Geometry_Engine::geometric_query( Geometric_Query_Interface* aGeometricQuery )
    {
        if ( aGeometricQuery->get_query_type() == Query_Type::INTERSECTION_NO_LOCATION )
        {
            mActiveGeometryIndex = aGeometricQuery->get_geomtric_index();
            return is_intersected( aGeometricQuery->get_query_entity_to_vertex_connectivity(),
                aGeometricQuery->get_query_indexed_coordinates() );
        }
        else if ( aGeometricQuery->get_query_type() == Query_Type::INTERSECTION_LOCATION )
        {
            // this preprocessing can be streamlined a lot
            mActiveGeometryIndex                                                                   = aGeometricQuery->get_geomtric_index();
            Matrix< IndexMat > const&                                         tEdgeToVertex        = aGeometricQuery->get_query_entity_to_vertex_connectivity();
            moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > >* tQueryIndexedCoords  = aGeometricQuery->get_query_indexed_coordinates();
            Matrix< IndexMat >                                                tParentEntityIndices = aGeometricQuery->get_query_parent_entity_connectivity();

            // annoying copy until I rewrite the queue intersection
            Cell< Matrix< DDRMat > > tParentEntityCoords( tParentEntityIndices.numel() );
            Matrix< DDUMat >         tParentEntityIndiceUINT( tParentEntityIndices.numel() );
            tParentEntityCoords.reserve( tParentEntityIndices.numel() * 3 );
            for ( moris::uint i = 0; i < tParentEntityIndices.numel(); i++ )
            {
                tParentEntityCoords( i )     = *( *tQueryIndexedCoords )( tParentEntityIndices( i ) );
                tParentEntityIndiceUINT( i ) = (uint)tParentEntityIndices( i );
            }

            return queue_intersection( tEdgeToVertex( 0 ),
                tEdgeToVertex( 1 ),
                aGeometricQuery->get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 0 ) ),
                aGeometricQuery->get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 1 ) ),
                *( *tQueryIndexedCoords )( tEdgeToVertex( 0 ) ),
                *( *tQueryIndexedCoords )( tEdgeToVertex( 1 ) ),
                tParentEntityIndiceUINT,
                tParentEntityCoords );
        }
        else
        {
            MORIS_ERROR( 0, "In progress" );
            return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::is_intersected(
        const Matrix< IndexMat >& aNodeIndices,
        const Matrix< DDRMat >&   aNodeCoordinates )
    {
        // Check input
        // MORIS_ASSERT(aNodeIndices.length() == aNodeCoordinates.n_rows(),
        //         "Geometry engine must be provided the same number of node indices as node coordinates for "
        //         "determining if an element is intersected or not.");
        // MORIS_ASSERT(aNodeIndices.length() > 0,
        //         "Geometry engine must be provided at least 1 node to determine if an element is intersected or not.");

        bool tIsIntersected = false;

        switch ( mIntersectionMode )
        {
        case Intersection_Mode::LEVEL_SET:
        {
            // Initialize by evaluating the first node
            real tMin = mGeometries( mActiveGeometryIndex )->get_field_value( 0, aNodeCoordinates.get_row( 0 ) );
            real tMax = tMin;

            // Evaluate the rest of the nodes
            for ( uint tNodeCount = 1; tNodeCount < aNodeIndices.length(); tNodeCount++ )
            {
                real tEval = mGeometries( mActiveGeometryIndex )->get_field_value( tNodeCount, aNodeCoordinates.get_row( tNodeCount ) );

                tMin = std::min( tMin, tEval );
                tMax = std::max( tMax, tEval );
            }

            bool tHasZeroLsvInbetween = ( tMax >= mIsocontourThreshold ) and ( tMin <= mIsocontourThreshold );
            bool tHasCloseToZeroLsvFromAbove = ( std::abs( tMax - mIsocontourThreshold ) < mIsocontourTolerance );
            bool tHasCloseToZeroLsvFromBelow = ( std::abs( tMin - mIsocontourThreshold ) < mIsocontourTolerance );
            tIsIntersected = tHasZeroLsvInbetween or tHasCloseToZeroLsvFromAbove or tHasCloseToZeroLsvFromBelow;

            break;
        }
        case Intersection_Mode::COLORING:
        {
            real tFieldValue = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( 0 ), aNodeCoordinates.get_row( aNodeIndices( 0 ) ) );

            // Evaluate the rest of the nodes
            for ( uint Ik = 0; Ik < aNodeIndices.length(); Ik++ )
            {
                real tEval = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( Ik ), aNodeCoordinates.get_row( aNodeIndices( Ik ) ) );

                if ( tFieldValue != tEval )
                {
                    tIsIntersected = true;
                    break;
                }
            }

            break;
        }
        default:
        {
            MORIS_ERROR( false, "Geometry_Engine::is_intersected(), unknown intersection type." );
        }
        }

        // Return result
        return tIsIntersected;
    }

    bool
    Geometry_Engine::is_intersected(
        const Matrix< IndexMat >&                                         aNodeIndices,
        moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > >* aNodeCoordinates )
    {
        // Check input
        // MORIS_ASSERT(aNodeIndices.length() == aNodeCoordinates.n_rows(),
        //         "Geometry engine must be provided the same number of node indices as node coordinates for "
        //         "determining if an element is intersected or not.");
        // MORIS_ASSERT(aNodeIndices.length() > 0,
        //         "Geometry engine must be provided at least 1 node to determine if an element is intersected or not.");

        bool tIsIntersected = false;

        switch ( mIntersectionMode )
        {
        case Intersection_Mode::LEVEL_SET:
        {
            // Initialize by evaluating the first node
            real tMin = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( 0 ), *( *aNodeCoordinates )( aNodeIndices( 0 ) ) );
            real tMax = tMin;

            // Evaluate the rest of the nodes
            for ( uint tNodeCount = 0; tNodeCount < aNodeIndices.length(); tNodeCount++ )
            {
                real tEval = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( tNodeCount ), *( *aNodeCoordinates )( aNodeIndices( tNodeCount ) ) );

                tMin = std::min( tMin, tEval );
                tMax = std::max( tMax, tEval );
            }

            tIsIntersected = ( tMax >= mIsocontourThreshold and tMin <= mIsocontourThreshold ) or ( std::abs( tMax - mIsocontourThreshold ) < mIsocontourTolerance ) or ( std::abs( tMin - mIsocontourThreshold ) < mIsocontourTolerance );

            break;
        }
        case Intersection_Mode::COLORING:
        {
            real tFieldValue = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( 0 ), *( *aNodeCoordinates )( aNodeIndices( 0 ) ) );

            // Evaluate the rest of the nodes
            for ( uint Ik = 0; Ik < aNodeIndices.length(); Ik++ )
            {
                real tEval = mGeometries( mActiveGeometryIndex )->get_field_value( aNodeIndices( Ik ), *( *aNodeCoordinates )( aNodeIndices( Ik ) ) );

                if ( tFieldValue != tEval )
                {
                    tIsIntersected = true;
                    break;
                }
            }

            break;
        }
        default:
        {
            MORIS_ERROR( false, "Geometry_Engine::is_intersected(), unknown intersection type." );
        }
        }

        // Return result
        return tIsIntersected;
    }


    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queue_intersection(
        uint                            aFirstNodeIndex,
        uint                            aSecondNodeIndex,
        const Matrix< DDRMat >&         aFirstNodeLocalCoordinates,
        const Matrix< DDRMat >&         aSecondNodeLocalCoordinates,
        const Matrix< DDRMat >&         aFirstNodeGlobalCoordinates,
        const Matrix< DDRMat >&         aSecondNodeGlobalCoordinates,
        const Matrix< DDUMat >&         aBackgroundElementNodeIndices,
        const Cell< Matrix< DDRMat > >& aBackgroundElementNodeCoordinates )
    {
        // Queue an intersection node
        switch ( mIntersectionMode )
        {
        case Intersection_Mode::LEVEL_SET:
        {
            switch ( mGeometries( mActiveGeometryIndex )->get_intersection_interpolation() )
            {
            case Intersection_Interpolation::LINEAR:
            {
                mQueuedIntersectionNode = std::make_shared< Intersection_Node_Linear >(
                    mPDVHostManager.get_intersection_node( aFirstNodeIndex ),
                    mPDVHostManager.get_intersection_node( aSecondNodeIndex ),
                    aFirstNodeIndex,
                    aSecondNodeIndex,
                    aFirstNodeGlobalCoordinates,
                    aSecondNodeGlobalCoordinates,
                    mGeometries( mActiveGeometryIndex ),
                    mIsocontourThreshold,
                    mIsocontourTolerance,
                    mIntersectionTolerance );
                break;
            }
            case Intersection_Interpolation::MULTILINEAR:
            {
                if ( mNumSpatialDimensions == 2 )
                {
                    mQueuedIntersectionNode = std::make_shared< Intersection_Node_Bilinear >(
                        mPDVHostManager.get_intersection_node( aFirstNodeIndex ),
                        mPDVHostManager.get_intersection_node( aSecondNodeIndex ),
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        aFirstNodeLocalCoordinates,
                        aSecondNodeLocalCoordinates,
                        aBackgroundElementNodeIndices,
                        aBackgroundElementNodeCoordinates,
                        mGeometries( mActiveGeometryIndex ),
                        mIsocontourThreshold,
                        mIsocontourTolerance,
                        mIntersectionTolerance );
                }
                else
                {
                    MORIS_ERROR( false, "Only bilinear intersections have been implemented." );
                }

                break;
            }
            default:
            {
                MORIS_ERROR( false, "Intersection interpolation type not implemented yet." );
            }
            }
            break;
        }
        case Intersection_Mode::COLORING:
        {
            // Determine if edge is intersected
            if ( mGeometries( mActiveGeometryIndex )->get_field_value( aFirstNodeIndex, aFirstNodeGlobalCoordinates ) != mGeometries( mActiveGeometryIndex )->get_field_value( aSecondNodeIndex, aSecondNodeGlobalCoordinates ) )
            {
                mQueuedIntersectionNode = std::make_shared< Intersection_Node_Linear >(
                    mPDVHostManager.get_intersection_node( aFirstNodeIndex ),
                    mPDVHostManager.get_intersection_node( aSecondNodeIndex ),
                    aFirstNodeIndex,
                    aSecondNodeIndex,
                    aFirstNodeGlobalCoordinates,
                    aSecondNodeGlobalCoordinates,
                    mGeometries( mActiveGeometryIndex ),
                    mIsocontourThreshold,
                    mIsocontourTolerance );
            }
            else
            {
                return false;
            }

            break;
        }
        default:
        {
            MORIS_ERROR( false, "Geometry_Engine::queue_intersection(), unknown intersection type." );
        }
        }

        return mQueuedIntersectionNode->parent_edge_is_intersected();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queued_intersection_first_parent_on_interface()
    {
        return mQueuedIntersectionNode->first_parent_on_interface();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Geometry_Engine::queued_intersection_second_parent_on_interface()
    {
        return mQueuedIntersectionNode->second_parent_on_interface();
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

    void
    Geometry_Engine::admit_queued_intersection( uint aNodeIndex )
    {
        // Assign as PDV host
        if ( mGeometries( mActiveGeometryIndex )->depends_on_advs() )
        {
            mPDVHostManager.set_intersection_node( aNodeIndex, mQueuedIntersectionNode );
        }
        else
        {
            mPDVHostManager.set_intersection_node( aNodeIndex, nullptr );
        }

        // Assign as child node
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            // tGeomProx.set_geometric_proximity();
            mGeometries( tGeometryIndex )->add_child_node( aNodeIndex, mQueuedIntersectionNode );
        }

        // admit the queued intersection geometric proximity
        this->admit_queued_intersection_geometric_proximity( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::update_queued_intersection(
        const moris_index& aNodeIndex,
        const moris_index& aNodeId,
        const moris_index& aNodeOwner )
    {
        mPDVHostManager.update_intersection_node( aNodeIndex, aNodeId, aNodeOwner );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::create_new_child_nodes(
        const Cell< moris_index >*                   aNewNodeIndices,
        Cell< moris::mtk::Cell* >*                   aNewNodeParentCell,
        Cell< std::shared_ptr< Matrix< DDRMat > > >* aParamCoordRelativeToParent,
        Cell< Matrix< DDRMat > >*                    aNodeCoordinates )
    {
        // resize proximities
        mVertexGeometricProximity.resize( mVertexGeometricProximity.size() + aNewNodeIndices->size(), Geometric_Proximity( mGeometries.size() ) );

        // Loop over nodes
        for ( uint tNode = 0; tNode < aNewNodeIndices->size(); tNode++ )
        {
            std::shared_ptr< Child_Node > tChildNode                                          = std::make_shared< Child_Node >( ( *aNewNodeParentCell )( tNode ), ( *aParamCoordRelativeToParent )( tNode ).get() );
            mVertexGeometricProximity( ( *aNewNodeIndices )( tNode ) ).mAssociatedVertexIndex = ( *aNewNodeIndices )( tNode );

            Matrix< DDRMat > const& tCoord = ( *aNodeCoordinates )( tNode );

            // Assign to geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                mGeometries( tGeometryIndex )->add_child_node( ( *aNewNodeIndices )( tNode ), tChildNode );

                real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( ( *aNewNodeIndices )( tNode ), tCoord );

                moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                if ( std::abs( tVertGeomVal - mIsocontourThreshold ) < mIsocontourTolerance )
                {

                    tGeomProxIndex = 1;
                }

                mVertexGeometricProximity( ( *aNewNodeIndices )( tNode ) ).set_geometric_proximity( tGeomProxIndex, tGeometryIndex );
            }
        }

        // Set max node index
        if ( aNewNodeIndices->size() > 0 )
        {
            mPDVHostManager.set_num_background_nodes( ( *aNewNodeIndices )( aNewNodeIndices->size() - 1 ) + 1 );
        }
    }

    void
    Geometry_Engine::create_new_child_nodes(
        const Cell< moris_index >&      aNewNodeIndices,
        const Cell< xtk::Topology* >&   aParentTopo,
        const Cell< Matrix< DDRMat > >& aParamCoordRelativeToParent,
        const Matrix< DDRMat >&         aGlobalNodeCoord )
    {
        // resize proximities
        mVertexGeometricProximity.resize( mVertexGeometricProximity.size() + aNewNodeIndices.size(), Geometric_Proximity( mGeometries.size() ) );

        // Loop over nodes
        for ( uint tNode = 0; tNode < aNewNodeIndices.size(); tNode++ )
        {
            Matrix< DDUMat >         tParentNodeIndices( aParentTopo( tNode )->get_node_indices().length(), 1 );
            Cell< Matrix< DDRMat > > tParentNodeCoordinates( tParentNodeIndices.length() );
            for ( uint tParentNode = 0; tParentNode < tParentNodeIndices.length(); tParentNode++ )
            {
                tParentNodeIndices( tParentNode )     = aParentTopo( tNode )->get_node_indices()( tParentNode );
                tParentNodeCoordinates( tParentNode ) = aGlobalNodeCoord.get_row( tParentNodeIndices( tParentNode ) );
            }
            std::shared_ptr< Child_Node > tChildNode = std::make_shared< Child_Node >(
                tParentNodeIndices, tParentNodeCoordinates, aParentTopo( tNode )->get_basis_function(), aParamCoordRelativeToParent( tNode ) );

            mVertexGeometricProximity( aNewNodeIndices( tNode ) ).mAssociatedVertexIndex = aNewNodeIndices( tNode );

            Matrix< DDRMat > tCoord = aGlobalNodeCoord.get_row( aNewNodeIndices( tNode ) );

            // Assign to geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                mGeometries( tGeometryIndex )->add_child_node( aNewNodeIndices( tNode ), tChildNode );

                real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( aNewNodeIndices( tNode ), tCoord );

                moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                if ( std::abs( tVertGeomVal - mIsocontourThreshold ) < mIsocontourTolerance )
                {

                    tGeomProxIndex = 1;
                }

                mVertexGeometricProximity( aNewNodeIndices( tNode ) ).set_geometric_proximity( tGeomProxIndex, tGeometryIndex );
            }
        }

        // Set max node index
        if ( aNewNodeIndices.size() > 0 )
        {
            mPDVHostManager.set_num_background_nodes( aNewNodeIndices( aNewNodeIndices.size() - 1 ) + 1 );
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
        const Matrix< DDRMat >& aCoordinates )
    {
        // Initialize bitset of geometry signs
        Geometry_Bitset tGeometrySigns( 0 );

        // Flip bits as needed
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            moris_index tProxIndex = mVertexGeometricProximity( aNodeIndex ).get_geometric_proximity( (moris_index)tGeometryIndex );
            tGeometrySigns.set( tGeometryIndex, tProxIndex == 2 );
        }

        return mPhaseTable.get_phase_index( tGeometrySigns );
    }

    moris_index
    Geometry_Engine::is_interface_vertex( moris_index aNodeIndex,
        moris_index                                   aGeometryIndex )
    {
        moris_index tProxIndex = mVertexGeometricProximity( aNodeIndex ).get_geometric_proximity( (moris_index)aGeometryIndex );

        if ( tProxIndex == 1 )
        {
            return true;
        }

        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Geometry_Engine::get_elem_phase_index( Matrix< IndexMat > const& aElemOnOff )
    {
        // FIXME
        Geometry_Bitset tGeometrySigns( 0 );
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            tGeometrySigns.set( tGeometryIndex, aElemOnOff( tGeometryIndex ) );
        }

        return mPhaseTable.get_phase_index( tGeometrySigns );
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_node_phase_index_wrt_a_geometry(
        uint aNodeIndex,
        uint aGeometryIndex )
    {
        moris_index tProxIndex = mVertexGeometricProximity( aNodeIndex ).get_geometric_proximity( aGeometryIndex );

        size_t tPhaseOnOff = 0;

        if ( tProxIndex == 2 )
        {
            tPhaseOnOff = 1;
        }

        return tPhaseOnOff;
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Geometry_Engine::get_num_geometries()
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

    const Matrix< DDSMat >&
    Geometry_Engine::get_num_refinements( uint aFieldIndex )
    {
        return mGeometries( aFieldIndex )->get_num_refinements();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDSMat >&
    Geometry_Engine::get_refinement_mesh_indices( uint aFieldIndex )
    {
        return mGeometries( aFieldIndex )->get_refinement_mesh_indices();
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > >
    Geometry_Engine::get_mtk_fields()
    {
        Cell< std::shared_ptr< mtk::Field > > tFields( mGeometries.size() + mProperties.size() );
        std::copy( mGeometries.begin(), mGeometries.end(), tFields.begin() );
        std::copy( mProperties.begin(), mProperties.end(), tFields.begin() + mGeometries.size() );
        return tFields;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Geometry_Engine::get_field_value(
        uint                    aFieldIndex,
        uint                    aNodeIndex,
        const Matrix< DDRMat >& aCoordinates )
    {
        // TODO can return property field too
        return mGeometries( aFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
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
    Geometry_Engine::create_pdvs( mtk::Mesh_Pair aMeshPair )
    {
        // Tracer
        Tracer tTracer( "GEN", "Create PDVs" );

        // Get meshes using first mesh on mesh manager: Lagrange mesh with numbered aura (default)
        mtk::Integration_Mesh*   tIntegrationMesh   = aMeshPair.get_integration_mesh();
        mtk::Interpolation_Mesh* tInterpolationMesh = aMeshPair.get_interpolation_mesh();

        // Initialize PDV type groups and mesh set info from integration mesh
        Cell< Cell< Cell< PDV_Type > > > tPdvTypes( tIntegrationMesh->get_num_sets() );
        Cell< PDV_Type >                 tPDVTypeGroup( 1 );
        Cell< Matrix< DDUMat > >         tMeshSetIndicesPerProperty( mProperties.size() );

        tic tTimer;

        // Loop over properties to create PDVs
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            // PDV type and mesh set names/indices from parameter list
            tPDVTypeGroup( 0 ) = mProperties( tPropertyIndex )->get_pdv_type();

            tMeshSetIndicesPerProperty( tPropertyIndex ) = mProperties( tPropertyIndex )->get_pdv_mesh_set_indices();
            Cell< std::string > tMeshSetNames            = mProperties( tPropertyIndex )->get_pdv_mesh_set_names();

            // Convert mesh set names to indices
            uint tNumSetIndices = tMeshSetIndicesPerProperty( tPropertyIndex ).length();
            tMeshSetIndicesPerProperty( tPropertyIndex ).resize( tNumSetIndices + tMeshSetNames.size(), 1 );

            // number of mesh sets for current property
            uint tNumberOfSets = tMeshSetIndicesPerProperty( tPropertyIndex ).length();

            // Set for each property index the list of mesh set indices
            for ( uint tIndex = tNumSetIndices; tIndex < tNumberOfSets; tIndex++ )
            {
                tMeshSetIndicesPerProperty( tPropertyIndex )( tIndex ) = tIntegrationMesh->get_set_index_by_name( tMeshSetNames( tIndex - tNumSetIndices ) );
            }

            // Assign PDV types to each mesh set
            for ( uint tIndex = 0; tIndex < tNumberOfSets; tIndex++ )
            {
                uint tMeshSetIndex = tMeshSetIndicesPerProperty( tPropertyIndex )( tIndex );

                tPdvTypes( tMeshSetIndex ).push_back( tPDVTypeGroup );// why not use mProperties(tPropertyIndex)->get_pdv_type()
            }
        }

        real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( " Assign PDV type to mesh sets on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );

        // start timer
        tic tTimer2;
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

        tElapsedTime = tTimer2.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( " Get global-to-local-maps from mtk on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );

        // start timer
        tic tTimer3;
        // Create PDV hosts
        this->create_interpolation_pdv_hosts(
            tInterpolationMesh,
            tIntegrationMesh,
            tPdvTypes );

        tElapsedTime = tTimer3.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( " Create interpolation PDV hosts on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );

        // Set integration PDV types
        if ( mShapeSensitivities )
        {
            // start timer
            tic tTimer4;

            this->set_integration_pdv_types( tIntegrationMesh );

            tElapsedTime = tTimer4.toc< moris::chronos::milliseconds >().wall;
            MORIS_LOG_INFO( " Set integration PDV types on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );
        }

        // start timer
        tic tTimer5;

        // Loop over properties to assign PDVs
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            // Assign PDVs
            if ( mProperties( tPropertyIndex )->is_interpolation_pdv() )
            {
                // PDV type and mesh set names/indices from parameter list
                tPDVTypeGroup( 0 ) = mProperties( tPropertyIndex )->get_pdv_type();

                mProperties( tPropertyIndex )->add_nodal_data( tInterpolationMesh );

                this->assign_property_to_pdv_hosts(
                    mProperties( tPropertyIndex ),
                    tPDVTypeGroup( 0 ),
                    tIntegrationMesh,
                    tMeshSetIndicesPerProperty( tPropertyIndex ) );
            }
            else
            {
                MORIS_ERROR( false, "Assignment of PDVs is only supported with an interpolation mesh right now." );
            }
        }

        tElapsedTime = tTimer5.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( " Assign property to pdv host on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );

        // start timer
        tic tTimer6;
        // Create PDV IDs
        mPDVHostManager.create_pdv_ids();

        tElapsedTime = tTimer6.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( " Compute and comunicate PDV IDs on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );
    }

    void
    Geometry_Engine::print_gen_vertices(
        std::string aFile,
        mtk::Mesh*  aMesh )
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

        for ( moris::uint iHeader = 0; iHeader < this->get_num_geometries(); iHeader++ )
        {
            tStringStream << "gval_" + std::to_string( iHeader ) << ",";
        }
        for ( moris::uint iHeader = 0; iHeader < this->get_num_geometries(); iHeader++ )
        {
            tStringStream << "gprox_" + std::to_string( iHeader );
            if ( iHeader != this->get_num_geometries() - 1 )
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
            for ( moris::uint iGeom = 0; iGeom < this->get_num_geometries(); iGeom++ )
            {
                tStringStream << mGeometries( iGeom )->get_field_value( tVertex.get_index(), tCoords ) << ",";
            }
            for ( moris::uint iGeom = 0; iGeom < this->get_num_geometries(); iGeom++ )
            {
                // 0 - G(x) < threshold
                // 1 - G(x) == threshold
                // 2 - G(x) > threshold

                moris_index tGeomProx = mVertexGeometricProximity( tVertex.get_index() ).get_geometric_proximity( iGeom );

                if ( tGeomProx == 0 )
                {
                    tStringStream << "-";
                }
                else if ( tGeomProx == 1 )
                {
                    tStringStream << "=";
                }
                else
                {
                    tStringStream << "+";
                }
                if ( iGeom != this->get_num_geometries() - 1 )
                {

                    tStringStream << ",";
                }
            }
            tStringStream << std::endl;
        }
        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }


    void
    Geometry_Engine::setup_diagnostics(
        bool               aDiagnostics,
        std::string const& aDiagnosticPath,
        std::string const& aDiagnosticLabel )
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

    std::string
    Geometry_Engine::get_diagnostic_file_name( std::string const& aLabel ) const
    {
        MORIS_ASSERT( mDiagnostics, "Only callable with diagnostics on" );
        return mDiagnosticPath + "/id_" + mDiagnosticId + "_p_" + std::to_string( moris::par_rank() ) + "_" + aLabel + ".diag";
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::initialize_pdv_type_list()
    {
        // Reserve of temporary pdv type list
        Cell< enum PDV_Type > tTemporaryPdvTypeList;

        tTemporaryPdvTypeList.reserve( static_cast< int >( PDV_Type::UNDEFINED ) + 1 );

        Matrix< DDUMat > tListToCheckIfEnumExist( ( static_cast< int >( PDV_Type::UNDEFINED ) + 1 ), 1, 0 );

        // PDV type map
        map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();

        // Loop over properties to build parallel consistent pdv list
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            // PDV type and mesh set names/indices from parameter list
            PDV_Type tPdvType = mProperties( tPropertyIndex )->get_pdv_type();

            if ( tListToCheckIfEnumExist( static_cast< int >( tPdvType ), 0 ) == 0 )
            {
                // Set 1 at position of the enum value
                tListToCheckIfEnumExist( static_cast< int >( tPdvType ), 0 ) = 1;

                tTemporaryPdvTypeList.push_back( tPdvType );
            }
        }

        // Shrink pdv type list to fit
        tTemporaryPdvTypeList.shrink_to_fit();

        // Communicate dof types so that all processors have the same unique list
        mPDVHostManager.communicate_dof_types( tTemporaryPdvTypeList );

        // Create a map
        mPDVHostManager.create_dv_type_map();
    }

    //--------------------------------------------------------------------------------------------------------------


        void Geometry_Engine::distribute_advs(
                mtk::Mesh_Pair aMeshPair,
                moris::Cell<std::shared_ptr< mtk::Field > > aFields )
        {
            // Tracer
            Tracer tTracer("GEN", "Distribute ADVs");

            // Gather all fields
            Cell<std::shared_ptr<Field>> tFields(mGeometries.size() + mProperties.size());
            std::copy(mGeometries.begin(), mGeometries.end(), tFields.begin());
            std::copy(mProperties.begin(), mProperties.end(), tFields.begin() + mGeometries.size());

            // Interpolation mesh
            mtk::Interpolation_Mesh* tMesh = aMeshPair.get_interpolation_mesh();

            //------------------------------------//
            // Determine owned and shared ADV IDs //
            //------------------------------------//
            clock_t tStart_Owned_Shared_ADVs = clock();

            // Set primitive IDs
            Matrix<DDSMat> tPrimitiveADVIds(mInitialPrimitiveADVs.length(), 1);
            for (uint tADVIndex = 0; tADVIndex < mInitialPrimitiveADVs.length(); tADVIndex++)
            {
                tPrimitiveADVIds(tADVIndex) = tADVIndex;
            }

            // Start with primitive IDs for owned IDs on processor 0
            Matrix<DDSMat> tOwnedADVIds(0, 0);
            if (par_rank() == 0)
            {
                tOwnedADVIds = tPrimitiveADVIds;
            }

            // Owned and shared ADVs per field
            Cell<Matrix<DDUMat>> tSharedCoefficientIndices(tFields.size());
            Cell<Matrix<DDSMat>> tSharedADVIds(tFields.size());
            Matrix<DDUMat> tAllOffsetIDs(tFields.size(), 1);

            // Get all node indices from the mesh (for now)
            uint tNumNodes = tMesh->get_num_nodes();
            Matrix<DDUMat> tNodeIndices(tNumNodes, 1);
            for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
            {
                tNodeIndices(tNodeIndex) = tNodeIndex;
            }

            moris::Cell< uint > tNumCoeff(tFields.size());
            // Loop over all geometries to get number of new ADVs
            sint tOffsetID = tPrimitiveADVIds.length();
            for (uint tFieldIndex = 0; tFieldIndex < tFields.size(); tFieldIndex++)
            {
                // Determine if level set will be created
                if (tFields(tFieldIndex)->intended_discretization())
                {
                    // Get discretization mesh index
                    uint tDiscretizationMeshIndex = tFields(tFieldIndex)->get_discretization_mesh_index();

                    uint tMaxNumberOfCoefficients = tMesh->get_max_num_coeffs_on_proc(tDiscretizationMeshIndex);

                    Matrix< IdMat > tAllCoefIds(tMaxNumberOfCoefficients,1,gNoID);
                    Matrix< IndexMat > tAllCoefIndices(tMaxNumberOfCoefficients,1,gNoIndex);
                    Matrix< IdMat > tAllCoefOwners(tMaxNumberOfCoefficients,1,gNoID);

                    for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
                    {
                        // check whether node has an underlying discretization on this processor
                        bool tNodeHasDiscretization = tMesh->
                                get_mtk_vertex(tNodeIndex).has_interpolation(tDiscretizationMeshIndex);

                        // process only nodes that have discretization
                        if ( tNodeHasDiscretization )
                        {
                            // get indices and IDs from mtk mesh - FIXME: should return const &
                            const Matrix< IndexMat > tCoefIndices = tMesh->get_coefficient_indices_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex);

                            const Matrix< IdMat > tCoefIds = tMesh->get_coefficient_IDs_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex);

                            const Matrix< IdMat > tCoefOwners = tMesh->get_coefficient_owners_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex);

                            // check that number of indices and ids are the same
                            MORIS_ASSERT( tCoefIds.numel() == tCoefIndices.numel(),
                                    "distribute_advs - numbers of coefficients and ids do not match.\n");

                            // get number of coefficients for current node
                            uint tNumCoefOfNode=tCoefIds.numel();

                            for (uint tCoefIndex=0;tCoefIndex<tNumCoefOfNode;++tCoefIndex)
                            {
                                // get coefficient index
                                moris_index tCurrentIndex = tCoefIndices(tCoefIndex);

                                // check whether mesh coefficient has already been set
                                if ( tAllCoefIds(tCurrentIndex) == -1 )
                                {
                                    // increase field coefficient count
                                    tNumCoeff( tFieldIndex )++;

                                    // populate mesh index to mesh coefficient id map
                                    tAllCoefIds(tCurrentIndex) = tCoefIds(tCoefIndex);

                                    tAllCoefOwners(tCurrentIndex) = tCoefOwners(tCoefIndex);
                                }
                                else
                                {
                                    // check for consistency
                                    MORIS_ASSERT( tAllCoefIds(tCurrentIndex) == tCoefIds(tCoefIndex),
                                            "distribute_advs - inconsistent index and ids.\n");
                                }
                            }
                        }
                    }

                    if(par_size() > 1)
                    {
                        this->communicate_missing_owned_coefficients(
                                aMeshPair,
                                tAllCoefIds,
                                tAllCoefOwners,
                                tNumCoeff,
                                tFieldIndex,
                                tDiscretizationMeshIndex);
                    }

                    uint tOwnedCounter = 0;
                    uint tSharedCounter = 0;

                    for( uint Ik = 0; Ik< tAllCoefIds.numel(); Ik++ )
                    {
                        if( tAllCoefIds( Ik ) != gNoID && tAllCoefOwners( Ik ) == par_rank() )
                        {
                            tOwnedCounter++;
                        }
                        else if( tAllCoefIds( Ik ) != gNoID )
                        {
                            tSharedCounter++;
                        }
                    }

                    Matrix<DDUMat> tOwnedCoefficients(tOwnedCounter,1);
                    Matrix<DDUMat> tSharedCoefficients( tSharedCounter,1);

                    tOwnedCounter = 0;
                    tSharedCounter = 0;

                    for( uint Ik = 0; Ik< tAllCoefIds.numel(); Ik++ )
                    {
                        if( tAllCoefIds( Ik ) != gNoID && tAllCoefOwners( Ik ) == par_rank() )
                        {
                            tOwnedCoefficients(tOwnedCounter++) = Ik;
                        }
                        else if( tAllCoefIds( Ik ) != gNoID )
                        {
                            tSharedCoefficients( tSharedCounter++)= Ik;
                        }
                    }

                    //                    // Get owned coefficients
                    //                    Matrix<DDUMat> tOwnedCoefficients = tMesh->get_owned_discretization_coefficient_indices(
                    //                            tNodeIndices,
                    //                            tDiscretizationMeshIndex);
                    //
                    //                    // Get shared coefficients
                    //                    Matrix<DDUMat> tSharedCoefficients = tMesh->get_shared_discretization_coefficient_indices(
                    //                            tNodeIndices,
                    //                            tDiscretizationMeshIndex);

                    // Sizes of ID vectors
                    uint tNumOwnedADVs = tOwnedADVIds.length();
                    uint tNumOwnedCoefficients = tOwnedCoefficients.numel();
                    uint tNumSharedCoefficients = tSharedCoefficients.numel();

                    // Resize ID lists and bounds
                    tOwnedADVIds.resize(tNumOwnedADVs + tNumOwnedCoefficients, 1);
                    mLowerBounds.resize(tNumOwnedADVs + tNumOwnedCoefficients, 1);
                    mUpperBounds.resize(tNumOwnedADVs + tNumOwnedCoefficients, 1);
                    tSharedADVIds(tFieldIndex).resize(tNumOwnedCoefficients + tNumSharedCoefficients, 1);
                    tSharedCoefficientIndices(tFieldIndex) = tOwnedCoefficients;

                    // Add owned coefficients to lists
                    for (uint tOwnedCoefficient = 0; tOwnedCoefficient < tNumOwnedCoefficients; tOwnedCoefficient++)
                    {
                        sint tADVId = tOffsetID + tMesh->get_glb_entity_id_from_entity_loc_index(
                                tOwnedCoefficients(tOwnedCoefficient),
                                EntityRank::BSPLINE,
                                tDiscretizationMeshIndex);

                        MORIS_ASSERT( tADVId-tOffsetID == tAllCoefIds( tOwnedCoefficients(tOwnedCoefficient) ),"check if this is a problem");

                        tOwnedADVIds(tNumOwnedADVs + tOwnedCoefficient) = tADVId;
                        mLowerBounds(tNumOwnedADVs + tOwnedCoefficient) = tFields(tFieldIndex)->get_discretization_lower_bound();
                        mUpperBounds(tNumOwnedADVs + tOwnedCoefficient) = tFields(tFieldIndex)->get_discretization_upper_bound();;
                        tSharedADVIds(tFieldIndex)(tOwnedCoefficient) = tADVId;
                    }

                    // Add shared coefficients to field-specific list
                    tSharedCoefficientIndices(tFieldIndex).resize(tNumOwnedCoefficients + tNumSharedCoefficients, 1);
                    for (uint tSharedCoefficient = 0; tSharedCoefficient < tNumSharedCoefficients; tSharedCoefficient++)
                    {
                        sint tADVId = tOffsetID + tMesh->get_glb_entity_id_from_entity_loc_index(
                                tSharedCoefficients(tSharedCoefficient),
                                EntityRank::BSPLINE,
                                tDiscretizationMeshIndex);

                        MORIS_ASSERT( tADVId-tOffsetID == tAllCoefIds( tSharedCoefficients(tSharedCoefficient) ),"check if this is a problem");

                        tSharedCoefficientIndices(tFieldIndex)(tNumOwnedCoefficients + tSharedCoefficient) = tSharedCoefficients(tSharedCoefficient);
                        tSharedADVIds(tFieldIndex)(tNumOwnedCoefficients + tSharedCoefficient) = tADVId;
                    }

                    // Update offset based on maximum ID
                    tAllOffsetIDs(tFieldIndex) = tOffsetID;
                    tOffsetID += tMesh->get_max_entity_id(EntityRank::BSPLINE, tDiscretizationMeshIndex);
                }
            }

            // Set owned ADV IDs
            mPDVHostManager.set_owned_adv_ids(tOwnedADVIds);

            MORIS_LOG_INFO("Time to collect owned and shared ADVs: %f sec",(moris::real) ( clock() - tStart_Owned_Shared_ADVs ) / CLOCKS_PER_SEC);

            //----------------------------------------//
            // Create owned ADV vector                //
            //----------------------------------------//
            clock_t tStart_Create_Owned_ADVs = clock();

            // Create factory for distributed ADV vector
            sol::Matrix_Vector_Factory tDistributedFactory;

            // Create map for distributed vectors
            sol::Dist_Map* tOwnedADVMap = tDistributedFactory.create_map(tOwnedADVIds);
            sol::Dist_Map* tPrimitiveADVMap = tDistributedFactory.create_map(tPrimitiveADVIds);

            // Create vectors
            sol::Dist_Vector* tNewOwnedADVs = tDistributedFactory.create_vector(tOwnedADVMap, 1, false, true);
            mPrimitiveADVs = tDistributedFactory.create_vector(tPrimitiveADVMap, 1, false, true);

            // Assign primitive ADVs
            if (par_rank() == 0)
            {
                tNewOwnedADVs->replace_global_values(tPrimitiveADVIds, mInitialPrimitiveADVs);
            }

            // Global assembly
            tNewOwnedADVs->vector_global_assembly();

            // Get primitive ADVs from owned vector
            mPrimitiveADVs->import_local_to_global(*tNewOwnedADVs);

            // Set field ADVs using distributed vector
            if (mInitialPrimitiveADVs.length() > 0)
            {
                for (uint tFieldIndex = 0; tFieldIndex < tFields.size(); tFieldIndex++)
                {
                    tFields(tFieldIndex)->set_advs(mPrimitiveADVs);
                }
            }

            MORIS_LOG_INFO("Time to create owned ADVs: %f sec",(moris::real) ( clock() - tStart_Create_Owned_ADVs ) / CLOCKS_PER_SEC);

            //----------------------------------------//
            // Convert to B-spline fields             //
            //----------------------------------------//

            clock_t tStart_Convert_to_Bspline_Fields = clock();

            //FIXME this hole section is super hacky and limiting. has to be rewritten from scratch.
            moris::map< std::string, uint > tFieldNameToIndexMap;
            for( uint Ik =0; Ik<aFields.size(); Ik++ )
            {
                std::string tLabel = aFields( Ik )->get_label();
                tFieldNameToIndexMap[ tLabel ] = Ik;
            }

            // Loop to find B-spline geometries
            for (uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                // Shape sensitivities logic
                mShapeSensitivities = (mShapeSensitivities or mGeometries(tGeometryIndex)->depends_on_advs());

                // Convert to B-spline field
                if (mGeometries(tGeometryIndex)->intended_discretization())
                {
                    // Always have shape sensitivities if B-spline field
                    mShapeSensitivities = true;

                    std::string tGeoName = mGeometries(tGeometryIndex)->get_name();

                    if( not tFieldNameToIndexMap.key_exists(tGeoName) )
                    {
                        // Create B-spline geometry FIXME Overwriting the given geometry is obviously wrong
                        mGeometries(tGeometryIndex) = std::make_shared<BSpline_Geometry>(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices(tGeometryIndex),
                                tSharedADVIds(tGeometryIndex),
                                tAllOffsetIDs(tGeometryIndex),
                                aMeshPair,
                                mGeometries(tGeometryIndex));
                    }
                    else
                    {
                        uint tMTKFieldIndex = tFieldNameToIndexMap.find( tGeoName );

                        // Create B-spline geometry
                        mGeometries(tGeometryIndex) = std::make_shared<BSpline_Geometry>(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices(tGeometryIndex),
                                tSharedADVIds(tGeometryIndex),
                                tAllOffsetIDs(tGeometryIndex),
                                aMeshPair,
                                mGeometries(tGeometryIndex),
                                aFields(tMTKFieldIndex));
                    }
                }
                // Store field values if needed. FIXME this is obviously wrong that GEN sets it's own mesh
                else if (mGeometries(tGeometryIndex)->intended_storage())
                {
                    // Create stored geometry FIXME this stored geometry stuff is kind of hacky
                    mGeometries(tGeometryIndex) = std::make_shared<Stored_Geometry>(
                            tMesh,
                            mGeometries(tGeometryIndex));

                    mGeometries(tGeometryIndex)->unlock_field();
                    mGeometries(tGeometryIndex)->set_mesh_pair( aMeshPair );
                }
                else
                {
                    // Every Field needs a mesh. FIXME setting the mesh here is to late
                    mGeometries(tGeometryIndex)->unlock_field();
                    mGeometries(tGeometryIndex)->set_mesh_pair( aMeshPair );
                    mGeometries(tGeometryIndex)->set_num_original_nodes(aMeshPair.get_interpolation_mesh()->get_num_nodes() );
                }
            }

            // Loop to find B-spline properties
            for (uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++)
            {
                // Convert to B-spline field
                if (mProperties(tPropertyIndex)->intended_discretization())
                {
                    // Create B-spline property
                    mProperties(tPropertyIndex) = std::make_shared<BSpline_Property>(
                            tNewOwnedADVs,
                            tSharedCoefficientIndices(mGeometries.size() + tPropertyIndex),
                            tSharedADVIds(mGeometries.size() + tPropertyIndex),
                            tAllOffsetIDs(mGeometries.size() + tPropertyIndex),
                            aMeshPair,
                            mProperties(tPropertyIndex));
                }
            }

            // Update dependencies
            std::copy(mGeometries.begin(), mGeometries.end(), tFields.begin());
            std::copy(mProperties.begin(), mProperties.end(), tFields.begin() + mGeometries.size());
            for (uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++)
            {
                mProperties(tPropertyIndex)->update_dependencies(tFields);
            }

            // Save new owned ADVs
            mOwnedADVs = tNewOwnedADVs;

            MORIS_LOG_INFO("Time to convert to Bspline fields: %f sec",(moris::real) ( clock() - tStart_Convert_to_Bspline_Fields ) / CLOCKS_PER_SEC);

            //----------------------------------------//
            // Communicate all ADV IDs to processor 0 //
            //----------------------------------------//

            clock_t tStart_Communicate_ADV_IDs = clock();

            // Sending mats
            Cell<Matrix<DDSMat>> tSendingIDs(0);
            Cell<Matrix<DDRMat>> tSendingLowerBounds(0);
            Cell<Matrix<DDRMat>> tSendingUpperBounds(0);

            // Receiving mats
            Cell<Matrix<DDSMat>> tReceivingIDs(0);
            Cell<Matrix<DDRMat>> tReceivingLowerBounds(0);
            Cell<Matrix<DDRMat>> tReceivingUpperBounds(0);

            // Set up communication list for communicating ADV IDs
            Matrix<IdMat> tCommunicationList(1, 1, 0);
            if (par_rank() == 0)
            {
                // Resize communication list and sending mats
                tCommunicationList.resize(par_size() - 1, 1);
                tSendingIDs.resize(par_size() - 1);
                tSendingLowerBounds.resize(par_size() - 1);
                tSendingUpperBounds.resize(par_size() - 1);

                // Assign communication list
                for (uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++)
                {
                    tCommunicationList(tProcessorIndex - 1) = tProcessorIndex;
                }
            }
            else
            {
                tSendingIDs = {tOwnedADVIds};
                tSendingLowerBounds = {mLowerBounds};
                tSendingUpperBounds = {mUpperBounds};
            }

            // Communicate mats
            communicate_mats(tCommunicationList, tSendingIDs, tReceivingIDs);
            communicate_mats(tCommunicationList, tSendingLowerBounds, tReceivingLowerBounds);
            communicate_mats(tCommunicationList, tSendingUpperBounds, tReceivingUpperBounds);

            MORIS_LOG_INFO("Time to communicate ADV IDs: %f sec",(moris::real) ( clock() - tStart_Communicate_ADV_IDs ) / CLOCKS_PER_SEC);

            // Assemble full ADVs/bounds
            clock_t tStart_ADV_Bounds = clock();

            if (par_rank() == 0)
            {
                // Start full IDs with owned IDs on processor 0
                mFullADVIds = tOwnedADVIds;

                // Assemble additional IDs/bounds from other processors
                for (uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++)
                {
                    // Get number of received ADVs
                    uint tFullADVsFilled = mFullADVIds.length();
                    uint tNumReceivedADVs = tReceivingIDs(tProcessorIndex - 1).length();

                    // Resize full ADV IDs and bounds
                    mFullADVIds.resize(tFullADVsFilled + tNumReceivedADVs, 1);
                    mLowerBounds.resize(tFullADVsFilled + tNumReceivedADVs, 1);
                    mUpperBounds.resize(tFullADVsFilled + tNumReceivedADVs, 1);

                    // Assign received ADV IDs
                    for (uint tADVIndex = 0; tADVIndex < tNumReceivedADVs; tADVIndex++)
                    {
                        mFullADVIds(tFullADVsFilled + tADVIndex) = tReceivingIDs(tProcessorIndex - 1)(tADVIndex);
                        mLowerBounds(tFullADVsFilled + tADVIndex) =
                                tReceivingLowerBounds(tProcessorIndex - 1)(tADVIndex);
                        mUpperBounds(tFullADVsFilled + tADVIndex) =
                                tReceivingUpperBounds(tProcessorIndex - 1)(tADVIndex);
                    }
                }
            }
            else
            {
                mLowerBounds.set_size(0, 0);
                mUpperBounds.set_size(0, 0);
            }

            MORIS_LOG_INFO("Time to assemble ADVs and bounds on Proc 0: %f sec",(moris::real) ( clock() - tStart_ADV_Bounds ) / CLOCKS_PER_SEC);

            // Reset mesh information
            clock_t tStart_Reset_Mesh_Info = clock();

            this->reset_mesh_information(tMesh);

            MORIS_LOG_INFO("Time to reset mesh information: %f sec",(moris::real) ( clock() - tStart_Reset_Mesh_Info ) / CLOCKS_PER_SEC);
        }
    void
    Geometry_Engine::reset_mesh_information( mtk::Interpolation_Mesh* aMesh )
    {
        // Register spatial dimension
        mNumSpatialDimensions = aMesh->get_spatial_dim();

        // Reset PDV host manager
        mPDVHostManager.reset();
        mPDVHostManager.set_num_background_nodes( aMesh->get_num_nodes() );

        // Reset info related to the mesh
        mActiveGeometryIndex = 0;
        for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
        {
            mGeometries( tGeometryIndex )->reset_nodal_data();
        }
        for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
        {
            mProperties( tPropertyIndex )->reset_nodal_data();
        }

        // Allocate proximity data
        this->setup_initial_geometric_proximities( aMesh );
    }


    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::communicate_missing_owned_coefficients(
        mtk::Mesh_Pair&  aMeshPair,
        Matrix< IdMat >& aAllCoefIds,
        Matrix< IdMat >& aAllCoefOwners,
        Cell< uint >&    aNumCoeff,
        uint             aFieldIndex,
        uint             aDiscretizationMeshIndex )
    {

        Matrix< IdMat > tCommTable = aMeshPair.get_interpolation_mesh()->get_communication_table();

        // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
        Matrix< DDSMat > tCommTableMap( tCommTable.max() + 1, 1, -1 );

        moris::uint tNumCommProcs = tCommTable.numel();

        // Loop over communication table to fill the communication table map
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            tCommTableMap( tCommTable( Ik ) ) = Ik;
        }

        moris::Cell< Matrix< IdMat > > tSharedCoeffsPosGlobal( tNumCommProcs );

        // Set Mat to store number of shared coeffs per processor
        Matrix< DDUMat > tNumSharedCoeffsPerProc( tNumCommProcs, 1, 0 );

        // Count number of coeffs per proc which have to be communicated
        for ( moris::uint Ib = 0; Ib < aAllCoefIds.numel(); Ib++ )
        {
            // Check if coeffs at this position is not NULL
            if ( aAllCoefIds( Ib ) != gNoID && aAllCoefOwners( Ib ) != par_rank() )
            {

                // get owning processor
                moris::moris_id tProcID = aAllCoefOwners( Ib );

                moris::sint tProcIdPos = tCommTableMap( tProcID );

                MORIS_ASSERT( tProcIdPos != gNoID,
                    "communicate_missing_owned_coefficients: Map returns proc rank -1. Check communication table" );

                // Add +1 to the processor number of shared coeffs per processor
                tNumSharedCoeffsPerProc( tProcIdPos )++;
            }
        }


        // Set size of the moris::Mats in the Cell
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            if ( tNumSharedCoeffsPerProc( Ik ) != 0 )
            {
                tSharedCoeffsPosGlobal( Ik ).set_size( tNumSharedCoeffsPerProc( Ik ), 1 );
            }
        }

        // Temporary Mat to add external coeffs ids at the next spot in the matrix which will be communicated
        Matrix< DDUMat > tShredCoeffPosPerProc( tNumCommProcs, 1, 0 );

        // Loop over coeffs per type
        for ( moris::uint Ia = 0; Ia < aAllCoefIds.numel(); Ia++ )
        {
            // Check if coeffs at this position is not NULL
            if ( aAllCoefIds( Ia ) != gNoID && aAllCoefOwners( Ia ) != par_rank() )
            {
                // Get owning processor
                moris::uint tProcID = aAllCoefOwners( Ia );

                moris::sint tProcIdPos = tCommTableMap( tProcID );

                // Add owning processor id to moris::Mat
                tSharedCoeffsPosGlobal( tProcIdPos )( tShredCoeffPosPerProc( tProcIdPos ) ) =
                    aAllCoefIds( Ia );

                tShredCoeffPosPerProc( tProcIdPos )++;
            }
        }

        // receiving list
        moris::Cell< Matrix< IdMat > > tMatsToReceive;

        barrier();

        // Communicate position of shared adofs to the owning processor
        communicate_mats(
            tCommTable,
            tSharedCoeffsPosGlobal,
            tMatsToReceive );

        map< moris_id, moris_index > tCoeffGlobaltoLocalMap;
        aMeshPair.get_interpolation_mesh()->get_adof_map(
            aDiscretizationMeshIndex,
            tCoeffGlobaltoLocalMap );

        // Loop over all Mats set dummy owned coeffs
        for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
        {
            for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
            {
                // Get owned coeff Index
                moris_id    tID            = tMatsToReceive( Ik )( Ii );
                moris_index tLocalCoeffInd = tCoeffGlobaltoLocalMap.find( tID );

                if ( aAllCoefIds( tLocalCoeffInd ) == gNoID )
                {
                    aAllCoefIds( tLocalCoeffInd )    = tID;
                    aAllCoefOwners( tLocalCoeffInd ) = par_rank();

                    aNumCoeff( aFieldIndex )++;
                }

                MORIS_ASSERT( aAllCoefIds( tLocalCoeffInd ) == tID,
                    "Field_Discrete::communicate_missing_owned_coefficients( ), coefficient IDs are nor parallel consistent" );
            }
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

        if ( aExodusFileName != "" )
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
            uint                tNumGeometries = mGeometries.size();
            uint                tNumProperties = mProperties.size();
            Cell< std::string > tFieldNames( tNumGeometries + tNumProperties );

            // Geometry field names
            for ( uint tGeometryIndex = 0; tGeometryIndex < tNumGeometries; tGeometryIndex++ )
            {
                tFieldNames( tGeometryIndex ) = mGeometries( tGeometryIndex )->get_name();
                if ( tFieldNames( tGeometryIndex ) == "" )
                {
                    tFieldNames( tGeometryIndex ) = "Geometry " + std::to_string( tGeometryIndex );
                }
            }

            // Property field names
            for ( uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++ )
            {
                tFieldNames( tNumGeometries + tPropertyIndex ) = mProperties( tPropertyIndex )->get_name();
                if ( tFieldNames( tNumGeometries + tPropertyIndex ) == "" )
                {
                    tFieldNames( tNumGeometries + tPropertyIndex ) = "Property " + std::to_string( tPropertyIndex );
                }
            }

            // write time to file
            tWriter.set_time( tTimeShift );

            // Set nodal fields based on field names
            tWriter.set_nodal_fields( tFieldNames );

            // Get all node coordinates
            Cell< Matrix< DDRMat > > tNodeCoordinates( aMesh->get_num_nodes() );
            for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
            {
                tNodeCoordinates( tNodeIndex ) = aMesh->get_node_coordinate( tNodeIndex );
            }

            // Loop over geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < tNumGeometries; tGeometryIndex++ )
            {
                // Create field vector
                Matrix< DDRMat > tFieldData( aMesh->get_num_nodes(), 1 );

                // Assign field to vector
                for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                {
                    tFieldData( tNodeIndex ) = mGeometries( tGeometryIndex )->get_field_value( tNodeIndex, tNodeCoordinates( tNodeIndex ) );
                }

                // Create field on mesh
                tWriter.write_nodal_field( tFieldNames( tGeometryIndex ), tFieldData );
            }

            // Loop over properties
            for ( uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++ )
            {
                // Create field vector
                Matrix< DDRMat > tFieldData( aMesh->get_num_nodes(), 1 );

                // Assign field to vector
                for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                {
                    tFieldData( tNodeIndex ) = mProperties( tPropertyIndex )->get_field_value( tNodeIndex, tNodeCoordinates( tNodeIndex ) );
                }

                // Create field on mesh
                tWriter.write_nodal_field( tFieldNames( tNumGeometries + tPropertyIndex ), tFieldData );
            }

            // Finalize
            tWriter.close_file( true );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::write_geometry_fields(
        mtk::Mesh*  aMesh,
        std::string aBaseFileName )
    {
        if ( aBaseFileName != "" )
        {
            // Get all node coordinates
            Cell< Matrix< DDRMat > > tNodeCoordinates( aMesh->get_num_nodes() );
            for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
            {
                tNodeCoordinates( tNodeIndex ) = aMesh->get_node_coordinate( tNodeIndex );
            }

            // Loop over geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                // Create file
                std::ofstream tOutFile( aBaseFileName + "_" + std::to_string( tGeometryIndex ) + ".txt" );

                // Write to file
                for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                {
                    // Coordinates
                    for ( uint tDimension = 0; tDimension < mNumSpatialDimensions; tDimension++ )
                    {
                        tOutFile << tNodeCoordinates( tNodeIndex )( tDimension ) << ", ";
                    }

                    // Fill unused dimensions with zeros
                    for ( uint tDimension = mNumSpatialDimensions; tDimension < 3; tDimension++ )
                    {
                        tOutFile << 0.0 << ", ";
                    }

                    // Level-set field
                    tOutFile << mGeometries( tGeometryIndex )->get_field_value( tNodeIndex, tNodeCoordinates( tNodeIndex ) ) << std::endl;
                }

                // Close file
                tOutFile.close();
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    // PRIVATE
    //--------------------------------------------------------------------------------------------------------------
    void
    Geometry_Engine::create_interpolation_pdv_hosts(
        mtk::Interpolation_Mesh*         aInterpolationMesh,
        mtk::Integration_Mesh*           aIntegrationMesh,
        Cell< Cell< Cell< PDV_Type > > > aPdvTypes )
    {
        // Get information from integration mesh
        uint tNumSets = aPdvTypes.size();

        Cell< Matrix< DDSMat > > tNodeIndicesPerSet( tNumSets );
        Cell< Matrix< DDSMat > > tNodeIdsPerSet( tNumSets );
        Cell< Matrix< DDSMat > > tNodeOwnersPerSet( tNumSets );
        Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( tNumSets );

        Cell< uint > tCounter( tNumSets, 0 );

        uint tNumSpatialDim = 0;

        // Determine if we need to do the below loop
        bool tDoJankLoop = false;
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            if ( aPdvTypes( tMeshSetIndex ).size() > 0 )
            {
                tDoJankLoop = true;
            }
        }

        // Determine the indices, IDs, and ownerships of base cell nodes in each set of integration mesh
        // List has redundant entries
        if ( tDoJankLoop )
        {
            // Loop through sets in integration mesh
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
            {
                mtk::Set* tSet = aIntegrationMesh->get_set_by_index( tMeshSetIndex );

                // Get number of clusters on set
                uint tNumberOfClusters = tSet->get_num_clusters_on_set();

                // Clusters per set
                for ( uint tClusterIndex = 0; tClusterIndex < tNumberOfClusters; tClusterIndex++ )
                {
                    // get pointer to cluster
                    const mtk::Cluster* tCluster = tSet->get_clusters_by_index( tClusterIndex );

                    // Indices, IDs, and ownership of base cell nodes in cluster
                    // FIXME this is really bad and slow. especially when building the pdvs; should return const &
                    mtk::Cell const* tBaseCell = tCluster->get_interpolation_cell( mtk::Master_Slave::MASTER ).get_base_cell();

                    Matrix< IndexMat > tNodeIndicesInCluster     = tBaseCell->get_vertex_inds();
                    Matrix< DDRMat >   tNodeCoordinatesInCluster = tBaseCell->get_vertex_coords();

                    // number of base nodes in cluster
                    uint tNumberOfBaseNodes = tNodeIndicesInCluster.length();

                    tCounter( tMeshSetIndex ) = tCounter( tMeshSetIndex ) + tNumberOfBaseNodes;

                    // number of spatial dimension
                    tNumSpatialDim = tNodeCoordinatesInCluster.n_cols();

                    // this is a hack. the hole function is super slow and should be rewritten
                    // hack was implemented intentionally to get something running. Redo completely when rewriting GEN
                    if ( tSet->get_set_type() == moris::SetType::DOUBLE_SIDED_SIDESET )
                    {
                        mtk::Cell const* tBaseCellSlave = tCluster->get_interpolation_cell( mtk::Master_Slave::SLAVE ).get_base_cell();

                        Matrix< IndexMat > tNodeIndicesInCluster = tBaseCellSlave->get_vertex_inds();

                        // number of base nodes in cluster
                        tNumberOfBaseNodes = tNodeIndicesInCluster.length();

                        tCounter( tMeshSetIndex ) = tCounter( tMeshSetIndex ) + tNumberOfBaseNodes;
                    }
                }
            }

            for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
            {
                tNodeIndicesPerSet( tMeshSetIndex ).resize( tCounter( tMeshSetIndex ), 1 );
                tNodeIdsPerSet( tMeshSetIndex ).resize( tCounter( tMeshSetIndex ), 1 );
                tNodeOwnersPerSet( tMeshSetIndex ).resize( tCounter( tMeshSetIndex ), 1 );
                tNodeCoordinatesPerSet( tMeshSetIndex ).resize( tCounter( tMeshSetIndex ), tNumSpatialDim );
            }

            // Loop through sets in integration mesh
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
            {
                uint      tCurrentNode = 0;
                mtk::Set* tSet         = aIntegrationMesh->get_set_by_index( tMeshSetIndex );

                // Get number of clusters on set
                uint tNumberOfClusters = tSet->get_num_clusters_on_set();

                // Clusters per set
                for ( uint tClusterIndex = 0; tClusterIndex < tNumberOfClusters; tClusterIndex++ )
                {
                    // get pointer to cluster
                    const mtk::Cluster* tCluster = tSet->get_clusters_by_index( tClusterIndex );

                    // Indices, IDs, and ownership of base cell nodes in cluster
                    // FIXME this is really bad and slow. especially when building the pdvs; should return const &
                    mtk::Cell const* tBaseCell = tCluster->get_interpolation_cell( mtk::Master_Slave::MASTER ).get_base_cell();

                    Matrix< IndexMat > tNodeIndicesInCluster     = tBaseCell->get_vertex_inds();
                    Matrix< IndexMat > tNodeIdsInCluster         = tBaseCell->get_vertex_ids();
                    Matrix< IndexMat > tNodeOwnersInCluster      = tBaseCell->get_vertex_owners();
                    Matrix< DDRMat >   tNodeCoordinatesInCluster = tBaseCell->get_vertex_coords();

                    // number of base nodes in cluster
                    uint tNumberOfBaseNodes = tNodeIndicesInCluster.length();

                    // check for consistency
                    MORIS_ASSERT( tNodeIdsInCluster.length() == tNumberOfBaseNodes && tNodeOwnersInCluster.length() == tNumberOfBaseNodes,
                        "Geometry_Engine::create_interpolation_pdv_hosts - inconsistent cluster information.\n" );

                    // FIXME we have nodes up to 8 times in this list in 3d
                    for ( uint tNodeInCluster = 0; tNodeInCluster < tNumberOfBaseNodes; tNodeInCluster++ )
                    {
                        tNodeIndicesPerSet( tMeshSetIndex )( tCurrentNode ) = tNodeIndicesInCluster( tNodeInCluster );
                        tNodeIdsPerSet( tMeshSetIndex )( tCurrentNode )     = tNodeIdsInCluster( tNodeInCluster );
                        tNodeOwnersPerSet( tMeshSetIndex )( tCurrentNode )  = tNodeOwnersInCluster( tNodeInCluster );

                        tNodeCoordinatesPerSet( tMeshSetIndex ).get_row( tCurrentNode ) =
                            tNodeCoordinatesInCluster.get_row( tNodeInCluster );

                        tCurrentNode++;
                    }

                    // this is a hack. the hole function is super slow and should be rewritten
                    // hack was implemented intentionally to get something running. Redo completely when rewriting GEN
                    if ( tSet->get_set_type() == moris::SetType::DOUBLE_SIDED_SIDESET )
                    {

                        mtk::Cell const* tBaseCellSlave = tCluster->get_interpolation_cell( mtk::Master_Slave::SLAVE ).get_base_cell();

                        Matrix< IndexMat > tNodeIndicesInCluster     = tBaseCellSlave->get_vertex_inds();
                        Matrix< IndexMat > tNodeIdsInCluster         = tBaseCellSlave->get_vertex_ids();
                        Matrix< IndexMat > tNodeOwnersInCluster      = tBaseCellSlave->get_vertex_owners();
                        Matrix< DDRMat >   tNodeCoordinatesInCluster = tBaseCellSlave->get_vertex_coords();

                        // number of base nodes in cluster
                        tNumberOfBaseNodes = tNodeIndicesInCluster.length();

                        // check for consistency
                        MORIS_ASSERT( tNodeIdsInCluster.length() == tNumberOfBaseNodes && tNodeOwnersInCluster.length() == tNumberOfBaseNodes,
                            "Geometry_Engine::create_interpolation_pdv_hosts - inconsistent cluster information.\n" );

                        // FIXME we have nodes up to 8 times in this list in 3d
                        for ( uint tNodeInCluster = 0; tNodeInCluster < tNumberOfBaseNodes; tNodeInCluster++ )
                        {
                            tNodeIndicesPerSet( tMeshSetIndex )( tCurrentNode ) = tNodeIndicesInCluster( tNodeInCluster );
                            tNodeIdsPerSet( tMeshSetIndex )( tCurrentNode )     = tNodeIdsInCluster( tNodeInCluster );
                            tNodeOwnersPerSet( tMeshSetIndex )( tCurrentNode )  = tNodeOwnersInCluster( tNodeInCluster );

                            tNodeCoordinatesPerSet( tMeshSetIndex ).get_row( tCurrentNode ) =
                                tNodeCoordinatesInCluster.get_row( tNodeInCluster );

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
            tNodeCoordinatesPerSet,
            aPdvTypes );
    }

    //--------------------------------------------------------------------------------------------------------------
    void
    Geometry_Engine::set_integration_pdv_types( mtk::Integration_Mesh* aIntegrationMesh )
    {
        // Get information from integration mesh
        uint tNumSets = aIntegrationMesh->get_num_sets();

        // Cell of IG PDV_Type types
        Cell< PDV_Type > tCoordinatePdvs( mNumSpatialDimensions );

        switch ( mNumSpatialDimensions )
        {
        case 2:
        {
            tCoordinatePdvs( 0 ) = PDV_Type::X_COORDINATE;
            tCoordinatePdvs( 1 ) = PDV_Type::Y_COORDINATE;
            break;
        }
        case 3:
        {
            tCoordinatePdvs( 0 ) = PDV_Type::X_COORDINATE;
            tCoordinatePdvs( 1 ) = PDV_Type::Y_COORDINATE;
            tCoordinatePdvs( 2 ) = PDV_Type::Z_COORDINATE;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Geometry Engine only works for 2D and 3D models." );
        }
        }

        // Loop through sets
        Cell< Cell< Cell< PDV_Type > > > tPdvTypes( tNumSets );
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            // PDV_Type types per set
            tPdvTypes( tMeshSetIndex ).resize( 1 );
            tPdvTypes( tMeshSetIndex )( 0 ) = tCoordinatePdvs;
        }

        // Set PDV types
        mPDVHostManager.set_integration_pdv_types( tPdvTypes );
        mPDVHostManager.set_requested_integration_pdv_types( tCoordinatePdvs );
    }

    void
    Geometry_Engine::assign_property_to_pdv_hosts(
        std::shared_ptr< Property > aPropertyPointer,
        PDV_Type                    aPdvType,
        mtk::Integration_Mesh*      aIntegrationMesh,
        Matrix< DDUMat >            aSetIndices )
    {
        for ( uint tSet = 0; tSet < aSetIndices.length(); tSet++ )
        {
            // get the mesh set from index
            mtk::Set* tSetPointer = aIntegrationMesh->get_set_by_index( aSetIndices( tSet ) );

            // get the list of cluster on mesh set
            Cell< mtk::Cluster const* > tClusterPointers = tSetPointer->get_clusters_on_set();

            // get number of clusters on mesh set
            uint tNumClusters = tClusterPointers.size();

            // loop over the clusters on mesh set
            for ( uint iClust = 0; iClust < tNumClusters; iClust++ )
            {
                // get the IP cell from cluster
                mtk::Cell const& tIPCell = tClusterPointers( iClust )->get_interpolation_cell( mtk::Master_Slave::MASTER );

                // get the vertices from IP cell
                Cell< mtk::Vertex* > tVertices = tIPCell.get_base_cell()->get_vertex_pointers();

                // get the number of vertices on IP cell
                uint tNumVerts = tVertices.size();

                // loop over vertices on IP cell
                for ( uint iVert = 0; iVert < tNumVerts; iVert++ )
                {
                    // get the vertex index
                    moris_index tVertIndex = tVertices( iVert )->get_index();

                    // ask pdv host manager to assign to vertex a pdv type and a property
                    mPDVHostManager.create_interpolation_pdv( uint( tVertIndex ), aPdvType, aPropertyPointer );
                }
            }

            // this is kind of a hack. recommending rewriting it properly when rewriting GEN
            if ( tSetPointer->get_set_type() == moris::SetType::DOUBLE_SIDED_SIDESET )
            {
                // loop over the clusters on mesh set
                for ( uint iClust = 0; iClust < tNumClusters; iClust++ )
                {
                    // get the IP cell from cluster
                    mtk::Cell const& tIPCell = tClusterPointers( iClust )->get_interpolation_cell( mtk::Master_Slave::SLAVE );

                    // get the vertices from IP cell
                    Cell< mtk::Vertex* > tVertices = tIPCell.get_base_cell()->get_vertex_pointers();

                    // get the number of vertices on IP cell
                    uint tNumVerts = tVertices.size();

                    // loop over vertices on IP cell
                    for ( uint iVert = 0; iVert < tNumVerts; iVert++ )
                    {
                        // get the vertex index
                        moris_index tVertIndex = tVertices( iVert )->get_index();

                        // ask pdv host manager to assign to vertex a pdv type and a property
                        mPDVHostManager.create_interpolation_pdv( uint( tVertIndex ), aPdvType, aPropertyPointer );
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::setup_initial_geometric_proximities( mtk::Interpolation_Mesh* aMesh )
    {
        mVertexGeometricProximity =
            Cell< Geometric_Proximity >( aMesh->get_num_nodes(), Geometric_Proximity( mGeometries.size() ) );

        // iterate through vertices then geometries
        for ( uint iV = 0; iV < aMesh->get_num_nodes(); iV++ )
        {
            Matrix< DDRMat > tCoords = aMesh->get_node_coordinate( moris_index( iV ) );

            mVertexGeometricProximity( iV ).mAssociatedVertexIndex = (moris_index)iV;

            for ( uint iGeometryIndex = 0; iGeometryIndex < mGeometries.size(); iGeometryIndex++ )
            {
                real tVertGeomVal = mGeometries( iGeometryIndex )->get_field_value( iV, tCoords );

                moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                mVertexGeometricProximity( iV ).set_geometric_proximity( tGeomProxIndex, iGeometryIndex );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Geometry_Engine::get_geometric_proximity_index( real const& aGeometricVal )
    {
        moris_index tGeometricProxIndex = MORIS_INDEX_MAX;
        if ( aGeometricVal - mIsocontourThreshold < -mIsocontourTolerance )
        {
            tGeometricProxIndex = 0;
        }
        else if ( std::abs( aGeometricVal - mIsocontourThreshold ) < mIsocontourTolerance )
        {
            tGeometricProxIndex = 1;
        }
        else if ( aGeometricVal - mIsocontourThreshold > mIsocontourTolerance )
        {
            tGeometricProxIndex = 2;
        }

        return tGeometricProxIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Geometry_Engine::get_queued_intersection_geometric_proximity_index( moris_index const& aGeomIndex )
    {
        // parent vertex
        moris_index tParentVertexIndex0 = mQueuedIntersectionNode->mAncestorNodeIndices( 0 );
        moris_index tParentVertexIndex1 = mQueuedIntersectionNode->mAncestorNodeIndices( 1 );

        // parent vertex proximity wrt aGeomIndex
        moris_index tParentProx0 = mVertexGeometricProximity( tParentVertexIndex0 ).get_geometric_proximity( aGeomIndex );
        moris_index tParentProx1 = mVertexGeometricProximity( tParentVertexIndex1 ).get_geometric_proximity( aGeomIndex );

        // 0 - G(x) < threshold
        // 1 - G(x) == threshold
        // 2 - G(x) > threshold
        // verify we dont transition across interface
        // MORIS_ERROR((tParentProx0 == 0 && tParentProx1 == 0) ||
        //             (tParentProx0 == 0 && tParentProx1 == 1) ||
        //             (tParentProx0 == 1 && tParentProx1 == 0) ||
        //             (tParentProx0 == 1 && tParentProx1 == 1) ||
        //             (tParentProx0 == 2 && tParentProx1 == 2) ||
        //             (tParentProx0 == 2 && tParentProx1 == 1) ||
        //             (tParentProx0 == 1 && tParentProx1 == 2),   "Invalid proximity data");
        //
        // add them together
        moris_index tSum = tParentProx0 + tParentProx1;

        // proximity value
        moris_index tProxIndex = MORIS_INDEX_MAX;
        if ( tSum == 0 )
        {
            tProxIndex = 0;
        }
        else if ( tSum == 1 )
        {
            tProxIndex = 0;
        }
        else if ( tSum == 2 )
        {
            tProxIndex = 1;
        }
        else if ( tSum == 3 )
        {
            tProxIndex = 2;
        }
        else if ( tSum == 4 )
        {
            tProxIndex = 2;
        }
        else
        {
            MORIS_ASSERT( 0, "Proximity determination failed." );
        }

        return tProxIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Geometry_Engine::admit_queued_intersection_geometric_proximity( uint aNodeIndex )
    {

        MORIS_ERROR( aNodeIndex == mVertexGeometricProximity.size(), "Index mismatch" );
        // initialize proximity data
        mVertexGeometricProximity.push_back( Geometric_Proximity( mGeometries.size() ) );

        // node index associated with this proximity
        mVertexGeometricProximity( aNodeIndex ).mAssociatedVertexIndex = aNodeIndex;

        // geometry iteration through previous geometries
        for ( uint tGeometryIndex = 0; tGeometryIndex < this->get_active_geometry_index(); tGeometryIndex++ )
        {
            moris_index tProxIndex = this->get_queued_intersection_geometric_proximity_index( tGeometryIndex );
            mVertexGeometricProximity( aNodeIndex ).set_geometric_proximity( tProxIndex, tGeometryIndex );
        }

        // place the current one on the interface
        mVertexGeometricProximity( aNodeIndex ).set_geometric_proximity( 1, this->get_active_geometry_index() );

        if ( this->get_active_geometry_index() != mGeometries.size() - 1 )
        {
            // iterate through following geometries (here we just compute the vertex value to determine proximity)
            for ( uint tGeometryIndex = this->get_active_geometry_index() + 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( aNodeIndex, mQueuedIntersectionNode->get_global_coordinates() );

                moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                mVertexGeometricProximity( aNodeIndex ).set_geometric_proximity( tGeomProxIndex, tGeometryIndex );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table
    Geometry_Engine::create_phase_table(
        Cell< Cell< ParameterList > > aParameterLists,
        std::shared_ptr< Library_IO > aLibrary )
    {
        // Get number of geometries
        uint tNumGeometries = aParameterLists( 1 ).size();

        // Recreate phase table via different methods if needed
        std::string tPhaseFunctionName = aParameterLists( 0 )( 0 ).get< std::string >( "phase_function_name" );
        if ( tPhaseFunctionName != "" )
        {
            // User-defined phase function
            return Phase_Table(
                aLibrary->load_function< PHASE_FUNCTION >( tPhaseFunctionName ),
                aParameterLists( 0 )( 0 ).get< sint >( "number_of_phases" ) );
        }
        else if ( aParameterLists( 0 )( 0 ).get< std::string >( "phase_table" ) != "" )
        {
            // User-defined bulk phases
            return Phase_Table( tNumGeometries, string_to_mat< DDUMat >( aParameterLists( 0 )( 0 ).get< std::string >( "phase_table" ) ) );
        }
        else
        {
            // Unique phase per geometry combination
            return Phase_Table( tNumGeometries );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table
    Geometry_Engine::create_phase_table(
        uint             aNumGeometries,
        Matrix< DDUMat > aBulkPhases,
        PHASE_FUNCTION   aPhaseFunction,
        uint             aNumPhases )
    {
        if ( aPhaseFunction )
        {
            return Phase_Table( aPhaseFunction, aNumPhases );
        }
        else if ( aBulkPhases.length() > 0 )
        {
            return Phase_Table( aNumGeometries, aBulkPhases );
        }
        else
        {
            return Phase_Table( aNumGeometries );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}// namespace ge
}// namespace moris
