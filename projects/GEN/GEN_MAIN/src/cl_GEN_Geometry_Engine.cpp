/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Engine.cpp
 *
 */

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
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"

// MTK
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Enums.hpp"

// XTK FIXME
#include "cl_XTK_Topology.hpp"

// SOL FIXME
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

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
                mtk::Mesh*                    aMesh )
                : mPhaseTable( create_phase_table( aParameterLists, aLibrary ) )
        {
            // Tracer
            Tracer tTracer( "GEN", "Create geometry engine" );

            // Level set options
            mEvaluateNewChildNodeAsLinear = aParameterLists( 0 )( 0 ).get< bool >( "evaluate_new_pts_as_linear" );

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

            // iterate through geometries if any are multilinear, we turn the linear flag on
            for ( uint iGeom = 0; iGeom < mGeometries.size(); iGeom++ )
            {
                if ( mGeometries( iGeom )->get_intersection_interpolation() == Intersection_Interpolation::MULTILINEAR )
                {
                    MORIS_LOG_INFO( "New Child Vertices will be evaluated as using linear background cells" );
                    mEvaluateNewChildNodeAsLinear = true;
                }
            }

            MORIS_ERROR( mGeometries.size() <= MAX_GEOMETRIES,
                    "Number of geometries exceeds MAX_GEOMETRIES, please change this in GEN_typedefs.hpp" );

            // Properties
            mProperties = create_properties(
                    aParameterLists( 2 ),
                    mInitialPrimitiveADVs,
                    mGeometries,
                    mLibrary );

            // Set requested PDVs
            Cell< std::string > tRequestedPdvNames = string_to_cell< std::string >( aParameterLists( 0 )( 0 ).get< std::string >( "PDV_types" ) );
            Cell< PDV_Type >    tRequestedPdvTypes( tRequestedPdvNames.size() );

            map< std::string, PDV_Type > tPdvTypeMap = get_pdv_type_map();

            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tRequestedPdvTypes.size(); tPdvTypeIndex++ )
            {
                tRequestedPdvTypes( tPdvTypeIndex ) = tPdvTypeMap[ tRequestedPdvNames( tPdvTypeIndex ) ];
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
                Geometry_Engine_Parameters aParameters )
                : mGeometries( aParameters.mGeometries )
                , mProperties( aParameters.mProperties )
                , mPhaseTable( create_phase_table( aParameters.mGeometries.size(), aParameters.mBulkPhases ) )
                , mInitialPrimitiveADVs( aParameters.mADVs )
                , mTimeOffset( aParameters.mTimeOffset )
        {
            // Tracer
            Tracer tTracer( "GEN", "Create geometry engine" );

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

        Pdv_Host_Manager*
        Geometry_Engine::get_pdv_host_manager()
        {
            return &mPDVHostManager;
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
                PHASE_FUNCTION      aPhaseFunction,
                uint                aNumPhases,
                Cell< std::string > aPhaseNames )
        {
            mPhaseTable.set_phase_function( aPhaseFunction, aNumPhases, aPhaseNames );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::set_dQIdp(
                moris::Cell< Matrix< DDRMat >* > adQIdp,
                Matrix< DDSMat >*                aMap )
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

        //--------------------------------------------------------------------------------------------------------------

        bool
        Geometry_Engine::geometric_query( Geometric_Query_Interface* aGeometricQuery )
        {
            if ( aGeometricQuery->get_query_type() == Query_Type::INTERSECTION_NO_LOCATION )
            {
                mActiveGeometryIndex = aGeometricQuery->get_geometric_index();
                return this->is_intersected( aGeometricQuery->get_query_entity_to_vertex_connectivity(),
                        aGeometricQuery->get_query_indexed_coordinates() );
            }
            else if ( aGeometricQuery->get_query_type() == Query_Type::INTERSECTION_LOCATION )
            {
                // this preprocessing can be streamlined a lot
                mActiveGeometryIndex = aGeometricQuery->get_geometric_index();

                // get the vertex indices connected to the queried edge
                Matrix< IndexMat > const & tEdgeToVertex = aGeometricQuery->get_query_entity_to_vertex_connectivity();

                // get the physical vertex coordinates in the parent BG cell
                moris::Cell< std::shared_ptr< Matrix< DDRMat > > >* tQueryIndexedCoords = aGeometricQuery->get_query_indexed_coordinates();

                // get the vertex indices in the parent BG cell
                Matrix< IndexMat > tParentEntityIndices = aGeometricQuery->get_query_parent_entity_connectivity();

                // Copy the information retrieved above into arrays compatible with the "queue_intersection" function below
                // TODO: the function below should be re-written to accept a different input format instead of this copying
                Cell< Matrix< DDRMat > > tParentEntityCoords( tParentEntityIndices.numel() );
                Matrix< DDUMat >         tParentEntityIndicesUINT( tParentEntityIndices.numel() );
                tParentEntityCoords.reserve( tParentEntityIndices.numel() * 3 );
                for ( uint i = 0; i < tParentEntityIndices.numel(); i++ )
                {
                    tParentEntityCoords( i )      = *( *tQueryIndexedCoords )( tParentEntityIndices( i ) );
                    tParentEntityIndicesUINT( i ) = (uint)tParentEntityIndices( i );
                }

                // compute intersection
                return queue_intersection(
                        tEdgeToVertex( 0 ),
                        tEdgeToVertex( 1 ),
                        aGeometricQuery->get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 0 ) ),
                        aGeometricQuery->get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 1 ) ),
                        *( *tQueryIndexedCoords )( tEdgeToVertex( 0 ) ),
                        *( *tQueryIndexedCoords )( tEdgeToVertex( 1 ) ),
                        tParentEntityIndicesUINT,
                        tParentEntityCoords );
            }
            else    // unsupported query type
            {
                MORIS_ERROR(
                        false,
                        "ge::Geometry_Engine::geometric_query() - "
                        "Only supports Query_Type::INTERSECTION_NO_LOCATION and Query_Type::INTERSECTION_LOCATION" );
                return false;
            }

        }    // end function: Geometry_Engine::geometric_query()

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

            switch ( mGeometries( mActiveGeometryIndex )->get_intersection_mode() )
            {
                case Intersection_Mode::LEVEL_SET:
                {
                    // get the current geometries level set parameters
                    real tIsocontourThreshold = mGeometries( mActiveGeometryIndex )->get_isocontour_threshold();
                    real tIsocontourTolerance = mGeometries( mActiveGeometryIndex )->get_isocontour_tolerance();

                    // Initialize by evaluating the first node
                    real tMin = mGeometries( mActiveGeometryIndex )->get_field_value( 0, aNodeCoordinates.get_row( 0 ) );
                    real tMax = tMin;

                    // Evaluate the rest of the nodes
                    for ( uint tNodeCount = 1; tNodeCount < aNodeIndices.length(); tNodeCount++ )
                    {
                        real tEval = mGeometries( mActiveGeometryIndex )->    //
                                     get_field_value( tNodeCount, aNodeCoordinates.get_row( tNodeCount ) );

                        tMin = std::min( tMin, tEval );
                        tMax = std::max( tMax, tEval );
                    }

                    tIsIntersected = ( tMax >= tIsocontourThreshold and tMin <= tIsocontourThreshold )
                                  or ( std::abs( tMax - tIsocontourThreshold ) < tIsocontourTolerance )
                                  or ( std::abs( tMin - tIsocontourThreshold ) < tIsocontourTolerance );

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
                case Intersection_Mode::SURFACE_MESH:
                {
                    // BRENDAN TODO
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry_Engine::is_intersected(), unknown intersection type." );
                }
            }

            // Return result
            return tIsIntersected;

        }    // end function: Geometry_Engine::is_intersected()

        //--------------------------------------------------------------------------------------------------------------

        bool
        Geometry_Engine::is_intersected(
                const Matrix< IndexMat >&                           aNodeIndices,
                moris::Cell< std::shared_ptr< Matrix< DDRMat > > >* aNodeCoordinates )
        {
            // Check input
            // MORIS_ASSERT(aNodeIndices.length() == aNodeCoordinates.n_rows(),
            //         "Geometry engine must be provided the same number of node indices as node coordinates for "
            //         "determining if an element is intersected or not.");
            // MORIS_ASSERT(aNodeIndices.length() > 0,
            //         "Geometry engine must be provided at least 1 node to determine if an element is intersected or not.");

            // get the current geometries intersection mode, isocontour threshold
            Intersection_Mode tIntersectionMode    = mGeometries( mActiveGeometryIndex )->get_intersection_mode();
            real              tIsocontourThreshold = mGeometries( mActiveGeometryIndex )->get_isocontour_threshold();
            real              tIsocontourTolerance = mGeometries( mActiveGeometryIndex )->get_isocontour_tolerance();

            bool tIsIntersected = false;

            switch ( tIntersectionMode )
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

                    tIsIntersected = ( tMax >= tIsocontourThreshold and tMin <= tIsocontourThreshold )
                                  or ( std::abs( tMax - tIsocontourThreshold ) < tIsocontourTolerance )
                                  or ( std::abs( tMin - tIsocontourThreshold ) < tIsocontourTolerance );

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
                case Intersection_Mode::SURFACE_MESH:
                {
                    // BRENDAN TODO
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
            // get access to the geometry
            std::shared_ptr< ge::Geometry > mGeometry = mGeometries( mActiveGeometryIndex );

            // Get the current geometries intersection mode
            Intersection_Mode tIntersectionMode = mGeometry->get_intersection_mode();

            // get the the intersection nodes
            std::shared_ptr< ge::Intersection_Node > tFirstIntersectionNode  = mPDVHostManager.get_intersection_node( aFirstNodeIndex );
            std::shared_ptr< ge::Intersection_Node > tSecondIntersectionNode = mPDVHostManager.get_intersection_node( aSecondNodeIndex );

            // Queue an intersection node
            switch ( tIntersectionMode )
            {
                case Intersection_Mode::LEVEL_SET:
                {
                    switch ( mGeometry->get_intersection_interpolation() )
                    {
                        case Intersection_Interpolation::LINEAR:
                        {
                            mQueuedIntersectionNode = std::make_shared< Intersection_Node_Linear >(
                                    tFirstIntersectionNode,
                                    tSecondIntersectionNode,
                                    aFirstNodeIndex,
                                    aSecondNodeIndex,
                                    aFirstNodeGlobalCoordinates,
                                    aSecondNodeGlobalCoordinates,
                                    mGeometry );
                            break;
                        }

                        case Intersection_Interpolation::MULTILINEAR:
                        {
                            Element_Intersection_Type tInterpolationType =
                                    mNumSpatialDimensions == 2 ? Element_Intersection_Type::Linear_2D : Element_Intersection_Type::Linear_3D;

                            mQueuedIntersectionNode = std::make_shared< Intersection_Node_Bilinear >(
                                    tFirstIntersectionNode,
                                    tSecondIntersectionNode,
                                    aFirstNodeIndex,
                                    aSecondNodeIndex,
                                    aFirstNodeLocalCoordinates,
                                    aSecondNodeLocalCoordinates,
                                    aBackgroundElementNodeIndices,
                                    aBackgroundElementNodeCoordinates,
                                    tInterpolationType,
                                    mGeometry );
                            break;
                        }

                        default:
                        {
                            MORIS_ERROR( false, "Intersection interpolation type not implemented yet." );
                        }

                    }    // end switch: Intersection Interpolation type

                    break;

                }    // end case: Intersection_Mode::LEVEL_SET

                case Intersection_Mode::SURFACE_MESH:
                {
                    // BRENDAN TODO
                    // aBackgroundElementNodeIndices and aBackgroundElementNodeCoordinates need to come from the geometry

                    mQueuedIntersectionNode = std::make_shared< Intersection_Node_Surface_Mesh >(
                        mPDVHostManager.get_intersection_node( aFirstNodeIndex ),
                        mPDVHostManager.get_intersection_node( aSecondNodeIndex ),
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        aFirstNodeLocalCoordinates,
                        aSecondNodeLocalCoordinates,
                        aBackgroundElementNodeIndices,
                        aBackgroundElementNodeCoordinates 
                    );
                }

                case Intersection_Mode::COLORING:
                {
                    // Determine if edge is intersected
                    real tFieldValue1 = mGeometry->get_field_value( aFirstNodeIndex, aFirstNodeGlobalCoordinates );
                    real tFieldValue2 = mGeometry->get_field_value( aSecondNodeIndex, aSecondNodeGlobalCoordinates );
                    if ( tFieldValue1 != tFieldValue2 )
                    {
                        mQueuedIntersectionNode = std::make_shared< Intersection_Node_Linear >(
                                tFirstIntersectionNode,
                                tSecondIntersectionNode,
                                aFirstNodeIndex,
                                aSecondNodeIndex,
                                aFirstNodeGlobalCoordinates,
                                aSecondNodeGlobalCoordinates,
                                mGeometry );
                    }
                    else
                    {
                        return false;
                    }

                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Geometry_Engine::queue_intersection() - Unknown intersection type." );
                }

            }    // end switch: Intersection Mode

            //
            bool tParentEdgeIsIntersected = mQueuedIntersectionNode->parent_edge_is_intersected();

            //
            return tParentEdgeIsIntersected;

        }    // end function: Geometry_Engine::queue_intersection()

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
            // get parent indices
            moris_index tFirstParentIndex  = mQueuedIntersectionNode->get_first_parent_node_index();
            moris_index tSecondParentIndex = mQueuedIntersectionNode->get_second_parent_node_index();

            // check if parent nodes are PDVs, i.e. are defined directly or indirectly on variable geometry
            bool tFirstParentIsPDV  = mPDVHostManager.get_intersection_node( tFirstParentIndex ) != nullptr;
            bool tSecondParentIsPDV = mPDVHostManager.get_intersection_node( tSecondParentIndex ) != nullptr;

            // Assign as PDV host if constructed on adv dependent geometry or parent nodes are adv dependent
            if ( mGeometries( mActiveGeometryIndex )->depends_on_advs() || tFirstParentIsPDV || tSecondParentIsPDV )
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
            // Tracer
            Tracer tTracer( "GEN", "Create new child nodes" );

            // resize proximities
            mVertexGeometricProximity.resize(
                    mVertexGeometricProximity.size() + aNewNodeIndices->size(),
                    Geometric_Proximity( mGeometries.size() ) );

            // Loop over nodes
            for ( uint tNode = 0; tNode < aNewNodeIndices->size(); tNode++ )
            {
                std::shared_ptr< Child_Node > tChildNode = std::make_shared< Child_Node >(
                        ( *aNewNodeParentCell )( tNode ),
                        ( *aParamCoordRelativeToParent )( tNode ).get(),
                        this->mEvaluateNewChildNodeAsLinear );

                mVertexGeometricProximity( ( *aNewNodeIndices )( tNode ) ).mAssociatedVertexIndex = ( *aNewNodeIndices )( tNode );

                // debug - needs to go
                Matrix< DDRMat > const & tCoord = ( *aNodeCoordinates )( tNode );

                // Assign to geometries
                for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
                {
                    mGeometries( tGeometryIndex )->add_child_node( ( *aNewNodeIndices )( tNode ), tChildNode );

                    // FIXME: need to get value from child element based on element interpolation
                    real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( ( *aNewNodeIndices )( tNode ), tCoord );

                    moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                    mVertexGeometricProximity( ( *aNewNodeIndices )( tNode ) ).set_geometric_proximity( tGeomProxIndex, tGeometryIndex );
                }
            }

            // Set max node index
            if ( aNewNodeIndices->size() > 0 )
            {
                mPDVHostManager.set_num_background_nodes( ( *aNewNodeIndices )( aNewNodeIndices->size() - 1 ) + 1 );
            }

        } // end function: Geometry_Engine::create_new_child_nodes()

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::create_new_child_nodes(
                const Cell< moris_index >&               aNewNodeIndices,
                const Cell< Element_Intersection_Type >& aParentIntersectionType,
                const Cell< Matrix< IndexMat > >&        tVertexIndices,
                const Cell< Matrix< DDRMat > >&          aParamCoordRelativeToParent,
                const Matrix< DDRMat >&                  aGlobalNodeCoord )
        {
            // Tracer
            Tracer tTracer( "GEN", "Create new child nodes" );

            // get current geometries level set info
            real tIsocontourThreshold = mGeometries( mActiveGeometryIndex )->get_isocontour_threshold();
            real tIsocontourTolerance = mGeometries( mActiveGeometryIndex )->get_isocontour_tolerance();

            // resize proximities
            mVertexGeometricProximity.resize(
                    mVertexGeometricProximity.size() + aNewNodeIndices.size(),
                    Geometric_Proximity( mGeometries.size() ) );

            // Loop over nodes
            for ( uint tNode = 0; tNode < aNewNodeIndices.size(); tNode++ )
            {
                Matrix< DDUMat >         tParentNodeIndices( tVertexIndices( tNode ).numel(), 1 );
                Cell< Matrix< DDRMat > > tParentNodeCoordinates( tParentNodeIndices.length() );

                for ( uint tParentNode = 0; tParentNode < tParentNodeIndices.length(); tParentNode++ )
                {
                    tParentNodeIndices( tParentNode )     = tVertexIndices( tNode )( tParentNode );
                    tParentNodeCoordinates( tParentNode ) = aGlobalNodeCoord.get_row( tParentNodeIndices( tParentNode ) );
                }

                std::shared_ptr< Child_Node > tChildNode = std::make_shared< Child_Node >(
                        tParentNodeIndices,
                        tParentNodeCoordinates,
                        aParentIntersectionType( tNode ),
                        aParamCoordRelativeToParent( tNode ) );

                mVertexGeometricProximity( aNewNodeIndices( tNode ) ).mAssociatedVertexIndex = aNewNodeIndices( tNode );

                Matrix< DDRMat > tCoord = aGlobalNodeCoord.get_row( aNewNodeIndices( tNode ) );

                // Assign to geometries
                for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
                {
                    mGeometries( tGeometryIndex )->add_child_node( aNewNodeIndices( tNode ), tChildNode );

                    real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( aNewNodeIndices( tNode ), tCoord );

                    moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                    if ( std::abs( tVertGeomVal - tIsocontourThreshold ) < tIsocontourTolerance )
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

        //--------------------------------------------------------------------------------------------------------------

        moris_index
        Geometry_Engine::is_interface_vertex(
                moris_index aNodeIndex,
                moris_index aGeometryIndex )
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
        Geometry_Engine::get_elem_phase_index( Matrix< IndexMat > const & aElemOnOff )
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

            // 0 - G(x) < threshold
            // 1 - G(x) == threshold
            // 2 - G(x) > threshold

            size_t tPhaseOnOff = 0;

            if ( tProxIndex == 2 )
            {
                tPhaseOnOff = 1;
            }

            return tPhaseOnOff;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index
        Geometry_Engine::get_node_proximity_wrt_a_geometry(
                uint aNodeIndex,
                uint aGeometryIndex )
        {
            return mVertexGeometricProximity( aNodeIndex ).get_geometric_proximity( aGeometryIndex );
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

            // Loop over properties to create PDVs
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                // Get PDV type from property
                tPDVTypeGroup( 0 ) = mProperties( tPropertyIndex )->get_pdv_type();

                // Get mesh set indices and names
                Matrix< DDUMat >    tMeshSetIndices = mProperties( tPropertyIndex )->get_pdv_mesh_set_indices();
                Cell< std::string > tMeshSetNames   = mProperties( tPropertyIndex )->get_pdv_mesh_set_names();

                // Convert mesh set names to indices
                uint tNumSetIndices = tMeshSetIndices.length();
                tMeshSetIndices.resize( tNumSetIndices + tMeshSetNames.size(), 1 );

                // number of mesh sets for current property
                uint tTotalNumberOfSets = tMeshSetIndices.length();

                // Set for each property index the list of mesh set indices TODO pass this to property to have it update its own mesh set indices
                for ( uint tSetIndexPosition = tNumSetIndices; tSetIndexPosition < tTotalNumberOfSets; tSetIndexPosition++ )
                {
                    tMeshSetIndices( tSetIndexPosition ) =
                            tIntegrationMesh->get_set_index_by_name( tMeshSetNames( tSetIndexPosition - tNumSetIndices ) );
                }

                // Assign PDV types to each mesh set
                for ( uint tSetIndexPosition = 0; tSetIndexPosition < tTotalNumberOfSets; tSetIndexPosition++ )
                {
                    uint tMeshSetIndex = tMeshSetIndices( tSetIndexPosition );
                    tPdvTypes( tMeshSetIndex ).push_back( tPDVTypeGroup );
                }

                // Add nodal data from the interpolation mesh
                mProperties( tPropertyIndex )->add_nodal_data( tInterpolationMesh );
            }

            // Set interpolation PDV types in host manager
            mPDVHostManager.set_interpolation_pdv_types( tPdvTypes );

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
                    tPdvTypes );

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
            for ( uint iHeader = 0; iHeader < aMesh->get_spatial_dim(); iHeader++ )
            {
                tStringStream << "Coords_" + std::to_string( iHeader ) << ",";
            }

            for ( uint iHeader = 0; iHeader < this->get_num_geometries(); iHeader++ )
            {
                tStringStream << "gval_" + std::to_string( iHeader ) << ",";
            }
            for ( uint iHeader = 0; iHeader < this->get_num_geometries(); iHeader++ )
            {
                tStringStream << "gprox_" + std::to_string( iHeader );
                if ( iHeader != this->get_num_geometries() - 1 )
                {
                    tStringStream << ",";
                }
            }
            tStringStream << std::endl;
            // iterate through vertices
            for ( uint iV = 0; iV < aMesh->get_num_nodes(); iV++ )
            {
                mtk::Vertex& tVertex = aMesh->get_mtk_vertex( (moris_index)iV );
                tStringStream << tVertex.get_id() << ",";
                tStringStream << tVertex.get_index() << ",";
                tStringStream << tVertex.get_owner() << ",";
                tStringStream << par_rank() << ",";
                Matrix< DDRMat > tCoords = tVertex.get_coords();

                for ( uint iSp = 0; iSp < aMesh->get_spatial_dim(); iSp++ )
                {
                    tStringStream << std::scientific << tCoords( iSp ) << ",";
                }
                for ( uint iGeom = 0; iGeom < this->get_num_geometries(); iGeom++ )
                {
                    tStringStream << mGeometries( iGeom )->get_field_value( tVertex.get_index(), tCoords ) << ",";
                }
                for ( uint iGeom = 0; iGeom < this->get_num_geometries(); iGeom++ )
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
        Geometry_Engine::induce_as_interface_vertex_on_active_geometry( moris_index aVertexIndex )
        {
            // to do this I change the geometric proximity to = for the given vertex
            mVertexGeometricProximity( aVertexIndex ).set_geometric_proximity( 1, this->get_active_geometry_index() );
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

        void
        Geometry_Engine::distribute_advs(
                mtk::Mesh_Pair                               aMeshPair,
                moris::Cell< std::shared_ptr< mtk::Field > > aFields,
                mtk::EntityRank                              aADVEntityRank )
        {
            // Tracer
            Tracer tTracer( "GEN", "Distribute ADVs" );

            // Gather all fields
            Cell< std::shared_ptr< Field > > tFields( mGeometries.size() + mProperties.size() );
            std::copy( mGeometries.begin(), mGeometries.end(), tFields.begin() );
            std::copy( mProperties.begin(), mProperties.end(), tFields.begin() + mGeometries.size() );

            // Interpolation mesh
            mtk::Interpolation_Mesh* tMesh = aMeshPair.get_interpolation_mesh();

            //------------------------------------//
            // Determine owned and shared ADV IDs //
            //------------------------------------//
            clock_t tStart_Owned_Shared_ADVs = clock();

            // Set primitive IDs
            Matrix< DDSMat > tPrimitiveADVIds( mInitialPrimitiveADVs.length(), 1 );
            for ( uint tADVIndex = 0; tADVIndex < mInitialPrimitiveADVs.length(); tADVIndex++ )
            {
                tPrimitiveADVIds( tADVIndex ) = tADVIndex;
            }

            // Start with primitive IDs for owned IDs on processor 0
            Matrix< DDSMat > tOwnedADVIds( 0, 0 );
            if ( par_rank() == 0 )
            {
                tOwnedADVIds = tPrimitiveADVIds;
            }

            // this is done to initialize primitive adv positions with gNoID
            mOwnedijklIds.set_size( tPrimitiveADVIds.numel(), 1, gNoID );

            // Owned and shared ADVs per field
            Cell< Matrix< DDUMat > > tSharedCoefficientIndices( tFields.size() );
            Cell< Matrix< DDSMat > > tSharedADVIds( tFields.size() );
            Matrix< DDUMat >         tAllOffsetIDs( tFields.size(), 1 );

            // Get all node indices from the mesh (for now)
            uint             tNumNodes = tMesh->get_num_nodes();
            Matrix< DDUMat > tNodeIndices( tNumNodes, 1 );
            for ( uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++ )
            {
                tNodeIndices( tNodeIndex ) = tNodeIndex;
            }

            moris::Cell< uint > tNumCoeff( tFields.size() );
            // Loop over all geometries to get number of new ADVs
            sint tOffsetID = tPrimitiveADVIds.length();
            for ( uint tFieldIndex = 0; tFieldIndex < tFields.size(); tFieldIndex++ )
            {
                // Determine if level set will be created
                if ( tFields( tFieldIndex )->intended_discretization() )
                {
                    // Get discretization mesh index
                    uint tDiscretizationMeshIndex = tFields( tFieldIndex )->get_discretization_mesh_index();

                    uint tMaxNumberOfCoefficients = tMesh->get_max_num_coeffs_on_proc( tDiscretizationMeshIndex );

                    Matrix< IdMat >    tAllCoefIds( tMaxNumberOfCoefficients, 1, gNoID );
                    Matrix< IndexMat > tAllCoefIndices( tMaxNumberOfCoefficients, 1, gNoIndex );
                    Matrix< IdMat >    tAllCoefOwners( tMaxNumberOfCoefficients, 1, gNoID );
                    Matrix< IdMat >    tAllCoefijklIDs( tMaxNumberOfCoefficients, 1, gNoID );

                    for ( uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++ )
                    {
                        // check whether node has an underlying discretization on this processor
                        bool tNodeHasDiscretization =
                                tMesh->get_mtk_vertex( tNodeIndex ).has_interpolation( tDiscretizationMeshIndex );

                        // process only nodes that have discretization
                        if ( tNodeHasDiscretization )
                        {
                            // get indices and IDs from mtk mesh - FIXME: should return const &
                            const Matrix< IndexMat > tCoefIndices = tMesh->get_coefficient_indices_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex );

                            const Matrix< IdMat > tCoefIds = tMesh->get_coefficient_IDs_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex );

                            const Matrix< IdMat > tCoefOwners = tMesh->get_coefficient_owners_of_node(
                                    tNodeIndex,
                                    tDiscretizationMeshIndex );

                            Matrix< IdMat > tCoeffijklIDs;

                            if ( mtk::MeshType::HMR == tMesh->get_mesh_type() )
                            {
                                tCoeffijklIDs = tMesh->get_coefficient_ijkl_IDs_of_node(
                                        tNodeIndex,
                                        tDiscretizationMeshIndex );
                            }

                            // check that number of indices and ids are the same
                            MORIS_ASSERT( tCoefIds.numel() == tCoefIndices.numel(),
                                    "distribute_advs - numbers of coefficients and ids do not match.\n" );

                            // get number of coefficients for current node
                            uint tNumCoefOfNode = tCoefIds.numel();

                            for ( uint tCoefIndex = 0; tCoefIndex < tNumCoefOfNode; ++tCoefIndex )
                            {
                                // get coefficient index
                                moris_index tCurrentIndex = tCoefIndices( tCoefIndex );

                                // check whether mesh coefficient has already been set
                                if ( tAllCoefIds( tCurrentIndex ) == -1 )
                                {
                                    // increase field coefficient count
                                    tNumCoeff( tFieldIndex )++;

                                    // populate mesh index to mesh coefficient id map
                                    tAllCoefIds( tCurrentIndex ) = tCoefIds( tCoefIndex );

                                    tAllCoefOwners( tCurrentIndex ) = tCoefOwners( tCoefIndex );

                                    if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
                                    {
                                        tAllCoefijklIDs( tCurrentIndex ) = tCoeffijklIDs( tCoefIndex );
                                    }
                                }
                                else
                                {
                                    // check for consistency
                                    MORIS_ASSERT( tAllCoefIds( tCurrentIndex ) == tCoefIds( tCoefIndex ),
                                            "distribute_advs - inconsistent index and ids.\n" );
                                }
                            }
                        }
                    }

                    if ( par_size() > 1 )
                    {
                        this->communicate_missing_owned_coefficients(
                                aMeshPair,
                                tAllCoefIds,
                                tAllCoefOwners,
                                tAllCoefijklIDs,
                                tNumCoeff,
                                tFieldIndex,
                                tDiscretizationMeshIndex,
                                tMesh->get_mesh_type() );
                    }

                    uint tOwnedCounter  = 0;
                    uint tSharedCounter = 0;

                    for ( uint Ik = 0; Ik < tAllCoefIds.numel(); Ik++ )
                    {
                        if ( tAllCoefIds( Ik ) != gNoID && tAllCoefOwners( Ik ) == par_rank() )
                        {
                            tOwnedCounter++;
                        }
                        else if ( tAllCoefIds( Ik ) != gNoID )
                        {
                            tSharedCounter++;
                        }
                    }

                    Matrix< DDUMat > tOwnedCoefficients( tOwnedCounter, 1 );
                    Matrix< DDUMat > tSharedCoefficients( tSharedCounter, 1 );

                    tOwnedCounter  = 0;
                    tSharedCounter = 0;

                    for ( uint Ik = 0; Ik < tAllCoefIds.numel(); Ik++ )
                    {
                        if ( tAllCoefIds( Ik ) != gNoID && tAllCoefOwners( Ik ) == par_rank() )
                        {
                            tOwnedCoefficients( tOwnedCounter++ ) = Ik;
                        }
                        else if ( tAllCoefIds( Ik ) != gNoID )
                        {
                            tSharedCoefficients( tSharedCounter++ ) = Ik;
                        }
                    }

                    // Sizes of ID vectors
                    uint tNumOwnedADVs          = tOwnedADVIds.length();
                    uint tNumOwnedCoefficients  = tOwnedCoefficients.numel();
                    uint tNumSharedCoefficients = tSharedCoefficients.numel();

                    // Resize ID lists and bounds
                    tOwnedADVIds.resize( tNumOwnedADVs + tNumOwnedCoefficients, 1 );
                    mLowerBounds.resize( tNumOwnedADVs + tNumOwnedCoefficients, 1 );
                    mUpperBounds.resize( tNumOwnedADVs + tNumOwnedCoefficients, 1 );
                    mOwnedijklIds.resize( tNumOwnedADVs + tNumOwnedCoefficients, 1 );
                    tSharedADVIds( tFieldIndex ).resize( tNumOwnedCoefficients + tNumSharedCoefficients, 1 );
                    tSharedCoefficientIndices( tFieldIndex ) = tOwnedCoefficients;

                    // Add owned coefficients to lists
                    for ( uint tOwnedCoefficient = 0; tOwnedCoefficient < tNumOwnedCoefficients; tOwnedCoefficient++ )
                    {
                        // HMR coeffs are not neccesarily consecutive. Therefore this is a really hacky implementation
                        sint tADVId = tOffsetID
                                    + tMesh->get_glb_entity_id_from_entity_loc_index(
                                            tOwnedCoefficients( tOwnedCoefficient ),
                                            aADVEntityRank,
                                            tDiscretizationMeshIndex );

                        MORIS_ASSERT( tADVId - tOffsetID == tAllCoefIds( tOwnedCoefficients( tOwnedCoefficient ) ), "check if this is a problem" );

                        tOwnedADVIds( tNumOwnedADVs + tOwnedCoefficient ) = tADVId;
                        mLowerBounds( tNumOwnedADVs + tOwnedCoefficient ) = tFields( tFieldIndex )->get_discretization_lower_bound();
                        mUpperBounds( tNumOwnedADVs + tOwnedCoefficient ) = tFields( tFieldIndex )->get_discretization_upper_bound();

                        if ( tMesh->get_mesh_type() == mtk::MeshType::HMR )
                        {
                            mOwnedijklIds( tNumOwnedADVs + tOwnedCoefficient ) = tAllCoefijklIDs( tOwnedCoefficients( tOwnedCoefficient ) );
                        }

                        tSharedADVIds( tFieldIndex )( tOwnedCoefficient ) = tADVId;
                    }

                    // Add shared coefficients to field-specific list
                    tSharedCoefficientIndices( tFieldIndex ).resize( tNumOwnedCoefficients + tNumSharedCoefficients, 1 );
                    for ( uint tSharedCoefficient = 0; tSharedCoefficient < tNumSharedCoefficients; tSharedCoefficient++ )
                    {
                        // HMR coeffs are not necessarily consecutive. Therefore this is a really hacky implementation
                        sint tADVId = tOffsetID
                                    + tMesh->get_glb_entity_id_from_entity_loc_index(
                                            tSharedCoefficients( tSharedCoefficient ),
                                            aADVEntityRank,
                                            tDiscretizationMeshIndex );

                        MORIS_ASSERT( tADVId - tOffsetID == tAllCoefIds( tSharedCoefficients( tSharedCoefficient ) ), "check if this is a problem" );

                        tSharedCoefficientIndices( tFieldIndex )( tNumOwnedCoefficients + tSharedCoefficient ) = tSharedCoefficients( tSharedCoefficient );
                        tSharedADVIds( tFieldIndex )( tNumOwnedCoefficients + tSharedCoefficient )             = tADVId;
                    }

                    // Update offset based on maximum ID
                    tAllOffsetIDs( tFieldIndex ) = tOffsetID;
                    tOffsetID += tMesh->get_max_entity_id( aADVEntityRank, tDiscretizationMeshIndex );
                }
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

            // Create map for distributed vectors
            sol::Dist_Map* tOwnedADVMap     = tDistributedFactory.create_map( tOwnedADVIds );
            sol::Dist_Map* tPrimitiveADVMap = tDistributedFactory.create_map( tPrimitiveADVIds );

            // Create vectors
            sol::Dist_Vector* tNewOwnedADVs = tDistributedFactory.create_vector( tOwnedADVMap, 1, false, true );
            mPrimitiveADVs                  = tDistributedFactory.create_vector( tPrimitiveADVMap, 1, false, true );

            // Assign primitive ADVs
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
                for ( uint tFieldIndex = 0; tFieldIndex < tFields.size(); tFieldIndex++ )
                {
                    tFields( tFieldIndex )->set_advs( mPrimitiveADVs );
                }
            }

            MORIS_LOG_INFO( "Time to create owned ADVs: %f sec", ( moris::real )( clock() - tStart_Create_Owned_ADVs ) / CLOCKS_PER_SEC );

            //----------------------------------------//
            // Convert to B-spline fields             //
            //----------------------------------------//

            clock_t tStart_Convert_to_Bspline_Fields = clock();

            // FIXME this hole section is super hacky and limiting. has to be rewritten from scratch.
            moris::map< std::string, uint > tFieldNameToIndexMap;
            for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                std::string tLabel             = aFields( Ik )->get_label();
                tFieldNameToIndexMap[ tLabel ] = Ik;
            }

            // Loop to find B-spline geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < mGeometries.size(); tGeometryIndex++ )
            {
                // Shape sensitivities logic
                mShapeSensitivities = ( mShapeSensitivities or mGeometries( tGeometryIndex )->depends_on_advs() );

                // Convert to B-spline field
                if ( mGeometries( tGeometryIndex )->intended_discretization() )
                {
                    // Always have shape sensitivities if B-spline field
                    mShapeSensitivities = true;

                    std::string tGeoName = mGeometries( tGeometryIndex )->get_name();

                    if ( not tFieldNameToIndexMap.key_exists( tGeoName ) )
                    {
                        // Create B-spline geometry FIXME Overwriting the given geometry is obviously wrong
                        mGeometries( tGeometryIndex ) = std::make_shared< BSpline_Geometry >(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices( tGeometryIndex ),
                                tSharedADVIds( tGeometryIndex ),
                                tAllOffsetIDs( tGeometryIndex ),
                                aMeshPair,
                                mGeometries( tGeometryIndex ) );
                    }
                    else
                    {
                        uint tMTKFieldIndex = tFieldNameToIndexMap.find( tGeoName );

                        // Create B-spline geometry
                        mGeometries( tGeometryIndex ) = std::make_shared< BSpline_Geometry >(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices( tGeometryIndex ),
                                tSharedADVIds( tGeometryIndex ),
                                tAllOffsetIDs( tGeometryIndex ),
                                aMeshPair,
                                mGeometries( tGeometryIndex ),
                                aFields( tMTKFieldIndex ) );
                    }
                }
                // Store field values if needed. FIXME this is obviously wrong that GEN sets it's own mesh
                else if ( mGeometries( tGeometryIndex )->intended_storage() )
                {
                    // Create stored geometry FIXME this stored geometry stuff is kind of hacky
                    mGeometries( tGeometryIndex ) = std::make_shared< Stored_Geometry >(
                            tMesh,
                            mGeometries( tGeometryIndex ) );

                    mGeometries( tGeometryIndex )->unlock_field();
                    mGeometries( tGeometryIndex )->set_mesh_pair( aMeshPair );
                }
                else
                {
                    // Every Field needs a mesh. FIXME setting the mesh here is to late
                    mGeometries( tGeometryIndex )->unlock_field();
                    mGeometries( tGeometryIndex )->set_mesh_pair( aMeshPair );
                    mGeometries( tGeometryIndex )->set_num_original_nodes( aMeshPair.get_interpolation_mesh()->get_num_nodes() );
                }
            }

            // Loop to find B-spline properties
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                // Convert to B-spline field
                if ( mProperties( tPropertyIndex )->intended_discretization() )
                {
                    // Always have shape sensitivities if B-spline field
                    mShapeSensitivities = true;

                    std::string tPropName = mProperties( tPropertyIndex )->get_name();

                    if ( not tFieldNameToIndexMap.key_exists( tPropName ) )
                    {
                        // Create B-spline property
                        mProperties( tPropertyIndex ) = std::make_shared< BSpline_Property >(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices( mGeometries.size() + tPropertyIndex ),
                                tSharedADVIds( mGeometries.size() + tPropertyIndex ),
                                tAllOffsetIDs( mGeometries.size() + tPropertyIndex ),
                                aMeshPair,
                                mProperties( tPropertyIndex ) );
                    }
                    else
                    {
                        uint tMTKFieldIndex = tFieldNameToIndexMap.find( tPropName );

                        // Create B-spline property
                        mProperties( tPropertyIndex ) = std::make_shared< BSpline_Property >(
                                tNewOwnedADVs,
                                tSharedCoefficientIndices( mGeometries.size() + tPropertyIndex ),
                                tSharedADVIds( mGeometries.size() + tPropertyIndex ),
                                tAllOffsetIDs( mGeometries.size() + tPropertyIndex ),
                                aMeshPair,
                                mProperties( tPropertyIndex ),
                                aFields( tMTKFieldIndex ) );
                    }
                }
            }

            // Update dependencies
            std::copy( mGeometries.begin(), mGeometries.end(), tFields.begin() );
            std::copy( mProperties.begin(), mProperties.end(), tFields.begin() + mGeometries.size() );
            for ( uint tPropertyIndex = 0; tPropertyIndex < mProperties.size(); tPropertyIndex++ )
            {
                mProperties( tPropertyIndex )->update_dependencies( tFields );
            }

            // Save new owned ADVs
            mOwnedADVs = tNewOwnedADVs;

            MORIS_LOG_INFO( "Time to convert to Bspline fields: %f sec", ( moris::real )( clock() - tStart_Convert_to_Bspline_Fields ) / CLOCKS_PER_SEC );

            //----------------------------------------//
            // Communicate all ADV IDs to processor 0 //
            //----------------------------------------//

            clock_t tStart_Communicate_ADV_IDs = clock();

            // Sending mats
            Cell< Matrix< DDSMat > > tSendingIDs( 0 );
            Cell< Matrix< DDRMat > > tSendingLowerBounds( 0 );
            Cell< Matrix< DDRMat > > tSendingUpperBounds( 0 );
            Cell< Matrix< DDSMat > > tSendingijklIDs( 0 );

            // Receiving mats
            Cell< Matrix< DDSMat > > tReceivingIDs( 0 );
            Cell< Matrix< DDRMat > > tReceivingLowerBounds( 0 );
            Cell< Matrix< DDRMat > > tReceivingUpperBounds( 0 );
            Cell< Matrix< DDSMat > > tReceivingjklIDs( 0 );

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
                tSendingijklIDs     = { mOwnedijklIds };
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
                    mFullijklIDs = mOwnedijklIds;
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

            // debug: this needs to go
            // Allocate proximity data
            this->setup_initial_geometric_proximities( aMesh );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::communicate_missing_owned_coefficients(
                mtk::Mesh_Pair&  aMeshPair,
                Matrix< IdMat >& aAllCoefIds,
                Matrix< IdMat >& aAllCoefOwners,
                Matrix< IdMat >& aAllCoefijklIds,
                Cell< uint >&    aNumCoeff,
                uint             aFieldIndex,
                uint             aDiscretizationMeshIndex,
                mtk::MeshType    aMeshType )
        {

            Matrix< IdMat > tCommTable = aMeshPair.get_interpolation_mesh()->get_communication_table();

            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap( tCommTable.max() + 1, 1, -1 );

            uint tNumCommProcs = tCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( tCommTable( Ik ) ) = Ik;
            }

            moris::Cell< Matrix< IdMat > > tSharedCoeffsPosGlobal( tNumCommProcs );
            moris::Cell< Matrix< IdMat > > tSharedCoeffsijklIdGlobal( tNumCommProcs );

            // Set Mat to store number of shared coeffs per processor
            Matrix< DDUMat > tNumSharedCoeffsPerProc( tNumCommProcs, 1, 0 );

            // Count number of coeffs per proc which have to be communicated
            for ( uint Ib = 0; Ib < aAllCoefIds.numel(); Ib++ )
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
            for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                if ( tNumSharedCoeffsPerProc( Ik ) != 0 )
                {
                    tSharedCoeffsPosGlobal( Ik ).set_size( tNumSharedCoeffsPerProc( Ik ), 1 );
                    tSharedCoeffsijklIdGlobal( Ik ).set_size( tNumSharedCoeffsPerProc( Ik ), 1 );
                }
            }

            // Temporary Mat to add external coeffs ids at the next spot in the matrix which will be communicated
            Matrix< DDUMat > tShredCoeffPosPerProc( tNumCommProcs, 1, 0 );

            // Loop over coeffs per type
            for ( uint Ia = 0; Ia < aAllCoefIds.numel(); Ia++ )
            {
                // Check if coeffs at this position is not NULL
                if ( aAllCoefIds( Ia ) != gNoID && aAllCoefOwners( Ia ) != par_rank() )
                {
                    // Get owning processor
                    uint tProcID = aAllCoefOwners( Ia );

                    moris::sint tProcIdPos = tCommTableMap( tProcID );

                    // Add owning processor id to moris::Mat
                    tSharedCoeffsPosGlobal( tProcIdPos )( tShredCoeffPosPerProc( tProcIdPos ) ) =
                            aAllCoefIds( Ia );

                    if ( aMeshType == mtk::MeshType::HMR )
                    {
                        tSharedCoeffsijklIdGlobal( tProcIdPos )( tShredCoeffPosPerProc( tProcIdPos ) ) =
                                aAllCoefijklIds( Ia );
                    }

                    tShredCoeffPosPerProc( tProcIdPos )++;
                }
            }

            // receiving list
            moris::Cell< Matrix< IdMat > > tMatsToReceive;
            moris::Cell< Matrix< IdMat > > tMatsToReceiveijklID;

            barrier();

            // Communicate position of shared adofs to the owning processor
            communicate_mats(
                    tCommTable,
                    tSharedCoeffsPosGlobal,
                    tMatsToReceive );

            barrier();

            if ( aMeshType == mtk::MeshType::HMR )
            {
                communicate_mats(
                        tCommTable,
                        tSharedCoeffsijklIdGlobal,
                        tMatsToReceiveijklID );

                MORIS_ASSERT( tMatsToReceiveijklID.size() == tMatsToReceive.size(), "size must be the same" );
            }

            map< moris_id, moris_index > tCoeffGlobalToLocalMap;
            aMeshPair.get_interpolation_mesh()->get_adof_map(
                    aDiscretizationMeshIndex,
                    tCoeffGlobalToLocalMap );

            // Loop over all Mats set dummy owned coeffs
            for ( uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                {
                    // Get owned coeff Index
                    moris_id    tID            = tMatsToReceive( Ik )( Ii );
                    moris_index tLocalCoeffInd = tCoeffGlobalToLocalMap.find( tID );

                    if ( aAllCoefIds( tLocalCoeffInd ) == gNoID )
                    {
                        aAllCoefIds( tLocalCoeffInd )    = tID;
                        aAllCoefOwners( tLocalCoeffInd ) = par_rank();

                        if ( aMeshType == mtk::MeshType::HMR )
                        {
                            aAllCoefijklIds( tLocalCoeffInd ) = tMatsToReceiveijklID( Ik )( Ii );
                        }

                        aNumCoeff( aFieldIndex )++;
                    }

                    MORIS_ASSERT( aAllCoefIds( tLocalCoeffInd ) == tID,
                            "communicate_missing_owned_coefficients( ), coefficient IDs are not parallel consistent" );

                    if ( aMeshType == mtk::MeshType::HMR )
                    {
                        MORIS_ASSERT( aAllCoefijklIds( tLocalCoeffInd ) == tMatsToReceiveijklID( Ik )( Ii ),
                                "communicate_missing_owned_coefficients( ), coefficient ijkl IDs are not parallel consistent" );
                    }
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
        Geometry_Engine::create_interpolation_pdvs(
                mtk::Interpolation_Mesh*         aInterpolationMesh,
                mtk::Integration_Mesh*           aIntegrationMesh,
                Cell< Cell< Cell< PDV_Type > > > aPdvTypes )
        {
            // Tracer
            Tracer tTracer( "GEN", "Create interpolation PDV hosts" );

            // Get information from integration mesh
            uint tNumSets = aPdvTypes.size();

            // Size node information cells
            Cell< Cell< uint > >     tNodeIndicesPerSet( tNumSets );
            Cell< Cell< sint > >     tNodeIdsPerSet( tNumSets );
            Cell< Cell< uint > >     tNodeOwnersPerSet( tNumSets );
            Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( tNumSets );

            // Get communication table and map
            Matrix< IdMat >  tCommTable             = aInterpolationMesh->get_communication_table();
            Cell< moris_id > tCommunicationTableMap = build_communication_table_map( tCommTable );

            // TODO change over to just use a cell to begin with
            Cell< moris_index > tCommunicationTable( tCommTable.length() );
            for ( uint iCommTableIndex = 0; iCommTableIndex < tCommunicationTable.size(); iCommTableIndex++ )
            {
                tCommunicationTable( iCommTableIndex ) = tCommTable( iCommTableIndex );
            }

            // Loop through sets in integration mesh
            for ( uint iMeshSetIndex = 0; iMeshSetIndex < tNumSets; iMeshSetIndex++ )
            {
                // Determine number of nodes if there are PDVs on this set
                if ( aPdvTypes( iMeshSetIndex ).size() > 0 )
                {
                    // Get set pointer
                    mtk::Set* tSet = aIntegrationMesh->get_set_by_index( iMeshSetIndex );

                    // Select sides of interpolation cells to get info from
                    Cell< mtk::Leader_Follower > tSetSides = mtk::get_leader_follower( tSet->get_set_type() );

                    // Get number of clusters on set
                    uint tNumberOfClusters = tSet->get_num_clusters_on_set();

                    // Number of nodes on this set
                    uint tNumberOfNodesInSet = 0;

                    // Number of shared nodes on this set per proc
                    Cell< uint > tNumSharedNodesPerProc( tCommunicationTable.size() );

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
                    Cell< uint > tNumOwnedNodesPerProc( tCommunicationTable.size() );
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
                        Cell< Cell< sint > > tSharedNodeIdsOnSet( tCommunicationTable.size() );
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
                        Cell< Cell< sint > > tOwnedNodeIdsOnSet( tCommunicationTable.size() );

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
            for ( uint iPropertyIndex = 0; iPropertyIndex < mProperties.size(); iPropertyIndex++ )
            {
                // Check if this is an interpolation PDV
                MORIS_ERROR( mProperties( iPropertyIndex )->is_interpolation_pdv(),
                        "Assignment of PDVs is only supported with an interpolation mesh right now." );

                // Loop through sets in integration mesh TODO this can be simplified once a property can set its own (total) mesh set indices, see TODO in create_pdvs()
                for ( uint iMeshSetIndex = 0; iMeshSetIndex < tNumSets; iMeshSetIndex++ )
                {
                    // Check with PDVs on this set
                    for ( uint iPdvTypeIndex = 0; iPdvTypeIndex < aPdvTypes( iMeshSetIndex ).size(); iPdvTypeIndex++ )
                    {
                        PDV_Type tPdvType = mProperties( iPropertyIndex )->get_pdv_type();
                        if ( tPdvType == aPdvTypes( iMeshSetIndex )( iPdvTypeIndex )( 0 ) )
                        {
                            for ( uint iNodeInSet = 0; iNodeInSet < tNodeIndicesPerSet( iMeshSetIndex ).size(); iNodeInSet++ )
                            {
                                mPDVHostManager.create_interpolation_pdv( tNodeIndicesPerSet( iMeshSetIndex )( iNodeInSet ), tPdvType, mProperties( iPropertyIndex ) );
                            }
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::set_integration_pdv_types( mtk::Integration_Mesh* aIntegrationMesh )
        {
            // Tracer
            Tracer tTracer( "GEN", "Set integration PDV types" );

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

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry_Engine::setup_initial_geometric_proximities( mtk::Interpolation_Mesh* aMesh )
        {
            // Tracer
            Tracer tTracer( "GEN", "Setup initial geometric proximities" );

            mVertexGeometricProximity =
                    Cell< Geometric_Proximity >( aMesh->get_num_nodes(), Geometric_Proximity( mGeometries.size() ) );

            // iterate through vertices then geometries
            for ( uint iVertex = 0; iVertex < aMesh->get_num_nodes(); iVertex++ )
            {
                Matrix< DDRMat > tCoords = aMesh->get_node_coordinate( moris_index( iVertex ) );

                mVertexGeometricProximity( iVertex ).mAssociatedVertexIndex = (moris_index)iVertex;

                for ( uint iGeometryIndex = 0; iGeometryIndex < mGeometries.size(); iGeometryIndex++ )
                {
                    real tVertGeomVal = mGeometries( iGeometryIndex )->get_field_value( iVertex, tCoords );

                    moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                    mVertexGeometricProximity( iVertex ).set_geometric_proximity( tGeomProxIndex, iGeometryIndex );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index
        Geometry_Engine::get_geometric_proximity_index( real const & aGeometricVal )
        {
            // get current geometries level set info
            real tIsocontourThreshold = mGeometries( mActiveGeometryIndex )->get_isocontour_threshold();
            real tIsocontourTolerance = mGeometries( mActiveGeometryIndex )->get_isocontour_tolerance();

            // initialize index to 1, i.e. vertex is on interface
            moris_index tGeometricProxIndex = 1;

            if ( aGeometricVal - tIsocontourThreshold < -tIsocontourTolerance )
            {
                tGeometricProxIndex = 0;
            }
            else if ( aGeometricVal - tIsocontourThreshold > tIsocontourTolerance )
            {
                tGeometricProxIndex = 2;
            }

            return tGeometricProxIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Geometry_Engine::check_queued_intersection_geometric_proximity_index(
                moris_index const & aProximIndex,
                moris_index const & aGeomIndex )
        {
            // parent vertex
            moris_index tParentVertexIndex0 = mQueuedIntersectionNode->get_first_parent_node_index();
            moris_index tParentVertexIndex1 = mQueuedIntersectionNode->get_second_parent_node_index();

            // parent vertex proximity wrt aGeomIndex
            moris_index tParentProx0 = mVertexGeometricProximity( tParentVertexIndex0 ).get_geometric_proximity( aGeomIndex );
            moris_index tParentProx1 = mVertexGeometricProximity( tParentVertexIndex1 ).get_geometric_proximity( aGeomIndex );

            // 0 - G(x) < threshold:  left of interface
            // 1 - G(x) == threshold: on interface
            // 2 - G(x) > threshold:  right of interface
            // add them together
            moris_index tSum = tParentProx0 + tParentProx1;

            // proximity value
            switch ( tSum )
            {
                case 0:    // both parents are left of interface -> child has to be left of interface
                {
                    return aProximIndex == 0;
                    break;
                }
                case 1:    // one parent is left of, one parent is on interface -> child has to be left of or on the interface
                {
                    return aProximIndex == 0 || aProximIndex == 1;
                    break;
                }
                case 2:    // one parent is left the other right of interface -> child position cannot be determined yet
                {
                    return true;
                    break;
                }
                case 3:    // one parent is on interface the other is right of interface -> child has to be right of or on the interface
                {
                    return aProximIndex == 2 || aProximIndex == 1;
                    break;
                }
                case 4:    // both right of interface -> child has to be  right of interface
                {
                    return aProximIndex == 2;
                    break;
                }
                default:
                {
                    MORIS_ASSERT( 0, "Proximity determination failed." );
                    return false;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index
        Geometry_Engine::get_queued_intersection_geometric_proximity_index( moris_index const & aGeomIndex )
        {
            // parent vertex
            moris_index tParentVertexIndex0 = mQueuedIntersectionNode->get_first_parent_node_index();
            moris_index tParentVertexIndex1 = mQueuedIntersectionNode->get_second_parent_node_index();

            // parent vertex proximity wrt aGeomIndex
            moris_index tParentProx0 = mVertexGeometricProximity( tParentVertexIndex0 ).get_geometric_proximity( aGeomIndex );
            moris_index tParentProx1 = mVertexGeometricProximity( tParentVertexIndex1 ).get_geometric_proximity( aGeomIndex );

            // 0 - G(x) < threshold:  left of interface
            // 1 - G(x) == threshold: on interface
            // 2 - G(x) > threshold:  right of interface
            // add them together
            moris_index tSum = tParentProx0 + tParentProx1;

            // proximity value
            switch ( tSum )
            {
                case 0:    // both parents are left of interface -> child is left of interface
                {
                    return 0;
                    break;
                }
                case 1:    // one parent is left of, one parent is on interface -> child is left of interface
                {
                    return 0;
                    break;
                }
                case 2:    // one parent is left the other right of interface -> child is on interface (correct?)
                {
                    return 1;
                    break;
                }
                case 3:    // one parent is on interface the other is right of interface -> child is right of interface
                {
                    return 2;
                    break;
                }
                case 4:    // both right of interface -> child is right of interface
                {
                    return 2;
                    break;
                }
                default:
                {
                    MORIS_ASSERT( 0, "Proximity determination failed." );
                    return MORIS_INDEX_MAX;
                }
            }
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
                    // FIXME: need to use level set value of child node
                    real tVertGeomVal = mGeometries( tGeometryIndex )->get_field_value( aNodeIndex, mQueuedIntersectionNode->get_global_coordinates() );

                    // compute proximity index
                    moris_index tGeomProxIndex = this->get_geometric_proximity_index( tVertGeomVal );

                    // check that tVertGeomVal is consistent with parent nodes
                    MORIS_ERROR( check_queued_intersection_geometric_proximity_index( tGeomProxIndex, tGeometryIndex ),
                            "Geometry_Engine::admit_queued_intersection_geometric_proximity - inconsistent proximity value." );

                    // save proximity index
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
            // Tracer
            Tracer tTracer( "GEN", "Create phase table" );

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

    }    // namespace ge
}    // namespace moris
