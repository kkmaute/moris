/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_STK_FEM.cpp
 *
 */

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_STK_FEM.hpp"

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Integration_Surface_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"
// #include "cl_GEN_Geometry_Engine.hpp" // brendan delete?
// #include "cl_XTK_Model.hpp" // brendan delete
#include "cl_MDL_Model.hpp"
#include "cl_MSI_QI_Manager_STK.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"

#include "fn_norm.hpp"

namespace moris::wrk
{
    //--------------------------------------------------------------------------------------------------------------

    Workflow_STK_FEM::Workflow_STK_FEM( wrk::Performer_Manager* aPerformerManager )
            : Workflow( aPerformerManager )
    {
        // log & trace this function
        Tracer tTracer( "WRK", "Workflow_STK_FEM", "Create" );
        MORIS_LOG_SPEC( "Par_Rank", par_rank() );
        MORIS_LOG_SPEC( "Par_Size", par_size() );

        // Performer set for this workflow
        mPerformerManager->mMTKPerformer.resize( 1 );
        mPerformerManager->mMDLPerformer.resize( 1 );

        // load the STK parameter list
        const Submodule_Parameter_Lists& tSTKParameterList = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::STK )( 0 );
        std::cout << "Attempting to get file name from stk: " << tSTKParameterList( 0 ).get< std::string >( "input_file" ) << "\n";
        mRequestedQIs = tSTKParameterList( 0 ).get_vector< std::string >( "requested_QIs" );

        // load the meshes
        mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();
        this->create_stk( tSTKParameterList );

        // create MDL performer
        mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary, 0 );

        // Set performer to MDL
        mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 0 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_STK_FEM::initialize(
            Vector< real >&  aADVs,
            Vector< real >&  aLowerBounds,
            Vector< real >&  aUpperBounds,
            Matrix< IdMat >& aIjklIDs )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real >
    Workflow_STK_FEM::perform( Vector< real >& aNewADVs )
    {
        // Stage 0: Create QI manager and register with MDL performer
        std::shared_ptr< MSI::QI_Manager_STK > tQIManager = std::make_shared< MSI::QI_Manager_STK >( mIgMesh, mRequestedQIs );

        mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface( tQIManager );

        // Stage 1: Compute SQIs (QIs on the mesh itself)
        this->compute_sqis( mPerformerManager->mLibrary );

        // Stage 2: MDL perform ---------------------------------------------------------------------
        mPerformerManager->mMDLPerformer( 0 )->initialize();

        // Build MDL components and solve
        mPerformerManager->mMDLPerformer( 0 )->perform();

        Vector< Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

        // Communicate IQIs
        for ( uint iIQIIndex = 0; iIQIIndex < tVal.size(); iIQIIndex++ )
        {
            tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
        }

        Vector< real > tCriteria( tVal.size(), 0.0 );

        for ( uint iCriteriaIndex = 0; iCriteriaIndex < tVal.size(); iCriteriaIndex++ )
        {
            tCriteria( iCriteriaIndex ) = tVal( iCriteriaIndex )( 0 );
        }

        return tCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Workflow_STK_FEM::compute_dcriteria_dadv()
    {
        return { {} };
    }

    void
    Workflow_STK_FEM::create_stk( const Submodule_Parameter_Lists& aParameterLists )
    {
        Tracer            tTracer( "STK", "Mesh", "InitializeMesh" );
        std::string       tMeshFile     = aParameterLists( 0 ).get< std::string >( "input_file" );
        mtk::MtkMeshData* tSuppMeshData = nullptr;

        mtk::Cell_Cluster_Input* tCellClusterData = nullptr;

        // construct the meshes
        mIpMesh = std::make_shared< mtk::Interpolation_Mesh_STK >( tMeshFile, tSuppMeshData, true );

        mIgMesh = std::make_shared< mtk::Integration_Mesh_STK >( *mIpMesh, tCellClusterData );

        mPerformerManager->mMTKPerformer( 0 )->register_mesh_pair( mIpMesh.get(), mIgMesh.get() );

        if ( aParameterLists( 0 ).get< bool >( "periodic_workspace" ) )
        {
            std::cout << " Periodic BCs " << '\n';
            mtk::Periodic_Boundary_Condition_Helper tPBCHelper( mPerformerManager->mMTKPerformer( 0 ), 0, aParameterLists( 0 ) );
            tPBCHelper.setup_periodic_boundary_conditions();

            // call some function in MTK to setup periodic boundary conditions
        }

        // output the mesh
        moris::mtk::Writer_Exodus tWriter( mIgMesh.get() );
        tWriter.write_mesh( "", "stk_fem_mesh.exo", "", "temp.exo" );
        tWriter.close_file();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Workflow_STK_FEM::compute_sqis( const std::shared_ptr< Library_IO > aLibrary )
    {
        Tracer tTracer( "STK", "Compute SQIs" );

        Submodule_Parameter_Lists tSQIParameterLists = aLibrary->get_parameters_for_module( Module_Type::STK )( 1 );

        // Get the design criteria manager from the geometry engine (PDV_Host_Manager/Design_Variable_Interface) brendan fix naming
        std::shared_ptr< MSI::Design_Variable_Interface > tDesignCriteriaManager = mPerformerManager->mMDLPerformer( 0 )->get_design_variable_interface();

        uint tNumSQIs = tSQIParameterLists.size();

        if ( tNumSQIs > 0 )    // todo brendan be more elegant
        {
            // Create a dist vector to store the sensitivities
            sol::Matrix_Vector_Factory tDistFactory;
            sol::Dist_Map*             tMap      = tDistFactory.create_map( tDesignCriteriaManager->get_my_local_global_map() );
            sol::Dist_Vector*          tdSQIdPDV = tDistFactory.create_vector( tMap, tNumSQIs );

            // Loop through requested SQIs and compute them
            for ( uint iSQI = 0; iSQI < tNumSQIs; iSQI++ )
            {
                // Get the SQI name and type
                std::string           tSQIName  = tSQIParameterLists( iSQI ).get< std::string >( "SQI_name" );
                mtk::QI_Type          tSQIType  = static_cast< mtk::QI_Type >( tSQIParameterLists( iSQI ).get< uint >( "SQI_type" ) );
                Vector< std::string > tMeshSets = tSQIParameterLists( iSQI ).get_vector< std::string >( "mesh_set_names" );

                // Create a surface mesh from the IG mesh for surface SQIs
                // FIXME brendan construct from phase names
                mtk::Integration_Surface_Mesh_Data tSurfaceMeshData( mIgMesh.get(), tMeshSets );
                mtk::Integration_Surface_Mesh      tSurfaceMesh( tSurfaceMeshData );

                tSurfaceMesh.write_to_file( "integ_mesh_iter_" + std::to_string( gLogger.get_opt_iteration() ) + ".obj" );

                // Get the IG to PDV ID map for this surface mesh
                Vector< Vector< moris_index > > tPDVIDs;
                Vector< gen::PDV_Type >         tPDVTypes = tSurfaceMesh.get_spatial_dimension() == 2 ? Vector< gen::PDV_Type >( { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE } )
                                                                                                      : Vector< gen::PDV_Type >( { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE, gen::PDV_Type::Z_COORDINATE } );
                tDesignCriteriaManager->get_ig_dv_ids_for_type_and_ind( tSurfaceMeshData.mLocalToGlobalVertexIndex, tPDVTypes, tPDVIDs );

                // Determine the SQI type and compute it
                real tSQIValue = MORIS_REAL_MAX;
                switch ( tSQIType )
                {
                    case ( mtk::QI_Type::VOLUME ):
                    {
                        tSQIValue = tSurfaceMesh.compute_volume();
                        tSurfaceMesh.compute_QI_sensitivities( tSQIType, tPDVIDs, tdSQIdPDV, iSQI );
                        break;
                    }
                    case ( mtk::QI_Type::RAYCAST_SHAPE_DIAMETER ):
                    {
                        real                                             tConeAngle      = tSQIParameterLists( iSQI ).get< real >( "cone_angle" );
                        uint                                             tNumPolarRays   = static_cast< uint >( tSQIParameterLists( iSQI ).get< moris_index >( "number_of_polar_rays" ) );
                        uint                                             tNumAzimuthRays = static_cast< uint >( tSQIParameterLists( iSQI ).get< moris_index >( "number_of_azimuth_rays" ) );
                        std::unique_ptr< mtk::Agglomeration_Parameters > tAgglom         = mtk::create_agglomeration_parameters( tSQIParameterLists( iSQI ), aLibrary );

                        tSQIValue = tSurfaceMesh.compute_global_shape_diameter_raycast(
                                *tAgglom,
                                tConeAngle,
                                tNumPolarRays,
                                tNumAzimuthRays );

                        tSurfaceMesh.compute_QI_sensitivities(
                                tSQIType,
                                tPDVIDs,
                                tdSQIdPDV,
                                iSQI );
                        break;
                    }
                    case ( mtk::QI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER ):
                    {
                        std::unique_ptr< mtk::Agglomeration_Parameters > tAgglom = mtk::create_agglomeration_parameters( tSQIParameterLists( iSQI ), aLibrary );

                        tSQIValue = tSurfaceMesh.compute_global_shape_diameter_inscribed_circle(
                                *tAgglom,
                                tSQIParameterLists( iSQI ).get< real >( "cone_angle" ),
                                tSQIParameterLists( iSQI ).get< real >( "minimum_relative_chord_length" ),
                                tSQIParameterLists( iSQI ).get< moris_index >( "number_of_samples" ) );

                        tSurfaceMesh.compute_QI_sensitivities(
                                tSQIType,
                                tPDVIDs,
                                tdSQIdPDV,
                                iSQI );
                        break;
                    }
                    case ( mtk::QI_Type::SHORTEST_DISTANCE_SHAPE_DIAMETER ):
                    {
                        std::unique_ptr< mtk::Agglomeration_Parameters > tAgglom = mtk::create_agglomeration_parameters( tSQIParameterLists( iSQI ), aLibrary );

                        tSQIValue = tSurfaceMesh.compute_global_shape_diameter_shortest_distance(
                                *tAgglom,
                                tSQIParameterLists( iSQI ).get< real >( "cone_angle" ),
                                tSQIParameterLists( iSQI ).get< moris_index >( "number_of_samples" ) );

                        tSurfaceMesh.compute_QI_sensitivities(
                                tSQIType,
                                tPDVIDs,
                                tdSQIdPDV,
                                iSQI );

                        break;
                    }
                    default:
                        MORIS_ERROR( false, "STK::Model::compute_SQIs - SQI Type not implemented." );
                }

                tDesignCriteriaManager->register_QI( tSQIName, Module_Type::STK, tSQIValue );
            }

            tdSQIdPDV->vector_global_assembly();
            tDesignCriteriaManager->update_QI_sensitivity( Module_Type::STK, tdSQIdPDV );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::wrk
