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
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MDL_Model.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"

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
        // mPerformerManager->mGENPerformer.resize( 1 );
        // mPerformerManager->mXTKPerformer.resize( 1 );
        mPerformerManager->mMTKPerformer.resize( 1 );
        mPerformerManager->mMDLPerformer.resize( 1 );

        // load the STK parameter list
        Module_Parameter_Lists tSTKParameterList = aPerformerManager->mLibrary->get_parameters_for_module( Parameter_List_Type::STK );

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

        // Stage 1: MDL perform ---------------------------------------------------------------------

        mPerformerManager->mMDLPerformer( 0 )->initialize();

        // Build MDL components and solve
        mPerformerManager->mMDLPerformer( 0 )->perform();

        Vector< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

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
    Workflow_STK_FEM::create_stk( Module_Parameter_Lists& aParameterLists )
    {
        Tracer            tTracer( "STK", "Mesh", "InitializeMesh" );
        std::string       tMeshFile     = aParameterLists( 0 )( 0 ).get< std::string >( "input_file" );
        mtk::MtkMeshData* tSuppMeshData = nullptr;

        mtk::Cell_Cluster_Input* tCellClusterData = nullptr;

        // construct the meshes
        mIpMesh = std::make_shared< mtk::Interpolation_Mesh_STK >( tMeshFile, tSuppMeshData, true );

        mIgMesh = std::make_shared< mtk::Integration_Mesh_STK >( *mIpMesh, tCellClusterData );

        mPerformerManager->mMTKPerformer( 0 )->register_mesh_pair( mIpMesh.get(), mIgMesh.get() );

        if ( aParameterLists( 0 )( 0 ).get< bool >( "periodic_workspace" ) )
        {
            std::cout << " Periodic BCs " << '\n';
            mtk::Periodic_Boundary_Condition_Helper tPBCHelper( mPerformerManager->mMTKPerformer( 0 ), 0, aParameterLists( 0 )( 0 ) );
            tPBCHelper.setup_periodic_boundary_conditions();

            // call some function in MTK to setup periodic boundary conditions
        }

        // output the mesh
        moris::mtk::Writer_Exodus tWriter( mIgMesh.get() );
        tWriter.write_mesh( "", "stk_fem_mesh.exo", "", "temp.exo" );
        tWriter.close_file();
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::wrk
