/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_Model.cpp
 *
 */

#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Map.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"           //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"    //MTK/src

#include "cl_MDL_Model.hpp"

#include <utility>

#include "cl_FEM_Enums.hpp"    //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"
#include "fn_MSI_get_mesh_index_for_dof_type.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_VIS_Output_Manager.hpp"

#include "fn_PRM_MSI_Parameters.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Library_IO.hpp"

// Parameter function
// typedef void ( *Parameter_Function ) ( Vector< Vector< moris::ParameterList > > & aParameterList );

namespace moris::mdl
{
    //------------------------------------------------------------------------------
    Model::Model(
            std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
            const uint                           aBSplineIndex,
            Vector< fem::Set_User_Info >&        aSetInfo,
            const moris_index                    aMeshPairIndex,
            const bool                           aUseMultigrid )
            : mMeshManager( std::move( aMeshManager ) )
            , mMeshPairIndex( aMeshPairIndex )
            , mBSplineIndex( aBSplineIndex )
            , mUseMultigrid( aUseMultigrid )
    {
        Tracer tTracer( "FEM", "Model", "Create" );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 0: initialize
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // get pointers to interpolation and integration mesh
        mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create the FEM model
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "Create Equation Model" );

            // build the FEM model from FEM parameter list
            mEquationModel = std::make_shared< fem::FEM_Model >(
                    mMeshManager,
                    mMeshPairIndex,
                    aSetInfo );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create the model solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "Create MSI" );

            // create the MSI parameter list
            moris::Parameter_List tMSIParameters = prm::create_msi_parameter_list();
            tMSIParameters.set( "multigrid", mUseMultigrid );

            if ( tInterpolationMesh->get_mesh_type() == mtk::MeshType::HMR )
            {
                tMSIParameters.set( "L2", (sint)mBSplineIndex );
                tMSIParameters.set( "TEMP", (sint)mBSplineIndex );
            }

            // build the model solver interface from the MSI parameter list
            mModelSolverInterface = new moris::MSI::Model_Solver_Interface(
                    tMSIParameters,
                    mEquationModel,
                    tInterpolationMesh );

            // finalize the model solver interface
            mModelSolverInterface->finalize();

            // finalize the fem sets
            mEquationModel->finalize_equation_sets( mModelSolverInterface );

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: create the solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "Create Solver Interface" );

            mSolverInterface = new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );
        }
    }

    //------------------------------------------------------------------------------
    Model::Model(
            std::shared_ptr< Library_IO > aLibrary,
            const uint                    aBSplineIndex,
            const moris_index             aMeshPairIndex )
            : mMeshPairIndex( aMeshPairIndex )
            , mBSplineIndex( aBSplineIndex )
            , mLibrary( std::move( aLibrary ) )
    {
        Tracer tTracer( "MDL", "Model", "Create" );

        // load the MSI parameter list
        mMSIParameterList = mLibrary->get_parameters_for_module( Parameter_List_Type::MSI );

        // load the SOL parameter list
        mSOLParameterList = mLibrary->get_parameters_for_module( Parameter_List_Type::SOL );

        // load the FEM parameter list
        mFEMParameterList = mLibrary->get_parameters_for_module( Parameter_List_Type::FEM );

        // load the VIS parameter list
        mVISParameterList = mLibrary->get_parameters_for_module( Parameter_List_Type::VIS );
    }

    //------------------------------------------------------------------------------
    Model::Model(
            std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
            const uint                           aBSplineIndex,
            Vector< fem::Set_User_Info >&        aSetInfo,
            MSI::Design_Variable_Interface*      aDesignVariableInterface,
            const moris_index                    aMeshPairIndex,
            const bool                           aUseMultigrid )
            : mMeshManager( std::move( aMeshManager ) )
            , mMeshPairIndex( aMeshPairIndex )
            , mBSplineIndex( aBSplineIndex )
            , mUseMultigrid( aUseMultigrid )
    {
        Tracer tTracer( "MDL", "Model", "Create" );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create the FEM model
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "Create Equation Model" );

            // build the FEM model from FEM parameter list
            mEquationModel = std::make_shared< fem::FEM_Model >(
                    mMeshManager,
                    mMeshPairIndex,
                    aSetInfo,
                    aDesignVariableInterface );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create the model solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // start timer
        {
            Tracer tTracer( "MDL", "Model", "Create MSI" );

            // Get pointers to interpolation and integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

            // create the MSI parameter list
            moris::Parameter_List tMSIParameters = prm::create_msi_parameter_list();
            tMSIParameters.set( "multigrid", mUseMultigrid );

            // create the model solver interface from the MSI parameter list
            mModelSolverInterface = new moris::MSI::Model_Solver_Interface(
                    tMSIParameters,
                    mEquationModel,
                    tInterpolationMesh );

            // finalize the model solver interface
            mModelSolverInterface->finalize();

            // finalize the fem sets
            // FIXME pass in the interface
            mEquationModel->finalize_equation_sets( mModelSolverInterface );

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create Solver Interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // create the solver interface
        {
            Tracer tTracer( "MDL", "Model", "Create Solver Interface" );

            mSolverInterface = new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );
        }
    }

    //------------------------------------------------------------------------------

    Model::~Model()
    {
        // delete SI
        delete mSolverInterface;

        // delete MSI
        delete mModelSolverInterface;

        // delete MSI
        if ( mOutputManager != nullptr && mOutputManagerOwned == true )
        {
            delete mOutputManager;

            mOutputManagerOwned = false;
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::free_memory()
    {
        // delete SI
        delete mSolverInterface;
        mSolverInterface = nullptr;

        // delete MSI
        delete mModelSolverInterface;
        mModelSolverInterface = nullptr;

        // delete MSI
        if ( mOutputManager != nullptr && mOutputManagerOwned == true )
        {
            delete mOutputManager;
            mOutputManager = nullptr;

            mOutputManagerOwned = false;
        }

        mEquationModel->free_memory();
    }

    //------------------------------------------------------------------------------
    void
    Model::set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
    {
        mMeshManager = std::move( aMTKPerformer );
    }

    //------------------------------------------------------------------------------

    void
    Model::initialize()
    {
        Tracer tTracer( "MDL", "Model", "Initialize" );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 0: delete model and interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if ( mModelSolverInterface != nullptr )
        {
            delete mModelSolverInterface;
        }
        if ( mSolverInterface != nullptr )
        {
            delete mSolverInterface;
        }
        if ( mOutputManager != nullptr )
        {
            delete mOutputManager;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create the FEM model
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "CreateEquationModel" );

            // build the FEM model from FEM parameter list
            mEquationModel = std::make_shared< fem::FEM_Model >(
                    mMeshManager,
                    mMeshPairIndex,
                    mFEMParameterList,
                    mDesignVariableInterface );
            // tDofTypeToBsplineMeshIndex );

            // build a map relating the dof types to their respective B-spline mesh indices and pass it to the FEM model
            std::unordered_map< MSI::Dof_Type, moris_index > tDofTypeToBsplineMeshIndex = this->build_dof_type_to_mesh_index();
            mEquationModel->set_dof_type_to_Bspline_mesh_index( tDofTypeToBsplineMeshIndex );

            // pass whether the new ghost sets are used
            mEquationModel->set_use_new_ghost_sets( mUseNewGhostSets );

            // initialize the Equation Model
            mEquationModel->initialize_from_inputfile( mLibrary );

            // set the equation model if using design variables
            if ( mDesignVariableInterface != nullptr )
            {
                mDesignVariableInterface->set_equation_model( mEquationModel );
            }
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create the model solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "CreateMSI" );

            // get pointers to interpolation and integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

            // Does not work with STK
            // MORIS_ERROR( tInterpolationMesh->get_mesh_type() != mtk::MeshType::STK, "Does not work for STK");

            // build the model solver interface

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface(
                    mMSIParameterList( 0 )( 0 ),
                    mEquationModel,
                    tInterpolationMesh );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create the solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        {
            Tracer tTracer( "MDL", "Model", "Create Solver Interface" );

            mSolverInterface = new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create the solver
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            {
                Tracer tTracer( "MDL", "Model", "Create Solver" );

                mSolverWarehouse = std::make_shared< sol::SOL_Warehouse >( mSolverInterface, mLibrary );

                mSolverWarehouse->set_parameterlist( mSOLParameterList );

                mSolverWarehouse->initialize();
            }

            // set the solver warehouse pointer
            mSolverInterface->set_solver_warehouse( mSolverWarehouse );

            // finalize the model solver interface
            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();
        }

        // finalize the fem sets
        mEquationModel->finalize_equation_sets( mModelSolverInterface );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 5: create the output manager
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        mOutputManager = new vis::Output_Manager( mVISParameterList( 0 ) );

        mOutputManagerOwned = true;
    }

    //------------------------------------------------------------------------------
    void
    Model::perform( uint aIndex )
    {
        if ( aIndex == 0 )
        {
            Tracer tTracer( "MDL", "Model", "Perform Forward Analysis" );

            this->perform_forward_analysis();
        }
        else if ( aIndex == 1 )
        {
            Tracer tTracer( "MDL", "Model", "Perform Sensitivity Analysis" );

            this->perform_sensitivity_analysis();
        }
    }

    //------------------------------------------------------------------------------
    void
    Model::perform_forward_analysis()
    {
        mEquationModel->set_is_forward_analysis();

        mEquationModel->initialize_IQIs();

        mSolverWarehouse->get_main_time_solver()->solve();
    }

    //------------------------------------------------------------------------------

    Vector< moris::Matrix< DDRMat > >
    Model::get_IQI_values()
    {
        return mEquationModel->get_IQI_values();
    }

    //------------------------------------------------------------------------------
    void
    Model::perform_sensitivity_analysis()
    {
        // get type of sensitivity analysis
        bool tIsAdjointSensitivityAnalysis = mSolverWarehouse->is_adjoint_sensitivity_analysis();

        // activate sensitivity analysis and set type of sensitivity analysis
        mEquationModel->set_sensitivity_analysis_type( tIsAdjointSensitivityAnalysis );

        // initialize computation of dQdp
        mEquationModel->initialize_explicit_and_implicit_dQIdp();

        // compute sensitivities
        mSolverWarehouse->get_main_time_solver()->solve_sensitivity();
    }

    //------------------------------------------------------------------------------
    void
    Model::set_design_variable_interface( MSI::Design_Variable_Interface* aDesignVariableInterface )
    {
        mDesignVariableInterface = aDesignVariableInterface;

        //            mEquationModel->set_design_variable_interface( mDesignVariableInterface );
        //
        //            mDesignVariableInterface->set_equation_model( mEquationModel );
    }

    //------------------------------------------------------------------------------
    void
    Model::set_weak_bcs( const Matrix< DDRMat >& aWeakBCs )
    {
        Vector< MSI::Equation_Object* > tFemClusters = mEquationModel->get_equation_objects();

        // set weak BCs
        for ( auto tElement : tFemClusters )
        {
            Matrix< DDRMat >& tNodalWeakBCs  = tElement->get_weak_bcs();
            uint              tNumberOfNodes = tElement->get_num_nodes();
            tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                // copy weak bc into element
                tNodalWeakBCs( k ) = aWeakBCs( tElement->get_node_index( k ) );
            }
        }
    }

    //------------------------------------------------------------------------------
    void
    Model::set_weak_bcs_from_nodal_field( moris_index aFieldIndex )
    {
        Vector< MSI::Equation_Object* > tFemClusters = mEquationModel->get_equation_objects();

        for ( auto tElement : tFemClusters )
        {
            Matrix< DDRMat >& tNodalWeakBCs  = tElement->get_weak_bcs();
            uint              tNumberOfNodes = tElement->get_num_nodes();
            tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                // copy weak bc into element
                tNodalWeakBCs( k ) = mMeshManager->get_interpolation_mesh( mMeshPairIndex )->get_value_of_scalar_field( aFieldIndex, mtk::EntityRank::NODE, tElement->get_node_index( k ) );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::set_dof_order( const uint aBSplineIndex )
    {
        mBSplineIndex = aBSplineIndex;
    }

    //------------------------------------------------------------------------------

    void
    Model::output_solution(
            const uint aVisMeshIndex,
            const real aTime,
            const bool aCloseFile )
    {
        Tracer tTracer( "MDL", "Model", "Post-Process and Output Solution" );

        // create vis mesh and setup output
        mOutputManager->setup_vis_mesh_for_output(
                aVisMeshIndex,
                mMeshManager,
                mMeshPairIndex,
                mEquationModel );

        // write requested fields
        mOutputManager->write_field(
                aVisMeshIndex,
                aTime,
                mEquationModel );

        if ( aCloseFile )
        {
            // end writing and delete vis mesh
            mOutputManager->end_writing( aVisMeshIndex );
        }

        // populate FEM fields if requested for output
        mEquationModel->populate_fields();
    }

    //------------------------------------------------------------------------------

    void Model::update_model()
    {
        mEquationModel->update_fields();
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > >
    Model::get_mtk_fields()
    {
        return mEquationModel->get_fields();
    }

    //------------------------------------------------------------------------------

    std::unordered_map< enum MSI::Dof_Type, moris_index >
    Model::build_dof_type_to_mesh_index()
    {
        // get the number of DoF-types there are
        uint tNumDofTypes = (uint)MSI::Dof_Type::END_ENUM;

        // get the enum to string map for the dof types
        map< enum MSI::Dof_Type, std::string > tDofEnumToStringMap = MSI::get_dof_type_name_map();

        // initialize map relating the MSI Dof Type to the mesh index
        std::unordered_map< MSI::Dof_Type, moris_index > tDofTypeToMeshIndex;

        // loop over list of dof types and retrieve the mesh index for each one
        for ( uint iDofType = 0; iDofType < tNumDofTypes; iDofType++ )
        {
            // get the Enum for the current DoF type
            enum MSI::Dof_Type tDofTypeEnum = (MSI::Dof_Type)iDofType;

            // making sure the DoF type exists in the list
            if ( tDofEnumToStringMap.key_exists( tDofTypeEnum ) )
            {
                // get the DoF type name

                // get the mesh index from the parameter list
                moris_index tMeshIndex = MSI::get_mesh_index_for_dof_type( tDofTypeEnum, mMSIParameterList( 0 )( 0 ) );

                // store mesh index for Dof type in map
                tDofTypeToMeshIndex[ tDofTypeEnum ] = tMeshIndex;
            }
        }

        // return the map built
        return tDofTypeToMeshIndex;
    }

    //------------------------------------------------------------------------------

    uint
    Model::get_spatial_dim()
    {
        mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

        return tInterpolationMesh->get_spatial_dim();
    }

}    // namespace moris::mdl
