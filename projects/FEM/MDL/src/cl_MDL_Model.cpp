
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_MDL_Model.hpp"

#include "cl_FEM_Enums.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_VIS_Output_Manager.hpp"

#include "fn_PRM_MSI_Parameters.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Library_IO.hpp"

namespace moris
{
    // Parameter function
    typedef void ( *Parameter_Function ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );

    namespace mdl
    {
        //------------------------------------------------------------------------------
        Model::Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const uint                          aBSplineIndex,
                moris::Cell< fem::Set_User_Info > & aSetInfo,
                const moris_index                   aMeshPairIndex,
                const bool                          aUseMultigrid )
        : mMeshManager( aMeshManager ),
          mMeshPairIndex( aMeshPairIndex ),
          mBSplineIndex ( aBSplineIndex ),
          mUseMultigrid( aUseMultigrid )
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
                Tracer tTracer( "FEM", "Model", "CreateEquationModel" );

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
                Tracer tTracer( "FEM", "Model", "CreateMSI" );

                // create the MSI parameter list
                moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
                tMSIParameters.set( "multigrid", mUseMultigrid );

                if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
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
                Tracer tTracer( "FEM", "Model", "CreateSolverInterface" );

                mSolverInterface = new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

                mSolverInterface->set_model( this );
            }
        }

        //------------------------------------------------------------------------------
        Model::Model(std::shared_ptr< Library_IO > aLibrary,
                const uint                         aBSplineIndex,
                const moris_index                  aMeshPairIndex )
        : mMeshPairIndex( aMeshPairIndex ),
          mBSplineIndex( aBSplineIndex ),
          mLibrary( aLibrary )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            // load the MSI parameter list
            std::string tMSIString = "MSIParameterList";
            Parameter_Function tMSIParameterListFunc = mLibrary->load_function<Parameter_Function>( tMSIString );
            tMSIParameterListFunc( mMSIParameterList );

            // load the SOL parameter list
            std::string tSOLString = "SOLParameterList";
            Parameter_Function tSOLParameterListFunc = mLibrary->load_function<Parameter_Function>( tSOLString );
            tSOLParameterListFunc( mSOLParameterList );

            // load the FEM parameter list
            std::string tFEMString = "FEMParameterList";
            Parameter_Function tFEMParameterListFunc = mLibrary->load_function<Parameter_Function>( tFEMString );
            tFEMParameterListFunc( mFEMParameterList );

            // load the VIS parameter list
            std::string tVISString = "VISParameterList";
            Parameter_Function tVISParameterListFunc = mLibrary->load_function<Parameter_Function>( tVISString );
            tVISParameterListFunc( mVISParameterList );
        }

        //------------------------------------------------------------------------------
        Model::Model(
                std::shared_ptr< mtk::Mesh_Manager >aMeshManager,
                const uint                          aBSplineIndex,
                moris::Cell< fem::Set_User_Info > & aSetInfo,
                MSI::Design_Variable_Interface    * aDesignVariableInterface,
                const moris_index                   aMeshPairIndex,
                const bool                          aUseMultigrid )
        : mMeshManager( aMeshManager ),
          mMeshPairIndex( aMeshPairIndex ),
          mBSplineIndex( aBSplineIndex ),
          mUseMultigrid( aUseMultigrid )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create the FEM model
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                Tracer tTracer( "FEM", "Model", "CreateEquationModel" );

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
                Tracer tTracer( "FEM", "Model", "CreateMSI" );

                // Get pointers to interpolation and integration mesh
                mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

                // create the MSI parameter list
                moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
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
                Tracer tTracer( "FEM", "Model", "CreateSolverInterface" );

                mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

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
            if( mOutputManager != nullptr && mOutputManagerOwned == true )
            {
                delete mOutputManager;

                mOutputManagerOwned = false;
            }
        }

        //------------------------------------------------------------------------------

        void Model::free_memory()
        {
            // delete SI
            delete mSolverInterface;
			mSolverInterface = nullptr;

            // delete MSI
            delete mModelSolverInterface;
			mModelSolverInterface = nullptr;

            // delete MSI
            if( mOutputManager != nullptr && mOutputManagerOwned == true )
            {
                delete mOutputManager;
				mOutputManager = nullptr;

                mOutputManagerOwned = false;
            }

            mEquationModel->free_memory();
        }

        //------------------------------------------------------------------------------
        void Model::set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
        {
            mMeshManager = aMTKPerformer;
        }

        //------------------------------------------------------------------------------

        void Model::initialize()
        {
            Tracer tTracer( "FEM", "Model", "Initialize" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: delete model and interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			if( mModelSolverInterface != nullptr)
			{
				delete mModelSolverInterface;
			}
			if( mSolverInterface != nullptr)
			{
				delete mSolverInterface;
			}
			if( mOutputManager != nullptr)
			{
				delete mOutputManager;
			}

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create the FEM model
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                Tracer tTracer( "FEM", "Model", "CreateEquationModel" );

                // build the FEM model from FEM parameter list
                mEquationModel = std::make_shared< fem::FEM_Model >(
                        mMeshManager,
                        mMeshPairIndex,
                        mFEMParameterList,
                        mLibrary,
                        mDesignVariableInterface );

                // set the equation model if using design variables
                if (mDesignVariableInterface != nullptr)
                {
                    mDesignVariableInterface->set_equation_model( mEquationModel );
                }
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create the model solver interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                Tracer tTracer( "FEM", "Model", "CreateMSI" );

                // get pointers to interpolation and integration mesh
                mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

                // Does not work with STK
                // MORIS_ERROR( tInterpolationMesh->get_mesh_type() != MeshType::STK, "Does not work for STK");

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
                Tracer tTracer( "FEM", "Model", "CreateSolverInterface" );

                mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

                mSolverInterface->set_model( this );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 4: create the solver
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // start timer
                {
                    Tracer tTracer( "FEM", "Model", "CreateSolver" );

                    mSolverWarehouse = std::make_shared< sol::SOL_Warehouse >( mSolverInterface, mLibrary );

                    mSolverWarehouse->set_parameterlist( mSOLParameterList );

                    mSolverWarehouse->initialize();
                }

                //set the solver warehouse pointer
                mSolverInterface->set_solver_warehouse(mSolverWarehouse);

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
        void Model::perform( uint aIndex )
        {
            if( aIndex == 0 )
            {
                Tracer tTracer( "FEM", "Model", "PerformForwardAnalysis" );

                this->perform_forward_analysis();
            }
            else if( aIndex == 1 )
            {
                Tracer tTracer( "FEM", "Model", "PerformSensitivityAnalysis" );

                this->perform_sensitivity_analysis();
            }
        }

        //------------------------------------------------------------------------------
        void Model::perform_forward_analysis()
        {
            mEquationModel->set_is_forward_analysis();

            mEquationModel->initialize_IQIs();

            mSolverWarehouse->get_main_time_solver()->solve();
        }

        //------------------------------------------------------------------------------

        moris::Cell< moris::Matrix< DDRMat > > Model::get_IQI_values()
        {
            return mEquationModel->get_IQI_values();
        }

        //------------------------------------------------------------------------------
        void Model::perform_sensitivity_analysis()
        {
            mEquationModel->set_is_sensitivity_analysis();

            mEquationModel->initialize_explicit_and_implicit_dQIdp();

            mSolverWarehouse->get_main_time_solver()->solve_sensitivity();
        }

        //------------------------------------------------------------------------------
        void Model::set_design_variable_interface( MSI::Design_Variable_Interface * aDesignVariableInterface )
        {
            mDesignVariableInterface = aDesignVariableInterface;

//            mEquationModel->set_design_variable_interface( mDesignVariableInterface );
//
//            mDesignVariableInterface->set_equation_model( mEquationModel );
        }

        //------------------------------------------------------------------------------
        void Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
            moris::Cell< MSI::Equation_Object * > tFemClusters = mEquationModel->get_equation_objects();

            // set weak BCs
            for( auto tElement : tFemClusters )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k < tNumberOfNodes; ++k )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( k ) = aWeakBCs( tElement->get_node_index( k ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Model::set_weak_bcs_from_nodal_field( moris_index aFieldIndex )
        {
            moris::Cell< MSI::Equation_Object * > tFemClusters = mEquationModel->get_equation_objects();

            for( auto tElement : tFemClusters )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( k ) = mMeshManager->get_interpolation_mesh( mMeshPairIndex ) ->
                            get_value_of_scalar_field(
                                    aFieldIndex,
                                    EntityRank::NODE,
                                    tElement->get_node_index( k ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void Model::set_dof_order( const uint aBSplineIndex )
        {
            mBSplineIndex = aBSplineIndex;
        }

        //------------------------------------------------------------------------------
        void Model::output_solution(
                const uint aVisMeshIndex,
                const real aTime,
                const bool aCloseFile )
        {
            Tracer tTracer( "FemModel", "PostProcessing", "OutputSolution" );

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

            if( aCloseFile )
            {
                // end writing and delete vis mesh
                mOutputManager->end_writing( aVisMeshIndex );
            }

            mEquationModel->populate_fields();
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< mtk::Field > > Model::get_mtk_fields()
        {
            return mEquationModel->get_fields();
        }

    } /* namespace mdl */
} /* namespace moris */
