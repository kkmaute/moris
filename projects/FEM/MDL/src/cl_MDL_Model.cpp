
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
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

#include "cl_PRM_MSI_Parameters.hpp"

#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------
    Model::Model(       mtk::Mesh_Manager                 * aMeshManager,
                  const uint                                aBSplineIndex,
                        moris::Cell< fem::Set_User_Info > & aSetInfo,
                  const moris_index                         aMeshPairIndex,
                  const bool                                aUseMultigrid )
    : mMeshManager( aMeshManager ),
      mMeshPairIndex( aMeshPairIndex ),
      mBSplineIndex ( aBSplineIndex ),
      mUseMultigrid( aUseMultigrid )
    {
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 0: initialize
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // get pointers to interpolation and integration mesh
        mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create the FEM model
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // start timer
        tic tTimer1;

        // build the FEM model from FEM parameter list
        mEquationModel = std::make_shared< fem::FEM_Model >( mMeshManager,
                                                             mMeshPairIndex,
                                                             aSetInfo );
        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( "Model: created FEM model in %5.3f seconds.",
                            ( double ) tElapsedTime / 1000 );
            MORIS_LOG_INFO( " " );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create the model solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // start timer
        tic tTimer2;

        // create the MSI parameter list
        moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
        tMSIParameters.set( "multigrid", mUseMultigrid );

        if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
        {
            tMSIParameters.set( "L2", (sint)mBSplineIndex );
            tMSIParameters.set( "TEMP", (sint)mBSplineIndex );
        }

        // build the model solver interface from the MSI parameter list
        mModelSolverInterface = new moris::MSI::Model_Solver_Interface( tMSIParameters,
                                                                        mEquationModel,
                                                                        tInterpolationMesh );

        // finalize the fem sets
        mEquationModel->finalize_equation_sets( mModelSolverInterface );

        // finalize the model solver interface
        mModelSolverInterface->finalize();

        // calculate AdofMap
        mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( "Model: created Model-Solver Interface in %5.3f seconds.",
                            ( double ) tElapsedTime / 1000 );
            MORIS_LOG_INFO( " " );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: create the solver interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        mSolverInterface = new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

        mSolverInterface->set_model( this );
    }

//------------------------------------------------------------------------------
    Model::Model(       std::shared_ptr< Library_IO > aLibrary,
                  const uint                          aBSplineIndex,
                  const moris_index                   aMeshPairIndex )
    : mMeshPairIndex( aMeshPairIndex ),
      mBSplineIndex( aBSplineIndex ),
	  mLibrary( aLibrary )
    {
        // start timer
        tic tTimer0;

        // load the MSI parameter list
        std::string tMSIString = "MSIParameterList";
        MORIS_PARAMETER_FUNCTION tMSIParameterListFunc = mLibrary->load_parameter_file( tMSIString );
        tMSIParameterListFunc( mMSIParameterList );

        // load the SOL parameter list
        std::string tSOLString = "SOLParameterList";
        MORIS_PARAMETER_FUNCTION tSOLParameterListFunc = mLibrary->load_parameter_file( tSOLString );
        tSOLParameterListFunc( mSOLParameterList );

        // load the FEM parameter list
        std::string tFEMString = "FEMParameterList";
        MORIS_PARAMETER_FUNCTION tFEMParameterListFunc = mLibrary->load_parameter_file( tFEMString );
        tFEMParameterListFunc( mFEMParameterList );

        // load the VIS parameter list
        std::string tVISString = "VISParameterList";
        MORIS_PARAMETER_FUNCTION tVISParameterListFunc = mLibrary->load_parameter_file( tVISString );
        tVISParameterListFunc( mVISParameterList );

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer0.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( "Model: unpack of the inputs in %5.3f seconds.",
                            ( double ) tElapsedTime / 1000 );
            MORIS_LOG_INFO( " " );
        }
    }

//------------------------------------------------------------------------------
    Model::Model
    (       mtk::Mesh_Manager*                  aMeshManager,
      const uint                                aBSplineIndex,
            moris::Cell< fem::Set_User_Info > & aSetInfo,
            MSI::Design_Variable_Interface    * aDesignVariableInterface,
      const moris_index                         aMeshPairIndex,
      const bool                                aUseMultigrid )
    : mMeshManager( aMeshManager ),
      mMeshPairIndex( aMeshPairIndex ),
      mBSplineIndex( aBSplineIndex ),
      mUseMultigrid( aUseMultigrid )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create the FEM model
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer1;

            // build the FEM model from FEM parameter list
            mEquationModel = std::make_shared< fem::FEM_Model >( mMeshManager,
                                                                 mMeshPairIndex,
                                                                 aSetInfo );
            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created FEM model in %5.3f seconds.",
                                ( double ) tElapsedTime / 1000 );
                MORIS_LOG_INFO( " " );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create the model solver interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer2;

            // Get pointers to interpolation and integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

            // create the MSI parameter list
            moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
            tMSIParameters.set( "multigrid", mUseMultigrid );

            // create the model solver interface from the MSI parameter list
            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( tMSIParameters,
                                                                            mEquationModel,
                                                                            tInterpolationMesh );

            // finalize the fem sets
            // FIXME pass in the interface
            mEquationModel->finalize_equation_sets( mModelSolverInterface );

            // finalize the model solver interface
            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created Model-Solver Interface in %5.3f seconds.",
                                ( double ) tElapsedTime / 1000 );
                MORIS_LOG_INFO( " " );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // create the solver interface
            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );
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
        void Model::set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
        {
            mMeshManager = aMTKPerformer.get();
        }

//------------------------------------------------------------------------------
        void Model::initialize()
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create the FEM model
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer1;

            // build the FEM model from FEM parameter list
            mEquationModel = std::make_shared< fem::FEM_Model >( mMeshManager,
                                                                 mMeshPairIndex,
                                                                 mFEMParameterList,
                                                                 mLibrary );

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created FEM model in %5.3f seconds.",
                                ( double ) tElapsedTime / 1000 );
                MORIS_LOG_INFO( " " );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create the model solver interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer2;

            // get pointers to interpolation and integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );

            // Does not work with STK
            MORIS_ERROR( tInterpolationMesh->get_mesh_type() != MeshType::STK, "Does not work for STK");

            // build the model solver interface
            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mMSIParameterList( 0 )( 0 ),
                                                                            mEquationModel,
                                                                            tInterpolationMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create the solver interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create the solver
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimerSolver;
            mSolverWarehouse = std::make_shared< sol::SOL_Warehouse >( mSolverInterface, mLibrary );

            mSolverWarehouse->set_parameterlist( mSOLParameterList );

            mSolverWarehouse->initialize();

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimerSolver.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created Solver Warehouse in %5.3f seconds.\n\n",
                                ( double ) tElapsedTime / 1000 );
            }

            // finalize the model solver interface
            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created Model-Solver Interface in %5.3f seconds.",
                                ( double ) tElapsedTime / 1000 );
                MORIS_LOG_INFO( " " );
            }

            // finalize the fem sets
            mEquationModel->finalize_equation_sets( mModelSolverInterface );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 5: create the output manager
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimerVisMesh;

            mOutputManager = new vis::Output_Manager( mVISParameterList( 0 )( 0 ) );

            mOutputManagerOwned = true;

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimerVisMesh.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created Vis Mesh in %5.3f seconds.\n\n",
                                ( double ) tElapsedTime / 1000 );
            }
        }

//------------------------------------------------------------------------------
        void Model::perform( uint aIndex )
        {
            if( aIndex == 0 )
            {
                this->perform_forward_analysis();
            }
            else if( aIndex == 1 )
            {
                this->perform_sensitivity_analysis();
            }
        }

//------------------------------------------------------------------------------
        void Model::perform_forward_analysis()
        {
            mEquationModel->set_is_forward_analysis();

            mSolverWarehouse->get_main_time_solver()->solve();
        }

//------------------------------------------------------------------------------

        moris::Cell< moris::Matrix< DDRMat > > Model::perform_post_processing()
        {
            return mEquationModel->compute_IQIs();
        }

//------------------------------------------------------------------------------
        void Model::perform_sensitivity_analysis()
        {
            mEquationModel->set_is_sensitivity_analysis();

            mSolverWarehouse->get_main_time_solver()->solve_sensitivity();

            mEquationModel->compute_implicit_dQIdp();

            mEquationModel->compute_explicit_dQIdp();
        }

//------------------------------------------------------------------------------
        void Model::set_design_variable_interface( MSI::Design_Variable_Interface * aDesignVariableInterface )
        {
            mDesignVariableInterface = aDesignVariableInterface;

            mEquationModel->set_design_variable_interface( mDesignVariableInterface );

            mDesignVariableInterface->set_equation_model( mEquationModel );
        }

//------------------------------------------------------------------------------
        void Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
            moris::Cell< MSI::Equation_Object * > tFemClusters
            = mEquationModel->get_equation_objects();

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
            moris::Cell< MSI::Equation_Object * > tFemClusters
            = mEquationModel->get_equation_objects();

            for( auto tElement : tFemClusters )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( k ) = mMeshManager->get_interpolation_mesh( mMeshPairIndex )
                                                     ->get_value_of_scalar_field( aFieldIndex,
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
        void Model::output_solution( const uint aVisMeshIndex,
                                     const real aTime,
                                     const bool aCloseFile )
        {
            // create vis mesh and setup output
            mOutputManager->setup_vis_mesh_for_output(aVisMeshIndex,
                                                      mMeshManager,
                                                      mMeshPairIndex,
                                                      mEquationModel );

            // write requested fields
            mOutputManager->write_field( aVisMeshIndex,
                                         aTime,
                                         mEquationModel );

            if( aCloseFile )
            {
                // end writing and delete vis mesh
                mOutputManager->end_writing( aVisMeshIndex );
            }
        }

    } /* namespace mdl */
} /* namespace moris */
