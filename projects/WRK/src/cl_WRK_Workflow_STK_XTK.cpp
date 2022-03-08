#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_STK_XTK.hpp"
#include "cl_Param_List.hpp"

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MDL_Model.hpp"


#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace wrk
    {
        //------------------------------------------------------------------------------

        // Parameter function
        typedef void ( *Parameter_Function ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );
        
        //--------------------------------------------------------------------------------------------------------------

        Workflow_STK_XTK::Workflow_STK_XTK( wrk::Performer_Manager * aPerformerManager )
        : Workflow( aPerformerManager )
        {
            Tracer tTracer( "Workflow", "STK_XTK", "Initialize" );
            MORIS_LOG_SPEC("Par_Rank",par_rank());
            MORIS_LOG_SPEC("Par_Size",par_size());

            // Performer set for this workflow
            mPerformerManager->mGENPerformer.resize( 1 );
            mPerformerManager->mXTKPerformer.resize( 1 );
            mPerformerManager->mMTKPerformer.resize( 2 );
            mPerformerManager->mMDLPerformer.resize( 1 );

            // load the STK parameter list
            std::string tSTKString = "STKParameterList";
            Parameter_Function tSTKParameterListFunc = mPerformerManager->mLibrary->load_function<Parameter_Function>( tSTKString );
            moris::Cell< moris::Cell< ParameterList > > tSTKParameterList;
            tSTKParameterListFunc( tSTKParameterList );

            // load the meshes
            mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();
            this->create_stk(tSTKParameterList);
            mPerformerManager->mMTKPerformer( 0 )->register_mesh_pair(mIpMesh.get(),mIgMesh.get());

            // moris::mtk::Cell tCell

            // // verify these meshes
            // mtk::Mesh_Checker tMeshChecker(0,mIpMesh.get(),mIgMesh.get());
            // tMeshChecker.perform();
            // tMeshChecker.print_diagnostics();

            // load gen parameter list
            std::string tGENString = "GENParameterList";
            Parameter_Function tGENParameterListFunc = mPerformerManager->mLibrary->load_function<Parameter_Function>( tGENString );
            moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
            tGENParameterListFunc( tGENParameterList );

            // Create GE performer
            mPerformerManager->mGENPerformer( 0 ) = std::make_shared< ge::Geometry_Engine >(
                    tGENParameterList,
                    mPerformerManager->mLibrary,
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0) );

            // create MTK performer - will be used for XTK mesh
            mPerformerManager->mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

            // create MDL performer
            mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary, 0 );

            // Set performer to MDL
            mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Workflow_STK_XTK::initialize(
                Matrix<DDRMat>& aADVs,
                Matrix<DDRMat>& aLowerBounds,
                Matrix<DDRMat>& aUpperBounds)
        {
            // Stage 2: Initialize Level set field in GEN -----------------------------------------------
            // Trace GEN
            Tracer tTracer( "GEN", "Levelset", "InitializeADVs" );

                mPerformerManager->mGENPerformer( 0 )->distribute_advs(
                        mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair(0),{} );

                // Get ADVs
                aADVs        = mPerformerManager->mGENPerformer( 0 )->get_advs();
                aLowerBounds = mPerformerManager->mGENPerformer( 0 )->get_lower_bounds();
                aUpperBounds = mPerformerManager->mGENPerformer( 0 )->get_upper_bounds();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow_STK_XTK::perform(const Matrix<DDRMat> & aNewADVs)
        {
            // Set new advs in GE
            mPerformerManager->mGENPerformer( 0 )->set_advs(aNewADVs);

            // Stage 1: HMR refinement

            // Stage 2: XTK -----------------------------------------------------------------------------
            this->create_xtk();

            // Compute level set data in GEN
            mPerformerManager->mGENPerformer( 0 )->reset_mesh_information(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ));

            // Output GEN fields, if requested
            mPerformerManager->mGENPerformer( 0 )->output_fields(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ));

            // XTK perform - decompose - enrich - ghost - multigrid
            mPerformerManager->mXTKPerformer( 0 )->perform_decomposition();
            mPerformerManager->mXTKPerformer( 0 )->perform_enrichment();

            // Assign PDVs
            mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0) );

            // Stage 3: MDL perform ---------------------------------------------------------------------

            mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                                mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

            mPerformerManager->mMDLPerformer( 0 )->initialize();

            mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

            // Build MDL components and solve
            mPerformerManager->mMDLPerformer( 0 )->perform();

            moris::Cell< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

            // Communicate IQIs
            for( uint iIQIIndex = 0; iIQIIndex < tVal.size(); iIQIIndex++ )
            {
                tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
            }

            moris::Matrix< DDRMat > tMat( tVal.size(), 1, 0.0 );

            for( uint Ik = 0; Ik < tVal.size(); Ik ++ )
            {
                tMat( Ik ) = tVal( Ik )( 0 );
            }

            return tMat;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow_STK_XTK::compute_dcriteria_dadv()
        {
            mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

            mPerformerManager->mMDLPerformer( 0 )->perform( 1 );

            Matrix<DDRMat> tDCriteriaDAdv = mPerformerManager->mGENPerformer( 0 )->get_dcriteria_dadv();

            if (par_rank() == 0)
            {
                MORIS_LOG_INFO ( "--------------------------------------------------------------------------------");
                MORIS_LOG_INFO ( "Gradients of design criteria wrt ADVs:");

                for (uint i=0;i<tDCriteriaDAdv.n_rows();++i)
                {
                    Matrix<DDRMat> tDIQIDAdv = tDCriteriaDAdv.get_row(i);

                    auto tItrMin = std::min_element(tDIQIDAdv.data(),tDIQIDAdv.data()+tDIQIDAdv.numel());
                    auto tIndMin = std::distance(tDIQIDAdv.data(),tItrMin);

                    auto tItrMax = std::max_element(tDIQIDAdv.data(),tDIQIDAdv.data()+tDIQIDAdv.numel());
                    auto tIndMax = std::distance(tDIQIDAdv.data(),tItrMax);

                    MORIS_LOG_INFO ( "Criteria(%i): norm = %e   min = %e  (index = %i)   max = %e  (index = %i)",
                                     i, norm(tDIQIDAdv),tDIQIDAdv.min(),tIndMin,tDIQIDAdv.max(),tIndMax);
                }

                MORIS_LOG_INFO ( "--------------------------------------------------------------------------------");
            }

            return tDCriteriaDAdv;
        }

        void Workflow_STK_XTK::create_xtk()
        {
            // Read parameter list from shared object
            Parameter_Function tXTKParameterListFunc = mPerformerManager->mLibrary->load_function<Parameter_Function>( "XTKParameterList" );
            moris::Cell< moris::Cell< ParameterList > > tXTKParameterList;
            tXTKParameterListFunc( tXTKParameterList );

            // Create XTK
            mPerformerManager->mXTKPerformer( 0 ) = std::make_shared< xtk::Model >( tXTKParameterList( 0 )( 0 ) );

            // Reset output MTK performer
            mPerformerManager->mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();
            mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

            // Set performers
            mPerformerManager->mXTKPerformer( 0 )->set_geometry_engine( mPerformerManager->mGENPerformer( 0 ).get() );
            mPerformerManager->mXTKPerformer( 0 )->set_input_performer( mPerformerManager->mMTKPerformer( 0 ) );
            mPerformerManager->mXTKPerformer( 0 )->set_output_performer( mPerformerManager->mMTKPerformer( 1 ) );
        }

        void
        Workflow_STK_XTK::create_stk(Cell< Cell<ParameterList> > & aParameterLists)
        {
            Tracer tTracer( "STK", "Mesh", "InitializeMesh" );
            std::string tMeshFile = aParameterLists(0)(0).get<std::string>("input_file");
            mtk::MtkMeshData* tSuppMeshData = nullptr;

            mtk::Cell_Cluster_Input* tCellClusterData = nullptr;

            // construct the meshes
            mIpMesh = std::make_shared<mtk::Interpolation_Mesh_STK>( tMeshFile, tSuppMeshData, true );
            mIgMesh = std::make_shared<mtk::Integration_Mesh_STK> ( *mIpMesh, tCellClusterData);
        }

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
