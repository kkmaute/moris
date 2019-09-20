#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

moris::real
LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 1.0;
    real mYn = 0.0;
    real mXc = 3.3;
    real mYc = 0.0;
    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc);
}

moris::real
CircleFunc(const moris::Matrix< moris::DDRMat > & aPoint )
{

    moris::real mXCenter = 0;
    moris::real mYCenter = 0;
    moris::real mRadius = 0.7;

    return    -((aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
            + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
            - (mRadius * mRadius));
}

Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
                                    fem::Geometry_Interpolator             * aGeometryInterpolator )
{
    return aCoeff( 0 );
}

TEST_CASE("HMR Interpolation XTK Cut Diffusion Model Lag Order 2 In 2D","[XTK_HMR_DIFF_2D]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {6}, {6}} );
         tParameters.set_domain_dimensions({ {4}, {4} });
         tParameters.set_domain_offset({ {-2}, {-2.0} });
         tParameters.set_bspline_truncation( true );
         tParameters.set_side_sets({ {2},{4} });

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {1} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {2} } );
         tParameters.set_bspline_patterns ( { {0} } );

         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );

         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2);

         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
         tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( CircleFunc );

        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
            tField->evaluate_scalar_function( CircleFunc );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_plane_interp_l2_b2_2D.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        xtk::Geom_Field tFieldAsGeom(tField);

        moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        size_t tModelDimension = 2;
        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine,tModelDimension);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = true;

        // Do the cutting
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        std::cout<<"tEnrIntegMesh = "<<tEnrIntegMesh.get_num_nodes()<<std::endl;

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create a list of IWG type
        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 3 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

        // number of groups of IWgs
        uint tNumSets = tIWGTypeList.size();

        // list of residual dof type
        moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResidualDofType( tNumSets );
        tResidualDofType( 0 ).resize( tIWGTypeList( 0 ).size(), { MSI::Dof_Type::TEMP } );
        tResidualDofType( 1 ).resize( tIWGTypeList( 1 ).size(), { MSI::Dof_Type::TEMP } );
        tResidualDofType( 2 ).resize( tIWGTypeList( 2 ).size(), { MSI::Dof_Type::TEMP } );
        tResidualDofType( 3 ).resize( tIWGTypeList( 3 ).size(), { MSI::Dof_Type::TEMP } );

        // list of IWG master dof dependencies
        moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > tMasterDofTypes( tNumSets );
        tMasterDofTypes( 0 ).resize( tIWGTypeList( 0 ).size(), {{ MSI::Dof_Type::TEMP }} );
        tMasterDofTypes( 1 ).resize( tIWGTypeList( 1 ).size(), {{ MSI::Dof_Type::TEMP }} );
        tMasterDofTypes( 2 ).resize( tIWGTypeList( 2 ).size(), {{ MSI::Dof_Type::TEMP }} );
        tMasterDofTypes( 3 ).resize( tIWGTypeList( 3 ).size(), {{ MSI::Dof_Type::TEMP }} );


        // list of IWG master property dependencies
        moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > tMasterPropTypes( tNumSets );
        tMasterPropTypes( 0 ).resize( tIWGTypeList( 0 ).size(), { fem::Property_Type::CONDUCTIVITY } );
        tMasterPropTypes( 1 ).resize( tIWGTypeList( 1 ).size(), { fem::Property_Type::CONDUCTIVITY } );
        tMasterPropTypes( 2 ).resize( tIWGTypeList( 2 ).size(), { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET } );
        tMasterPropTypes( 3 ).resize( tIWGTypeList( 3 ).size(), { fem::Property_Type::TEMP_NEUMANN } );


        // build an IWG user defined info
        fem::IWG_User_Defined_Info tIWGUserDefinedInfo( tIWGTypeList,
                                                        tResidualDofType,
                                                        tMasterDofTypes,
                                                        tMasterPropTypes );

        // list of property type
        Cell< fem::Property_Type > tPropertyTypeList = {{ fem::Property_Type::CONDUCTIVITY   },
                                                        { fem::Property_Type::TEMP_DIRICHLET },
                                                        { fem::Property_Type::TEMP_NEUMANN   }};

        // list of property dependencies
        Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 3 );

        // list of the property coefficients
        Cell< Cell< Matrix< DDRMat > > > tCoeffList( 3 );
        tCoeffList( 0 ).resize( 1 );
        tCoeffList( 0 )( 0 )= {{ 1.0 }};
        tCoeffList( 1 ).resize( 1 );
        tCoeffList( 1 )( 0 )= {{ 5.0 }};
        tCoeffList( 2 ).resize( 1 );
        tCoeffList( 2 )( 0 )= {{ 20.0 }};

        // cast free function into std::function
        fem::PropertyFunc tValFunction0 = tConstValFunction;

        // create the list with function pointers for the value
        Cell< fem::PropertyFunc > tValFuncList( 3, tValFunction0 );

        // create the list with cell of function pointers for the derivatives
        Cell< Cell< fem::PropertyFunc > > tDerFuncList( 3 );

        // collect properties info
        fem::Property_User_Defined_Info tPropertyUserDefinedInfo( tPropertyTypeList,
                                                                  tPropertyDofList,
                                                                  tCoeffList,
                                                                  tValFuncList,
                                                                  tDerFuncList );

        // interface side name
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);

        // create a list of active block-sets
        moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_side_set_index(tInterfaceSideSetName),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0")};

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };


        // create a list of BC type for the side-sets
        moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                           fem::BC_Type::NEUMANN};

        // create a list of active double side-sets
        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              tBSplineMeshIndex,
                                              &tIWGUserDefinedInfo,
                                              tSetList, tSetTypeList,
                                              &tPropertyUserDefinedInfo,
                                              0,
                                              false);

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        //        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
        //        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        //        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        NLA::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();


        // TODO: add gold solution data for this problem

        //        print_fancy(tFullSol,"Full Solution");

        // verify solution
        //        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);

        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

         // Declare the fields related to enrichment strategy in output options
         Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );



        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);


        print_fancy(tSTKIntegSol,"tSTKIntegSol");
        //        print_fancy(tIntegSol,"tIntegSol");


        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);


        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}

//TEST_CASE("HMR Interpolation XTK Cut Diffusion Model  With STK Lag Order 2 In 2D","[XTK_HMR_STK_DIFF_2D]")
//{
//    if(par_size() == 1)
//    {
//        std::string tFieldName = "Cylinder";
//
//         moris::uint tLagrangeMeshIndex = 0;
//         moris::uint tBSplineMeshIndex = 0;
//
//         moris::hmr::Parameters tParameters;
//
//         tParameters.set_number_of_elements_per_dimension( { {6}, {6}} );
//         tParameters.set_domain_dimensions({ {4}, {4} });
//         tParameters.set_domain_offset({ {-2}, {-2.0} });
//         tParameters.set_bspline_truncation( true );
//         tParameters.set_side_sets({ {2},{4} });
//
//         tParameters.set_output_meshes( { {0} } );
//
//         tParameters.set_lagrange_orders  ( { {1} });
//         tParameters.set_lagrange_patterns({ {0} });
//
//         tParameters.set_bspline_orders   ( { {1} } );
//         tParameters.set_bspline_patterns ( { {0} } );
//
//         tParameters.set_union_pattern( 2 );
//         tParameters.set_working_pattern( 3 );
//
//         tParameters.set_refinement_buffer( 2 );
//         tParameters.set_staircase_buffer( 2);
//
//         Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
//         tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        hmr::HMR tHMR( tParameters );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        // create field
//        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//        tField->evaluate_scalar_function( CircleFunc );
//
//        for( uint k=0; k<9; ++k )
//        {
//            tHMR.flag_surface_elements( tField );
//            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
//            tHMR.update_refinement_pattern( 0 );
//
//            tField->evaluate_scalar_function( CircleFunc );
//        }
//
//        tHMR.finalize();
//
//        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_plane_interp_l2_b2_2D.e" );
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        xtk::Geom_Field tFieldAsGeom(tField);
//
//        moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        // Tell the geometry engine about the discrete field mesh and how to interpret phases
//        size_t tModelDimension = 2;
//        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
//
//        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine,tModelDimension);
//        tXTKModel.mSameMesh = true;
//        tXTKModel.mVerbose = true;
//
//        // Do the cutting
//        tXTKModel.decompose(tDecompositionMethods);
//
//        std::cout<<"mModel->mBackgroundMesh.get_num_entities(EntityRank::NODE) = "<<tXTKModel.get_background_mesh().get_num_entities(EntityRank::NODE) <<std::endl;
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = true;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//
//        std::cout<<"tIntegMesh Nodes = "<<tIntegMesh1->get_num_entities(EntityRank::NODE)<<std::endl;
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);
//
//        // create a list of IWG type
//        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
//        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
//        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );
//
//
//        moris::Cell< moris_index >  tBlocksetList = { tIntegMesh1->get_block_set_index("child_0"),
//                                                      tIntegMesh1->get_block_set_index("parent_0") };
//
//
//        // create a list of active side-sets
//        moris::Cell< moris_index >  tSidesetList = { tIntegMesh1->get_side_set_index("SideSet_1"),
//                                                     tIntegMesh1->get_side_set_index("SideSet_2" )};
//
//
//        // create a list of BC type for the side-sets
//        moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
//                                                           fem::BC_Type::NEUMANN};
//
//        // create a list of active double side-sets
//        moris::Cell< moris_index >  tDoubleSidesetList = {  };
//
//        // create model
//        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex, tIWGTypeList,
//                                              tBlocksetList, tSidesetList,
//                                              tSidesetBCTypeList,
//                                              tDoubleSidesetList );
//
//        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 1: create linear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        dla::Solver_Factory  tSolFactory;
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//
//        dla::Linear_Solver tLinSolver;
//
//        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 2: create nonlinear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        //        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        //        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        //        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        //        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//        NLA::Nonlinear_Solver tNonlinearSolver;
//
//        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 3: create time Solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        tsa::Time_Solver_Factory tTimeSolverFactory;
//        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//        tsa::Time_Solver tTimeSolver;
//
//        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//        NLA::SOL_Warehouse tSolverWarehouse;
//
//        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
//        tTimeSolver.set_dof_type_list( tDofTypes1 );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 4: Solve and check
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        tTimeSolver.solve();
//
//
//        // TODO: add gold solution data for this problem
//
//        //        print_fancy(tFullSol,"Full Solution");
//
//        // verify solution
//        //        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);
//
//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//        //        print_fancy(tIntegSol,"tIntegSol");
//
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);
//
//
//        Matrix<DDRMat> tFullSol;
//        tTimeSolver.get_full_solution(tFullSol);
//
//        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d.e";
//        tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//        tModel->output_solution( "Circle" );
//        tField->put_scalar_values_on_field( tModel->get_mSolHMR() );
//        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l2_b2_2D.e" );
//
//        //    delete tInterpMesh1;
//        delete tModel;
//        delete tIntegMesh1;
//    }
//}
