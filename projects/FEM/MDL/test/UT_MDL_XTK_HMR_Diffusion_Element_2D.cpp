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
#include "cl_FEM_Constitutive_User_Defined_Info.hpp"      //FEM/INT/src

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
    real mXc = 2.577;
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

         tParameters.set_number_of_elements_per_dimension( { {4}, {1}} );
         tParameters.set_domain_dimensions({ {4}, {1} });
         tParameters.set_domain_offset({ {0}, {0} });
         tParameters.set_bspline_truncation( true );
         tParameters.set_side_sets({ {2},{4} });

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {1} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {1} } );
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

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<0; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
            tField->evaluate_scalar_function( LevelSetPlaneFunction );
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

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create IWG user defined info
        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
        tIWGUserDefinedInfo( 0 ).resize( 1 );
        tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 1 ).resize( 1 );
        tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 2 ).resize( 1 );
        tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );
        tIWGUserDefinedInfo( 3 ).resize( 1 );
        tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::TEMP_NEUMANN },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );

        // create property user defined info
        Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 4 );
        tPropertyUserDefinedInfo( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 1 ).resize( 1 );
        tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 2 ).resize( 2 );
        tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 2 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 5.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 3 ).resize( 1 );
        tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 20.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );


        // create constitutive user defined info
        Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 4 );
        tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      2 );
        tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      2 );

        // create a list of active block-sets
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
        moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0"),
                                                  tEnrIntegMesh.get_side_set_index(tInterfaceSideSetName)};


        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              tBSplineMeshIndex,
                                              tIWGUserDefinedInfo,
                                              tSetList, tSetTypeList,
                                              tPropertyUserDefinedInfo,
                                              tConstitutiveUserDefinedInfo,
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


        // verify solution
        moris::Matrix<DDRMat> tGoldSolution =
          {{ 5.000000e+00},
           { 2.500000e+01},
           { 4.500000e+01},
           { 6.500000e+01},
           { 5.000000e+00},
           { 2.500000e+01},
           { 4.500000e+01},
           { 6.500000e+01}};

        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);

        CHECK(norm(tFullSol - tGoldSolution)<1e-08);

        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

         // Declare the fields related to enrichment strategy in output options
        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);

            tSTKIntegSol(i) = tIntegSol(tMyIndex);
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}

TEST_CASE("HMR Interpolation XTK Cut Diffusion Model  With STK Lag Order 1 In 2D","[XTK_HMR_STK_DIFF_2D]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;

         tParameters.set_number_of_elements_per_dimension( { {4}, {1}} );
         tParameters.set_domain_dimensions({ {4}, {1} });
         tParameters.set_domain_offset({ {0}, {0} });
         tParameters.set_bspline_truncation( true );
         tParameters.set_side_sets({ {2},{4} });

         tParameters.set_output_meshes( { {0} } );

         tParameters.set_lagrange_orders  ( { {1} });
         tParameters.set_lagrange_patterns({ {0} });

         tParameters.set_bspline_orders   ( { {1} } );
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

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<0; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
            tField->evaluate_scalar_function( LevelSetPlaneFunction );
        }

        tHMR.finalize();


        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l1_b1_2D.e" );

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

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = true;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);

        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );


        // create IWG user defined info
        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
        tIWGUserDefinedInfo( 0 ).resize( 1 );
        tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 1 ).resize( 1 );
        tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 2 ).resize( 1 );
        tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );
        tIWGUserDefinedInfo( 3 ).resize( 1 );
        tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, 2, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::TEMP_NEUMANN },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );

        // create property user defined info
        Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 4 );
        tPropertyUserDefinedInfo( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 1 ).resize( 1 );
        tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 2 ).resize( 2 );
        tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 2 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 5.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 3 ).resize( 1 );
        tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 20.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );


        // create constitutive user defined info
        Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 4 );
        tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      2 );
        tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      2 );

        // create a list of active block-sets
        moris::Cell< moris_index >  tSetList = {  tIntegMesh1->get_block_set_index("child_0"),
                                                  tIntegMesh1->get_block_set_index("parent_0"),
                                                  tIntegMesh1->get_side_set_index("SideSet_2"),
                                                  tIntegMesh1->get_side_set_index("iside_0")};



        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              tBSplineMeshIndex,
                                              tIWGUserDefinedInfo,
                                              tSetList, tSetTypeList,
                                              tPropertyUserDefinedInfo,
                                              tConstitutiveUserDefinedInfo,
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


        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);

        // verify solution
        moris::Matrix<DDRMat> tGoldSolution =
          {{ 5.000000e+00},
           { 2.500000e+01},
           { 4.500000e+01},
           { 6.500000e+01},
           { 5.000000e+00},
           { 2.500000e+01},
           { 4.500000e+01},
           { 6.500000e+01}};

        CHECK(norm(tFullSol - tGoldSolution)<1e-08);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);

        tMeshOutputFile = "./mdl_exo/xtk_hmr_stk_bar_plane_integ_2d.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        tModel->output_solution( "Circle" );
        tField->put_scalar_values_on_field( tModel->get_mSolHMR() );
        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l1_b1_2D.e" );

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}
