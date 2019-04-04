
#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "fn_norm.hpp"

moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 0.5;
}

namespace moris
{
    namespace mdl
    {
        TEST_CASE( "Diffusion_2x2x2", "[moris],[mdl],[Diffusion_2x2x2]" )
        {
        if(par_size() == 1 )
        {
            // Create a 3D mesh of HEX8 using MTK ------------------------------------------
            std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
            //------------------------------------------------------------------------------

            std::string tPrefix = std::getenv("MORISROOT");
            std::string tMeshFileName = tPrefix + "projects/FEM/INT/test/data/Cube_with_side_sets.g";
            std::cout<<"Mesh input name = "<<tMeshFileName<<std::endl;

            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName1 = "Temp_Field";
            tNodeField1.set_field_name( tFieldName1 );
            tNodeField1.set_field_entity_rank( EntityRank::NODE );

            // Initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;

            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

            // Declare some supplementary fields
            mtk::MtkMeshData tMeshData;
            tMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Mesh* tMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tMeshData );

            // create a list of IWG type
            Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
            tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
            tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // create a list of active sidesets
            moris::Cell< moris_index >  tSidesetList = { 3, 5 };

            // create a list of BC type for the sidesets
            moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                               fem::BC_Type::NEUMANN };

            // create model
            mdl::Model * tModel = new mdl::Model( tMesh, 1, tIWGTypeList,
                                                  tSidesetList, tSidesetBCTypeList );

            //solve
            moris::Matrix< DDRMat > tSolution11;
            tModel->solve( tSolution11 );

            // checking the solution--------------------------------------------------------
            //------------------------------------------------------------------------------
            // Expected solution
            Matrix< DDRMat > tExpectedSolution = {{ 25.0, 25.0, 25.0,
                                                    25.0,  5.0, 25.0,
                                                    45.0, 25.0,  5.0,
                                                    25.0, 45.0, 25.0,
                                                     5.0, 25.0, 45.0,
                                                     5.0, 45.0,  5.0,
                                                    45.0,  5.0, 45.0,
                                                     5.0, 45.0,  5.0,
                                                    45.0,  5.0, 45.0 }};

            // define an epsilon environment
            double tEpsilon = 1E-12;

            // define a bool for solution check
            bool tCheckNodalSolution = true;

            // number of mesh nodes
            uint tNumOfNodes = tMesh->get_num_nodes();

            // loop over the node and chyeck solution
            for ( uint i = 0; i < tNumOfNodes; i++ )
            {
                // check solution
                tCheckNodalSolution = tCheckNodalSolution
                                   && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
            }
            // check bool is true
            REQUIRE( tCheckNodalSolution );

            tModel->output_solution( tFieldName1 );

        }/* if( par_size() */
    }

        TEST_CASE( "Element_Diffusion_3", "[moris],[mdl],[Diffusion_block_7x8x9]" )
        {
        if(par_size() == 1 )
        {
            // Create a 3D mesh of HEX8 using MTK ------------------------------------------
            std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
            //------------------------------------------------------------------------------

            std::string tPrefix = std::getenv("MORISROOT");
            std::string tMeshFileName = tPrefix + "projects/FEM/MDL/test/data/Block_7x8x9.g";

            std::cout<<"Mesh input name = "<< tMeshFileName<<std::endl;

            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName1 = "Temp_Field";
            tNodeField1.set_field_name(tFieldName1);
            tNodeField1.set_field_entity_rank(EntityRank::NODE);

            // Initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;

            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

            // Declare some supplementary fields
            mtk::MtkMeshData tMeshData;
            tMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Mesh* tMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tMeshData );

            //1) Create the fem nodes ------------------------------------------------------
            std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------
            Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
            tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
            tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // create a list of active sidesets
            moris::Cell< moris_index >  tSidesetList = { 3, 5 };

            // create a list of BC type for the sidesets
            moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                               fem::BC_Type::NEUMANN };

            // create model
            mdl::Model * tModel = new mdl::Model( tMesh, 1, tIWGTypeList,
                                                  tSidesetList, tSidesetBCTypeList );

            //solve
            moris::Matrix< DDRMat > tSolution11;
            tModel->solve( tSolution11 );


            tModel->output_solution( tFieldName1 );

            // print(tSolution11,"Solution");

            // Expected solution
            Matrix< DDRMat > tExpectedSolution =
                 {{+2.500000000875184e+01, +2.500000000987847e+01, +2.500000001157717e+01,
                   +2.500000001368304e+01, +2.500000001479958e+01, +2.500000001683809e+01,
                   +2.500000002126988e+01, +2.500000002173021e+01, +2.500000000560846e+01,
                   +2.500000000793639e+01, +2.500000001182473e+01, +2.500000001665577e+01,
                   +2.500000002234817e+01, +2.500000002784368e+01, +2.500000003077465e+01,
                   +2.500000002928050e+01, +2.499999999946052e+01, +2.500000000375664e+01,
                   +2.500000001180321e+01, +2.500000001954013e+01, +2.500000002360115e+01,
                   +2.500000003273075e+01, +2.500000003627261e+01, +2.500000002781392e+01,
                   +2.499999998871123e+01 }};

            // define an epsilon environment
            double tEpsilon = 1E-12;

            // define a bool for solution check
            bool tCheckNodalSolution = true;

            // loop over the node and chyeck solution
            for ( uint i = 0; i < 25; i++ )
            {
            	// check solution
            	tCheckNodalSolution = tCheckNodalSolution
            			&& ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
            }
            // check bool is true
            REQUIRE( tCheckNodalSolution );

        }/* if( par_size() */
    }
	
     TEST_CASE( "Diffusion_hmr_10x4x4", "[moris],[mdl],[Diffusion_hmr_10x4x4]" )
     {
        if(par_size() == 1 )
        {
            // Create a 3D mesh of HEX8 using MTK ------------------------------------------
            std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
            //------------------------------------------------------------------------------

            moris::uint tBplineOrder = 1;
            moris::uint tLagrangeOrder = 1;
            moris::uint tMyCoeff = 1;

            hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "10, 4, 4" );
            tParameters.set( "domain_dimensions", "10, 4, 4" );
            tParameters.set( "domain_offset", "-9.0, -2.0, -2.0" );
            tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");
            tParameters.set( "verbose", 0 );
            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "lagrange_orders", "1" );

            tParameters.set( "use_multigrid", 0 );

            tParameters.set( "refinement_buffer", 2 );
            tParameters.set( "staircase_buffer", 1 );

             hmr::HMR tHMR( tParameters );

             std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

             // create field
             std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeOrder );

             for( uint k=0; k<3; ++k )
             {
                 tField->evaluate_scalar_function( LevelSetFunction );
                 tHMR.flag_surface_elements( tField );
                 tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
                 tHMR.update_refinement_pattern();
             }

             tHMR.finalize();

             // evaluate node values
//             tField->evaluate_scalar_function( LevelSetFunction );
//
//             tHMR.save_to_exodus( "Circle_diff.exo" );

            //1) Create the fem nodes ------------------------------------------------------
            std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------
            Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
            tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
            tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // create a list of active sidesets
            Cell< moris_index >  tSidesetList = { 3, 5 };

            // create a list of BC type for the sidesets
            Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                        fem::BC_Type::NEUMANN };

            // create model
            mdl::Model * tModel = new mdl::Model( tMesh.get(), tBplineOrder, tIWGTypeList,
                                                  tSidesetList, tSidesetBCTypeList );

            //solve
            moris::Matrix< DDRMat > tSolution11;
            tModel->solve( tSolution11 );

            tModel->output_solution( "Circle" );

            tField->put_scalar_values_on_field( tModel->get_mSolHMR() );

            tHMR.save_to_exodus( "Circle_diff_temp.exo" );

            // Expected solution
            Matrix< DDRMat > tExpectedSolution = {{ 5.000000000439168e+00,    2.499999999484336e+01,    4.499999998931914e+01,
                                                    6.499999998398192e+01,    8.499999997909634e+01,    1.049999999750879e+02,
                                                    1.249999999726532e+02,    1.349999999677999e+02,    1.349999999669917e+02,
                                                    1.349999999688249e+02,    1.349999999678056e+02,    1.449999999632745e+02,
                                                    1.549999999576012e+02,    1.449999999607157e+02,    1.549999999559840e+02,
                                                    1.449999999644306e+02,    1.549999999607020e+02,    1.449999999638932e+02,
                                                    1.549999999578126e+02,    1.649999999541855e+02,    1.749999999460777e+02,
                                                    1.649999999497775e+02,    1.749999999406942e+02,    1.649999999575286e+02,
                                                    1.749999999499552e+02,    1.649999999538873e+02,    1.749999999470990e+02,
                                                    1.849999999366943e+02,    1.949999999299542e+02,    1.849999999312525e+02,
                                                    1.949999999239139e+02,    1.849999999408975e+02,    1.949999999326252e+02,
                                                    1.849999999358517e+02,    1.949999999274455e+02,    2.049999999278416e+02,
                                                    2.049999999213750e+02,    2.049999999302681e+02,    2.049999999247313e+02,
                                                    5.000000000440735e+00,    2.499999999481117e+01,    4.499999998922714e+01,
                                                    6.499999998376413e+01,    8.499999997856834e+01,    1.049999999739319e+02 }};

            // define an epsilon environment
            double tEpsilon = 1E-8;

            // define a bool for solution check
            bool tCheckNodalSolution = true;

            // loop over the node and chyeck solution
            for ( uint i = 0; i < 45; i++ )
            {
                // check solution
                tCheckNodalSolution = tCheckNodalSolution
                                   && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
            }

            print(tSolution11, "Solution");

            // check bool is true
            REQUIRE( tCheckNodalSolution );
        }/* if( par_size() */
    }

    TEST_CASE( "Diffusion_hmr2_10x4x4", "[moris],[mdl],[Diffusion_hmr2_10x4x4]" )
    {
       if(par_size() == 1 )
       {
           // Create a 3D mesh of HEX8 using MTK ------------------------------------------
           std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
           //------------------------------------------------------------------------------

           moris::uint tBplineOrder = 2;
           moris::uint tLagrangeOrder = 1;
           moris::uint tMyCoeff = 1;

           hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

           tParameters.set( "number_of_elements_per_dimension", "2, 2, 2" );
           tParameters.set( "domain_dimensions", "2, 2, 2" );
           tParameters.set( "domain_offset", "0, 0.0, 0.0" );
           tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");
           tParameters.set( "verbose", 0 );
           tParameters.set( "truncate_bsplines", 1 );
           tParameters.set( "bspline_orders", "2" );
           tParameters.set( "lagrange_orders", "1" );

           tParameters.set( "use_multigrid", 0 );

           tParameters.set( "refinement_buffer", 2 );
           tParameters.set( "staircase_buffer", 1 );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeOrder );

            for( uint k=0; k<0; ++k )
            {
                tField->evaluate_scalar_function( LevelSetFunction );
                tHMR.flag_surface_elements( tField );
                tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
                tHMR.update_refinement_pattern();
            }

            tHMR.finalize();

           // evaluate node values
           tField->evaluate_scalar_function( LevelSetFunction );

           tHMR.save_to_exodus( 2,"Circle_diff.exo" );

           tHMR.save_bsplines_to_vtk("DLA_BSplines.vtk");

           //1) Create the fem nodes ------------------------------------------------------
           std::cout<<" Create the fem nodes "<<std::endl;
           //------------------------------------------------------------------------------
           Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
           tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
           tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
           tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

           // create a list of active sidesets
           Cell< moris_index >  tSidesetList = { 3, 5 };

           // create a list of BC type for the sidesets
           Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                       fem::BC_Type::NEUMANN };

           // create model
           mdl::Model * tModel = new mdl::Model( tMesh.get(), tBplineOrder, tIWGTypeList,
                                                 tSidesetList, tSidesetBCTypeList );

           //solve
           moris::Matrix< DDRMat > tSolution11;
           tModel->solve( tSolution11 );

           print(tSolution11,"tSolution11");

           tModel->output_solution( "Circle" );

           tField->put_scalar_values_on_field( tModel->get_mSolHMR() );

           tHMR.save_to_exodus( 2,"Circle_diff_temp.exo" );

           // Expected solution
           Matrix< DDRMat > tExpectedSolution = {{ 5.000000000439168e+00,    2.499999999484336e+01,    4.499999998931914e+01,
                                                   6.499999998398192e+01,    8.499999997909634e+01,    1.049999999750879e+02,
                                                   1.249999999726532e+02,    1.349999999677999e+02,    1.349999999669917e+02,
                                                   1.349999999688249e+02,    1.349999999678056e+02,    1.449999999632745e+02,
                                                   1.549999999576012e+02,    1.449999999607157e+02,    1.549999999559840e+02,
                                                   1.449999999644306e+02,    1.549999999607020e+02,    1.449999999638932e+02,
                                                   1.549999999578126e+02,    1.649999999541855e+02,    1.749999999460777e+02,
                                                   1.649999999497775e+02,    1.749999999406942e+02,    1.649999999575286e+02,
                                                   1.749999999499552e+02,    1.649999999538873e+02,    1.749999999470990e+02,
                                                   1.849999999366943e+02,    1.949999999299542e+02,    1.849999999312525e+02,
                                                   1.949999999239139e+02,    1.849999999408975e+02,    1.949999999326252e+02,
                                                   1.849999999358517e+02,    1.949999999274455e+02,    2.049999999278416e+02,
                                                   2.049999999213750e+02,    2.049999999302681e+02,    2.049999999247313e+02,
                                                   5.000000000440735e+00,    2.499999999481117e+01,    4.499999998922714e+01,
                                                   6.499999998376413e+01,    8.499999997856834e+01,    1.049999999739319e+02 }};

           // define an epsilon environment
           double tEpsilon = 1E-12;

           // define a bool for solution check
           bool tCheckNodalSolution = true;

           // loop over the node and chyeck solution
           for ( uint i = 0; i < 45; i++ )
           {
               // check solution
               tCheckNodalSolution = tCheckNodalSolution
                                  && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
           }
           // check bool is true
           REQUIRE( tCheckNodalSolution );
       }/* if( par_size() */
   }

    TEST_CASE( "Diffusion_hmr3_10x4x4", "[moris],[mdl],[Diffusion_hmr3_10x4x4]" )
        {
           if(par_size() == 1 )
           {
               // Create a 3D mesh of HEX8 using MTK ------------------------------------------
               std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
               //------------------------------------------------------------------------------

               moris::uint tBplineOrder = 2;
               moris::uint tLagrangeOrder = 2;
               moris::uint tMyCoeff = 1;

               hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

               tParameters.set( "number_of_elements_per_dimension", "2, 2, 2" );
               tParameters.set( "domain_dimensions", "2, 2, 2" );
               tParameters.set( "domain_offset", "0, 0.0, 0.0" );
               tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");
               tParameters.set( "verbose", 0 );
               tParameters.set( "truncate_bsplines", 1 );
               tParameters.set( "bspline_orders", "2" );
               tParameters.set( "lagrange_orders", "2" );

               tParameters.set( "use_multigrid", 0 );

               tParameters.set( "refinement_buffer", 2 );
               tParameters.set( "staircase_buffer", 1 );

                hmr::HMR tHMR( tParameters );

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

                // create field
                std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeOrder );

                for( uint k=0; k<0; ++k )
                {
                    tField->evaluate_scalar_function( LevelSetFunction );
                    tHMR.flag_surface_elements( tField );
                    tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
                    tHMR.update_refinement_pattern();
                }

                tHMR.finalize();

               // evaluate node values
               tField->evaluate_scalar_function( LevelSetFunction );

               tHMR.save_to_exodus( 2,"Circle_diff.exo" );


               //1) Create the fem nodes ------------------------------------------------------
               std::cout<<" Create the fem nodes "<<std::endl;
               //------------------------------------------------------------------------------
               Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
               tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
               tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
               tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

               // create a list of active sidesets
               Cell< moris_index >  tSidesetList = { 3, 5 };

               // create a list of BC type for the sidesets
               Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                           fem::BC_Type::NEUMANN };

               // create model
               mdl::Model * tModel = new mdl::Model( tMesh.get(), tBplineOrder, tIWGTypeList,
                                                     tSidesetList, tSidesetBCTypeList );

               //solve
               moris::Matrix< DDRMat > tSolution11;
               tModel->solve( tSolution11 );

               print(tSolution11,"tSolution11");

               tModel->output_solution( "Circle" );

               tField->put_scalar_values_on_field( tModel->get_mSolHMR() );

               tHMR.save_to_exodus( 2,"Circle_diff_temp.exo" );

               // Expected solution
               Matrix< DDRMat > tExpectedSolution = {{ 5.000000000439168e+00,    2.499999999484336e+01,    4.499999998931914e+01,
                                                       6.499999998398192e+01,    8.499999997909634e+01,    1.049999999750879e+02,
                                                       1.249999999726532e+02,    1.349999999677999e+02,    1.349999999669917e+02,
                                                       1.349999999688249e+02,    1.349999999678056e+02,    1.449999999632745e+02,
                                                       1.549999999576012e+02,    1.449999999607157e+02,    1.549999999559840e+02,
                                                       1.449999999644306e+02,    1.549999999607020e+02,    1.449999999638932e+02,
                                                       1.549999999578126e+02,    1.649999999541855e+02,    1.749999999460777e+02,
                                                       1.649999999497775e+02,    1.749999999406942e+02,    1.649999999575286e+02,
                                                       1.749999999499552e+02,    1.649999999538873e+02,    1.749999999470990e+02,
                                                       1.849999999366943e+02,    1.949999999299542e+02,    1.849999999312525e+02,
                                                       1.949999999239139e+02,    1.849999999408975e+02,    1.949999999326252e+02,
                                                       1.849999999358517e+02,    1.949999999274455e+02,    2.049999999278416e+02,
                                                       2.049999999213750e+02,    2.049999999302681e+02,    2.049999999247313e+02,
                                                       5.000000000440735e+00,    2.499999999481117e+01,    4.499999998922714e+01,
                                                       6.499999998376413e+01,    8.499999997856834e+01,    1.049999999739319e+02 }};

               // define an epsilon environment
               double tEpsilon = 1E-12;

               // define a bool for solution check
               bool tCheckNodalSolution = true;

               // loop over the node and chyeck solution
               for ( uint i = 0; i < 45; i++ )
               {
                   // check solution
                   tCheckNodalSolution = tCheckNodalSolution
                                      && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
               }
               // check bool is true
               REQUIRE( tCheckNodalSolution );
           }/* if( par_size() */
       }

    }/* namespace fem */
}/* namespace moris */
