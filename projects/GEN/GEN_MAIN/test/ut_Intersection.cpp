#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Level_Set.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_XTK_Edge_Topology.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_PRM_HMR_Parameters.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Interface sensitivity test", "[GEN], [interface], [sensitivity], [interface sensitivity]")
        {
            if (par_size() == 1)
            {
                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string("2, 2"));
                tParameters.set( "domain_dimensions", std::string("2, 2") );
                tParameters.set( "domain_offset", std::string("-1.0, -1.0") );
                tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
                tParameters.set( "lagrange_output_meshes", std::string("0") );

                tParameters.set( "lagrange_orders", std::string("2") );
                tParameters.set( "lagrange_pattern", std::string("0") );
                tParameters.set( "bspline_orders", std::string("2") );
                tParameters.set( "bspline_pattern", std::string("0") );

                tParameters.set( "lagrange_to_bspline", std::string("0") );

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 3 );
                tParameters.set( "staircase_buffer", 3 );
                tParameters.set( "initial_refinement", 0 );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement( 0 );
                tHMR.finalize();

                hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(0);

                // Create geometry
                real tRadius = 0.25;
                Matrix<DDRMat> tADVs = {{0.0, 0.0, tRadius,
                        -1.0, -1.0, -1.0, -1.0,
                        -1.0, -1.0, 1.0, 0.0,
                        -0.5, 0.5, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0}};
                Cell<std::shared_ptr<Geometry>> tGeometry(2);
                tGeometry(0) = std::make_shared<Circle>(tADVs,
                                                        Matrix<DDUMat>({{0, 1, 2}}),
                                                        Matrix<DDUMat>({{0, 1, 2}}),
                                                        Matrix<DDRMat>(0, 0));

                tGeometry(1) = std::make_shared<Level_Set>(tADVs,
                    Matrix<DDUMat>({{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}}),
                    Matrix<DDUMat>({{3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}}),
                    Matrix<DDRMat>(0, 0),
                    tInterpolationMesh);

                Phase_Table tPhaseTable (2, Phase_Table_Structure::EXP_BASE_2);
                Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable, tInterpolationMesh, tADVs);

                xtk::Model tXTKModel(2, tInterpolationMesh, &tGeometryEngine);
                tXTKModel.mVerbose = false;

                //Specify decomposition Method and Cut Mesh ---------------------------------------
                Cell<Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
                tXTKModel.decompose(tDecompositionMethods);

                tXTKModel.perform_basis_enrichment( EntityRank::NODE, 0 );

                xtk::Enriched_Interpolation_Mesh &tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh &tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

                // Write mesh
//                mtk::Writer_Exodus writer( &tEnrIntegMesh );
//                writer.write_mesh("", "./xtk_temp.exo");
//                writer.close_file();

                // place the pair in mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();
                tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);

                // Create PDVs on integration mesh
                tGeometryEngine.create_pdvs(tMeshManager);

                // Make sure PDVs contain sensitivity information
                Pdv_Host_Manager *tPdvHostManager = dynamic_cast<Pdv_Host_Manager*>(tGeometryEngine.get_design_variable_interface());
                Matrix<DDRMat> tTempSolution = 
                        {{1.207106781186547, 1.207106781186547, 0.0, 0.0, -1.345092053441875e-01, -3.220092053441875e-01,
                          -2.494445008120904e-01, -6.385090209472533e-01, -7.997163847591288e-01, -4.959036601081719e-01,
                          -5.813368245917696e-01, -2.062115710797114, -1.930816000249107, -2.852889095279686e-01,
                          -4.247331928546001e-01, -1.541033041547084, -1.312946098940674, -1.113944547639843e-01,
                          -9.284086907492073e-02}};
                CHECK(norm(tPdvHostManager->compute_diqi_dadv() - tTempSolution) < 1E-8);

                delete tInterpolationMesh;

            }
        }
    }
}
