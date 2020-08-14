#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_GEN_create_properties.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "fn_GEN_create_simple_mesh.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        // Dummy values so I don't need to create a model for the sensitivity test
        uint tNumADVs = 36;
        Matrix<DDRMat> tDiqiDpdv(1, tNumADVs, 1.0);
        Matrix<DDRMat> Pdv_Host_Manager::compute_diqi_dadv()
        {
            return tDiqiDpdv * this->compute_dpdv_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Interpolation PDV creation", "[gen], [pdv], [interpolation pdv]")
        {
            // Create PDV_Type host manager
            Pdv_Host_Manager tPdvHostManager;
            
            // ----------------- Interpolation PDVs ---------------------- //
            // Node indices per set
            Cell<Matrix<DDSMat>> tIpNodeIndicesPerSet(2);
            tIpNodeIndicesPerSet(0).resize(4, 1);
            tIpNodeIndicesPerSet(1).resize(4, 1);
            tIpNodeIndicesPerSet(0) = {{0, 1, 2, 3}};
            tIpNodeIndicesPerSet(1) = {{2, 3, 4, 5}};

            // PDV_Type types per set
            Cell<Cell<Cell<PDV_Type>>> tIpPdvTypes(2);
            tIpPdvTypes(0).resize(2);
            tIpPdvTypes(1).resize(2);
            tIpPdvTypes(0)(0).resize(1);
            tIpPdvTypes(0)(1).resize(1);
            tIpPdvTypes(1)(0).resize(1);
            tIpPdvTypes(1)(1).resize(1);
            tIpPdvTypes(0)(0)(0) = PDV_Type::DENSITY;
            tIpPdvTypes(0)(1)(0) = PDV_Type::TEMPERATURE;
            tIpPdvTypes(1)(0)(0) = PDV_Type::TEMPERATURE;
            tIpPdvTypes(1)(1)(0) = PDV_Type::ELASTIC_MODULUS;

            // Create PDV_Type hosts
            tPdvHostManager.create_ip_pdv_hosts(tIpNodeIndicesPerSet, Cell<Matrix<F31RMat>>(6), tIpPdvTypes);

            // Set PDVs
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++)
            {
                for (uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++)
                    {
                        tPdvHostManager.create_ip_pdv(
                                (uint)tIpNodeIndicesPerSet(tMeshSetIndex)(tNodeIndex),
                                tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0),
                                (real)tMeshSetIndex);
                    }
                }
            }

            // Check PDVs
            Cell<Matrix<DDRMat>> tPdvValues;
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++)
            {
                for (uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++)
                {
                    tPdvValues.clear();
                    tPdvHostManager.get_ip_pdv_value(tIpNodeIndicesPerSet(tMeshSetIndex), tIpPdvTypes(tMeshSetIndex)(tPdvIndex), tPdvValues);
                    for (uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++)
                    {
                        CHECK(tPdvValues(0)(tNodeIndex) == tMeshSetIndex + (tMeshSetIndex == 0) * (tNodeIndex > 1)
                        * (tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0) == PDV_Type::TEMPERATURE));
                    }
                }
            }

            // ------------------- Check global map ----------------------- //
            const Matrix<DDSMat> & tLocalGlobalMap = tPdvHostManager.get_my_local_global_map();

            REQUIRE(tLocalGlobalMap.length() == 14);
            for (int tGlobalPdvIndex = 0; tGlobalPdvIndex < 14; tGlobalPdvIndex++)
            {
                CHECK(tLocalGlobalMap(tGlobalPdvIndex) == tGlobalPdvIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Intersection PDV creation", "[gen], [pdv], [intersection pdv]")
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh();

            // Set up geometry
            Cell<std::shared_ptr<Geometry>> tGeometries(2);
            Matrix<DDRMat> tADVs(0, 0);

            // Circle
            real tRadius = 0.5;
            ParameterList tCircleParameterList = prm::create_geometry_parameter_list();
            tCircleParameterList.set("type", "circle");
            tCircleParameterList.set("constant_parameters", "0.0, 0.0, " + std::to_string(tRadius));
            tGeometries(0) = create_geometry(tCircleParameterList, tADVs);

            // Plane
            ParameterList tPlaneParameterList = prm::create_geometry_parameter_list();
            tPlaneParameterList.set("type", "plane");
            tPlaneParameterList.set("constant_parameters", "0.25, 0.0, 1.0, 0.0");
            tGeometries(1) = create_geometry(tPlaneParameterList, tADVs);

            // Create geometry engine
            Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine tGeometriesEngine(tGeometries, tPhaseTable, tMesh);

            // Determine if intersected
            for (uint tElementIndex = 0; tElementIndex < tMesh->get_num_elems(); tElementIndex++)
            {
                CHECK(tGeometriesEngine.is_intersected(
                        tMesh->get_mtk_cell(tElementIndex).get_vertex_inds(),
                        tMesh->get_mtk_cell(tElementIndex).get_vertex_coords()));
            }

            // Clean up
            delete tMesh;
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("PDV sensitivities test", "[gen], [pdv], [sensitivity], [pdv sensitivity]")
        {
            // Create PDV_Type host manager
            Pdv_Host_Manager tPdvHostManager;
            tPdvHostManager.set_num_advs(tNumADVs);

            // Create discrete property
            ParameterList tParameterList = moris::prm::create_gen_property_parameter_list();;
            tParameterList.set("type", "discrete");
            tParameterList.set("property_variable_indices", "all");
            tParameterList.set("adv_indices", "all");
            tParameterList.set("pdv_type", "DENSITY");

            // Create property
            Matrix<DDRMat> tADVs(tNumADVs, 1);
            std::shared_ptr<Property> tProperty = create_property(tParameterList, tADVs, Cell<std::shared_ptr<moris::ge::Property>>(0));

            // Node indices per set
            Cell<Matrix<DDSMat>> tIpNodeIndicesPerSet(1);
            tIpNodeIndicesPerSet(0).resize(tNumADVs, 1);
            for (uint tNodeIndex = 0; tNodeIndex < tNumADVs; tNodeIndex++)
            {
                tIpNodeIndicesPerSet(0)(tNodeIndex) = tNodeIndex;
            }

            // PDV_Type types per set
            Cell<Cell<Cell<PDV_Type>>> tIpPdvTypes(1);
            tIpPdvTypes(0).resize(1);
            tIpPdvTypes(0)(0).resize(1);
            tIpPdvTypes(0)(0)(0) = PDV_Type::DENSITY;

            // Create PDV_Type hosts
            tPdvHostManager.create_ip_pdv_hosts(tIpNodeIndicesPerSet, Cell<Matrix<F31RMat>>(tNumADVs), tIpPdvTypes);

            // Set PDVs
            for (uint tNodeIndex = 0; tNodeIndex < tNumADVs; tNodeIndex++)
            {
                tPdvHostManager.create_ip_pdv(uint(tIpNodeIndicesPerSet(0)(tNodeIndex)), tIpPdvTypes(0)(0)(0), tProperty);
            }

            // Check sensitivities
            CHECK(norm(tPdvHostManager.compute_diqi_dadv() - tDiqiDpdv) <= 1E-12);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
