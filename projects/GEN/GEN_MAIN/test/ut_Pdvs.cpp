#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"

namespace moris
{
    namespace ge
    {
        TEST_CASE("PDV_Type creation through host manager", "[GE], [PDV_Type]")
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
            tPdvHostManager.create_ip_pdv_hosts(6, tIpNodeIndicesPerSet, tIpPdvTypes);

            // Set PDVs
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++)
            {
                for (uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++)
                    {
                        tPdvHostManager.create_ip_pdv(uint(tIpNodeIndicesPerSet(tMeshSetIndex)(tNodeIndex)), tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0), real(tMeshSetIndex));
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
                        CHECK(tPdvValues(0)(tNodeIndex) == tMeshSetIndex + (tMeshSetIndex == 0) * (tNodeIndex > 1) * (tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0) == PDV_Type::TEMPERATURE));
                    }
                }
            }

            // ----------------- Integration PDVs ---------------------- //
            // Node indices per set
            Cell<Matrix<DDSMat>> tIgNodeIndicesPerSet(3);
            tIgNodeIndicesPerSet(0).resize(4, 1);
            tIgNodeIndicesPerSet(1).resize(3, 1);
            tIgNodeIndicesPerSet(2).resize(3, 1);
            tIgNodeIndicesPerSet(0) = {{0, 1, 2, 3}};
            tIgNodeIndicesPerSet(1) = {{0, 4, 5}};
            tIgNodeIndicesPerSet(2) = {{1, 2, 6, 7}};

            // IG PDV_Type types
            Cell<PDV_Type> tCoordinatePdvs(3);
            tCoordinatePdvs(0) = PDV_Type::X_COORDINATE;
            tCoordinatePdvs(1) = PDV_Type::Y_COORDINATE;
            tCoordinatePdvs(2) = PDV_Type::Z_COORDINATE;

            // PDV_Type types per set
            Cell<Cell<Cell<PDV_Type>>> tIgPdvTypes(3);
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 3; tMeshSetIndex++)
            {
                tIgPdvTypes(tMeshSetIndex).resize(1);
                tIgPdvTypes(tMeshSetIndex)(0) = tCoordinatePdvs;
            }

            // Create PDV_Type hosts
            tPdvHostManager.create_ig_pdv_hosts(8, tIgNodeIndicesPerSet, tIgPdvTypes);

            // Set PDVs
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 3; tMeshSetIndex++)
            {
                for (uint tNodeIndex = 0; tNodeIndex < tIgNodeIndicesPerSet(tMeshSetIndex).length(); tNodeIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < 3; tPdvIndex++)
                    {
                        tPdvHostManager.create_ig_pdv(uint(tIgNodeIndicesPerSet(tMeshSetIndex)(tNodeIndex)), tIgPdvTypes(tMeshSetIndex)(0)(tPdvIndex), real(tIgPdvTypes(tMeshSetIndex)(0)(tPdvIndex)));
                    }
                }
            }

            // Check PDVs
            for (uint tMeshSetIndex = 0; tMeshSetIndex < 3; tMeshSetIndex++)
            {
                tPdvValues.clear();
                tPdvHostManager.get_ig_pdv_value(tIgNodeIndicesPerSet(tMeshSetIndex), tIgPdvTypes(tMeshSetIndex)(0), tPdvValues);
                for (uint tNodeIndex = 0; tNodeIndex < tIgNodeIndicesPerSet(tMeshSetIndex).length(); tNodeIndex++)
                {
                    for (uint tDimension = 0; tDimension < 3; tDimension++)
                    {
                        CHECK(tPdvValues(tDimension)(tNodeIndex) == real(tDimension));
                    }
                }
            }

            // ------------------- Check global map ----------------------- //
            const Matrix<DDSMat> & tLocalGlobalMap = tPdvHostManager.get_my_local_global_map();

            CHECK(tLocalGlobalMap.length() == 38);
            for (int tGlobalPdvIndex = 0; tGlobalPdvIndex < 38; tGlobalPdvIndex++)
            {
                CHECK(tLocalGlobalMap(tGlobalPdvIndex) == tGlobalPdvIndex);
            }
        }
    }
}
