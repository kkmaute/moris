#include "catch.hpp"
#include "cl_Matrix.hpp"

#include "fn_GEN_create_geometries.hpp"
#include "fn_GEN_create_properties.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

#define protected public
#define private   public
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Intersection_Node.hpp"
#undef protected
#undef private

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

        TEST_CASE("Interpolation PDV creation test", "[gen], [pdv], [interpolation pdv]")
        {
            if( par_size() == 1)
            {
                // Create PDV_Type host manager
                Pdv_Host_Manager tPdvHostManager;

                tPdvHostManager.mPdvTypeList = {PDV_Type::DENSITY, PDV_Type::TEMPERATURE, PDV_Type::ELASTIC_MODULUS};
                tPdvHostManager.mPdvTypeMap.set_size(10, 1, -1);
                tPdvHostManager.mPdvTypeMap( 3 ) = 0;
                tPdvHostManager.mPdvTypeMap( 4 ) = 1;
                tPdvHostManager.mPdvTypeMap( 5 ) = 2;

                // ----------------- Interpolation PDVs ---------------------- //
                // Node indices per set
                Cell<Matrix<DDSMat>> tIpNodeIndicesPerSet(2);
                tIpNodeIndicesPerSet(0).resize(4, 1);
                tIpNodeIndicesPerSet(1).resize(4, 1);
                tIpNodeIndicesPerSet(0) = {{0, 1, 2, 3}};
                tIpNodeIndicesPerSet(1) = {{2, 3, 4, 5}};

                Cell<Matrix<DDSMat>> tIpNodeIdsPerSet(2);
                tIpNodeIdsPerSet(0).resize(4, 1);
                tIpNodeIdsPerSet(1).resize(4, 1);
                tIpNodeIdsPerSet(0) = {{0, 1, 2, 3}};
                tIpNodeIdsPerSet(1) = {{2, 3, 4, 5}};

                Cell<Matrix<DDSMat>> tIpNodeOwnersPerSet(2);
                tIpNodeOwnersPerSet(0).resize(4, 1);
                tIpNodeOwnersPerSet(1).resize(4, 1);
                tIpNodeOwnersPerSet(0) = {{0, 0, 0, 0}};
                tIpNodeOwnersPerSet(1) = {{0, 0, 0, 0}};

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
                tPdvHostManager.create_interpolation_pdv_hosts(
                        tIpNodeIndicesPerSet,
                        tIpNodeIdsPerSet,
                        tIpNodeOwnersPerSet,
                        Cell<Matrix<F31RMat>>(6),
                        tIpPdvTypes);

                // Set PDVs
                for (uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++)
                {
                    for (uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++)
                    {
                        for (uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++)
                        {
                            tPdvHostManager.create_interpolation_pdv(
                                    (uint)tIpNodeIndicesPerSet(tMeshSetIndex)(tNodeIndex),
                                    tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0),
                                    (real)tMeshSetIndex);
                        }
                    }
                }

                tPdvHostManager.create_pdv_ids();

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
                            CHECK(tPdvValues(0)(tNodeIndex) == tMeshSetIndex +
                                    (tMeshSetIndex == 0) * (tNodeIndex > 1) *
                                    (tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0) == PDV_Type::TEMPERATURE));
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
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Interpolation PDV creation test parallel", "[gen], [pdv], [interpolation pdv parallel]")
        {
            if( par_size() == 2)
            {
                // Create PDV_Type host manager
                Pdv_Host_Manager tPdvHostManager;

                tPdvHostManager.mPdvTypeList = {PDV_Type::DENSITY, PDV_Type::TEMPERATURE};
                tPdvHostManager.mPdvTypeMap.set_size(10, 1, -1);
                tPdvHostManager.mPdvTypeMap( 3 ) = 0;
                tPdvHostManager.mPdvTypeMap( 4 ) = 1;

                // ----------------- Interpolation PDVs ---------------------- //

                Cell<Matrix<DDSMat>> tIpNodeIndicesPerSet(1);
                tIpNodeIndicesPerSet(0).resize(4, 1);

                Cell<Matrix<DDSMat>> tIpNodeIdsPerSet(1);
                tIpNodeIdsPerSet(0).resize(4, 1);

                Cell<Matrix<DDSMat>> tIpNodeOwnersPerSet(1);
                tIpNodeOwnersPerSet(0).resize(4, 1);

                // PDV_Type types per set
                Cell<Cell<Cell<PDV_Type>>> tIpPdvTypes(1);
                tIpPdvTypes(0).resize(2);
                tIpPdvTypes(0)(0).resize(1);
                tIpPdvTypes(0)(1).resize(1);
                tIpPdvTypes(0)(0)(0) = PDV_Type::DENSITY;
                tIpPdvTypes(0)(1)(0) = PDV_Type::TEMPERATURE;

                if( par_rank() == 0)
                {
                    // Node indices per set
                    tIpNodeIndicesPerSet(0) = {{0, 1, 2, 3}};

                    tIpNodeIdsPerSet(0) = {{0, 1, 2, 3}};

                    tIpNodeOwnersPerSet(0) = {{0, 0, 0, 1}};

                    tPdvHostManager.mCommTable.set_size( 2, 1, 0);
                    tPdvHostManager.mCommTable( 1, 0 ) = 1;

                    tPdvHostManager.mIPVertexIdtoIndMap[ 2 ] = 2;
                }
                else if( par_rank() == 1)
                {
                    // Node indices per set
                    tIpNodeIndicesPerSet(0) = {{0, 1, 2, 3}};

                    tIpNodeIdsPerSet(0) = {{2, 3, 4, 5}};

                    tIpNodeOwnersPerSet(0) = {{0, 1, 1, 1}};

                    tPdvHostManager.mCommTable.set_size( 2, 1, 1);
                    tPdvHostManager.mCommTable( 1, 0 ) = 0;

                    tPdvHostManager.mIPVertexIdtoIndMap[ 3 ] = 1;
                }

                // Create PDV_Type hosts
                tPdvHostManager.create_interpolation_pdv_hosts(
                        tIpNodeIndicesPerSet,
                        tIpNodeIdsPerSet,
                        tIpNodeOwnersPerSet,
                        Cell<Matrix<F31RMat>>(4),
                        tIpPdvTypes);

                // Set PDVs
                for (uint tMeshSetIndex = 0; tMeshSetIndex < 1; tMeshSetIndex++)
                {
                    for (uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++)
                    {
                        for (uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++)
                        {
                            tPdvHostManager.create_interpolation_pdv(
                                    (uint)tIpNodeIndicesPerSet(tMeshSetIndex)(tNodeIndex),
                                    tIpPdvTypes(tMeshSetIndex)(tPdvIndex)(0),
                                    (real)tMeshSetIndex);
                        }
                    }
                }

                tPdvHostManager.create_pdv_ids();

                // ------------------- Check global map ----------------------- //
                const Matrix<DDSMat> & tLocalGlobalMap = tPdvHostManager.get_my_local_global_map();
                const Matrix<DDSMat> & tLocalGlobalOSMap = tPdvHostManager.get_my_local_global_overlapping_map();

                //print( tLocalGlobalMap, "tLocalGlobalMap");
                //print( tLocalGlobalOSMap, "tLocalGlobalOSMap");

                REQUIRE(tLocalGlobalMap.length() == 6);
                REQUIRE(tLocalGlobalOSMap.length() == 8);

                if( par_rank() == 0)
                {
                    CHECK(tLocalGlobalMap(0) == 0);                    CHECK(tLocalGlobalMap(1) == 1);
                    CHECK(tLocalGlobalMap(2) == 2);                    CHECK(tLocalGlobalMap(3) == 3);
                    CHECK(tLocalGlobalMap(4) == 4);                    CHECK(tLocalGlobalMap(5) == 5);

                    CHECK(tLocalGlobalOSMap(0) == 0);                  CHECK(tLocalGlobalOSMap(1) == 1);
                    CHECK(tLocalGlobalOSMap(2) == 2);                  CHECK(tLocalGlobalOSMap(3) == 6);
                    CHECK(tLocalGlobalOSMap(4) == 3);                  CHECK(tLocalGlobalOSMap(5) == 4);
                    CHECK(tLocalGlobalOSMap(6) == 5);                  CHECK(tLocalGlobalOSMap(7) == 9);
                }
                if( par_rank() == 1)
                {
                    CHECK(tLocalGlobalMap(0) == 6);                    CHECK(tLocalGlobalMap(1) == 7);
                    CHECK(tLocalGlobalMap(2) == 8);                    CHECK(tLocalGlobalMap(3) == 9);
                    CHECK(tLocalGlobalMap(4) == 10);                   CHECK(tLocalGlobalMap(5) == 11);

                    CHECK(tLocalGlobalOSMap(0) == 2);                  CHECK(tLocalGlobalOSMap(1) == 6);
                    CHECK(tLocalGlobalOSMap(2) == 7);                  CHECK(tLocalGlobalOSMap(3) == 8);
                    CHECK(tLocalGlobalOSMap(4) == 5);                  CHECK(tLocalGlobalOSMap(5) == 9);
                    CHECK(tLocalGlobalOSMap(6) == 10);                 CHECK(tLocalGlobalOSMap(7) == 11);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("PDV sensitivities test", "[gen], [pdv], [sensitivity], [pdv sensitivity]")
        {
            // Create PDV_Type host manager
            Pdv_Host_Manager tPdvHostManager;
            tPdvHostManager.set_num_advs(tNumADVs);

            tPdvHostManager.mPdvTypeList = {PDV_Type::DENSITY };
            tPdvHostManager.mPdvTypeMap.set_size(10, 1, -1);
            tPdvHostManager.mPdvTypeMap( 3 ) = 0;

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
            Cell<Matrix<DDSMat>> tIpNodeIdsPerSet(1);
            Cell<Matrix<DDSMat>> tIpNodeOWnersPerSet(1);
            tIpNodeIndicesPerSet(0).set_size(tNumADVs, 1);
            tIpNodeIdsPerSet(0).set_size(tNumADVs, 1);
            tIpNodeOWnersPerSet(0).set_size(tNumADVs, 1, 0);
            for (uint tNodeIndex = 0; tNodeIndex < tNumADVs; tNodeIndex++)
            {
                tIpNodeIndicesPerSet(0)(tNodeIndex) = tNodeIndex;
                tIpNodeIdsPerSet(0)(tNodeIndex) = tNodeIndex;
            }

            // PDV_Type types per set
            Cell<Cell<Cell<PDV_Type>>> tIpPdvTypes(1);
            tIpPdvTypes(0).resize(1);
            tIpPdvTypes(0)(0).resize(1);
            tIpPdvTypes(0)(0)(0) = PDV_Type::DENSITY;

            // Create PDV_Type hosts
            tPdvHostManager.create_interpolation_pdv_hosts(
                    tIpNodeIndicesPerSet,
                    tIpNodeIdsPerSet,
                    tIpNodeOWnersPerSet,
                    Cell<Matrix<F31RMat>>(tNumADVs),
                    tIpPdvTypes);

            // Set PDVs
            for (uint tNodeIndex = 0; tNodeIndex < tNumADVs; tNodeIndex++)
            {
                tPdvHostManager.create_interpolation_pdv(uint(tIpNodeIndicesPerSet(0)(tNodeIndex)), tIpPdvTypes(0)(0)(0), tProperty);
            }

            tPdvHostManager.create_pdv_ids();

            // Check sensitivities
            CHECK(norm(tPdvHostManager.compute_diqi_dadv() - tDiqiDpdv) <= 1E-12);
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Intersection PDV creation test parallel", "[gen], [pdv], [intersection pdv parallel]")
                {
                    if( par_size() == 2)
                    {
                        // Create PDV_Type host manager
                        Pdv_Host_Manager tPdvHostManager;

                        Matrix<IdMat> tIpNodeIdsPerSet(4,1);

                        Matrix<DDSMat> tIpNodeOwnersPerSet(4,1);

                        Cell< std::shared_ptr<Intersection_Node> > tIntersectionNodes(4);

                        if( par_rank() == 0)
                        {
                            tIpNodeIdsPerSet = {{0}, {1}, {2}, {3}};

                            tIpNodeOwnersPerSet = {{0}, {0}, {0}, {1}};

                            tPdvHostManager.mCommTable.set_size( 2, 1, 0);
                            tPdvHostManager.mCommTable( 1, 0 ) = 1;

                            tPdvHostManager.mIGVertexIdtoIndMap[ 2 ] = 2;
                        }
                        else if( par_rank() == 1)
                        {
                            tIpNodeIdsPerSet = {{2}, {3}, {4}, {5}};

                            tIpNodeOwnersPerSet = {{0}, {1}, {1}, {1}};

                            tPdvHostManager.mCommTable.set_size( 2, 1, 1);
                            tPdvHostManager.mCommTable( 1, 0 ) = 0;

                            tPdvHostManager.mIGVertexIdtoIndMap[ 3 ] = 1;
                        }


                        for( sint Ik = 0; Ik < 4; Ik++ )
                        {
                            tIntersectionNodes( Ik ) = std::make_shared<Intersection_Node>();

                            tIntersectionNodes( Ik )->mGlobalCoordinates.set_size( 2, 1 );
                            tPdvHostManager.set_intersection_node( Ik, tIntersectionNodes( Ik ) );
                            tPdvHostManager.update_intersection_node( Ik, tIpNodeIdsPerSet( Ik ), tIpNodeOwnersPerSet( Ik ));
                        }

                        tPdvHostManager.create_pdv_ids();

                        // ------------------- Check global map ----------------------- //
                        const Matrix<DDSMat> & tLocalGlobalMap = tPdvHostManager.get_my_local_global_map();
                        const Matrix<DDSMat> & tLocalGlobalOSMap = tPdvHostManager.get_my_local_global_overlapping_map();

                        //print( tLocalGlobalMap, "tLocalGlobalMap");
                        //print( tLocalGlobalOSMap, "tLocalGlobalOSMap");

                        REQUIRE(tLocalGlobalMap.length() == 6);
                        REQUIRE(tLocalGlobalOSMap.length() == 8);

                        if( par_rank() == 0)
                        {
                            CHECK(tLocalGlobalMap(0) == 0);                    CHECK(tLocalGlobalMap(1) == 1);
                            CHECK(tLocalGlobalMap(2) == 2);                    CHECK(tLocalGlobalMap(3) == 3);
                            CHECK(tLocalGlobalMap(4) == 4);                    CHECK(tLocalGlobalMap(5) == 5);

                            CHECK(tLocalGlobalOSMap(0) == 0);                  CHECK(tLocalGlobalOSMap(1) == 1);
                            CHECK(tLocalGlobalOSMap(2) == 2);                  CHECK(tLocalGlobalOSMap(3) == 3);
                            CHECK(tLocalGlobalOSMap(4) == 4);                  CHECK(tLocalGlobalOSMap(5) == 5);
                            CHECK(tLocalGlobalOSMap(6) == 6);                  CHECK(tLocalGlobalOSMap(7) == 7);
                        }
                        if( par_rank() == 1)
                        {
                            CHECK(tLocalGlobalMap(0) == 6);                    CHECK(tLocalGlobalMap(1) == 7);
                            CHECK(tLocalGlobalMap(2) == 8);                    CHECK(tLocalGlobalMap(3) == 9);
                            CHECK(tLocalGlobalMap(4) == 10);                   CHECK(tLocalGlobalMap(5) == 11);

                            CHECK(tLocalGlobalOSMap(0) == 4);                  CHECK(tLocalGlobalOSMap(1) == 5);
                            CHECK(tLocalGlobalOSMap(2) == 6);                  CHECK(tLocalGlobalOSMap(3) == 7);
                            CHECK(tLocalGlobalOSMap(4) == 8);                  CHECK(tLocalGlobalOSMap(5) == 9);
                            CHECK(tLocalGlobalOSMap(6) == 10);                 CHECK(tLocalGlobalOSMap(7) == 11);
                        }
                    }
                }

    }
}
