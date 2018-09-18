/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#include <memory>
#include <mpi.h>
#include "catch.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

// XTKL: Geometry  Include
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"

#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/fn_write_element_ownership_as_field.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "geometry/cl_Composite_Fiber.hpp"
#include "geometry/cl_Gyroid.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "geomeng/fn_Triangle_Geometry.hpp" // For surface normals


namespace xtk
{

TEST_CASE("XFEM TOOLKIT CORE TESTING SERIAL","[XTK][SERIAL]")
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize==1)
    {


        SECTION("Regular Subdivision Method","[XTK] [REGULAR_SUBDIVISION]"){
            // Geometry Engine Setup -----------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius = 0.25;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0;
            Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

            Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

            // Create Mesh ---------------------------------
            std::string tMeshFileName = "generated:1x1x1";
            Cell<std::string> tScalarFields(0);
            mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
            std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);

            // Setup XTK Model -----------------------------
            size_t tModelDimension = 3;
            Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

            //Specify your decomposition methods and start cutting
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
            tXTKModel.decompose(tDecompositionMethods);

            // Access the decomposed XTK Mesh
            Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

            // Do some testing
            size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
            size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

            // 1 element was subdivided
            CHECK(tNumNodesAfterDecompositionXTK == 15);
            CHECK(tNumElementsAfterDecompositionXTK == 24);

            moris::Matrix< Default_Matrix_Real > tNodeCoordinates = tMeshData->get_all_node_coordinates_loc_inds();
            moris::Matrix< Default_Matrix_Real > tExpectedNodeCoordinates(
            {{0, 0, 0},
                {1, 0, 0},
                {0, 1, 0},
                {1, 1, 0},
                {0, 0, 1},
                {1, 0, 1},
                {0, 1, 1},
                {1, 1, 1},
                {0.5, 0, 0.5},
                {1, 0.5, 0.5},
                {0.5, 1, 0.5},
                {0, 0.5, 0.5},
                {0.5, 0.5, 0},
                {0.5, 0.5, 1},
                {0.5, 0.5, 0.5}});
            CHECK(equal_to(tNodeCoordinates,tExpectedNodeCoordinates));

            //TODO: MORE TESTING.

            /*
             * Tests to finish
             *  - Signed volume (check for inverted jacobians)
             *  - Signed surface area
             *  - Make sure paraview or trelis can compute volume
             */
            std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);


            std::string tPrefix = std::getenv("XTKOUTPUT");
            std::string tMeshOutputFile = tPrefix + "/xtk_test_output_regular_subdivision.e";

            tCutMeshData->write_output_mesh(tMeshOutputFile);



        }

        SECTION("Regular Subdivision and Nodal Hierarchy Subdivision","[XTK][CONFORMAL]"){
            // Geometry Engine Setup ---------------------------------------------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius  = 0.25;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0.0;
            Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
            Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:1x1x2";
            Cell<std::string> tScalarFields(0);
            mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
            std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            // Access the Cut Mesh-------------------------------------------------------------
            Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

            // Do some testing
            size_t tNumNodesAfterDecompositionSTK = tMeshData->get_num_entities(EntityRank::NODE);
            size_t tNumElementsAfterDecompositionSTK = tMeshData->get_num_entities(EntityRank::ELEMENT);
            size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
            size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

            CHECK(tNumNodesAfterDecompositionXTK == 22 );
            CHECK(tNumElementsAfterDecompositionXTK == 42);

            // XTK does not add any additional Elements to the stk mesh
            CHECK(tNumElementsAfterDecompositionSTK == 2);
            // But it does add nodes because this is where the coordinates are stored
            // There are 4 more nodes in the STK mesh because only 1/2 elements are decomposed
            CHECK(tNumNodesAfterDecompositionSTK == 26);

            moris::Matrix< Default_Matrix_Real > tNodeCoordinates = tMeshData->get_all_node_coordinates_loc_inds();
            moris::Matrix< Default_Matrix_Real > tExpectedNodeCoordinates(
           {{+0.000000000000000e+00,  +0.000000000000000e+00,  +0.000000000000000e+00},
            {+1.000000000000000e+00,  +0.000000000000000e+00,  +0.000000000000000e+00},
            {+0.000000000000000e+00,  +1.000000000000000e+00,  +0.000000000000000e+00},
            {+1.000000000000000e+00,  +1.000000000000000e+00,  +0.000000000000000e+00},
            {+0.000000000000000e+00,  +0.000000000000000e+00,  +1.000000000000000e+00},
            {+1.000000000000000e+00,  +0.000000000000000e+00,  +1.000000000000000e+00},
            {+0.000000000000000e+00,  +1.000000000000000e+00,  +1.000000000000000e+00},
            {+1.000000000000000e+00,  +1.000000000000000e+00,  +1.000000000000000e+00},
            {+0.000000000000000e+00,  +0.000000000000000e+00,  +2.000000000000000e+00},
            {+1.000000000000000e+00,  +0.000000000000000e+00,  +2.000000000000000e+00},
            {+0.000000000000000e+00,  +1.000000000000000e+00,  +2.000000000000000e+00},
            {+1.000000000000000e+00,  +1.000000000000000e+00,  +2.000000000000000e+00},
            {+5.000000000000000e-01,  +0.000000000000000e+00,  +5.000000000000000e-01},
            {+1.000000000000000e+00,  +5.000000000000000e-01,  +5.000000000000000e-01},
            {+5.000000000000000e-01,  +1.000000000000000e+00,  +5.000000000000000e-01},
            {+0.000000000000000e+00,  +5.000000000000000e-01,  +5.000000000000000e-01},
            {+5.000000000000000e-01,  +5.000000000000000e-01,  +0.000000000000000e+00},
            {+5.000000000000000e-01,  +5.000000000000000e-01,  +1.000000000000000e+00},
            {+5.000000000000000e-01,  +5.000000000000000e-01,  +5.000000000000000e-01},
            {+1.000000000000000e+00,  +9.375000000000000e-01,  +0.000000000000000e+00},
            {+1.000000000000000e+00,  +1.000000000000000e+00,  +6.250000000000000e-02},
            {+9.375000000000000e-01,  +1.000000000000000e+00,  +0.000000000000000e+00},
            {+1.000000000000000e+00,  +9.687500000000000e-01,  +3.125000000000000e-02},
            {+9.687500000000000e-01,  +1.000000000000000e+00,  +3.125000000000000e-02},
            {+9.687500000000000e-01,  +9.687500000000000e-01,  +0.000000000000000e+00},
            {+9.791666666666666e-01,  +9.791666666666666e-01,  +2.083333333333334e-02}}

            );

            CHECK(equal_to(tNodeCoordinates,tExpectedNodeCoordinates));

            moris::Matrix< Default_Matrix_Real > tSingleCoord(1,3);
            moris::Matrix< Default_Matrix_Real > tLevelSetValues(1,tNumNodesAfterDecompositionXTK);
            for(size_t i = 0;  i<tNumNodesAfterDecompositionXTK; i++)
            {
                (tLevelSetValues)(0,i) = tLevelsetSphere.evaluate_field_value_with_coordinate(i,tNodeCoordinates);
            }

            moris::Matrix< Default_Matrix_Real > tExpectedLevelSetValues({{1.9375, 0.9375, 0.9375, -0.0625, 2.9375, 1.9375, 1.9375, 0.9375, 5.9375, 4.9375, 4.9375, 3.9375, 1.4375, 0.4375, 0.4375, 1.4375, 0.4375, 1.4375, 0.6875, -0.05859375, -0.05859375, -0.05859375}});
            CHECK(xtk::equal_to(tLevelSetValues,tExpectedLevelSetValues));


            /*
             * Tests to finish
             *  - Signed volume (check for inverted jacobians)
             *  - Signed surface area
             */

            std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);
            std::string tPrefix = std::getenv("XTKOUTPUT");
            std::string tMeshOutputFile = tPrefix + "/unit_test_conformal.e";
            tCutMeshData->write_output_mesh(tMeshOutputFile);

        }
    }
}

TEST_CASE("XFEM TOOLKIT CORE TESTING PARALLEL","[XTK][PARALLEL]")
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize!=1)
    {
    SECTION("Regular Subdivision Method Parallel","[XTK][Parallel][n2]")
    {
        //     Geometry Engine Setup -----------------------
        //     Using a Levelset Sphere as the Geometry
        real tRadius = 0.25;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 0;
        Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
        Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:1x1x6";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);
        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model<real, size_t,Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
        tXTKModel.decompose(tDecompositionMethods);

        Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

        CHECK(tCutMesh.get_num_simple_meshes() == 1);

        XTK_INFO<<"Number of Nodes: "<< tCutMesh.get_num_entities(EntityRank::NODE)<<std::endl;
        XTK_INFO<<"Number of Elements: "<< tCutMesh.get_num_entities(EntityRank::ELEMENT)<<std::endl;
        //TODO: MORE TESTING

    }

    SECTION("Regular Subdivision and Node Hierarchy Method Parallel","[XTK][Parallel][n2]"){
        if(tProcSize == 2)
        {
        // Geometry Engine Setup -----------------------
        // Using a Levelset Sphere as the Geometry
        real tRadius =  0.6;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 1.0;
        Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
        Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:1x1x2";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

        CHECK(tCutMesh.get_num_simple_meshes() == 2);

        //TODO: MORE TESTING
        // -signed volume
        // -signed surface area
        XTK_INFO<<"Number of Nodes both: "<< tCutMesh.get_num_entities(EntityRank::NODE);
        XTK_INFO<<"Number of Elements both: " << tCutMesh.get_num_entities(EntityRank::ELEMENT);
        }
        }
}
}

TEST_CASE("XFEM TOOLKIT DISCRETE LEVELSET FIELD","[DISCRETE_SPHERE]"){
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;

    // Geometry Engine Setup -----------------------
    // Using a Level Set Sphere as the Geometry
    real tRadius  = 2.4;
    real tXCenter = 2.5;
    real tYCenter = 2.5;
    real tZCenter = 2.5;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    std::string tLevelSetMeshFileName = "generated:5x5x5";
    xtk::Cell<std::string> tScalarFieldNames = {"LEVEL_SET_SPHERE"};
    xtk::Cell<xtk::Geometry<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer>*> tLevelSetFunctions = {&tLevelSetSphere};
    xtk::Discrete_Level_Set<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tLevelSetFunctions,tLevelSetMeshFileName,tScalarFieldNames,tMeshBuilder);
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetMesh,tPhaseTable);
    tGeometryEngine.mComputeDxDp = true;


    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGeometryEngine);

    tXTKModel.mSameMesh = true;

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                           Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    // Access the decomposed XTK Mesh
    Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

    // Do some testing
    size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
    size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

    CHECK(tNumNodesAfterDecompositionXTK == 2804);
    CHECK(tNumElementsAfterDecompositionXTK == 6576);

    /*
     * Get the output mesh and write to exodus file
     */

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Set the sensitivity field names
    tOutputOptions.mPackageDxDpSparsely = true;
    tOutputOptions.mDxDpName = "dxdp_"; // base of the dxdp data (appended with a norm)
    tOutputOptions.mDxDpIndicesName = "dxdp_inds_"; // base of the dxdp indices (appended with a number)
    tOutputOptions.mDxDpNumIndicesName = "dxdp_ninds"; // number of indices for a given node


    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);
    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/discrete_sphere.e";


    tCutMeshData->write_output_mesh(tMeshOutputFile);
}



TEST_CASE("Simple Side Set", "[SIMPLE_SIDE_SET]")
{

    /*
     * Initialize Matrix Factory
     */

    size_t tRadIts = 10;
    size_t tXIts = 1;
    size_t tYIts = 1;
    size_t tZIts = 1;

//    real tRadStart = 0.6633333333333333;
    real tRadStart = 0.33;
    real tRadEnd = 0.99;
    real tXStart = 0.0 ;
    real tXEnd = 1;
    real tYStart = 0.0;
    real tYEnd = 1;
    real tZStart = 0.0;
    real tZEnd = 1;

    //    Rad: 0.654 xc: 0.394 yc: -0.01 zc: 0.192

    // Increment size
    real tRadInc = (tRadEnd-tRadStart)/tRadIts;
    real tXInc = (tXEnd-tXStart)/tXIts;
    real tYInc = (tYEnd-tYStart)/tYIts;
    real tZInc = (tZEnd-tZStart)/tZIts;

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = false;

    tOutputOptions.change_phases_to_output(2,Cell<size_t>(1,1));;


    for(size_t iRad = 0; iRad< tRadIts; iRad++)
    {
        real tRadius  = tRadStart + tRadInc*iRad;

        for(size_t iX = 0; iX<tXIts; iX++)
        {
            real tXCenter = tXStart+tXInc*iX;

            for(size_t iY = 0; iY<tYIts; iY++)
            {
                real tYCenter = tYStart + tYInc*iY;
                for(size_t iZ = 0; iZ<tZIts; iZ++)
                {
                    real tZCenter = tZStart+tZInc*iZ;

                    std::cout<< "Rad: "<< tRadius << " xc: "<< tXCenter << " yc: "<< tYCenter<< " zc: "<< tZCenter<<std::endl;

                    /*
                     * Load Mesh which is a unit cube with 2 faces belonging to a side set
                     */
                    std::string tMeshFileName = "generated:1x1x1|sideset:xXyY";
                    mesh::Mesh_Builder_Stk<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
                    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName,{},true);

                    moris::Matrix< Default_Matrix_Integer > tElementFaces = tMeshData->get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT, EntityRank::FACE);
                    moris::Matrix< Default_Matrix_Integer > tElementEdges = tMeshData->get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT, EntityRank::EDGE);

                    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius,tXCenter,tYCenter,tZCenter);
                    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
                    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetSphere,tPhaseTable);

                    /*
                     * Setup XTK Model and tell it how to cut
                     */
                    size_t tModelDimension = 3;
                    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
                    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

                    /*
                     * Decompose
                     */
                    tXTKModel.decompose(tDecompositionMethods);


                    /*
                     * Get the output mesh and write to exodus file
                     */

                    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tOutputMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);
                    std::string tPrefix = std::getenv("XTKOUTPUT");
                    std::string tMeshOutputFile = tPrefix + "/unit_simple_side_set.exo";

                    tOutputMeshData->write_output_mesh(tMeshOutputFile);

                    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCutMesh =  tXTKModel.get_cut_mesh();
//                    tCutMesh.print_node_to_entity_connectivity_with_ancestry(0);
                }
            }
        }
    }

}

TEST_CASE("Propagate Mesh Sets","[SET_PROPOGATION]")
{
    /*
     * Loads an exodus file with a Block Set and Side Set already populated
     * Performs regular subdivison method and then checks to see if the
     * children XTK mesh entities have the same as their parents
     */

    /*
     * Set up:
     * Geometry,
     * Geometry Engine,
     * Mesh
     */

    real tRadius = 5.1;
    real tXCenter = 0.0;
    real tYCenter = 0.0;
    real tZCenter = 0.0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius,tXCenter,tYCenter,tZCenter);
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetSphere,tPhaseTable);

    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp = false;

    /*
     * Load Mesh which has 3 block sets. These blocks are named:
     *  - top_bread
     *  - meat
     *  - bottom_bread
     *
     * Side Sets will eventually be named
     *  - top_crust
     *  - bottom_crust
     */
    std::string tPrefix;
    tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";
    xtk::Cell<std::string> tFieldNames;
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName,tFieldNames,true);

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};

    std::clock_t start;
    start = std::clock();
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);

//    tXTKModel.convert_mesh_tet4_to_tet10();

    // Do the enrichment with a graph based method
//     Enrichment<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tEnrichment(2);
//     tEnrichment.perform_enrichment(tXTKModel.get_cut_mesh(), tXTKModel.get_xtk_mesh());


    /*
     * Get the output mesh and write to exodus file
     */

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Set the sensitivity field names
    tOutputOptions.mPackageDxDpSparsely = true;
    tOutputOptions.mDxDpName = "dxdp_"; // base of the dxdp data (appended with a norm)
    tOutputOptions.mDxDpIndicesName = "dxdp_inds_"; // base of the dxdp indices (appended with a number)
    tOutputOptions.mDxDpNumIndicesName = "dxdp_ninds"; // number of indices for a given node


    // Specify there are 2 possible phases
    size_t tNumPhases = 2;

    // Say I only want to output phase 1
    Cell<size_t> tPhasesToOutput = {0,1};
    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    // Add field for enrichment
    tOutputOptions.mIntElementExternalFieldNames = {"owner"};
    tOutputOptions.mInternalUseFlag = true;

    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tOutputMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);
    tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/unit_sandwich_sphere_05_threshold.e";


    // Add element owner as field
    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCutMesh = tXTKModel.get_cut_mesh();
    XTK_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tXTKMesh = tXTKModel.get_xtk_mesh();
    write_element_ownership_as_field(tOutputOptions.mIntElementExternalFieldNames(0),tXTKMesh,tCutMesh,*tOutputMeshData);

   tOutputMeshData->write_output_mesh(tMeshOutputFile,{},{},tOutputOptions.mIntElementExternalFieldNames,{},{});

}


TEST_CASE("Propagate Mesh Sets Discrete","[SET_PROPOGATION_DISCRETE]"){
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;

    // Geometry Engine Setup -----------------------
    // Using a Level Set Sphere as the Geometry
    real tRadius = 5.1;
    real tXCenter = 5.0;
    real tYCenter = 5.0;
    real tZCenter = 5.0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    std::string tLevelSetMeshFileName = "generated:10x10x10";
    xtk::Cell<std::string> tScalarFieldNames = {"LEVEL_SET_SPHERE"};
    xtk::Cell<xtk::Geometry<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer>*> tLevelSetFunctions = {&tLevelSetSphere};
    xtk::Discrete_Level_Set<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tLevelSetFunctions,tLevelSetMeshFileName,tScalarFieldNames,tMeshBuilder);
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetMesh,tPhaseTable);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGeometryEngine);

    tXTKModel.mSameMesh = true;

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                           Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    // Access the decomposed XTK Mesh
    Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();

    // Do some testing
    size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
    size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;
    tOutputOptions.change_phases_to_output(2,Cell<size_t>(1,1));

    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/unit_sandwich_sphere_05_threshold_discrete.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);


    //TODO: MORE TESTING.
}

}
