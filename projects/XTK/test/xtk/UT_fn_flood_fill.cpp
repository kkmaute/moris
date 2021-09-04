/*
 * fn_flood_fill.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: ktdoble
 */

#include "xtk/cl_XTK_Child_Mesh.hpp"
#include "catch.hpp"


#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "cl_Mesh_Enums.hpp"

// XTKL: Geometry
#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_Sphere.hpp"
#include "geometry/cl_Mesh_Field_Geometry.hpp"


// XTKL: Linear Algebra Includes
#include "cl_XTK_Matrix_Base_Utilities.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace xtk
{

TEST_CASE("Generic Floodfill Consecutive Subdomain" ,"[GEN_FLOOD_FILL]")
        {
    // This test case colors a 6 element mesh with 2 possible phase indices (0,1)
    // (element id, phase index)
    // *-------*-------*-------*
    // |       |       |       |
    // |  0,1  |  1,0  |  2,1  |
    // |       |       |       |
    // *-------*-------*-------* *-------*
    // |       |       |       | |       |
    // |  3,1  |  4,0  |  5,0  | |  6,0  |
    // |       |       |       | |       |
    // *-------*-------*-------* *-------*
    //
    // Element 6 shows up in the element to element connectivity but I am not interested in considering it
    // This considers a part of the subdomain but the one outside of the subdomain is at the end. Next test will consider
    //
    // Number of phases
    size_t tNumPhases = 2;

    // Element to Element Connectivity
    size_t tMax = std::numeric_limits<size_t>::max();
    moris::Matrix< moris::DDSTMat >tElementToElement(7,3);
    tElementToElement.fill(tMax);
    (tElementToElement)(0,0) = 1;    (tElementToElement)(0,1) = 3;
    (tElementToElement)(1,0) = 0;    (tElementToElement)(1,1) = 2;    (tElementToElement)(1,2) = 4;
    (tElementToElement)(2,0) = 1;    (tElementToElement)(2,1) = 5;
    (tElementToElement)(3,0) = 0;    (tElementToElement)(3,1) = 4;
    (tElementToElement)(4,0) = 1;    (tElementToElement)(4,1) = 3;    (tElementToElement)(4,2) = 5;
    (tElementToElement)(5,0) = 2;    (tElementToElement)(5,1) = 4;    (tElementToElement)(5,1) = 6;
    (tElementToElement)(6,0) = 5;

    // Element Phase Indices
    moris::Matrix< moris::DDSTMat >tElementPhase(1,7);

    (tElementPhase)(0,0) = 1;
    (tElementPhase)(0,1) = 0;
    (tElementPhase)(0,2) = 1;
    (tElementPhase)(0,3) = 1;
    (tElementPhase)(0,4) = 0;
    (tElementPhase)(0,5) = 0;
    (tElementPhase)(0,6) = 0;


    SECTION("Include all elements")
    {
        // Active Element Indices (all of them in this case)
        moris::Matrix< moris::DDSTMat > tActiveElements(1,7);
        (tActiveElements)(0,0) = 0;
        (tActiveElements)(0,1) = 1;
        (tActiveElements)(0,2) = 2;
        (tActiveElements)(0,3) = 3;
        (tActiveElements)(0,4) = 4;
        (tActiveElements)(0,5) = 5;
        (tActiveElements)(0,6) = 6;

        moris::Matrix< moris::DDSTMat > tIncludedElementMarker(1,7);
        tIncludedElementMarker.fill(1); // Mark all elements as included


        // Run flood fill Algorithm

        moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                                                                       tElementPhase,
                                                                       tActiveElements,
                                                                       tIncludedElementMarker,
                                                                       tNumPhases,
                                                                       tMax,
                                                                       true);

        moris::Matrix< moris::DDSTMat >tExpElementSubphase(1,7);
        (tExpElementSubphase)(0,0) = 0;    (tExpElementSubphase)(0,1) = 1;    (tExpElementSubphase)(0,2) = 2;
        (tExpElementSubphase)(0,3) = 0;    (tExpElementSubphase)(0,4) = 1;    (tExpElementSubphase)(0,5) = 1;
        (tExpElementSubphase)(0,6) = 1;


        CHECK(equal_to(tElementSubphase,tExpElementSubphase));

    }
    SECTION("Excluding element 6")
    {
        // Active Element Indices (all of them in this case)
        moris::Matrix< moris::DDSTMat > tActiveElements(1,6);
        (tActiveElements)(0,0) = 0;
        (tActiveElements)(0,1) = 1;
        (tActiveElements)(0,2) = 2;
        (tActiveElements)(0,3) = 3;
        (tActiveElements)(0,4) = 4;
        (tActiveElements)(0,5) = 5;

        moris::Matrix< moris::DDSTMat > tIncludedElementMarker(1,7);
        tIncludedElementMarker.fill(1); // Mark all elements as included
        (tIncludedElementMarker)(0,6) = 0; // overwrite element 6 and say don't include it


        // Run flood fill Algorithm
        moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                                                                       tElementPhase,
                                                                       tActiveElements,
                                                                       tIncludedElementMarker,
                                                                       tNumPhases,
                                                                       tMax);


        moris::Matrix< moris::DDSTMat >tExpElementSubphase(1,6);
        (tExpElementSubphase)(0,0) = 0;    (tExpElementSubphase)(0,1) = 1;    (tExpElementSubphase)(0,2) = 2;
        (tExpElementSubphase)(0,3) = 0;    (tExpElementSubphase)(0,4) = 1;    (tExpElementSubphase)(0,5) = 1;

        CHECK(equal_to(tElementSubphase,tExpElementSubphase));

    }
        }




TEST_CASE("Generic Floodfill Nonconsectives Subdomain" ,"[NC_FLOOD_FILL]")
{
    // This test case colors a 6 element mesh with 2 possible phase indices (0,1)
    // (element id, phase index)
    // *-------*-------*-------*
    // |       |       |       |
    // |  0,1  |  6,0  |  2,1  |
    // |       |       |       |
    // *-------*-------*-------*-------*
    // |       |  (i)  |       |       |
    // |  3,1  |  4,1  |  5,1  |  1,0  |
    // |       |       |       |       |
    // *-------*-------*-------*-------*
    //
    // Ignoring element 4 (i)

    // Number of phases
    size_t tNumPhases = 2;

    // Element to Element Connectivity
    size_t tMax = std::numeric_limits<size_t>::max();
    moris::Matrix< moris::DDSTMat >tElementToElement(7,3);
    tElementToElement.fill(tMax);
    (tElementToElement)(0,0) = 3;    (tElementToElement)(0,1) = 6;
    (tElementToElement)(1,0) = 5;
    (tElementToElement)(2,0) = 6;    (tElementToElement)(2,1) = 5;
    (tElementToElement)(3,0) = 4;    (tElementToElement)(3,1) = 0;
    (tElementToElement)(4,0) = 3;    (tElementToElement)(4,1) = 5;    (tElementToElement)(4,2) = 6;
    (tElementToElement)(5,0) = 4;    (tElementToElement)(5,1) = 2;
    (tElementToElement)(6,0) = 0;    (tElementToElement)(6,1) = 2;    (tElementToElement)(6,2) = 4;

    // Element Phase Indices
    moris::Matrix< moris::DDSTMat >tElementPhase(1,7);

    (tElementPhase)(0,0) = 1;
    (tElementPhase)(0,1) = 0;
    (tElementPhase)(0,2) = 1;
    (tElementPhase)(0,3) = 1;
    (tElementPhase)(0,4) = 1;
    (tElementPhase)(0,5) = 1;
    (tElementPhase)(0,6) = 0;

    // Active Element Indices (all of them in this case)
    // Note in this case they are not sequential
    moris::Matrix< moris::DDSTMat > tActiveElements(1,6);
    (tActiveElements)(0,0) = 0;
    (tActiveElements)(0,1) = 6;
    (tActiveElements)(0,2) = 2;
    (tActiveElements)(0,3) = 3;
    (tActiveElements)(0,4) = 5;
    (tActiveElements)(0,5) = 1;

    moris::Matrix< moris::DDSTMat > tIncludedElementMarker(1,7);
    tIncludedElementMarker.fill(1); // Mark all elements as included
    (tIncludedElementMarker)(0,4) = 0; // overwrite element 4 and say don't include it

    // Run flood fill Algorithm
    moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                                                                   tElementPhase,
                                                                   tActiveElements,
                                                                   tIncludedElementMarker,
                                                                   tNumPhases,
                                                                   tMax);

    moris::Matrix< moris::DDSTMat > tExpElementSubphase(1,6);
    tExpElementSubphase(0,0) = 0;    tExpElementSubphase(0,1) = 1;    tExpElementSubphase(0,2) = 2;
    tExpElementSubphase(0,3) = 0;    tExpElementSubphase(0,4) = 2;    tExpElementSubphase(0,5) = 3;

    CHECK(equal_to(tElementSubphase,tExpElementSubphase));

}

TEST_CASE("Flood Fill on Child Mesh","[FLOOD_FILL]")
{
    if(par_size() == 1 || par_size() ==2)
    {
        // Specify a Sphere to Use
        real tRadius =  0.99;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 0;
        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
        Geometry_Engine tGeometryEngine(tLevelsetSphere);

        // Create Mesh --------------------------------------------------------------------
        std::string tMeshFileName = "generated:1x1x1";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);


        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        // Access the Cut Mesh-------------------------------------------------------------
        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();

        // At this point there is one child element in the cut mesh, and every element has its own phase value.
        // Now lets pass this child mesh to the enrichment flood fill
        Child_Mesh<real,size_t, moris::DDRMat,moris::DDSTMat> & tChildMesh = tCutMesh.get_child_mesh(0);

        // Perform floodfill
        size_t tChildMeshIndex = 0;
        moris::Matrix< moris::DDSTMat > tElementSubphase = local_child_mesh_flood_fill(tChildMesh);

        moris::Matrix< moris::DDSTMat > tExpElementSubphase({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});


        // Specify and enrichment level field on the output mesh
        Output_Options<size_t> tOutputOptions;
        tOutputOptions.mIntElementExternalFieldNames = {"enrich","phase_val"};
        tOutputOptions.mInternalUseFlag = true;

        // Tell the model
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tOutputMesh = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

        // Add element subphase to output
        //Add enrichment values to mesh
        size_t     tElementIndex    = 0;
        size_t     tElementPhase    = 0;
        size_t     tElementSubPhase = 0;
        size_t     tNumElements     = tChildMesh.get_num_entities(EntityRank::ELEMENT);
        Cell<real> tPhaseVal(tNumElements);
        Cell<real> tPhaseIndex(tNumElements);
        moris::Matrix< moris::DDSTMat > const & tElementId = tChildMesh.get_element_ids();
        for(size_t i = 0; i<tNumElements; i++)
        {
            tElementPhase = tChildMesh.get_element_phase_index(i)*10;
            tElementSubPhase = (tElementSubphase)(0,i);
            tElementIndex = tOutputMesh->get_loc_entity_index_from_entity_glb_id(tElementId(0,i), EntityRank::ELEMENT);
            tPhaseIndex(tElementIndex) = (real)(tElementPhase + tElementSubPhase);
            tPhaseVal(tElementIndex) = (real)tElementPhase;
        }

        tOutputMesh->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(0), EntityRank::ELEMENT, tPhaseIndex);
        tOutputMesh->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(1), EntityRank::ELEMENT, tPhaseVal);

        // Output the Mesh to Exodus File
        std::string tPrefix = std::getenv("XTKOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/flood_fill_cm.e";
        tOutputMesh->write_output_mesh(tMeshOutputFile,{},{},tOutputOptions.mIntElementExternalFieldNames,{},{});

    }
}



TEST_CASE("3 Subphase","[3_subphase_ff]")
{
    if(par_size() == 1 || par_size() ==2)
    {
        // Create Mesh --------------------------------------------------------------------
        std::string tMeshFileName = "generated:1x1x1";
        Cell<std::string> tScalarFields({"lsf"});
        mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

        // Add node level set values
        Cell<real> tLSVal({ 0.75,
            -1.00,
            -1.00,
            0.75,
            -1.00,
            0.75,
            0.75,
            -1.00});

        tMeshData->add_mesh_field_data_loc_indices(tScalarFields(0), EntityRank::NODE, tLSVal);

        // Output the Mesh to Exodus File
        std::string tPrefix = std::getenv("XTKOUTPUT");
        std::string tMeshInputFile = tPrefix + "/flood_fill_multisubphase_input.e";
        tMeshData->write_output_mesh(tMeshInputFile,tScalarFields,{},{},{},{});

        Mesh_Field_Geometry tDiscreteMesh(tMeshData,tScalarFields);

        Geometry_Engine tGeometryEngine(tDiscreteMesh);

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);


        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);



        // Access the Cut Mesh-------------------------------------------------------------
        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();

        // At this point there is one child element in the cut mesh, and every element has its own phase value.
        // Now lets pass this child mesh to the enrichment flood fill
        Child_Mesh<real,size_t, moris::DDRMat,moris::DDSTMat> & tChildMesh = tCutMesh.get_child_mesh(0);

        // Perform floodfill
        size_t tChildMeshIndex = 0;
        moris::Matrix< moris::DDSTMat > tElementSubphase = local_child_mesh_flood_fill(tChildMesh);

        moris::Matrix< moris::DDSTMat > tExpElementSubphase({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});


        // Specify and enrichment level field on the output mesh
        Output_Options<size_t> tOutputOptions;
        tOutputOptions.mRealNodeExternalFieldNames = {"lsf"};
        tOutputOptions.mIntElementExternalFieldNames = {"enrich","phase_val"};
        tOutputOptions.mInternalUseFlag = true;

        // Tell the model
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tOutputMesh = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

        // Add element subphase to output
        //Add enrichment values to mesh
        size_t     tNodeIndex    = 0;
        size_t     tElementPhase    = 0;
        size_t     tElementSubPhase = 0;
        size_t     tNumElements     = tChildMesh.get_num_entities(EntityRank::ELEMENT);
        Cell<real> tPhaseVal(tNumElements);
        Cell<real> tPhaseIndex(tNumElements);
        moris::Matrix< moris::DDSTMat > const & tElementId = tChildMesh.get_element_ids();
        for(size_t i = 0; i<tNumElements; i++)
        {
            tElementPhase = tChildMesh.get_element_phase_index(i)    ;
            tElementSubPhase = (tElementSubphase)(0,i);
            tNodeIndex = tOutputMesh->get_loc_entity_index_from_entity_glb_id(tElementId(0,i), EntityRank::ELEMENT);
            tPhaseIndex(tNodeIndex) = (real)(tElementPhase + tElementSubPhase);
            tPhaseVal(tNodeIndex) = (real)tElementPhase;
        }

        tOutputMesh->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(0), EntityRank::ELEMENT, tPhaseIndex);
        tOutputMesh->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(1), EntityRank::ELEMENT, tPhaseVal);

        size_t tNumNodes = tMeshData->get_num_entities(EntityRank::NODE);
        size_t tNodeId = 0;
        Cell<real> tNodeLSV(tNumNodes);
        for(size_t i = 0; i<tNumNodes; i++)
        {
            tNodeId = tMeshData->get_glb_entity_id_from_entity_loc_index(i, EntityRank::NODE);
            tNodeIndex = tOutputMesh->get_loc_entity_index_from_entity_glb_id(tNodeId,EntityRank::NODE);

            real tLSV = tXTKModel.get_geom_engine().get_entity_phase_val(i,0);

            size_t tNodeId = tMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);

            size_t tNodeIndex = tOutputMesh->get_loc_entity_index_from_entity_glb_id(tNodeId,EntityRank::NODE);

            tNodeLSV(tNodeIndex) = tLSV;
        }

        tOutputMesh->add_mesh_field_data_loc_indices(tOutputOptions.mRealNodeExternalFieldNames(0), EntityRank::NODE, tNodeLSV);



        // Output the Mesh to Exodus File
        std::string tMeshOutputFile = tPrefix + "/flood_fill_multisubphase_output.e";
        tOutputMesh->write_output_mesh(tMeshOutputFile,tOutputOptions.mRealNodeExternalFieldNames,{},tOutputOptions.mIntElementExternalFieldNames,{},{});
    }
}


}




