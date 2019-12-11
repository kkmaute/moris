/*
 * UT_FRF_Matching.hpp
 *
 *  Created on: Aug 9, 2019
 *      Author: doble
 */

#include "catch.hpp"
#include "cl_Multi_Geometry_KS_AD.hpp"
#include "cl_SphereBox_AD.hpp"
#include "cl_XTK_Model.hpp"
#include <Sacado.hpp>

using namespace moris;
namespace xtk
{

TEST_CASE("Box problem","[Box_FRF]")
        {
    if(par_size() == 1)
    {
        // initial parameters
        Sacado::Fad::DFad<moris::real> Width  = 6.00;
        Sacado::Fad::DFad<moris::real> Depth  = 6.00;
        Sacado::Fad::DFad<moris::real> Height = 6.00;
        Sacado::Fad::DFad<moris::real> H_thk  = 1.10;
        Sacado::Fad::DFad<moris::real> V_thk  = 1.10;
        Sacado::Fad::DFad<moris::real> slot_w = 0.75;

        // declare as independent variables
        Width.diff(0,5);
        Depth.diff(1,5);
        Height.diff(2,5);
        H_thk.diff(3,6);
        V_thk.diff(4,6);
        slot_w.diff(5,6);


        // bottom  back left corner off set
        moris::real tXOff = 2.04;
        moris::real tYOff = 2.05;
        moris::real tZOff = 4.02;

        // box sphere exponent
        moris::real tNExp = 10;

        // beta of KS
        moris::real tBeta = 4;

        // Bottom Plate
        Sacado::Fad::DFad<moris::real> tSx = Depth/2;
        Sacado::Fad::DFad<moris::real> tSy = Width/2;
        Sacado::Fad::DFad<moris::real> tSz = V_thk/2;
        Sacado::Fad::DFad<moris::real> tCx = Depth/2 + tXOff;
        Sacado::Fad::DFad<moris::real> tCy = Width/2 + tYOff;
        Sacado::Fad::DFad<moris::real> tCz = V_thk/2 + tZOff;
        Sphere_Box_AD tBottomPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // right plate (pos y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = Width - H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        Sphere_Box_AD tRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // Left plate (neg y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        Sphere_Box_AD tLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top left plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4 ;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        Sphere_Box_AD tTopLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top right plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width + slot_w)/2 +  (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        Sphere_Box_AD tTopRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);


        // place in geometry
        moris::Cell<Geometry_AD*> tGeomVect = {&tBottomPlate,&tRightPlate,&tLeftPlate,&tTopLeftPlate, &tTopRightPlate};
        //    moris::Cell<Geometry*> tGeomVect = {&tBottomPlate,&tRightPlate};
        Multi_Geometry_KS_AD tSubAssembly(tGeomVect,tBeta,6);


        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tSubAssembly,tPhaseTable);

        tGeometryEngine.mThresholdValue = 0.0;

//        std::string tPrefix = std::getenv("MORISROOT");
//        std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/frf_background.exo";
        std::string tMeshFileName = "../../../frf_background.exo";

        // Declare scalar node field
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName1 = "gd_0";
        tNodeField1.set_field_name(tFieldName1);
        tNodeField1.set_field_entity_rank(EntityRank::NODE);

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Place the node field into the field info container
        add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

        // Declare some supplementary fields
        mtk::MtkMeshData tInputMeshData;
        tInputMeshData.FieldsInfo = &tFieldsInfo;

        // fill in the parallel fields
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tInputMeshData  );


        Matrix<DDRMat> tBGGeomField(1,tMeshData->get_num_entities(EntityRank::NODE));
        for(moris::uint i = 0; i < tBGGeomField.numel(); i++)
        {
            Matrix<DDRMat> tNodeCoord = tMeshData->get_node_coordinate((moris_index)i);

            tBGGeomField(i) = tSubAssembly.evaluate_field_value_with_coordinate(0,tNodeCoord);
        }

        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tFieldName1,EntityRank::NODE,tBGGeomField);
        std::string tBGMeshOutputFile = "./xtk_exo/frf_box_bg.e";
        tMeshData->create_output_mesh(tBGMeshOutputFile);


        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};

        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose  =  false;
        tXTKModel.decompose(tDecompositionMethods);
        tXTKModel.compute_sensitivity();


        Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = true;
        tOutputOptions.mAddSideSets = true;

        // Add field for enrichment
        tOutputOptions.mInternalUseFlag = false;
        tOutputOptions.mPackageDxDpSparsely = false;
        tOutputOptions.mPackageDxDpDensely = false;



        moris::mtk::Mesh* tOutputMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tMeshOutputFile = "./xtk_exo/frf_box.e";
        tOutputMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tOutputMeshData;
    }
        }

TEST_CASE("Box problem Hex","[Box_FRF_HEX]")
        {
    if(par_size() == 1)
    {
        // initial parameters
        Sacado::Fad::DFad<moris::real> Width  = 6.00;
        Sacado::Fad::DFad<moris::real> Depth  = 6.00;
        Sacado::Fad::DFad<moris::real> Height = 6.00;
        Sacado::Fad::DFad<moris::real> H_thk  = 1.10;
        Sacado::Fad::DFad<moris::real> V_thk  = 1.10;
        Sacado::Fad::DFad<moris::real> slot_w = 0.75;

        // declare as independent variables
        Width.diff(0,5);
        Depth.diff(1,5);
        Height.diff(2,5);
        H_thk.diff(3,6);
        V_thk.diff(4,6);
        slot_w.diff(5,6);


        // bottom  back left corner off set
        moris::real tXOff = 2.04;
        moris::real tYOff = 2.05;
        moris::real tZOff = 4.02;

        // box sphere exponent
        moris::real tNExp = 100;

        // beta of KS
        moris::real tBeta = 4;

        // Bottom Plate
        Sacado::Fad::DFad<moris::real> tSx = Depth/2;
        Sacado::Fad::DFad<moris::real> tSy = Width/2;
        Sacado::Fad::DFad<moris::real> tSz = V_thk/2;
        Sacado::Fad::DFad<moris::real> tCx = Depth/2 + tXOff;
        Sacado::Fad::DFad<moris::real> tCy = Width/2 + tYOff;
        Sacado::Fad::DFad<moris::real> tCz = V_thk/2 + tZOff;
        Sphere_Box_AD tBottomPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // right plate (pos y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = Width - H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        Sphere_Box_AD tRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // Left plate (neg y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        Sphere_Box_AD tLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top left plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4 ;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        Sphere_Box_AD tTopLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top right plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width + slot_w)/2 +  (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        Sphere_Box_AD tTopRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);


//
//        std::cout<<"\nBottom Plate"<<std::endl;
//        tBottomPlate.print();
//        std::cout<<"\nRight Plate"<<std::endl;
//        tRightPlate.print();
//        std::cout<<"\nLeft Plate"<<std::endl;
//        tLeftPlate.print();
//        std::cout<<"\nTop Left Plate"<<std::endl;
//        tTopLeftPlate.print();
//        std::cout<<"\nTop Right Plate"<<std::endl;
//        tTopRightPlate.print();


        // place in geometry
        moris::Cell<Geometry_AD*> tGeomVect = {&tBottomPlate,&tRightPlate,&tLeftPlate,&tTopLeftPlate, &tTopRightPlate};
        //    moris::Cell<Geometry*> tGeomVect = {&tBottomPlate,&tRightPlate};
        Multi_Geometry_KS_AD tSubAssembly(tGeomVect,tBeta,6);


        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tSubAssembly,tPhaseTable);

        tGeometryEngine.mThresholdValue = 0.0;

//        std::string tPrefix = std::getenv("MORISROOT");
//        std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/frf_background.exo";
        std::string tMeshFileName = "generated:40x40x40|bbox:0,0,0,10,10,10";

        // Declare scalar node field
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName1 = "gd_0";
        tNodeField1.set_field_name(tFieldName1);
        tNodeField1.set_field_entity_rank(EntityRank::NODE);

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Place the node field into the field info container
        add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

        // Declare some supplementary fields
        mtk::MtkMeshData tInputMeshData;
        tInputMeshData.FieldsInfo = &tFieldsInfo;

        // fill in the parallel fields
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tInputMeshData  );


        Matrix<DDRMat> tBGGeomField(1,tMeshData->get_num_entities(EntityRank::NODE));
        for(moris::uint i = 0; i < tBGGeomField.numel(); i++)
        {
            Matrix<DDRMat> tNodeCoord = tMeshData->get_node_coordinate((moris_index)i);

            tBGGeomField(i) = tSubAssembly.evaluate_field_value_with_coordinate(0,tNodeCoord);
        }

        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tFieldName1,EntityRank::NODE,tBGGeomField);
        std::string tBGMeshOutputFile = "./xtk_exo/frf_box_hex_bg.e";
        tMeshData->create_output_mesh(tBGMeshOutputFile);


        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};

        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose  =  false;
        tXTKModel.decompose(tDecompositionMethods);
        tXTKModel.compute_sensitivity();


        Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = true;
        tOutputOptions.mAddSideSets = true;

        // Add field for enrichment
        tOutputOptions.mInternalUseFlag = false;
        tOutputOptions.mPackageDxDpSparsely = false;
        tOutputOptions.mPackageDxDpDensely = false;



        moris::mtk::Mesh* tOutputMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tMeshOutputFile = "./xtk_exo/frf_box_hex.e";
        tOutputMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tOutputMeshData;
    }
        }
}


