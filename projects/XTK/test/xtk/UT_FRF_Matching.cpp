/*
 * UT_FRF_Matching.hpp
 *
 *  Created on: Aug 9, 2019
 *      Author: doble
 */

#include "catch.hpp"
//#include "cl_Multi_Geometry.hpp"
//#include "cl_SphereBox.hpp"
#include "cl_XTK_Model.hpp"

#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/ripped/geometry/cl_GEN_Multi_Geometry.hpp"
#include "../projects/GEN/src/ripped/geometry/cl_GEN_Sphere_Box.hpp"

using namespace moris;
namespace xtk
{

TEST_CASE("Box problem","[Box_FRF]")
        {
    if(par_size() == 1)
    {
        // initial parameters
        moris::real Width  = 6.00;
        moris::real Depth  = 6.00;
        moris::real Height = 6.00;
        moris::real H_thk  = 0.50;
        moris::real V_thk  = 0.50;
        moris::real slot_w = 0.75;

        // bottom  back left corner off set
        moris::real tXOff = 2.0;
        moris::real tYOff = 2.0;
        moris::real tZOff = 4.02;

        // box sphere exponent
        moris::real tNExp = 10;

        // Bottom Plate
        moris::real tSx = Depth/2;
        moris::real tSy = Width/2;
        moris::real tSz = V_thk/2;
        moris::real tCx = Depth/2 + tXOff;
        moris::real tCy = Width/2 + tYOff;
        moris::real tCz = V_thk/2 + tZOff;
        moris::ge::Sphere_Box tBottomPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // right plate (pos y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = Width - H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        moris::ge::Sphere_Box tRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // Left plate (neg y face)
        tSx = Depth/2;
        tSy = H_thk/2;
        tSz = Height/2-V_thk*0.75;
        tCx = Depth/2 + tXOff;
        tCy = H_thk/2 + tYOff;
        tCz = Height/2 + tZOff;
        moris::ge::Sphere_Box tLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top left plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4 ;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        moris::ge::Sphere_Box tTopLeftPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);

        // top right plate
        tSx = Depth/2;
        tSy = (Width - slot_w)/4;
        tSz = V_thk/2;
        tCx = Depth/2 + tXOff;
        tCy = (Width + slot_w)/2 +  (Width - slot_w)/4 + tYOff;
        tCz = (Height) - V_thk/2 + tZOff;
        moris::ge::Sphere_Box tTopRightPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);


        // place in geometry
        moris::Cell<moris::ge::GEN_Geometry*> tGeomVect = {&tBottomPlate,&tRightPlate,&tLeftPlate,&tTopLeftPlate, &tTopRightPlate};
        //    moris::Cell<Geometry*> tGeomVect = {&tBottomPlate,&tRightPlate};
        moris::ge::Multi_Geometry tSubAssembly(tGeomVect);


        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tSubAssembly,tPhaseTable);

        tGeometryEngine.mThresholdValue = 0.0;
        tGeometryEngine.mComputeDxDp = false;

        std::string tPrefix = std::getenv("MORISROOT");
        std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/frf_background.exo";


        // fill in the parallel fields
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, nullptr  );


        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};

        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;
        tXTKModel.decompose(tDecompositionMethods);


        Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = true;
        tOutputOptions.mAddSideSets = true;

        // Add field for enrichment
        tOutputOptions.mInternalUseFlag = false;


        moris::mtk::Mesh* tOutputMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tMeshOutputFile = "./xtk_exo/frf_box.e";
        tOutputMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tOutputMeshData;
    }
        }
}


