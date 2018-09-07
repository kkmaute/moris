// XTKL: Linear Algebra Includes
#include "catch.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"


#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Mesh.hpp"

#include "containers/cl_XTK_Cell.hpp"

#include "geometry/cl_Multi_Cylinder.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"

namespace xtk
{
    TEST_CASE("Finite element of a pressure vessel","[FE_CYLINDER]")
    {
      real tCordLength = 10;
      Cell<Cell<real>> tCenter = {{0,0,0}};
      Cell<real> tRadius = {tCordLength/4};
      Cell<real> tLength = {1.1*tCordLength};
      Cell<Cell<real>> tAxis   = {{1,0,0}};

      Multi_Cylinder<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMultiCylinder(tCenter,tRadius,tLength, tAxis);
      Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
      Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGEIn(tMultiCylinder,tPhaseTable);
      tGEIn.mComputeDxDp = false;

      // Load the mesh
      std::string tPrefix = std::getenv("XTKROOT");
      std::string tMeshFileName = tPrefix + "/TestExoFiles/tet_cube_mesh.e";
      xtk::Cell<std::string> tFieldNames;
      mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
      std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName,tFieldNames,true);

      // Setup XTK Model -----------------------------
      size_t tModelDimension = 3;
      Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGEIn);
      tXTKModel.mSameMesh = true;

      //Specify your decomposition methods and start cutting
      Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
      tXTKModel.decompose(tDecompositionMethods);

      Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tGEOut = tXTKModel.get_geom_engine();
      XTK_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tXTKMesh = tXTKModel.get_xtk_mesh();

      Output_Options<size_t> tOutputOptions;
      tOutputOptions.mAddNodeSets = true;
      tOutputOptions.mAddSideSets = true;
      Cell<size_t> tPhasesToOutput = {0};

      tOutputOptions.change_phases_to_output(2,tPhasesToOutput);

      std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

      tPrefix = std::getenv("XTKOUTPUT");
      std::string tMeshOutputFile = tPrefix + "/fem_cylinder.e";

      tCutMeshData->write_output_mesh(tMeshOutputFile);

    }
}
