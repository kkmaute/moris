/*
 * cl_XTK_Sensitivity.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: ktdoble
 */


#include <memory>
#include <mpi.h>
#include "catch.hpp"
#include "ios/cl_Logger.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"



#include "xtk/cl_XTK_Sensitivity.hpp"
#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Mesh.hpp"

#include "containers/cl_XTK_Cell.hpp"

#include "geometry/cl_Discrete_Level_Set.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"





namespace xtk
{
  TEST_CASE("Sensitivity Data Structure", "[SENSITIVITY]")
        {
    size_t tNumDesignVars = 4;
    size_t tNumNodesWithSensitivity = 5;
    std::string tDxDpBase = "dxdp_";
    std::string tDxDpADVIndBase = "dxdp_ind";
    std::string tDxDpNADVIndBase = "dxdp_nind";

    xtk::Sensitivity<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tSensitivity(tNumDesignVars,
                                                                                          tNumNodesWithSensitivity,
                                                                                          tDxDpBase,
                                                                                          tDxDpADVIndBase,
                                                                                          tDxDpNADVIndBase);

    // Node 1 has sensitivity related to adv 1 and adv 3
    moris::Matrix<size_t,Default_Matrix_Integer> tADVIndicesN1({{1,3}});

    // Sensitivity Data related to node 1
    moris::Matrix<real, Default_Matrix_Real> tDxDpN1({{10.0,33.0,4.4},{3.0,1.1,2.4}});
    moris::Matrix<real, Default_Matrix_Real> tDxDpN2({{4.0,3.0,2.4},{1.1,33,40}});

    // Add the data
    tSensitivity.add_node_sensitivity(1,tADVIndicesN1,tDxDpN1);
    tSensitivity.add_node_sensitivity(1,tADVIndicesN1,tDxDpN1);

    tSensitivity.add_node_sensitivity(2,tADVIndicesN1,tDxDpN2);

    CHECK(tSensitivity.get_sensitivity_data(0).n_rows() == tNumNodesWithSensitivity);
    CHECK(tSensitivity.get_sensitivity_data(0).n_cols() == 3);

    CHECK(tSensitivity.get_sensitivity_data(1).n_rows() == tNumNodesWithSensitivity);
    CHECK(tSensitivity.get_sensitivity_data(1).n_cols() == 3);

    CHECK(tSensitivity.get_sensitivity_data(2).n_rows() == tNumNodesWithSensitivity);
    CHECK(tSensitivity.get_sensitivity_data(2).n_cols() == 3);

    tSensitivity.commit_sensitivities();

    CHECK(tSensitivity.get_sensitivity_data(0).n_rows() == 0);
    CHECK(tSensitivity.get_sensitivity_data(0).n_cols() == 3);

    CHECK(tSensitivity.get_sensitivity_data(1).n_rows() == 2);
    CHECK(tSensitivity.get_sensitivity_data(1).n_cols() == 3);

    CHECK(tSensitivity.get_sensitivity_data(3).n_rows() == 2);
    CHECK(tSensitivity.get_sensitivity_data(3).n_cols() == 3);


    moris::Matrix<size_t,Default_Matrix_Integer> tExpN1Map({{1, 3}, {0, 0}});
    moris::Matrix<size_t,Default_Matrix_Integer> tExpN2Map({{1, 3}, {1, 1}});

    print(tSensitivity.get_node_dxdp_map(2),"out");
    print(tExpN2Map,"exp");

    CHECK(equal_to(tSensitivity.get_node_dxdp_map(1),tExpN1Map));
    CHECK(equal_to(tSensitivity.get_node_dxdp_map(2),tExpN2Map));
        }

   TEST_CASE("Finite differencing check on a single tet mesh","[FD_TRI]")
   {
     // This is the interface node with dxdp information
     size_t tInterfaceNodeInd = 5;

     // background node level set value to perturb
     size_t tNodeToPerturbInd = 0;

     // Amount to perturb (for finite difference check)
     real   tPerturbVal = 1e-6;

     // Tolerance in the 2 norm
     real tTol = 1e-8;

     // Level Set Values
     Cell<Cell<real>> tLSF({{ 0.4, 0.3, 0.5, -0.5},  // base values
                             {0.4, 0.3+tPerturbVal, 0.5, -0.5},  // perturb up
                             {0.4, 0.3-tPerturbVal, 0.5, -0.5}});// perturb down

    // interface node location for each decomposition
    Cell<moris::Matrix<real,Default_Matrix_Real>> tInterfaceNodeLocation(2);

    // computed sensitivity
    moris::Matrix<real,Default_Matrix_Real> tDxDp(1,1);
    moris::Matrix<size_t,Default_Matrix_Integer> tDxDpInds(1,1);

    // Random threshold value between 0.3 and -0.5
    real tMax =  0.3;
    real tMin = -0.5;
    real tF = (real)rand() / RAND_MAX;
    real tThreshold = (real) (tMin + tF * (tMax - tMin));

    moris::Matrix<real,Default_Matrix_Real> tDxDpComp(3,1,10.0);
    for(size_t i = 0; i<3; i++)
    {
       std::string tPrefix;
       tPrefix = std::getenv("XTKROOT");
       std::string tMeshFileName = tPrefix + "/TestExoFiles/single_tet_mesh.e";
       Cell<std::string> tScalarFieldNames = {"lsf"};

       mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
       std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFieldNames, true);

       tMeshData->add_mesh_field_data_loc_indices(tScalarFieldNames(0),EntityRank::NODE,tLSF(i));


       xtk::Discrete_Level_Set<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tMeshData,tScalarFieldNames);
       Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
       Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGEIn(tLevelSetMesh,tPhaseTable);
       tGEIn.mComputeDxDp = true;
       tGEIn.mThresholdValue = tThreshold;
       // Setup XTK Model -----------------------------
       size_t tModelDimension = 3;
       Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGEIn);
       tXTKModel.mSameMesh = true;

       //Specify your decomposition methods and start cutting
       Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
       tXTKModel.decompose(tDecompositionMethods);

       Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tGEOut = tXTKModel.get_geom_engine();
       XTK_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tXTKMesh = tXTKModel.get_xtk_mesh();

       // store the xtk computed derivative
        if( i == 0)
        {
         tDxDp     = tGEOut.get_node_dx_dp(tInterfaceNodeInd);
         tDxDpInds = tGEOut.get_node_adv_indices(tInterfaceNodeInd);
        }
        else
        {
          tInterfaceNodeLocation(i-1) = tXTKMesh.get_mesh_data().get_selected_node_coordinates_loc_inds({{tInterfaceNodeInd}});
        }
     }

    Default_Matrix_Real tDxDpFdMat = (tInterfaceNodeLocation(1).matrix_data()-tInterfaceNodeLocation(0).matrix_data())/(2*tPerturbVal);
    moris::Matrix<real,Default_Matrix_Real> tDxDpFD(tDxDpFdMat);

    moris::Matrix<real,Default_Matrix_Real> tDxDpRow = tDxDp.get_row(0);
//    real t2Norm = (tDxDpFD.matrix_data()-tDxDpRow.matrix_data()).squaredNorm();
//    CHECK(t2Norm < tTol);

   }

   TEST_CASE("Finite differencing check on a single hex mesh","[FD_HEX]")
   {
     // This is the interface node with dxdp information
     size_t tInterfaceNodeInd = 15;

     // background node level set value to perturb
     size_t tNodeToPerturbInd = 0;

     // Amount to perturb (for finite difference check)
     real   tPerturbVal = 1e-6;

     // Tolerance in the 2 norm
     real tTol = 1e-8;

     // Level Set Values
     Cell<Cell<real>> tLSF({{ 0.4, 0.3, 0.5, -0.5,                0.2,  0.1,  0.3, 0.4},  // base values
                             {0.4, 0.3, 0.5, -0.5 + tPerturbVal,  0.2,  0.1, 0.3, 0.4},  // perturb up
                             {0.4, 0.3, 0.5, -0.5 - tPerturbVal,  0.2,  0.1, 0.3, 0.4  }});// perturb down

    // interface node location for each decomposition
    Cell<moris::Matrix<real,Default_Matrix_Real>> tInterfaceNodeLocation(2);

    // computed sensitivity
    moris::Matrix<real,Default_Matrix_Real> tDxDp(1,1);
    moris::Matrix<size_t,Default_Matrix_Integer> tDxDpInds(1,1);

    // Random threshold value between 0.3 and -0.5
    real tMax =  0.1;
    real tMin = -0.5;
    real tF = (real)rand() / RAND_MAX;
    real tThreshold = (real) (tMin + tF * (tMax - tMin));
    // real tThreshold = 0.0;

    moris::Matrix<real,Default_Matrix_Real> tDxDpComp(3,1,10.0);
    for(size_t i = 0; i<3; i++)
    {
       std::string tMeshFileName = "generated:1x1x1";
       Cell<std::string> tScalarFieldNames = {"lsf"};

       mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
       std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFieldNames, true);

       tMeshData->add_mesh_field_data_loc_indices(tScalarFieldNames(0),EntityRank::NODE,tLSF(i));


       xtk::Discrete_Level_Set<xtk::real, xtk::size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tMeshData,tScalarFieldNames);
       Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
       Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGEIn(tLevelSetMesh,tPhaseTable);
       tGEIn.mComputeDxDp = true;
       tGEIn.mThresholdValue = tThreshold;
       // Setup XTK Model -----------------------------
       size_t tModelDimension = 3;
       Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGEIn);
       tXTKModel.mSameMesh = true;

       //Specify your decomposition methods and start cutting
       Cell<enum Subdivision_Method> tDecompositionMethods =
          {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
           Subdivision_Method::C_HIERARCHY_TET4};

       tXTKModel.decompose(tDecompositionMethods);

       Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tGEOut = tXTKModel.get_geom_engine();
       XTK_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tXTKMesh = tXTKModel.get_xtk_mesh();

       // store the xtk computed derivative
        if( i == 0)
        {
         tDxDp     = tGEOut.get_node_dx_dp(tInterfaceNodeInd);
         tDxDpInds = tGEOut.get_node_adv_indices(tInterfaceNodeInd);
        }
        else
        {
          tInterfaceNodeLocation(i-1) = tXTKMesh.get_mesh_data().get_selected_node_coordinates_loc_inds({{tInterfaceNodeInd}});
        }
     }

    Default_Matrix_Real tDxDpFdMat = (tInterfaceNodeLocation(1).matrix_data()-tInterfaceNodeLocation(0).matrix_data())/(2*tPerturbVal);
    moris::Matrix<real,Default_Matrix_Real> tDxDpFD(tDxDpFdMat);

    moris::Matrix<real,Default_Matrix_Real> tDxDpRow = tDxDp.get_row(1);
//    real t2Norm = (tDxDpFD.matrix_data()-tDxDpRow.matrix_data()).squaredNorm();
//    CHECK(t2Norm < tTol);

   }


}
