/*
 * fn_generate_sphere_model.hpp
 *
 *  Created on: May 21, 2018
 *      Author: ktdoble
 */

#ifndef UNIT_TEST_SRC_XTK_FN_GENERATE_SPHERE_MODEL_HPP_
#define UNIT_TEST_SRC_XTK_FN_GENERATE_SPHERE_MODEL_HPP_
#include <memory>
#include <mpi.h>
#include "catch.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"


// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"

#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"

#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

namespace xtk
{
/*
 * Performs XFEM model generation with a spherical inclusion, and a mesh generated from string
 * @param[in] aRadius - sphere's radius
 * @param[in] aXCenter - sphere's x center
 * @param[in] aYCenter - sphere's y center
 * @param[in] aZCenter - sphere's z center
 * @param[in] aNumEleX - Number of elements in x direction
 * @param[in] aNumEleY - Number of elements in y direction
 * @param[in] aNumEleZ - Number of elements in z direction
 * @param[in] aHasSides - indicates the mesh is going to have side sets
 * @param[in] aLXSIde - Add a side set on left x side
 * @param[in] aLXSIde - Add a side set on right x side
 * @param[in] aBYSide - Add a side set on back y side
 * @param[in] aFYSide - Add a side set on front y side
 * @param[in] aBZSide - Add a side set on bottom z side
 * @param[in] aTZSide - Add a side set on top z side
 */
Model<real,size_t,Default_Matrix_Real,Default_Matrix_Integer>
generate_model_from_geometry_and_string(Geometry<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>  & aGeometry,
                                        std::string & aMeshFileName)
{
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(aGeometry,tPhaseTable);

    // Create Mesh ---------------------------------
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string(aMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    return tXTKModel;
}
}

#endif /* UNIT_TEST_SRC_XTK_FN_GENERATE_SPHERE_MODEL_HPP_ */
