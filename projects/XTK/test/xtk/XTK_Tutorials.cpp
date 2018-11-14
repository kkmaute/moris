/*
 * XTK_Tutorials.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: doble
 */

// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src

moris::Comm_Manager gMorisComm;


/*
 * XTK at a minimum needs:
 * 1.)An XTK model
 *
 */
#include "xtk/cl_XTK_Model.hpp"

/*
 * 4.)The Matrix class
 */
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

/*
 * 3.) The mesh
 */
#include "mesh/cl_Mesh_Data.hpp" // Mesh Data API
#include "mesh/cl_Mesh_Builder_Stk.hpp" // Method to build a mesh
#include "mesh/cl_Mesh_Enums.hpp" // convenient enums for mesh (i.e. entity ranks)

/*
 * 5.) A geometry and geometry engine
 */
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"



// Set namespaces to use
using namespace mesh;
using namespace moris;
using namespace xtk;

int
main( int    argc,
      char * argv[] )
{
    // Initialize the communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);


    /*!
     *  Specify the type of mesh building implementation to use, in this case,
     *  a STK mesh builder is used meaning a STK mesh will be the implementation
     *  for the background mesh
     * This allocates a 3x3 matrix but provides no values
     * \code{.cpp}
     * mesh::Mesh_Builder_Stk<real, size_t, DDRMat, DDSTMat> tMeshBuilder;
     * \endcode
     */
    mesh::Mesh_Builder_Stk<real, size_t, DDRMat, DDSTMat> tMeshBuilder;

    /*!
     * Next a file name (absolute path) needs to be specified)
     * XTKROOT is an environment variable specifying a directory where
     * some commonly used meshes are stored
     * \code{.cpp}
     * std::string tPrefix = std::getenv("XTKROOT");
     * std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";
     * \endcode
     */
    std::string tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";

    /*!
     * Specify the nodal fields to declare on mesh.
     * note: these do note need to be on the exodus file.
     *
     * \code{.cpp}
     * xtk::Cell<std::string> tNodalFieldNames = {"NODEFIELD1"} ;
     * \endcode
     */
    xtk::Cell<std::string> tNodalFieldNames = {"NODEFIELD1"} ;

    /*!
     * Specify whether or not to create faces and edges
     * Note: XTK needs faces and edges, so in general
     * this should be set to true
     * \code{.cpp}
     * bool tCreateFaces = true;
     * \endcode
     */
    bool tCreateFaces = true;

    /*!
     * Load the mesh into the Mesh_Data
     * \code{.cpp}
     *  std::shared_ptr<Mesh_Data<real, size_t, DDRMat, DDSTMat>> tBackgroundMesh
     *      = tMeshBuilder.build_mesh_from_string( tMeshFileName,
     *                                             tNodalFieldNames,
     *                                             tCreateFaces);
     * \endcode
     */
    std::shared_ptr<Mesh_Data<real, size_t, DDRMat, DDSTMat>> tBackgroundMesh
    = tMeshBuilder.build_mesh_from_string( tMeshFileName,
                                           tNodalFieldNames,
                                           tCreateFaces);


    /*!
     * Setup the geometries and geometry engine
     *
     * In this case, a problem with two geometries is set-up
     * \code{.cpp}
     *  size_t tNumGeometries = 2;
     * \endcode
     */
    size_t tNumGeometries = 2;

    /*!
     * The first sphere parameters
     *
     *
     * \code{.cpp}
     * real tRadius1  = 5.1;
     * real tXCenter1 = 0.0;
     * real tYCenter1 = 0.0;
     * real tZCenter1 = 0.0;
     * Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere1
     *     = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius1,tXCenter1,tYCenter1,tZCenter1);
     * \endcode
     *
     */
    real tRadius1  = 5.1;
    real tXCenter1 = 0.0;
    real tYCenter1 = 0.0;
    real tZCenter1 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere1
    = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius1,tXCenter1,tYCenter1,tZCenter1);

    /*!
     * The second sphere parameters
     *\code{.cpp}
     *     real tRadius2 = 3.1;
     *     real tXCenter2 = 0.0;
     *     real tYCenter2 = 0.0;
     *     real tZCenter2 = 0.0;
     *     Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere2
     *          = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius2,tXCenter2,tYCenter2,tZCenter2);
     *\endcode
     */
    real tRadius2 = 3.1;
    real tXCenter2 = 0.0;
    real tYCenter2 = 0.0;
    real tZCenter2 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere2
    = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius2,tXCenter2,tYCenter2,tZCenter2);

    /*!
     * Place pointers to geometries in a vector
     *
     * \code{.cpp}
     *  xtk::Cell<Geometry<real,size_t,DDRMat,DDSTMat>*> tGeometryVector = {&tLevelSetSphere1, &tLevelSetSphere2};
     * \end{code}
     */
    xtk::Cell<Geometry<real,size_t,DDRMat,DDSTMat>*> tGeometryVector = {&tLevelSetSphere1, &tLevelSetSphere2};

    /*!
     * The phase table is used to interpret an locations,
     * inside vs outside with respect to a given geometry
     * and interpret that as a material phase. This uses a
     * 2^n material phase table
     * \code{.cpp}
     *  Phase_Table<size_t, DDSTMat> tPhaseTable (tNumGeometries,  Phase_Table_Structure::EXP_BASE_2);
     * \end{code}
     *
     * Setup the geometry engine with the geomtry vector and phase table
     * \code{.cpp}
     * Geometry_Engine<real, size_t, DDRMat, DDSTMat> tGeometryEngine =
     *    \Geometry_Engine<real, size_t, DDRMat, DDSTMat> (tGeometryVector,tPhaseTable);
     * \endcode
     */

    Phase_Table<size_t, DDSTMat> tPhaseTable (tNumGeometries,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, DDRMat, DDSTMat> tGeometryEngine
    = Geometry_Engine<real, size_t, DDRMat, DDSTMat> (tGeometryVector,tPhaseTable);

    /*!
     * Specify the isocontour of the level set field.
     * and whether or not to compute the sensitivity field
     * \code{.cpp}
     *     tGeometryEngine.mThresholdValue = 0.0;
     *     tGeometryEngine.mComputeDxDp = true;
     * \endcode
     */
    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp = true;


    /*!
     * Tell the XTK model about:
     *  the model dimension,
     *  the background mesh
     *  the geometry engine
     *  \code{.cpp}
     *      size_t tModelDimension = 3;
     *      Model<real, size_t, DDRMat, DDSTMat> tXTKModel(tModelDimension, tBackgroundMesh, tGeometryEngine);
     *  \endcode
     */

    size_t tModelDimension = 3;
    Model<real, size_t, DDRMat, DDSTMat> tXTKModel(tModelDimension, tBackgroundMesh, tGeometryEngine);

    /*!
     * Specify the method to use for decomposing the mesh
     * Note: only the following two methods are implemented
     *
     * A Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 method
     * means the hex8 background xtk::Cells intersected by the geometries
     * will be regularly subdivided into 24 TET4s
     *
     * A Subdivision_Method::C_HIERARCHY_TET4 method
     * takes the regularly subdivided mesh from method 1
     * and generates a conformal interface mesh using a node
     * ids as a metric for subdivision
     *
     * \code{.cpp}
     *     xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
     *                                                            Subdivision_Method::C_HIERARCHY_TET4};
     * \endcode
     */
    xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                           Subdivision_Method::C_HIERARCHY_TET4};

    /*!
     * Tell the model to decompose the mesh with the specified methods
     * \code{.cpp}
     * tXTKModel.decompose(tDecompositionMethods);
     * \endcode
     */
    tXTKModel.decompose(tDecompositionMethods);



}

