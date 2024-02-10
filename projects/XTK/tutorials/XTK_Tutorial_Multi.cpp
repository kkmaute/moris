/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * XTK_Tutorial_Multi.cpp
 *
 */

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src
#include "paths.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

/*
 * XTK at a minimum needs:
 * 1.)An XTK model
 */
#include "xtk/cl_XTK_Model.hpp"

/*
 * 4.)The Matrix class
 */
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

/*
 * 3.) The mesh (via MTK)
 */
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

/*
 * 5.) A geometry and geometry engine
 */
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

// Set namespaces to use
using namespace mesh;
using namespace moris;
using namespace xtk;

/* @brief
 * This tutorial shows how to create an XTK model with
 * multiple immersed geometries.
 */
int main( int argc,
        char* argv[] )
{

    // Initialize the communication manager
    gMorisComm.initialize( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    /*!
     * \section XFEM Model With Multiple Geometries
     * The XTK supports iteration over an arbitrary number of geometries allowing for many-material and
     * many-physic model generation. XTK is agnostic to physics and operates on a vector of
     * geometries and a discrete domain. To generate an XFEM model from the multi-geometry problem,
     * the following things are needed:
     *  -# A \ref XTKBackgroundMesh "background mesh", a discrete domain or mesh,\f$\Omega\f$, in which to immerse the geometry vector.
     *  -# A \ref XTKGeometry "geometry vector", the geometries to immerse in \f$\Omega\f$
     *  -# A \ref XTKGeometry "phase table", to interpret a cell's relationship to the geometry vector
     *  -# A \ref XTKGeometry "geometry engine", to determine interface information in a cell
     *  -# A \ref XTKModel "XTK model", to facilitate the model generation process and become the resulting model.
     *
     *
     * \section  Problem Setup
     * @subsection background_mesh_setup 1.) Background Mesh Setup
     * The mesh file name (absolute path) needs to be specified.
     * XTKROOT is a convenient  environment variable specifying a directory where
     * some commonly used meshes are stored. It does not need to be used in general.
     * \code{.cpp}
     * std::string tPrefix = moris::get_base_moris_dir();
     * std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/sandwich.e";
     * \endcode
     *
     */
    std::string tPrefix       = moris::get_base_moris_dir();
    std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/sandwich.e";

    std::cout << "Mesh input name = " << tMeshFileName << std::endl;

    /*!
     * Load the mesh into the MTK library which provides an API to mesh functions. Allowing, for
     * many libraries to become the background mesh for XTK.
     * \code{.cpp}
     *  moris::mtk::Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName );
     * \endcode
     *
     * In this example a mesh file called sandwich.e which has multiple block sets, side sets and node sets is used.
     * @image html ./figures/sandwich_base_mesh.png "Sandwich Background Mesh"
     */
    std::shared_ptr< moris::mtk::Interpolation_Mesh > tBackgroundMesh = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName );

    /*!
     * \subsection geometry_vector_sectup 2.) Geometry Vector Setup
     * Each geometry is constructed independently, then placed into a geometry vector which is then passed
     * to the geometry engine.
     *
     * First specify the number of geometries which will be used.
     * \code{.cpp}
     *  size_t tNumGeometries = 2;
     * \endcode
     */
    size_t tNumGeometries = 2;

    /*!
     * In order for a geometry to be compatible with XTK, it needs to inherit from the Geometry base class. This allows
     * for use of a wide range of different geometries, i.e. a discrete field on the mesh or a geometric primitive.
     * The geometry only needs to provide a "distance" to an isocontour of interest.
     *
     * In this example, two spheres are going to be used. The level set field of a sphere is defined as \n
     * \f$ (x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2 - r^2 = \phi  \f$
     *
     * We are interested in where an isosurface,\f$ \phi = \phi_0 \f$, of the geometry intersects the provided mesh.
     * In general, \f$\phi_0 = 0\f$ is used but this can be changed in the geometry engine which is discussed later.
     *
     * The first geometry is user specified as a sphere with \f$r = 5.1\f$ and centered
     * at the origin.
     *
     * \code{.cpp}
     * real r1  = 5.3;
     * real xc1 = 0.0;
     * real yc1 = 0.0;
     * real zc1 = 0.0;
     * Sphere tLevelSetSphere1(r1,xc1,yc1,zc1);
     * \endcode
     *
     */
    real   r1  = 5.3;
    real   xc1 = 0.0;
    real   yc1 = 0.0;
    real   zc1 = 0.0;
    Sphere tLevelSetSphere1( r1, xc1, yc1, zc1 );

    /*!
     * Similarly, the second sphere is defined with \f$r = 3.1\f$ and centered at the
     * origin.
     *\code{.cpp}
     *     real r2  = 3.1;
     *     real xc2 = 0.0;
     *     real yc2 = 0.0;
     *     real zc2 = 0.0;
     *     Sphere tLevelSetSphere2(r2,xc2,yc2,zc2);
     *\endcode
     */
    real   r2  = 3.1;
    real   xc2 = 0.0;
    real   yc2 = 0.0;
    real   zc2 = 0.0;
    Sphere tLevelSetSphere2( r2, xc2, yc2, zc2 );

    /*!
     * The spheres are placed into a geometry vector, this works because the sphere
     * class inherits the geometry base class and therefore is type compatible with the geometry
     * class.
     *
     * \code{.cpp}
     *  Vector<Geometry*> tGeometryVector =
     *    {&tLevelSetSphere1, &tLevelSetSphere2};
     * \endcode
     */
    Vector< Geometry* > tGeometryVector = { &tLevelSetSphere1, &tLevelSetSphere2 };

    /*!
     * \subsection geom_eng_setup  4.) Geometry Engine Setup
     *
     * The geometry engine receives the geometry vector and phase table which have already been created.
     *
     * \code{.cpp}
     *  Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);
     * \endcode
     */
    Geometry_Engine tGeometryEngine( tGeometryVector );

    /*!
     * Specify the isocontour, \f$\phi_0\f$ of the level set field.
     * and whether or not to compute the sensitivity field.
     * \code{.cpp}
     *     tGeometryEngine.mThresholdValue = 0.0;
     *     tGeometryEngine.mComputeDxDp = true;
     * \endcode
     */
    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp    = true;

    /*!
     * \subsection xtk_setup 5.) XTK Model Setup
     *  The XTK model receives the model dimension, the background mesh and the geometry engine.
     *  The XTK model verbose can be set to output timing data, if the verbose flag is set to false then there will be
     *  no console outputs from the XTK model directly. This completes initialization of
     *  the XTK model.  Now we are ready to use some of the core capabilities of the library.
     *  \code{.cpp}
     *      size_t tModelDimension = 3;
     *      Model tXTKModel(tModelDimension, tBackgroundMesh, tGeometryEngine);
     *      tXTKModel.mVerbose  =  false;
     *  \endcode
     */
    size_t tModelDimension = 3;
    Model  tXTKModel( tModelDimension, tBackgroundMesh.get(), tGeometryEngine );
    tXTKModel.mVerbose = false;

    /*!
     * \section XTKCapabilities Core Capabilities
     * XTK offers a few core capabilities which are demonstrated in this tutorial. These core capabilities are:
     *  -# \ref XTKDecomp "Geometry aware domain decomposition"
     *  -# \ref XTKSensitivity "Interface design sensitivity computation"
     *  -# \ref XTKUnzip  "Interface unzipping"
     *  -# \ref XTKEnrich "Basis enrichment"
     *  -# \ref XTKOutput "Exporting the Model"
     *
     * \subsection decomp  1.) Geometry Aware Domain Decomposition
     *
     * Specify the method to use for decomposing the mesh. These decomposition methods will be iterated through.
     * Note: only the following two methods are implemented
     *
     * \code{.cpp}
     *     xtk::Cell<enum Subdivision_Method> tDecompositionMethods
     *                = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
     *                   Subdivision_Method::C_HIERARCHY_TET4};
     * \endcode
     *
     * The Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 method
     * takes linear 8 node hexahedron background cell subdivides the cell into 24 linear tetrahedrons
     * if the cell is intersected by \f$\phi_0\f$ of a geometry. This method is shown in the following figure.
     * @image html ./figures/regular_subdivision.png "Regular subdivision template colored by level set field"
     *
     * A Subdivision_Method::C_HIERARCHY_TET4 method starts from a linear 4-node tetrahedron and generates
     * an interface along a \f$\phi_0\f$. The result is a mesh which approximates the isosurface at \f$\phi_0\f$.
     * Using the mesh output from the regular subdivision as the input for this method,
     * we receive the following mesh where contour along \f$\phi_0\f$ isocontour is
     * highlighted with a white line.
     *  @image html ./figures/conformal_subdivision.png "Conformal subdivision colored by level set field with contour at threshold$"
     *
     */
    Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
        Subdivision_Method::C_HIERARCHY_TET4 };

    /*!
     * With the decomposition methods specified, calling decompose through the XTK model generates the geometry aware domain decomposition.
     * \code{.cpp}
     * tXTKModel.decompose(tDecompositionMethods);
     * \endcode
     */
    tXTKModel.decompose( tDecompositionMethods );

    /*!
     * Beyond the geometry aware discretization created, XTK has produced and inheritance of all entities to the background mesh and all child
     * element are aware in which part of the domain they belong. To access the discretized mesh, the user is referred to the Cut_Mesh and
     * Background_Mesh classes.
     * \code{.cpp}
     * Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();
     * Background_Mesh const & tXTKMesh = tXTKModel.get_background_mesh();
     * \endcode
     *
     * Other tutorials will go into more detail about using these data structures.
     */

    /*!
     * \subsection unzip 4.) Interface unzipping
     */

    //      tXTKModel.unzip_interface();

    /*!
     * \subsection enrichment 4.) Basis Enrichment
     * For XFEM, a basis enrichment strategy is needed. Enrichment is performed by XTK
     * through the call to perform_basis_enrichment.
     * \code{.cpp}
     *  tXTKModel.perform_basis_enrichment();
     * \endcode
     *
     */
    //      tXTKModel.perform_basis_enrichment();

    /*!
     *
     * \subsection export 5.) Exporting the Model
     * The model can be exported to an MTK mesh via a call to get_output_mesh(). This function call
     * creates the output database and returns a pointer to it. In the future XTK may provide a standalone MTK interface.
     * After the call to get_output_mesh(), the mtk mesh can be exported to an exodus file with the create_output_mesh() call.
     * \code{.cpp}
     *  moris::mtk::Mesh* tOutputMesh = tXTKModel.get_output_mesh();
     *
     *  std::string tMeshOutputFile = tPrefix + "/XTK_Tutorial_Multi.e";
     *  tOutputMesh->create_output_mesh(tMeshOutputFile);
     * \endcode
     *
     * @image html ./figures/output_mesh.png "Output mesh. Red line - isocontour surface of Sphere 1. Blue line - isocontour surface of Sphere 2$
     */
    moris::mtk::Mesh* tOutputMesh     = tXTKModel.get_output_mesh();
    std::string       tMeshOutputFile = "./xtk_exo/XTK_Tutorial_Multi.e";
    tOutputMesh->create_output_mesh( tMeshOutputFile );

    delete tOutputMesh;

    gMorisComm.finalize();
}
