/*
 * XTK_Tutorial_1.cpp
 *
 *  Created on: Oct 11, 2018
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
 * 3.) The mesh (via MTK)
 */
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"


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
int
main( int    argc,
      char * argv[] )
{

    // Initialize the communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);


    /*!
     * \section XFEM Model With Multiple Geometries
     * The XTK supports iterating of an arbitrary number of geometries allowing for many-material and
     * many-physic model generation. XTK is agnostic to physics and operates only on a collection of
     * geometries and a domain. The model generated has many block sets which are determined based on
     * a cell's relationship to the many geometries. Applying materials and physics to these block
     * sets is done external to XTK.
     * To generate an XFEM model from the multi-geometry problem, the following things are needed:
     *  - 1.) A discrete domain (a mesh),\f$\Omega\f$, in which to immerse a geometry vector
     *  - 2.) A geometry vector, the geometries to immerse in \f$\Omega\f$
     *  - 3.) A phase table in order to interpret a cell's relationship to the geometry vector
     *  - 4.) A geometry engine, to determine interface location in a cell
     *  - 5.) A decomposition method vector, to specify how to immerse the geometry vector in \f$\Omega\f$.
     *  - 6.) A XTK model, facilitates the decomposition process and becomes the resulting model.
     *
     *
     * The mesh file name (absolute path) needs to be specified.
     * XTKROOT is an environment variable specifying a directory where
     * some commonly used meshes are stored. It does not need to be used in general
     * but provides a convenient relative path to the desired exodus file.
     * \code{.cpp}
     * std::string tPrefix = std::getenv("XTKROOT");
     * std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";
     * \endcode
     *
     */
    std::string tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";


    /*!
     * Load the mesh into the Mesh_Data
     * \code{.cpp}
     * moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );
     * \endcode
     *
     * The mesh sandwich.e has multiple block sets, side sets and node sets.
     * @image html ./figures/sandwich_base_mesh.png "Sandwich Background Mesh"
     */
    moris::mtk::Mesh* tBackgroundMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );



    /*!
     * Setup the geometries and geometry engine. Each geometry is constructed
     * independently, then placed into a geometry vector and finally passed
     * to the geometry engine.
     *
     * For this example a two spheres is used. Before setting up the two sphere we specify the number of geometries that
     * will be used.
     * \code{.cpp}
     *  size_t tNumGeometries = 2;
     * \endcode
     */
    size_t tNumGeometries = 2;

    /*!
     * An isosurface of a sphere is defined as \n
     * \f$ (x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2 - r^2 = \phi  \f$
     *
     * We are interested in where an isosurface of the geometry intersects the provided mesh.
     * In general, the isosurface corresponding to \f$\phi = 0\f$ is used but this can be changed
     * in the geometry engine which is discussed later.
     *
     * The first geometry is user specified as a sphere with \f$r = 5.1\f$ and centered
     * at the origin.
     *
     * \code{.cpp}
     * real r1  = 5.1;
     * real xc1 = 0.0;
     * real yc1 = 0.0;
     * real zc1 = 0.0;
     * Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere1(r1,xc1,yc1,zc1);
     * \endcode
     *
     */
    real r1  = 5.1;
    real xc1 = 0.0;
    real yc1 = 0.0;
    real zc1 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere1(r1,xc1,yc1,zc1);

    /*!
     * Similarly, the second sphere is defined with \f$r = 3.1\f$ and centered at the
     * origin.
     *\code{.cpp}
     *     real r2  = 3.1;
     *     real xc2 = 0.0;
     *     real yc2 = 0.0;
     *     real zc2 = 0.0;
     *     Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere2(r2,xc2,yc2,zc2);
     *\endcode
     */
    real r2  = 3.1;
    real xc2 = 0.0;
    real yc2 = 0.0;
    real zc2 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere2(r2,xc2,yc2,zc2);

    /*!
     * The spheres are placed into a geometry vector, this works because the sphere
     * class inherits the geometry base class and therefore is type compatible with the geometry
     * class.
     *
     * \code{.cpp}
     *  xtk::Cell<Geometry<real,size_t,DDRMat,DDSTMat>*> tGeometryVector =
     *    {&tLevelSetSphere1, &tLevelSetSphere2};
     * \endcode
     */
    xtk::Cell<Geometry<real,size_t,DDRMat,DDSTMat>*> tGeometryVector =
           {&tLevelSetSphere1, &tLevelSetSphere2};

    /*!
     * The phase table is used to interpret if an element inside or outside with respect
     * to a given geometry. The collection of inside/outside information is then
     * condensed into a single index indicated the elemental phase index. A more detail ed
     * discussion can be found [here](doc/multi_phase.pdf).
     * TODO: LINK TO DISCUSSION ON PHASE TABLE
     *
     * For this example a \f$2^n\f$ material phase table.
     * \code{.cpp}
     *  Phase_Table<size_t, DDSTMat> tPhaseTable (tNumGeometries,  Phase_Table_Structure::EXP_BASE_2);
     * \endcode
     *
     * Setup the geometry engine with the geometry vector and phase table
     * \code{.cpp}
     * Geometry_Engine<real, size_t, DDRMat, DDSTMat> tGeometryEngine =
     *    Geometry_Engine<real, size_t, DDRMat, DDSTMat> (tGeometryVector,tPhaseTable);
     * \endcode
     */

    Phase_Table<size_t, DDSTMat> tPhaseTable (tNumGeometries,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, DDRMat, DDSTMat> tGeometryEngine
    = Geometry_Engine<real, size_t, DDRMat, DDSTMat> (tGeometryVector,tPhaseTable);

    /*!
     * Specify the isocontour, \f$\phi\f$ of the level set field.
     * and whether or not to compute the sensitivity field. The sensitivity field
     * is discussed in
     * TODO: ADD discussion of sensitivity link
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
     * \code{.cpp}
     *     xtk::Cell<enum Subdivision_Method> tDecompositionMethods
     *                = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
     *                   Subdivision_Method::C_HIERARCHY_TET4};
     * \endcode
     *
     * The Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 method
     * means the linear 8 node hexahedron background cell will be subdivided into 24 linear tetrahedrons
     * if it is interesected by a the isocontour of the allocated spheres. This regular subdivision
     * is shown in the following figure.
     * @image html ./figures/regular_subdivision.png "Regular subdivision template colored by level set field"
     *
     * A Subdivision_Method::C_HIERARCHY_TET4 method starts from a linear 4-node tetrahedron and generates
     * an interface along a isocontour. The result is a conformal mesh. For example, using the same example as shown
     * for the regular subdivision we receive the mesh in the following figure. The contour along the interface is
     * highlighted with a white line. More detail is provided in:
     *  @image html ./figures/conformal_subdivision.png "Conformal subdivision colored by level set field with contour at \f$\phi = 0\f$
     *
     */
    xtk::Cell<enum Subdivision_Method> tDecompositionMethods
            = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
               Subdivision_Method::C_HIERARCHY_TET4};

    /*!
     * Tell the model to decompose the mesh with the specified methods
     * \code{.cpp}
     * tXTKModel.decompose(tDecompositionMethods);
     * \endcode
     */
    tXTKModel.decompose(tDecompositionMethods);

    /*!
     * At this stage, the conformal model has been created. During the process, XTK
     * has produced an inheritance between elements in the decomposed model and the element
     * in the background mesh (input mesh), interface sensitivity for nodes on. This information
     * and more can be accessed via the Cut Mesh and XTK Meshes contained in the model. References
     * to these data structures can be accessed by
     * \code{.cpp}
     * Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();
     * XTK_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tXTKMesh = tXTKModel.get_xtk_mesh();
     * \endcode
     *
     * TODO: GO INTO MORE DETAIL ABOUT WORKING WITH THE CUT/XTK Mesh and the type of information in
     * TODO: each. This could be on its own page.
     *
     *
     */

    /*!
     * The model can be exported to MTK via STK via a call to get_output_mesh(). This function not only returns the output
     * mesh but it constructs it. In the future XTK may provide a standalone MTK interface. The call to XTK creates the MTK mesh
     * and the calls to create_output_mesh write the mesh to an exodus file.
     * \code{.cpp}
     *  moris::mtk::Mesh* tOutputMesh = tXTKModel.get_output_mesh();
     *
     *  std::string tMeshOutputFile = tPrefix + "/XTK_Tutorial_Multi.e";
     *  tOutputMesh->create_output_mesh(tMeshOutputFile);
     * \endcode
     *
     * @image html ./figures/output_mesh.png "Output mesh. Red line - isocontour surface of Sphere 1. Blue line - isocontour surface of Sphere 2$
     */
    moris::mtk::Mesh* tOutputMesh = tXTKModel.get_output_mesh();
    tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/XTK_Tutorial_Multi.e";
    tOutputMesh->create_output_mesh(tMeshOutputFile);

    delete tOutputMesh;
    delete tBackgroundMesh;

}


