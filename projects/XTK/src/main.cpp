/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */
//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_HMR_Field.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
// XTK
#include "xtk/cl_XTK_Model.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "cl_Geom_Field.hpp"
#include "xtk_typedefs.hpp"
#include "geometry/cl_Geom_Field.hpp"

//------------------------------------------------------------------------------



// select namespaces
using namespace moris;
using namespace hmr;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
/*!
 * \section Tutorial_4: Using the MTK Interface
 *
 * In this example, we will refine a mesh according to the surface of
 * a sphere. When we are done, we will create an MTK interface object.
 *
 */
real
SphereFunction( const Matrix< DDRMat > & aPoint )
{
//    return norm( aPoint ) - 1.2;
    /*
     * Composite Fiber Parameters
     */
    real mFiberRadius = 0.12;    // radius of fibers
    real mFiberFrq    = 2.0;     // frequency of weave
    real mFiberExpo   = 2.0;     // cuboid exponent
    real mFiberDelX   =  1.0;
    real mFiberDelY   =  2.5;
    real mFiberDelZ   = 10.0;
    real mFiberXctr   =  0.00;
    real mFiberYctr   = -1.25;
    real mFiberZctr   = -5.00;
    size_t mNumSpheres  = 300;   // number of spheres
    real mFiberAmp  = mFiberDelX/2.0-1.1*mFiberRadius;
    real mFiberKmax = std::floor(mFiberDelY/(6.0*mFiberRadius))+1;
    real MATH_PI = 3.14159265359;
    real xcoord  = aPoint(0);
    real ycoord  = aPoint(1);
    real zcoord  = aPoint(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/mFiberRadius, mFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<mNumSpheres;i++)
    {
        // wavy fiber
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres);
            yci=mFiberYctr+k*6.0*mFiberRadius;
            zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        // straight fibers
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*(2.0*k-0.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*((2.0*k-1.0)-0.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        // wavy fiber
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres+2.0/3.0*MATH_PI);
            yci=mFiberYctr+2.0*mFiberRadius+k*6.0*mFiberRadius;
            zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }


        // straight fibers
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*(2.0*k-2.0/3.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));

            xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*((2.0*k-1.0)-2.0/3.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        // wavy fiber
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres+4.0/3.0*MATH_PI);
            yci=mFiberYctr+4.0*mFiberRadius+k*6.0*mFiberRadius;
            zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        // straight fibers
        for (uint k=0;k<mFiberKmax;k++)
        {
            xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*(2.0*k-4.0/3.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));

            xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
            yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
            zci=mFiberZctr+mFiberDelZ*((2.0*k-1)-4.0/3.0)/mFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }
    }


    return LSval;

}


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

    /*!
     * <b> Step 1: Create a Parameter list </b>
     */

    /*!
     * We create a parameterlist and build a mesh with 4x4x4 elements.
     * We want it to be 2x2x2 units long, and the origin at (0,0,0)
     * \code{.cpp}
     * ParameterList tParameters = create_parameter_list();
     * tParameters.set( "number_of_elements_per_dimension", "4, 4, 4" );
     * tParameters.set( "domain_dimensions",                "2, 2, 2" );
     * tParameters.set( "domain_offset",                    "0, 0, 0" );
     * \endcode
     */
    ParameterList tParameters = create_parameter_list();
    tParameters.set( "number_of_elements_per_dimension", "2, 2, 2" );
    tParameters.set( "domain_dimensions",                "2, 2, 2" );
    tParameters.set( "domain_offset",                    "0.0, 0.0, 0.0" );


    /*!
     * We want to refine the mesh four times. We set the surface refinement
     * level to 4
     * \code{.cpp}
     * tParameters.set( "max_surface_refinement_level", 4 );
     * \endcode
     */
    tParameters.set( "max_surface_refinement_level", 4 );


    /*!
     * Make HMR talkative
     * \code{.cpp}
     * tParameters.set( "verbose", 1 );
     * \endcode
     */
    tParameters.set( "verbose", 1 );
//------------------------------------------------------------------------------

    /*!
     * <b> Step 2: Create HMR object and refine it to the sphere </b>
     */

    /*!
     * Create a field with respect to the Parameterlist.
     * \code{.cpp}
     * HMR tHMR( tParameters );
     * \endcode
     */
    HMR tHMR( tParameters );

    /*!
     * The following two lines create an MTK mesh and a field
     * that is linked to this mesh.
     * \code{.cpp}
     * auto tMesh = tHMR.create_mesh();
     * auto tField = tMesh->create_field( "Sphere" );
     * \endcode
     */
    auto tMesh = tHMR.create_mesh();
    auto tField = tMesh->create_field( "Sphere" );


    /*!
     * No we evaluate the field, flag all elements on the surface
     * and perform the refinement.
     * The command update_refinement_pattern() is required at this
     * state. It will be removed in the future.
     * We repeat this process four times.
     *
     * \code{.cpp}
     * for( uint k=0; k<4; ++k )
     * {
     *   // evaluate field
     *   tField->evaluate_scalar_function( SphereFunction );
     *   tHMR.flag_surface_elements( tField );
     *   tHMR.perform_refinement();
     *   tHMR.update_refinement_pattern();
     * }
     * \endcode
     */

    for( uint k=0; k<3; ++k )
    {
        tField->evaluate_scalar_function( SphereFunction );
        tHMR.flag_surface_elements( tField );
        tHMR.perform_refinement();
        tHMR.update_refinement_pattern();
    }

    /*!
     * We call finalize in order to make the T-Matrices, Surfaces, Edges
     * and IDs available to MTK.
     *
     * \code{.cpp}
     * tHMR.finalize();
     * \endcode
     */
    tHMR.finalize();

    /**
     * Evaluate the sphere function one final time. After the mesh is finalized.
     * \code{.cpp}
     * tField->evaluate_scalar_function( SphereFunction );
     * \endcode
     */
    tField->evaluate_scalar_function( SphereFunction );

    /*!
     * Let HMR write the sphere as Exodus file.
     * \code{.cpp}
     * tHMR.save_to_exodus( "Sphere.exo" );
     * \endcode
     */
    tHMR.save_to_exodus( "Sphere.exo" );

//------------------------------------------------------------------------------
    // Using a the field as the geometry
    xtk::Geom_Field<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tFieldAsGeom(tField);
    xtk::Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    xtk::Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tFieldAsGeom,tPhaseTable);

    /**
     * Setup xtk model with HMR MTK mesh
     */
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMesh.get(),tGeometryEngine);

    //Specify your decomposition methods and start cutting
    xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                Subdivision_Method::C_HIERARCHY_TET4};

    // Decompose the mesh
    tXTKModel.decompose(tDecompositionMethods);

    // Get the XTK mesh as an MTK mesh
    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();


    /*!
      * Write the output mesh
      * \code{.cpp}
      * std::string tPrefix = std::getenv("XTKOUTPUT");
      * std::string tMeshOutputFile = tPrefix + "/hmr_to_xtk_intersected.e";
      * tCutMeshData->create_output_mesh(tMeshOutputFile);
      * \endcode
      */
    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/hmr_to_xtk_intersected.e";

    tCutMeshData->create_output_mesh(tMeshOutputFile);

    /*!
     * <b> Step 3:Using the MTK API</b>
     *
     * Now we can use the MTK API to ask for entity connectivities etc.
     *
     * some examples:
     *
     * \code{.cpp}
     * std::cout << "Number of Elements on the Mesh :"    << tMesh->get_num_elems()  << std::endl;
     * std::cout << "Number of Faces on the Mesh    :"    << tMesh->get_num_faces()  << std::endl;
     * std::cout << "Number of Edges on the Mesh    :"    << tMesh->get_num_edges()  << std::endl;
     * std::cout << "Number of Nodes on the Mesh    :"    << tMesh->get_num_nodes()  << std::endl;
     * std::cout << "Number of DOFs on the Mesh     :"     << tMesh->get_num_coeffs() << std::endl;
     * \endcode
     */

    std::cout << "Number of Elements on the Mesh :"    << tMesh->get_num_elems()  << std::endl;
    std::cout << "Number of Faces on the Mesh    :"    << tMesh->get_num_faces()  << std::endl;
    std::cout << "Number of Edges on the Mesh    :"    << tMesh->get_num_edges()  << std::endl;
    std::cout << "Number of Nodes on the Mesh    :"    << tMesh->get_num_nodes()  << std::endl;
    std::cout << "Number of DOFs on the Mesh     :"    << tMesh->get_num_coeffs() << std::endl;


//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
