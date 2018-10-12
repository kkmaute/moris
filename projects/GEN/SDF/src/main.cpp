/*
 * main.cpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

// MORIS header files.


#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Matrix.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_SDF_Triangle_Mesh.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Core.hpp"
#include "fn_print.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace sdf;

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // print welcome banner and system information
    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------

    sdf::Triangle_Mesh tObject( "Part_1.obj" );


    // load mesh from file
    mtk::Mesh * tInput = mtk::create_mesh( MeshType::STK, "TensorMesh.exo" , nullptr);

    // create SDF Wrapper
    Mesh tMesh( tInput );

    Data tData( tObject );

    Core tCore( tMesh, tData );

    Matrix< DDRMat > tSDF;

    tCore.calculate_raycast();


    tCore.save_to_vtk("sdf.vtk");

    //Core tCore ( )
    std::cout << tMesh.get_num_nodes() << std::endl;

    delete tInput;

    return 0;

}
