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

#include "cl_Cell.hpp"


#include "banner.hpp" // COR/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_Mesh_Factory.hpp"

#include "SDF_Tools.hpp"

#include "cl_SDF_State.hpp"
#include "cl_SDF_Arguments.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Object.hpp"
#include "cl_SDF_Data.hpp"
#include "cl_SDF_Core.hpp"
#include "cl_SDF_STK.hpp"
#include "cl_SDF_Field.hpp"
moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace sdf;

//------------------------------------------------------------------------------
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );


    Object tObject( "/home/messe/bird3.d/stl.d/group_2.stl");

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
