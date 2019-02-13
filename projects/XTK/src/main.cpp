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
#include "cl_Cell.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// XTK
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_Sphere.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src
//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
