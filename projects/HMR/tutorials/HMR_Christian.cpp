/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Christian.cpp
 *
 */

#include <memory>
#include <string>

// dynamik linker function
#include "dlfcn.h"

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"

//------------------------------------------------------------------------------
// from LINALG
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "HDF5_Tools.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_sum.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

//------------------------------------------------------------------------------

// geometry engine
//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR

#include "cl_HMR.hpp"
#include "cl_HMR_BSpline.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR_Paramfile.hpp"
#include "cl_Logger.hpp"
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
//------------------------------------------------------------------------------

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

//------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;

}

