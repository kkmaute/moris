/*
 * GE_Tutorial_Intersection.hpp
 *
 *  Created on: May 23, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_TUTORIALS_GE_TUTORIAL_INTERSECTION_HPP_
#define PROJECTS_GEN_TUTORIALS_GE_TUTORIAL_INTERSECTION_HPP_

//------------------------------------------------------------------------------
// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

// necessary includes for tutorial
#include "catch.hpp"
// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"
// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

//------------------------------------------------------------------------------
using namespace moris;
using namespace ge;

moris::Comm_Manager gMorisComm;

/*!
 *  (1) define a hollow sphere using two analytic sphere functions
 *  (2) create a line intersection object
 *  (3) compute intersection and other information
 */

int
main( int    argc,
      char * argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    //------------------------------------------------------------------------------
    const std::string tFileName2 = "generated:6x6x6";



    //------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;
}


#endif /* PROJECTS_GEN_TUTORIALS_GE_TUTORIAL_INTERSECTION_HPP_ */
