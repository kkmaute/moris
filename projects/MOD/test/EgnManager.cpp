/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * EgnManager.cpp
 *
 */

#include <catch.hpp>
#include "algorithms.hpp"
#include "cl_Parameter_List.hpp"    // CON/src
#include "cl_EqnManager.hpp" // MOD/src

TEST_CASE("moris::model::EqnManager",
        "[moris],[model],[EqnManager]")
        {

    // Creating the EqnManager, and checking, if the ArgList can be read in the build_global_system
    // One element is created in get_eqn_obj and transmitted with the ArgList of get_eqn_obj to the build_global_system

    // Create Eqn Manager Object
    moris::EqnManager myManager;

    // Set number of equation objects
    myManager.set_NumberOfEqnObjects(1);

    // Define initial GlbSolVec
    moris::Mat< moris::real > GlbSolVec = {{2},{1}};

    // Initialize the global system
    //mySolverPtr.initialize_global_system(); //PSEUDO CODE

    // Clear global system
    //mySolverPtr.clear_global_system(); //PSEUDO CODE

    // Construct argument list for build_global_system
    moris::ArgListEqnMgr_build_global_system aL = myManager.mArgListEqnMgr_build_global_system;

    // Fill GlbSolVec in build_global_system argument list
    aL.GlbSolVecU = GlbSolVec;

    // Create global system
    myManager.build_global_system(aL);

    // Solve global system
    //DeltaGlbSolVec =   mySolverPtr.solve_global_system(); //PSEUDO CODE

    // Print outputs for testing purpose
//    aL.GlbSolVecU.print("aL.GlbSolVecU");
        }

