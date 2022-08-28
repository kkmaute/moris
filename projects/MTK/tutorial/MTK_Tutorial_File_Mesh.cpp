/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * MTK_Tutorial_File_Mesh.cpp
 *
 */

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src
#include "paths.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

/*
 * Include the necessary files from MTK
 */
#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

// Set namespaces to use
using namespace moris;
using namespace mtk;

/* @brief
 * This tutorial shows how to create an MTK mesh
 * from a file  using the mtk mesh data structure
 * in serial with multiple element types
 */
int
main( int    argc,
      char * argv[] )
{

    // Initialize the communication manager
    gMorisComm.initialize(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // Assert this is a serial tutorial
    MORIS_ASSERT(par_size() == 1,"MTK_Tutorial_Data_Mesh.cpp is only a serial test");

    /* @brief
     * The goal of this tutorial is to load an exodus file with all field
     * into MTK and declare some fields not on the file.
     */

    /*!
     * First an absolute path to the exodus file needs to be specified,
     * For this example we use the MORISROOT environment variable that
     * is declared in .cshrc or .bashrc file. If this is an mpirun
     * then there needs to be a decomposed exodus file with *.e.np.pr
     * where np is the number of procs in the mpi run and pr is the proc
     * rank.
     *
     * In this example a mesh file called sandwich.e which has multiple block
     * sets, side sets and node sets is used.
     * @image html ./figures/sandwich_base_mesh.png "Sandwich Exodus Mesh"
     *
     * \code{.cpp}
     *   std::string tPrefix = moris::get_base_moris_dir();
     *   std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";
     * \endcode
     */
    std::string tPrefix = moris::get_base_moris_dir();
    std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";

    std::cout<<"Mesh input name = "<< tMeshFileName<<std::endl;

    /*
     * If there are fields that will be populated that are not members of the
     * exodus file to be loaded in, they can be declared as follows
     */

    /*
     * Load the file and the mesh fields to declare into the MTK data structure
     * using the mesh factory
     */
    moris::mtk::Mesh* tBackgroundMesh = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );

}

