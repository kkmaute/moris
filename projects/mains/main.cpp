/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */
#include<iostream>
#include<vector>
#include<string>
#include<fn_print.hpp>


#include <cstdio>		// nicer than streams in some respects
// C system files
#include <unistd.h>
// C++ system files
#include <stdio.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
//// TPL header files
//#ifdef PARALLEL
//#include "mpi.h"
//#endif
// MORIS header files.
#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
// other header files
//#include <catch.hpp>
//#include "fn_equal_to.hpp" //ALG
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Stopwatch.hpp" //CHR/src
#include "op_move.hpp"

#include "cl_Matrix.hpp"

#include "cl_Profiler.hpp" //profiler header
#include <Sacado_No_Kokkos.hpp>		// for FAD and RAD
#include <Eigen/Dense>
//#include <armadillo>

moris::Comm_Manager gMorisComm;

//#include "gperftools/profiler.h"

using namespace moris;
//---------------------------------------------------------------


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );
//    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
//  The main executable is just for developing and testing.
//  Please do not push this file to git.
//------------------------------------------------------------------------------
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> aMat(3,3);
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> bMat(3,3);
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cMat(3,3);

    for(int l=0; l<3; l++)
    {
    	for(int k=0; k<3; k++)
    	{
    		aMat(l,k)=l+k;
    		bMat(l,k)=2*(l+k);
    	}
    }

    Profiler tProf("/home/sonne/Desktop/temp_profile");

    	//code to profile
    for(int a=0; a<1000; a++)
    {
    	cMat=aMat*bMat;
    	std::cout<<cMat<<std::endl;
    }



    tProf.stop();





//	kcachegrind /tmp/gprofmoris.callgrind;


    // finalize MORIS global communication manager
    gMorisComm.finalize();
    return 0;
}
