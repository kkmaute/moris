/*
 * main.cpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#include <mpi.h>

#include "cl_Communication_Manager.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include "cl_Matrix.hpp"
#include "fn_comp_abs.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "fn_trans.hpp"
#include "fn_det.hpp"
#include "fn_print.hpp"

#include <iostream>
#include <ctime>


#include "linalg_typedefs.hpp"
// Global variables
moris::Comm_Manager      gMorisComm;


using namespace moris;

int
main( int    argc,
      char * argv[] )
{

    // initialize MORIS global communication manager
    gMorisComm = Comm_Manager(&argc, &argv);

    size_t its = 100000;
    std::clock_t startf;

    startf = std::clock();

    Matrix< real, DDRMat > tOut;
    for( size_t i = 0; i<its; i++)
    {
        Matrix< real, DDRMat> tMat(3,3);
        tMat(0,0) = 10.0;
        tMat(0,1) = 11.0;
        tMat(0,2) = 13.0;
        tMat(1,0) = 10.0;
        tMat(1,1) = 3.0;
        tMat(1,2) = -10.0;
        tMat(2,0) = 10.0;
        tMat(2,1) = 10.0;
        tMat(2,2) = 14.0;


        tOut = tMat*(trans(tMat)+tMat);
        tOut(0,0) = 0;
    }

    tOut(0,0)=0;

    std::cout << "Time Matrix: " << (std::clock() - startf) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

    std::clock_t    start3;

    start3 = std::clock();
    DDRMat tOutTPL;
    for( size_t i = 0; i<its; i++)
    {
        DDRMat tMatTPL(3,3);


        tMatTPL(0,0) = 10.0;
        tMatTPL(0,1) = 11.0;
        tMatTPL(0,2) = 13.0;
        tMatTPL(1,0) = -10.0;
        tMatTPL(1,1) = 3.0;
        tMatTPL(1,2) = 10.0;
        tMatTPL(2,0) = 10.0;
        tMatTPL(2,1) = 10.0;
        tMatTPL(2,2) = 14.0;

        tOutTPL(0,0) = 0;
//        tOutTPL = tMatTPL*(tMatTPL.transpose()+tMatTPL);
//        tOutTPL(0,0) = 0;
    }
    std::cout << "Time Direct TPL: " << (std::clock() - start3) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

    gMorisComm.finalize();

    return 0;
}
