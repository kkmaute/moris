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
#include "op_times.hpp"

#include <iostream>
#include <ctime>


//#include <armadillo>

//#include "Arma_Impl/cl_Matrix_Arma_Dynamic.hpp"
#include "Eigen_Impl/cl_Matrix_Eigen_3x3.hpp"
#include "Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp"
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

    size_t its = 1000000;
    std::clock_t    startf;

    startf = std::clock();

    Mat_New< real, DDRMat> tOut;
    for( size_t i = 0; i<its; i++)
    {
        Mat_New< real, DDRMat> tMatFixed(3,3);
        tMatFixed(0,0) = 10.0;
        tMatFixed(0,1) = 11.0;
        tMatFixed(0,2) = 13.0;
        tMatFixed(1,0) = 10.0;
        tMatFixed(1,1) = 3.0;
        tMatFixed(1,2) = 10.0;
        tMatFixed(2,0) = 10.0;
        tMatFixed(2,1) = 10.0;
        tMatFixed(2,2) = 14.0;

        tOut = tMatFixed*tMatFixed;
        tOut(0,0) = 0;
    }

    tOut(0,0)=0;
    std::cout << "Time Mat_New: " << (std::clock() - startf) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;


    std::clock_t    start3;

    start3 = std::clock();
    for( size_t i = 0; i<its; i++)
    {
        DDRMat tMatEig(3,3);

        tMatEig(0,0) = 10.0;
        tMatEig(0,1) = 10.0;
        tMatEig(0,2) = 10.0;
        tMatEig(1,0) = 10.0;
        tMatEig(1,1) = 10.0;
        tMatEig(1,2) = 10.0;
        tMatEig(2,0) = 10.0;
        tMatEig(2,1) = 10.0;
        tMatEig(2,2) = 10.0;

        DDRMat tOut = tMatEig*tMatEig;
        tOut(0,0) = 0;
    }
    std::cout << "Time Direct TPL: " << (std::clock() - start3) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;



    gMorisComm.finalize();

    return 0;
}
