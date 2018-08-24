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

#include "cl_LINALG_Parameters.hpp"
#include "cl_LINALG_Enums.hpp"
#include "cl_Matrix.hpp"

//#include <armadillo>
//#include <Eigen/Dense>
#include <iostream>


// Global variables
moris::Comm_Manager      gMorisComm;

//namespace moris
//{
//Linalg_Parameters gParameters;
//}


using namespace moris;

int
main( int    argc,
      char * argv[] )
{

    // initialize MORIS global communication manager
    gMorisComm = Comm_Manager(&argc, &argv);

//    // Set up linear algebra parameters
//    moris::gParameters = moris::Linalg_Parameters(Backend_Dense_Matrix::ARMA_DYNAMIC);
//
//    Mat_New<double,arma::Mat<double>> tMat(4,4);
//
//    std::cout<<"tMat(0,0) = " << tMat(0,0)<<std::endl;
//    tMat(0,0) = 10.0;
//    std::cout<<"tMat(0,0) = " << tMat(0,0)<<std::endl;
//
//
//    moris::gParameters = moris::Linalg_Parameters(Backend_Dense_Matrix::EIGEN_DYNAMIC);
//    Mat_New<double,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> tMat2(4,4);
//
//    std::cout<<"tMat(0,0) = " << tMat2(0,0)<<std::endl;
//    tMat2(0,0) = 10.0;
//    std::cout<<"tMat(0,0) = " << tMat2(0,0)<<std::endl;

    gMorisComm.finalize();

    return 0;
}
