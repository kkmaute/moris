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

//#include "cl_Profiler.hpp" //profiler header
#include <Sacado_No_Kokkos.hpp>		// for FAD and RAD
//#include <Eigen/Dense>
//#include <armadillo>
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;


//#include "gperftools/profiler.h"

using namespace moris;
//---------------------------------------------------------------
moris::Matrix< DDRMat >
lCoords_4nQuad(int a)
{
	// gives the xi and eta values, respectively, for the node determined by input 'a' of a
	// 4 node quad, following convention of first node at bottom left and continuing CCW
	moris::Matrix< DDRMat > lCoord(2,4);
	lCoord(0,0)=-1.0; lCoord(0,1)=1.0; lCoord(0,2)=1.0; lCoord(0,3)=-1.0; //xi coords, node 1,2,3,4
	lCoord(1,0)=-1.0; lCoord(1,1)=-1.0; lCoord(1,2)=1.0; lCoord(1,3)=1.0; //eta coords, node 1,2,3,4

	moris::Matrix< DDRMat > nodePoint(1,2);
	nodePoint(0,0)=lCoord(0,a-1); nodePoint(0,1)=lCoord(1,a-1); // (xi,eta)
	return nodePoint;
}

moris::Matrix< DDRMat >
shapeFuncs_4nQuad(moris::Matrix< DDRMat > & nodePoint)
{
	// gives the shape functions of a 4 node quad, following convention
	// of first node at bottom left and continuing CCW
	moris::Matrix< DDRMat > Ns(1,4);
	    Ns(0,0)=0.25*(1.0-nodePoint(0,0))*(1.0-nodePoint(0,1)); //N1
	    Ns(0,1)=0.25*(1.0+nodePoint(0,0))*(1.0-nodePoint(0,1)); //N2
	    Ns(0,2)=0.25*(1.0+nodePoint(0,0))*(1.0+nodePoint(0,1)); //N3
	    Ns(0,3)=0.25*(1.0-nodePoint(0,0))*(1.0+nodePoint(0,1)); //N4
	return Ns;
}

moris::Matrix< DDRMat >
xiDerv_4nQuad(moris::Matrix< DDRMat > & nodePoint)
{
	// gives the xi partial derivatives of a 4 node quad, node 1,2,3,4
	    moris::Matrix< DDRMat >  dN_xi(1,4);
	    dN_xi(0,0)=-0.25*(1.0-nodePoint(0,1)); dN_xi(0,1)= 0.25*(1.0-nodePoint(0,1));
	    dN_xi(0,2)= 0.25*(1.0+nodePoint(0,1)); dN_xi(0,3)=-0.25*(1.0+nodePoint(0,1));
	    return dN_xi;
}

moris::Matrix< DDRMat >
etaDerv_4nQuad(moris::Matrix< DDRMat > & nodePoint)
{
	// gives the eta partial derivatives of a 4 node quad, node 1,2,3,4
	    moris::Matrix< DDRMat >  dN_eta(1,4);
	    dN_eta(0,0)=-0.25*(1.0-nodePoint(0,0)); dN_eta(0,1)=-0.25*(1.0+nodePoint(0,0));
	    dN_eta(0,2)= 0.25*(1.0+nodePoint(0,0)); dN_eta(0,3)= 0.25*(1.0-nodePoint(0,0));
	    return dN_eta;
}

moris::Matrix< DDRMat >
makeJacobian_4nQuad(moris::Matrix< DDRMat > & dN_xi, moris::Matrix< DDRMat > & dN_eta,
		moris::Matrix< DDRMat > & xs, moris::Matrix< DDRMat > & ys)
{
	// input: xi partial matrix, eta partial matrix, x coord matrix, y coord matrix
	// returns the 2x2 jacobian matrix
			moris::Matrix< DDRMat > Jac(2,2);
			for(int i=0; i<4; i++)
			{
				Jac(0,0)+=dN_xi(0,i)*xs(i,0);
			}
			for(int i=0; i<4; i++)
			{
				Jac(0,1)+=dN_xi(0,i)*ys(i,0);
			}
			for(int i=0; i<4; i++)
			{
				Jac(1,0)+=dN_eta(0,i)*xs(i,0);
			}
			for(int i=0; i<4; i++)
			{
				Jac(1,1)+=dN_eta(0,i)*ys(i,0);
			}
		return Jac;
}

double
detJ_4nQuad(moris::Matrix< DDRMat > & Jac)
{
	double detJ;
	detJ=(Jac(0,0)*Jac(1,1) - Jac(1,0)*Jac(0,1));
	return detJ;
}

moris::Matrix< DDRMat >
xDerv_4nQuad(double & detJ, moris::Matrix< DDRMat > & dN_xi,
		moris::Matrix< DDRMat > & dN_eta, moris::Matrix< DDRMat > & Jac)
{
	// input: det(Jacobian), xi partial matrix, eta partial matrix, jacobian matrix
	// returns partials of shape functions with respect to x
		moris::Matrix< DDRMat > dNx(1,4);
	    dNx(0,0)=(1/detJ)*((dN_xi(0,0)*Jac(1,1)) - (dN_eta(0,0)*Jac(0,1)));
	    dNx(0,1)=(1/detJ)*((dN_xi(0,1)*Jac(1,1)) - (dN_eta(0,1)*Jac(0,1)));
	    dNx(0,2)=(1/detJ)*((dN_xi(0,2)*Jac(1,1)) - (dN_eta(0,2)*Jac(0,1)));
	    dNx(0,3)=(1/detJ)*((dN_xi(0,3)*Jac(1,1)) - (dN_eta(0,3)*Jac(0,1)));
	    return dNx;
}

moris::Matrix< DDRMat >
yDerv_4nQuad(double & detJ, moris::Matrix< DDRMat > & dN_xi,
		moris::Matrix< DDRMat > & dN_eta, moris::Matrix< DDRMat > & Jac)
{
	// input: det(Jacobian), xi partial matrix, eta partial matrix, jacobian matrix
	// returns partials of shape functions with respect to y
		moris::Matrix< DDRMat > dNy(1,4);
		dNy(0,0)=(1/detJ)*((dN_xi(0,0)*-Jac(1,0)) + (dN_eta(0,0)*Jac(0,0)));
		dNy(0,1)=(1/detJ)*((dN_xi(0,1)*-Jac(1,0)) + (dN_eta(0,1)*Jac(0,0)));
		dNy(0,2)=(1/detJ)*((dN_xi(0,2)*-Jac(1,0)) + (dN_eta(0,2)*Jac(0,0)));
		dNy(0,3)=(1/detJ)*((dN_xi(0,3)*-Jac(1,0)) + (dN_eta(0,3)*Jac(0,0)));
	    return dNy;
}

moris::Matrix< DDRMat >
quad4_sF_partial_x(moris::Matrix< DDRMat > & xs, moris::Matrix< DDRMat > & ys, int & a)
{
// input: matrix of x coordinates, matrix of y coordinates, node to evaluate
    moris::Matrix< DDRMat > nodePoint(1,2);
    nodePoint=lCoords_4nQuad(a);

    moris::Matrix< DDRMat > sF(1,4);
    sF=shapeFuncs_4nQuad(nodePoint);

    moris::Matrix< DDRMat > dN_xi(1,4);
    dN_xi=xiDerv_4nQuad(nodePoint);

    moris::Matrix< DDRMat > dN_eta(1,4);
    dN_eta=etaDerv_4nQuad(nodePoint);

    moris::Matrix< DDRMat > Jac(2,2);
    Jac=makeJacobian_4nQuad(dN_xi, dN_eta, xs, ys);

    double detJ;
    detJ=detJ_4nQuad(Jac);

    moris::Matrix< DDRMat > dNx(1,4);
    dNx=xDerv_4nQuad(detJ, dN_xi, dN_eta, Jac);

    return dNx;
}

moris::Matrix< DDRMat >
quad4_sF_partial_y(moris::Matrix< DDRMat > & xs, moris::Matrix< DDRMat > & ys, int & a)
{
// input: matrix of x coordinates, matrix of y coordinates, node to evaluate
    moris::Matrix< DDRMat > nodePoint(1,2);
    nodePoint=lCoords_4nQuad(a);

    moris::Matrix< DDRMat > sF(1,4);
    sF=shapeFuncs_4nQuad(nodePoint);

    moris::Matrix< DDRMat > dN_xi(1,4);
    dN_xi=xiDerv_4nQuad(nodePoint);

    moris::Matrix< DDRMat > dN_eta(1,4);
    dN_eta=etaDerv_4nQuad(nodePoint);

    moris::Matrix< DDRMat > Jac(2,2);
    Jac=makeJacobian_4nQuad(dN_xi, dN_eta, xs, ys);

    double detJ;
    detJ=detJ_4nQuad(Jac);

    moris::Matrix< DDRMat > dNy(1,4);
    dNy=yDerv_4nQuad(detJ, dN_xi, dN_eta, Jac);

    return dNy;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );
//    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
//  The main executable is just for developing and testing.
//  Please do not push this file to git.
//------------------------------------------------------------------------------

//build quad 4 shape functions and derivatives:
    moris::Matrix< DDRMat > xs(4,1);
    xs(0,0)=0; xs(1,0)=2; xs(2,0)=2; xs(3,0)=0;
    moris::Matrix< DDRMat > ys(4,1);
    ys(0,0)=0; ys(1,0)=0; ys(2,0)=1; ys(3,0)=1;
    print(xs,"x vals"); print(ys,"y vals");

    for(int i=1; i<5; i++)
    {
    moris::Matrix< DDRMat > dNx(1,4);
    dNx=quad4_sF_partial_x(xs, ys, i);
    print(dNx, "SF x partials ");

    moris::Matrix< DDRMat > dNy(1,4);
    dNy=quad4_sF_partial_y(xs, ys, i);
    print(dNy, "SF y partials ");
    }



    //by hand:

    //with AD:



    // finalize MORIS global communication manager
    gMorisComm.finalize();
    return 0;
}
