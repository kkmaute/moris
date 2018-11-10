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
#include <armadillo>

#include <Sacado_No_Kokkos.hpp>		// for FAD and RAD
//#include "Eigen/Dense"
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
// TPL header files
#ifdef PARALLEL
#include "mpi.h"
#endif
// MORIS header files.
#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
// other header files
//#include <catch.hpp>
//#include "fn_equal_to.hpp" //ALG
//#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Stopwatch.hpp" //CHR/src
moris::Comm_Manager gMorisComm;
using namespace moris;
//---------------------------------------------------------------
typedef Sacado::Fad::DFad<double>   F;
typedef Sacado::Fad::SFad<double,2> F2; // FAD with # of ind. vars fixed at 2
typedef Sacado::Rad::ADvar<double>  R;  // for RAD

template <typename T>
const T func2( T &a, T &b) 	// sample function of 2 variables
{

	return sqrt(a*a + 2*b*b);
}
template <typename T>
const T func(int n, T *x)	// sample function of n variables
				// == func2 when n == 2
{
	int i;
	T t = 0;
	for(i = 1; i < n; i++)
		t += i*x[i]*x[i];
	return sqrt(t);
}
//---------------------------------------------------------
//---------------------------------------------------------
//void
//timeTest1( Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat,
//		Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat1 )
//{
//	fmat1 = fmat*fmat*fmat*fmat*fmat;
//}
//void
//timeTest2( Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat,
//		Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat1 )
//{
//	rmat1 = rmat*rmat*rmat*rmat*rmat;
//}
//void
//timeTest3( Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat,
//		Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat1 )
//{
//	Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat2 (fmat);
//	Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat3 (fmat);
//	Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat4 (fmat);
//
//	fmat1 = fmat*fmat;
//	fmat2 = fmat1*fmat;
//	fmat3 = fmat2*fmat;
//	fmat4 = fmat3*fmat;
//}
//void
//timeTest4( Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat,
//		Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat1 )
//{
//	Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat2 (rmat);
//	Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat3 (rmat);
//	Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat4 (rmat);
//
//	rmat1 = rmat*rmat;
//	rmat2 = rmat1*rmat;
//	rmat3 = rmat2*rmat;
//	rmat4 = rmat3*rmat;
//}
//void
//timeTest5(Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat,
//		Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fmat1)
//{
//	fmat1 = fmat.cwiseProduct(fmat.cwiseProduct(fmat.cwiseProduct(fmat.cwiseProduct(fmat))));
//}
//void
//timeTest6(Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat,
//		Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> rmat1 )
//{
//	rmat1 = rmat.cwiseProduct(rmat.cwiseProduct(rmat.cwiseProduct(rmat.cwiseProduct(rmat))));
//}
void
timeTestFAD_array( F a[], F & b)
{
	for (int t = 0; t < 100; t++)
	{
		b *= a[t];
	}
}
void
timeTestReal_array( double a[], double & b)
{
	for (int t = 0; t < 100; t++)
	{
		b *= a[t];
	}
}
void
timeTestFAD_array1( F b[10][10], F result[10][10] )
{
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			for (int k=0; k<10; k++)
			{
				result[i][j] += b[i][k]*b[k][j]; //matrix multiplication by hand
			}
		}
	}
}
void
timeTestReal_array1( double b[10][10], double result[10][10] )
{
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			for (int k=0; k<10; k++)
			{
				result[i][j] += b[i][k]*b[k][j]; //matrix multiplication by hand
			}
		}
	}
}

//---------------------------------------------------------
//---------------------------------------------------------

// demo of FAD (forward-mode AD), with general number of ind. vars
void
Fad_demo()
{
	F a, b, x[5], y;
	int i, n;

	printf("Fad_demo...\n\n");

	// first try n == 2
	a = 1.;
	b = 2.;
	// indicate the independent variables, and initialize their partials to 1:
	a.diff(0,2);	// 0 ==> this is the first independent var. of 2
	b.diff(1,2);	// 1 ==> this is the second ind. var. of 2

	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	// When we know the result was not constant (i.e., did involve ind. vars)
	// or when hasFastAccess() is true, we access partials more quickly
	// by using member function fastAccessDx rather than dx

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
				y.fastAccessDx(0), y.fastAccessDx(1));

	// Similar exercise with general n, in this case n == 5
	n = 5;
	for(i = 0; i < n; i++) {
		x[i] = i;
		x[i].diff(i, n);
	}
	y = func(n, x);
	printf("\nfunc(5,x) for x = (0,1,2,3,4) = %g\n", y.val());
	for(i = 0; i < n; i++)
		printf("d func / d x[%d] = %g == %g\n", i, y.dx(i), y.fastAccessDx(i));
}
// Fad_demo2 == repeat first part of Fad_Demo with type F2 instead of F
// i.e., with fixed-size allocations
void
Fad2_demo()
{
	F2 a, b, y;

	printf("\n\nFad2_demo...\n\n");

	a = 1.;
	b = 2.;
	// indicate the independent variables, and initialize their partials to 1:
	a.diff(0,2);	// 0 ==> this is the first independent var., of 2
	b.diff(1,2);	// 1 ==> this is the second ind. var.

	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
				y.fastAccessDx(0), y.fastAccessDx(1));
}
// Fad_demo3 == repeat of Fad_Demo2 with a different constructor, one that
// indicates the independent variables and their initial values
// and removes the need to invoke .diff()
void
Fad3_demo()
{
	F2 a(2,0,1.), b(2,1,2.), y;

	printf("\n\nFad3_demo...\n\n");

	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
				y.fastAccessDx(0), y.fastAccessDx(1));
}
// Rad_demo == repeat of Fad_Demo with type R instead of F,
// i.e., with reverse-mode rather than forward-mode AD
void
Rad_demo()
{
	R a, b, x[5], y;
	int i, n;

	printf("\n\nRad_demo...\n\n");

	// first try n == 2
	a = 1.;
	b = 2.;

	y = func2(a,b);

	R::Gradcomp();	// do the reverse sweep

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", a.adj(), b.adj());

	// Similar exercise with general n, in this case n == 5
	n = 5;
	for(i = 0; i < n; i++)
		x[i] = i;
	y = func(n, x);
	printf("\nfunc(5,x) for x = (0,1,2,3,4) = %g\n", y.val());

	// the .val() values are always available; we must call Gradcomp
	// before accessing the adjoints, i.e., the .adj() values...

	R::Gradcomp();

	for(i = 0; i < n; i++)
		printf("d func / d x[%d] = %g\n", i, x[i].adj());
}

//Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>
//EigenMoveFunc( Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> & temp ){
//	Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> moved;
//	moved = std::move(temp);
//	return moved;
//}
arma::Mat< real >
ArmaMoveFunc( arma::Mat< real > & temp ){
	arma::Mat< real > moved;
	moved = std::move(temp);
	return moved;
}
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
//-------------------initialize vectors/matrices-------------------//
//    F  fad[100], newfad, fad1[10][10], newfad1[10][10];
//    for (int k = 0; k < 100; k++)
//    {
//    	fad[k] = k+1;
////    	fad[k].diff(k,100);
//    }
//    printf("\nfad element [26] = %g\n", fad[25].val());
//    printf("initiated partial of fad element [24] as:  %g\n\n", fad[25].dx(25));
//    Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fadmat(10,10);
//    int n = 0;
//    for (int i = 0; i < 10; i++)
//    {
//    	for (int j = 0; j < 10; j++)
//    	{
//    		fadmat(i,j) = fad[n];
//    		n++;
//    	}
//    }
//    Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> fadmat1(10,10);
//    fadmat1.fill( 0 );
//    int m = 0;
//    for (int a=0; a<10; a++)
//    {
//    	for (int b=0; b<10; b++)
//    	{
//    		fad1[a][b] = fad[m];
//    		m++;
//    	}
//    }
////---end fad initialization---
//    double stan[100], stan1[10][10], newstan1[10][10];
//    for (int kk = 0; kk < 100; kk++)
//    {
//    	stan[kk] = kk+1;
//    }
//    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> stanmat(10,10);
//    int nn = 0;
//    for (int i = 0; i < 10; i++)
//    {
//    	for (int j = 0; j < 10; j++)
//    	{
//    		stanmat(i,j) = stan[nn];
//    		nn++;
//    	}
//    }
//    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> stanmat1(10,10);
//    stanmat1.fill( 0 );
//    int mm = 0;
//    for (int i=0; i<10; i++)
//    {
//    	for (int j=0; j<10; j++)
//    	{
//    		stan1[i][j] = stan[mm];
//    		mm++;
//    	}
//    }
//
//    double newstan;
////---end stan initialization---
//
////-------------------end initialization-------------------//
////-------------------first
//    moris::tic tTimer0;
//    for (int loop = 0; loop < 10000; loop++)
//    {
//    	timeTestFAD_array( fad, newfad );
//    }
//    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//
//    tElapsedTime0 /= 10000;
//
//    std::cout<< "FAD data type, array element multiplication:  "
//    		<< tElapsedTime0 << std::endl << std::endl;
////-------------------second
//    moris::tic tTimer00;
//    for (int loop = 0; loop < 10000; loop++)
//    {
//    	timeTestReal_array( stan, newstan );
//    }
//    real tElapsedTime00 = tTimer00.toc<moris::chronos::milliseconds>().wall;
//
//    tElapsedTime00 /= 10000;
//
//    std::cout<< "stan data type, array element multiplication:  "
//    		<< tElapsedTime00 << std::endl << std::endl;
////-------------------third
//    moris::tic tTimer000;
//    for (int loop = 0; loop < 10000; loop++)
//    {
//    	timeTestFAD_array1( fad1, newfad1 );
//    }
//    real tElapsedTime000 = tTimer000.toc<moris::chronos::milliseconds>().wall;
//
//    tElapsedTime000 /= 10000;
//
//    std::cout<< "FAD data type, 10x10 matrix multiplication by hand:  "
//    		<< tElapsedTime000 << std::endl << std::endl;
////-------------------fourth
//    moris::tic tTimer0000;
//    for (int loop = 0; loop < 10000; loop++)
//    {
//    	timeTestReal_array1( stan1, newstan1 );
//    }
//    real tElapsedTime0000 = tTimer0000.toc<moris::chronos::milliseconds>().wall;
//
//    tElapsedTime0000 /= 10000;
//
//    std::cout<< "stan data type, 10x10 matrix multiplication by hand:  "
//    		<< tElapsedTime0000 << std::endl << std::endl;
//
////------------------------------------------------------------------------------

    class MYCLASS{
    public:
    	std::string str = "my string";
    	double a = 10;
    };
    MYCLASS myclass;
    std::cout<< "before move str = "<< myclass.str << std::endl;
    std::cout<< "and a = "<< myclass.a << std::endl;

    std::string s = std::move(myclass.str);
    double b = std::move(myclass.a);

    std::cout<< "after move str = "<< myclass.str << std::endl;
    std::cout<< "and a = "<< myclass.a << std::endl;

    std::cout<< "after move, b = "<< b <<std::endl;
    std::cout<< "after move, s = "<< s <<std::endl;

//    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> temp(10,10);
//    double n=0.;
//    for (int i=0; i<10; i++)
//    {
//    	for (int j=0; j<10; j++)
//    	{
//    		temp(i,j) = n;
//    		n++;
//    	}
//    }
//    std::cout<<"before move, temp = \n"<<temp<<std::endl;
//    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> moved(10,10);
//    moved=EigenMoveFunc(temp);
//
//    std::cout<<"after move, temp = \n"<<temp<<std::endl;
//    std::cout<<"after move, moved = \n"<<moved<<std::endl;

    arma::Mat< real > temp(10,10);
    double n=0.;
    for (int i=0; i<10; i++)
    {
    	for (int j=0; j<10; j++)
    	{
    		temp(i,j) = n;
    		n++;
    	}
    }
    std::cout<<"before move, temp = \n"<<temp<<std::endl;
    arma::Mat< real > moved;
    moved=ArmaMoveFunc(temp);
    std::cout<<"after move, temp = \n"<<temp<<std::endl;
    std::cout<<"after move, moved = \n"<<moved<<std::endl;


    // finalize MORIS global communication manager
    gMorisComm.finalize();
    return 0;
}
