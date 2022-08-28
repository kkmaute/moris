/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bspline.cpp
 *
 */

#include "cl_Bspline.hpp" // MOD/src
using namespace moris;

Mat<real>
Bspline::build_spline_1d(real const & aXi,
        Mat<real> & aKnotVector,
        uint const & aPolynomialDegree)
{
    uint n = aKnotVector.length()-aPolynomialDegree-1;   // Number of basis functions
    Mat<real> Bspline_all(aKnotVector.length()-1,aPolynomialDegree+1,0); // All basis functions, which are needed to create the specific polynomial degree
    Mat<real> Bspline(1,n,0);
    uint count=0; // Counter is needed for a loop
    real a,b;

    if(aXi==aKnotVector(aKnotVector.length()-1)) // If aXi is equal to the last knot
    {
        for (uint  i = aKnotVector.length()-1; i>=0;i--) //count number of repeated last knot
        {
            if(aKnotVector(i-1)<aKnotVector(i))
            {
                break;
            }
            count++;
        }
        Bspline_all(n-count+aPolynomialDegree,0)=1; // Raise the basis to 1 at last knot
    }

    for (uint  i = 0; i<aKnotVector.length()-1;i++) //Loop over constant basis functions p = 0
    {
        if(aXi>=aKnotVector(i) && aXi<aKnotVector(i+1))
        {
            Bspline_all(i,0) = 1;
        }
    }
    for (uint  p = 1; p<=aPolynomialDegree;p++) //Loop over degree of basis functions
    {
        for (uint  i = 0; i<aKnotVector.length()-p-1;i++) //Loop over degree of basis functions
        {
            if(Bspline_all(i,p-1)!=0 && (aKnotVector(i+p)-aKnotVector(i))!=0)
            {
                a = (aXi-aKnotVector(i))/(aKnotVector(i+p)-aKnotVector(i)); // first part of the Cox-de Boor recursion formula
            }
            else
            {
                a = 0;
            }
            if(Bspline_all(i+1,p-1)!=0 && (aKnotVector(i+p+1)-aKnotVector(i+1))!=0)
            {
                b = (aKnotVector(i+p+1)-aXi)/(aKnotVector(i+p+1)-aKnotVector(i+1));  // second part of the Cox-de Boor recursion formula
            }
            else
            {
                b = 0;
            }
            Bspline_all(i,p) = a*Bspline_all(i,p-1) + b*Bspline_all(i+1,p-1); // Sum of the Cox-de Boor recursion formula
        }
    }
    for (uint  i = 0; i<n;i++) //Loop over degree of basis functions
    {
        Bspline(i) = Bspline_all(i,aPolynomialDegree); // Extract only the basis functions of the polynomial degree
    }
    return Bspline;
}

Mat<real>
Bspline::build_spline_deriv_1d(real const & aXi,
        Mat<real> & aKnotVector,
        uint const & aPolynomialDegree)
{
    uint tLengthKnotvector = aKnotVector.length();
    Mat<real> tBsplineDeriv(tLengthKnotvector-aPolynomialDegree-1,1,0);
    real a,b;
    Mat<real> tBspline = Bspline::build_spline_1d(aXi,aKnotVector,aPolynomialDegree-1);  //B-Splines of a lower polynomial degree are required and computed for p = p-1
    for (uint  i = 0; i<tLengthKnotvector-aPolynomialDegree-1;i++) //Loop over degree of basis functions
    {
        if((aKnotVector(i+(real)aPolynomialDegree)-aKnotVector(i))!=0.0)
        {
            a = (aPolynomialDegree)/(aKnotVector(i+aPolynomialDegree)-aKnotVector(i)); // first part of the Cox-de Boor recursion formula
        }
        else
        {
            a = 0;
        }
        if((aKnotVector(i+(real)aPolynomialDegree+1)-aKnotVector(i+1))!=0.0)
        {
            b = (aPolynomialDegree)/(aKnotVector(i+aPolynomialDegree+1)-aKnotVector(i+1));  // second part of the Cox-de Boor recursion formula
        }
        else
        {
            b = 0;
        }
        tBsplineDeriv(i) = a*tBspline(i) - b*tBspline(i+1); // Sum of the Cox-de Boor recursion formula
    }
    return tBsplineDeriv;
}

Mat<real>
Bspline::build_spline_nd(Mat<real> const & aXi,
        Mat<real> & aKnotVector1,
        Mat<real> & aKnotVector2,
        Mat<real> & aKnotVector3,
        Mat<uint> const & aPolynomialDegree,
        Mat<uint> const & aElement_flag,
        uint const & aDim)
{
    Mat<real> tbspline;
    Mat<real> tbspline1;
    Mat<real> tbspline2;
    Mat<real> tbspline3;
    uint count = 0; // Counter is needed for a loop
    if(aDim>=2 && aDim<=3)
    {
        tbspline = build_spline_1d(aXi(0),aKnotVector1,aPolynomialDegree(0));  //B-Splines for the direction Xi_1
        tbspline1 = tbspline.cols( aElement_flag(0), aElement_flag(0)+aPolynomialDegree(0));
        tbspline = build_spline_1d(aXi(1),aKnotVector2,aPolynomialDegree(1));  //B-Splines for the direction Xi_2
        tbspline2 = tbspline.cols( aElement_flag(1), aElement_flag(1)+aPolynomialDegree(1));
        if(aDim==3)
        {
            tbspline = build_spline_1d(aXi(2),aKnotVector3,aPolynomialDegree(2));  //B-Splines for the direction Xi_3
            tbspline3 = tbspline.cols( aElement_flag(2), aElement_flag(2)+aPolynomialDegree(2));
        }
    }
    else
    {
        MORIS_LOG_ERROR << "Dimension is out of range 1 < dim < 4 ";
    }
    Mat<real> bspline;
    if( aDim == 2)
    {
        bspline.set_size(1,(aPolynomialDegree(0)+1)*(aPolynomialDegree(1)+1),0); // Initialize the size of the vector
        for (uint  j = 0; j<aPolynomialDegree(1)+1;j++) //Loop over basis functions in Xi_2 direction
        {
            for (uint  k = 0; k<aPolynomialDegree(0)+1;k++) //Loop over basis functions in Xi_1 direction
            {
                bspline(count) = tbspline1(k)*tbspline2(j);
                count++;
            }
        }
    }
    else if( aDim == 3)
    {
        bspline.set_size(1,(aPolynomialDegree(0)+1)*(aPolynomialDegree(1)+1)*(aPolynomialDegree(2)+1),0); // Initialize the size of the vector
        for (uint  i = 0; i<aPolynomialDegree(2)+1;i++) //Loop over basis functions in Xi_3 direction
        {
            for (uint  j = 0; j<aPolynomialDegree(1)+1;j++) //Loop over basis functions in Xi_2 direction
            {
                for (uint  k = 0; k<aPolynomialDegree(0)+1;k++) //Loop over basis functions in Xi_1 direction
                {
                    bspline(count) = tbspline1(k)*tbspline2(j)*tbspline3(i);
                    count++;
                }
            }
        }
    }
    return bspline;
}

Mat<real>
Bspline::build_spline_deriv_nd(Mat<real> const & aXi,
        Mat<real> & aKnotVector1,
        Mat<real> & aKnotVector2,
        Mat<real> & aKnotVector3,
        Mat<uint> const & aPolynomialDegree,
        Mat<uint> const & aElement_flag,
        uint const & aDim)
{
    Mat<real> tbspline;
    Mat<real> tbspline1;
    Mat<real> tbspline2;
    Mat<real> tbspline3;
    Mat<real> ttBsplineDeriv;
    Mat<real> ttBsplineDeriv1;
    Mat<real> ttBsplineDeriv2;
    Mat<real> ttBsplineDeriv3;
    uint count = 0; // Counter is needed for a loop
    if(aDim>=2 && aDim<=3)
    {
        tbspline = build_spline_1d(aXi(0),aKnotVector1,aPolynomialDegree(0));  //B-Splines for the direction Xi_1
        tbspline1 = tbspline.cols( aElement_flag(0), aElement_flag(0)+aPolynomialDegree(0));
        ttBsplineDeriv = build_spline_deriv_1d(aXi(0),aKnotVector1,aPolynomialDegree(0));  //Derivative of the B-Splines for the direction Xi_1
        ttBsplineDeriv1 = ttBsplineDeriv.cols( aElement_flag(0), aElement_flag(0)+aPolynomialDegree(0));
        tbspline = build_spline_1d(aXi(1),aKnotVector2,aPolynomialDegree(1));  //B-Splines for the direction Xi_2
        tbspline2 = tbspline.cols( aElement_flag(1), aElement_flag(1)+aPolynomialDegree(1));
        ttBsplineDeriv = build_spline_deriv_1d(aXi(1),aKnotVector2,aPolynomialDegree(1));  //Derivative of the B-Splines for the direction Xi_2
        ttBsplineDeriv2 = ttBsplineDeriv.cols( aElement_flag(1), aElement_flag(1)+aPolynomialDegree(1));
        if(aDim==3)
        {
            tbspline = build_spline_1d(aXi(2),aKnotVector3,aPolynomialDegree(2));  //B-Splines for the direction Xi_3
            tbspline3 = tbspline.cols( aElement_flag(2), aElement_flag(2)+aPolynomialDegree(2) );
            ttBsplineDeriv = build_spline_deriv_1d(aXi(2),aKnotVector3,aPolynomialDegree(2));  //Derivative of the B-Splines for the direction Xi_2
            ttBsplineDeriv3 = ttBsplineDeriv.cols( aElement_flag(2), aElement_flag(2)+aPolynomialDegree(2));
        }
    }
    else
    {
        MORIS_LOG_ERROR << "Dimension is out of range 1 < dim < 4 ";
    }

    Mat<real> tBsplineDeriv;
    if( aDim == 2)
    {
        tBsplineDeriv.set_size(aDim,(aPolynomialDegree(0)+1)*(aPolynomialDegree(1)+1),0); // Initialize the size of the vector
        for (uint  j = 0; j<aPolynomialDegree(1)+1;j++) //Loop over basis functions in Xi_2 direction
        {
            for (uint  k = 0; k<aPolynomialDegree(0)+1;k++) //Loop over basis functions in Xi_1 direction
            {
                tBsplineDeriv(0,count) = ttBsplineDeriv1(k)*tbspline2(j);
                tBsplineDeriv(1,count) = tbspline1(k)*ttBsplineDeriv2(j);
                count++;
            }
        }
    }
    else if( aDim == 3)
    {
        tBsplineDeriv.set_size(aDim,(aPolynomialDegree(0)+1)*(aPolynomialDegree(1)+1)*(aPolynomialDegree(2)+1),0); // Initialize the size of the vector
        for (uint  i = 0; i<aPolynomialDegree(2)+1;i++) //Loop over basis functions in Xi_3 direction
        {
            for (uint  j = 0; j<aPolynomialDegree(1)+1;j++) //Loop over basis functions in Xi_2 direction
            {
                for (uint  k = 0; k<aPolynomialDegree(0)+1;k++) //Loop over basis functions in Xi_1 direction
                {
                    tBsplineDeriv(0,count) = ttBsplineDeriv1(k)*tbspline2(j)*tbspline3(i);
                    tBsplineDeriv(1,count) = tbspline1(k)*ttBsplineDeriv2(j)*tbspline3(i);
                    tBsplineDeriv(2,count) = tbspline1(k)*tbspline2(j)*ttBsplineDeriv3(i);
                    count++;
                }
            }
        }
    }
    return tBsplineDeriv;
}

Mat<real>
Bspline::build_spline_uniform_nd(Mat<real> const& aXi,
        Mat<real> & aKnotVector,
        uint const & aPolynomialDegree,
        Mat<uint> const & aElement_flag,
        uint const & aDim)
{
    Mat<real> bspline_uniform(pow(aPolynomialDegree+1,aDim),1,0); // Initialize the size of the vector
    Mat<real> tbspline;
    Mat<real> tbspline1;
    Mat<real> tbspline2;
    Mat<real> tbspline3;
    uint count = 0; // Counter is needed for a loop
    if(aDim>=2 && aDim<=3)
    {
        tbspline = build_spline_1d(aXi(0),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_1
        tbspline1 = tbspline.cols( aElement_flag(0), aElement_flag(0)+aPolynomialDegree);
        tbspline = build_spline_1d(aXi(1),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_2
        tbspline2 = tbspline.cols( aElement_flag(1), aElement_flag(1)+aPolynomialDegree);
        if(aDim==3)
        {
            tbspline = build_spline_1d(aXi(2),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_3
            tbspline3 = tbspline.cols( aElement_flag(2), aElement_flag(2)+aPolynomialDegree);
        }
    }
    else
    {
        MORIS_LOG_ERROR << "Dimension is out of range 1 < dim < 4 ";
    }

    if(aDim==2) // Creates 2D basis functions
    {
        for (uint  j = 0; j<aPolynomialDegree+1;j++) //Loop over basis functions in Xi_2 direction
        {
            for (uint  k = 0; k<aPolynomialDegree+1;k++) //Loop over basis functions in Xi_1 direction
            {
                bspline_uniform(count) = tbspline1(k)*tbspline2(j);
                count++;
            }
        }
    }
    else if(aDim==3) // Creates 3D basis functions
    {
        for (uint  i = 0; i<aPolynomialDegree+1;i++) //Loop over basis functions in Xi_3 direction
        {
            for (uint  j = 0; j<aPolynomialDegree+1;j++) //Loop over basis functions in Xi_2 direction
            {
                for (uint  k = 0; k<aPolynomialDegree+1;k++) //Loop over basis functions in Xi_1 direction
                {
                    bspline_uniform(count) = tbspline1(k)*tbspline2(j)*tbspline3(i);
                    count++;
                }
            }
        }
    }
    return bspline_uniform;
}

Mat<real>
Bspline::build_spline_uniform_1d(
        real const & aXi,
        uint const & aPolynomialDegree)
{
    Mat<real> tBspline(aPolynomialDegree+1,1,0);
    if( aPolynomialDegree == 1 )
    {
        tBspline(0) = 1 - aXi;
        tBspline(1) = aXi;
    }
    else if( aPolynomialDegree == 2 )
    {
        tBspline(0) = (1-aXi)*(1-aXi)*0.5;
        tBspline(1) = -0.5*(2*aXi*aXi-2*aXi-1);
        tBspline(2) = 0.5*aXi*aXi;
    }
    else if( aPolynomialDegree == 3 )
    {
        tBspline(0) = 1.0/6.0*pow(1-aXi,3);
        tBspline(1) = 1.0/6.0*pow(1-aXi,2)*(aXi+2.0)-1.0/6.0*(2.0-aXi)*(2.0*aXi*aXi-2.0*aXi-1.0);
        tBspline(2) = 1.0/6.0*(3.0-aXi)*aXi*aXi-1.0/6.0*(aXi+1.0)*(2.0*aXi*aXi-2*aXi-1.0);
        tBspline(3) = aXi*aXi*aXi/6.0;
    }
    else
    {
        Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
        for(uint i = 0; i < tKnotVector.length(); i++)
            tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
        tBspline = build_spline_1d(aXi,tKnotVector,aPolynomialDegree);
    }
    return tBspline;
}

Mat<real>
Bspline::build_spline_uniform_nd(
        Mat<real> const& aXi,
        uint const & aPolynomialDegree,
        uint const & aDim)
{
    Mat<real> tBspline(pow(aPolynomialDegree+1,aDim),1,0); // Initialize the size of the vector
    if( aDim == 2)
    {
        if( aPolynomialDegree == 1 )
        {
            tBspline(0) = (1 - aXi(0)) * (1 - aXi(1));
            tBspline(1) = aXi(0)       * (1 - aXi(1));
            tBspline(2) = (1 - aXi(0)) * aXi(1);
            tBspline(3) = aXi(0)       * aXi(1);
        }
        else if( aPolynomialDegree == 2 )
        {
            tBspline(0) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5);
            tBspline(1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5);
            tBspline(2) = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5);
            tBspline(3) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1)));
            tBspline(4) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1)));
            tBspline(5) = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1)));
            tBspline(6) = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1));
            tBspline(7) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1));
            tBspline(8) = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1));
        }
        else if( aPolynomialDegree == 3 )
        {
            tBspline(0) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * 1.0/6.0*pow(1.0-aXi(1),3);
            tBspline(1) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * 1.0/6.0*pow(1.0-aXi(1),3);
            tBspline(2) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * 1.0/6.0*pow(1.0-aXi(1),3);
            tBspline(3) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * 1.0/6.0*pow(1.0-aXi(1),3);
            tBspline(4) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));
            tBspline(5) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));
            tBspline(6) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));
            tBspline(7) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));
            tBspline(8) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));
            tBspline(9) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));
            tBspline(10) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))      * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));
            tBspline(11) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                        * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));
            tBspline(12) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                       * aXi(1)*aXi(1)*aXi(1)/6.0;
            tBspline(13) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0))* aXi(1)*aXi(1)*aXi(1)/6.0;
            tBspline(14) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))      * aXi(1)*aXi(1)*aXi(1)/6.0;
            tBspline(15) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                        * aXi(1)*aXi(1)*aXi(1)/6.0;
        }
        else
        {
            Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
            for(uint i = 0; i < tKnotVector.length(); i++)
                tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
            Mat<uint> tElement_flag(aDim,1,0);
            tBspline = build_spline_uniform_nd(aXi,tKnotVector,aPolynomialDegree,tElement_flag,aDim);
        }
    }
    else if( aDim == 3)
    {
        if( aPolynomialDegree == 1 )
        {
            tBspline(0) = (1 - aXi(0)) * (1 - aXi(1)) * (1 - aXi(2));
            tBspline(1) = aXi(0)       * (1 - aXi(1)) * (1 - aXi(2));
            tBspline(2) = (1 - aXi(0)) * aXi(1)       * (1 - aXi(2));
            tBspline(3) = aXi(0)       * aXi(1)       * (1 - aXi(2));
            tBspline(4) = (1 - aXi(0)) * (1 - aXi(1)) * aXi(2);
            tBspline(5) = aXi(0)       * (1 - aXi(1)) * aXi(2);
            tBspline(6) = (1 - aXi(0)) * aXi(1)       * aXi(2);
            tBspline(7) = aXi(0)       * aXi(1)       * aXi(2);
        }
        else if( aPolynomialDegree == 2 )
        {
            tBspline(0)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(1)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(2)  = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(3)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(4)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(5)  = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(6)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(7)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(8)  = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);
            tBspline(9)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(10) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(11) = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(12) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(13) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(14) = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(15) = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(16) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(17) = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(18) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));
            tBspline(19) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));
            tBspline(20) = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));
            tBspline(21) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));
            tBspline(22) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));
            tBspline(23) = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));
            tBspline(24) = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));
            tBspline(25) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));
            tBspline(26) = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));
        }
        else if( aPolynomialDegree == 3 )
        {
            tBspline(0)  = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(1)  = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(2)  = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(3)  = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(4)  = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(5)  = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(6)  = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(7)  = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(8)  = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(9)  = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(10) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(11) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(12) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(13) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(14) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(15) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);
            tBspline(16) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(17) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(18) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(19) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(20) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(21) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(22) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(23) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(24) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(25) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(26) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(27) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(28) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(29) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(30) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(31) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(32) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(33) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(34) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(35) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(36) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(37) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(38) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(39) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(40) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(41) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(42) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(43) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(44) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(45) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(46) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(47) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));
            tBspline(48) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(49) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(50) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(51) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(52) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(53) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(54) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(55) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(56) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(57) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(58) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(59) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(60) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(61) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(62) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;
            tBspline(63) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;
        }
        else
        {
            Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
            for(uint i = 0; i < tKnotVector.length(); i++)
                tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
            Mat<uint> tElement_flag(aDim,1,0);
            tBspline = build_spline_uniform_nd(aXi,tKnotVector,aPolynomialDegree,tElement_flag,aDim);
        }
    }
    return tBspline;
}

Mat<real>
Bspline::build_spline_deriv_uniform_1d(
        real const & aXi,
        uint const & aPolynomialDegree)
{
    Mat<real> tBsplineDeriv(aPolynomialDegree+1,1,0);
    if( aPolynomialDegree == 1 )
    {
        tBsplineDeriv(0) = -1.0;
        tBsplineDeriv(1) = 1.0;
    }
    else if( aPolynomialDegree == 2 )
    {
        tBsplineDeriv(0) = aXi-1.0;
        tBsplineDeriv(1) = 1.0-2.0*aXi;
        tBsplineDeriv(2) = aXi;
    }
    else if( aPolynomialDegree == 3 )
    {
        tBsplineDeriv(0) = -0.5*pow(aXi-1.0,2);
        tBsplineDeriv(1) = 1.5*aXi*aXi-2.0*aXi;
        tBsplineDeriv(2) = -1.5*aXi*aXi+aXi+0.5;
        tBsplineDeriv(3) = aXi*aXi/2.0;
    }
    else
    {
        Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
        for(uint i = 0; i < tKnotVector.length(); i++)
            tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
        tBsplineDeriv = build_spline_deriv_1d(aXi,tKnotVector,aPolynomialDegree);
    }
    return tBsplineDeriv;
}

Mat<real>
Bspline::build_spline_deriv_uniform_nd(
        Mat<real> const & aXi,
        uint const & aPolynomialDegree,
        uint const & aDim)
{
    Mat<real> tBsplineDeriv(pow(aPolynomialDegree+1,aDim),aDim,0); // Initialize the size of the vector
    if( aDim == 2)
    {
        if( aPolynomialDegree == 1 )
        {
            tBsplineDeriv(0,0) = -1.0 * (1 - aXi(1));  tBsplineDeriv(0,1) = (1 - aXi(0)) * (-1.0);
            tBsplineDeriv(1,0) = 1.0  * (1 - aXi(1));  tBsplineDeriv(1,1) = aXi(0)       * (-1.0);
            tBsplineDeriv(2,0) = -1.0 * aXi(1);        tBsplineDeriv(2,1) = (1 - aXi(0)) * 1.0;
            tBsplineDeriv(3,0) = 1.0  * aXi(1);        tBsplineDeriv(3,1) = aXi(0)       * 1.0;
        }
        else if( aPolynomialDegree == 2 )
        {
            tBsplineDeriv(0,0) = (aXi(0)-1.0)     * (1.0-aXi(1))*(1.0-aXi(1))*0.5;            tBsplineDeriv(0,1) = ((1.0-aXi(0))*(1.0-aXi(0))*0.5)         * (aXi(1)-1.0)     ;
            tBsplineDeriv(1,0) = (1.0-2.0*aXi(0)) * (1.0-aXi(1))*(1.0-aXi(1))*0.5;            tBsplineDeriv(1,1) = ((-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1.0))) * (aXi(1)-1.0);
            tBsplineDeriv(2,0) = aXi(0)           * (1.0-aXi(1))*(1.0-aXi(1))*0.5;            tBsplineDeriv(2,1) = (0.5*aXi(0)*aXi(0))                     * (aXi(1)-1.0)         ;
            tBsplineDeriv(3,0) = (aXi(0)-1.0)     * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1.0)));  tBsplineDeriv(3,1) = ((1.0-aXi(0))*(1.0-aXi(0))*0.5)         * (1.0-2.0*aXi(1))   ;
            tBsplineDeriv(4,0) = (1.0-2.0*aXi(0)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1.0)));  tBsplineDeriv(4,1) = ((-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1.0))) * (1.0-2.0*aXi(1));
            tBsplineDeriv(5,0) = aXi(0)           * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1.0)));  tBsplineDeriv(5,1) = (0.5*aXi(0)*aXi(0))                     * (1.0-2.0*aXi(1))        ;
            tBsplineDeriv(6,0) = (aXi(0)-1.0)     * (0.5*aXi(1)*aXi(1));                      tBsplineDeriv(6,1) = ((1.0-aXi(0))*(1.0-aXi(0))*0.5)         * aXi(1)    ;
            tBsplineDeriv(7,0) = (1.0-2.0*aXi(0)) * (0.5*aXi(1)*aXi(1));                      tBsplineDeriv(7,1) = ((-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1.0))) * aXi(1) ;
            tBsplineDeriv(8,0) = aXi(0)           * (0.5*aXi(1)*aXi(1));                      tBsplineDeriv(8,1) = (0.5*aXi(0)*aXi(0))                     * aXi(1)        ;
        }
        else if( aPolynomialDegree == 3 )
        {
            tBsplineDeriv(0 ,0) = -0.5*pow(aXi(0)-1.0,2)          * 1.0/6.0*pow(1.0-aXi(1),3);                                                                          tBsplineDeriv(0 ,1) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * (-0.5*pow(aXi(1)-1.0,2));
            tBsplineDeriv(1 ,0) = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))  * 1.0/6.0*pow(1.0-aXi(1),3);                                                                          tBsplineDeriv(1 ,1) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * (-0.5*pow(aXi(1)-1.0,2));
            tBsplineDeriv(2 ,0) = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5) * 1.0/6.0*pow(1.0-aXi(1),3);                                                                          tBsplineDeriv(2 ,1) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * (-0.5*pow(aXi(1)-1.0,2));
            tBsplineDeriv(3 ,0) = aXi(0)*aXi(0)/2.0               * 1.0/6.0*pow(1.0-aXi(1),3);                                                                          tBsplineDeriv(3 ,1) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * (-0.5*pow(aXi(1)-1.0,2));
            tBsplineDeriv(4 ,0) = -0.5*pow(aXi(0)-1.0,2)          * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));   tBsplineDeriv(4 ,1) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * (1.5*aXi(1)*aXi(1)-2.0*aXi(1));
            tBsplineDeriv(5 ,0) = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))  * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));   tBsplineDeriv(5 ,1) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * (1.5*aXi(1)*aXi(1)-2.0*aXi(1));
            tBsplineDeriv(6 ,0) = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5) * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));   tBsplineDeriv(6 ,1) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * (1.5*aXi(1)*aXi(1)-2.0*aXi(1));
            tBsplineDeriv(7 ,0) = aXi(0)*aXi(0)/2.0               * (1.0/6.0*pow(1.0-aXi(1),2)*(aXi(1)+2.0)-1.0/6.0*(2.0-aXi(1))*(2.0*aXi(1)*aXi(1)-2.0*aXi(1)-1.0));   tBsplineDeriv(7 ,1) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * (1.5*aXi(1)*aXi(1)-2.0*aXi(1));
            tBsplineDeriv(8 ,0) = -0.5*pow(aXi(0)-1.0,2)          * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));         tBsplineDeriv(8 ,1) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5);
            tBsplineDeriv(9 ,0) = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))  * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));         tBsplineDeriv(9 ,1) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5);
            tBsplineDeriv(10,0) = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5) * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));         tBsplineDeriv(10,1) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5);
            tBsplineDeriv(11,0) = aXi(0)*aXi(0)/2.0               * (1.0/6.0*(3.0-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1.0)*(2.0*aXi(1)*aXi(1)-2*aXi(1)-1.0));         tBsplineDeriv(11,1) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5);
            tBsplineDeriv(12,0) = -0.5*pow(aXi(0)-1.0,2)          * aXi(1)*aXi(1)*aXi(1)/6.0;                                                                           tBsplineDeriv(12,1) = 1.0/6.0*pow(1.0-aXi(0),3)                                                                        * aXi(1)*aXi(1)/2.0;
            tBsplineDeriv(13,0) = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))  * aXi(1)*aXi(1)*aXi(1)/6.0;                                                                           tBsplineDeriv(13,1) = (1.0/6.0*pow(1.0-aXi(0),2)*(aXi(0)+2.0)-1.0/6.0*(2.0-aXi(0))*(2.0*aXi(0)*aXi(0)-2.0*aXi(0)-1.0)) * aXi(1)*aXi(1)/2.0;
            tBsplineDeriv(14,0) = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5) * aXi(1)*aXi(1)*aXi(1)/6.0;                                                                           tBsplineDeriv(14,1) = (1.0/6.0*(3.0-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1.0)*(2.0*aXi(0)*aXi(0)-2*aXi(0)-1.0))       * aXi(1)*aXi(1)/2.0;
            tBsplineDeriv(15,0) = aXi(0)*aXi(0)/2.0               * aXi(1)*aXi(1)*aXi(1)/6.0;                                                                           tBsplineDeriv(15,1) = aXi(0)*aXi(0)*aXi(0)/6.0                                                                         * aXi(1)*aXi(1)/2.0;
        }
        else
        {
            Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
            for(uint i = 0; i < tKnotVector.length(); i++)
                tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
            Mat<uint> tElement_flag(aDim,1,0);
            tBsplineDeriv = build_spline_deriv_uniform_nd(aXi,tKnotVector,aPolynomialDegree,tElement_flag,aDim);
        }
    }
    else if( aDim == 3)
    {
        if( aPolynomialDegree == 1 )
        {
            tBsplineDeriv(0,0) = -1.0 * (1 - aXi(1)) * (1 - aXi(2));  tBsplineDeriv(0,1) = (1 - aXi(0)) * (-1.0) * (1 - aXi(2));  tBsplineDeriv(0,2) = (1 - aXi(0)) * (1 - aXi(1)) * (-1.0);
            tBsplineDeriv(1,0) = 1.0  * (1 - aXi(1)) * (1 - aXi(2));  tBsplineDeriv(1,1) = aXi(0)       * (-1.0) * (1 - aXi(2));  tBsplineDeriv(1,2) = aXi(0)       * (1 - aXi(1)) * (-1.0);
            tBsplineDeriv(2,0) = -1.0 * aXi(1)       * (1 - aXi(2));  tBsplineDeriv(2,1) = (1 - aXi(0)) * 1.0    * (1 - aXi(2));  tBsplineDeriv(2,2) = (1 - aXi(0)) * aXi(1)       * (-1.0);
            tBsplineDeriv(3,0) = 1.0  * aXi(1)       * (1 - aXi(2));  tBsplineDeriv(3,1) = aXi(0)       * 1.0    * (1 - aXi(2));  tBsplineDeriv(3,2) = aXi(0)       * aXi(1)       * (-1.0);
            tBsplineDeriv(4,0) = -1.0 * (1 - aXi(1)) * aXi(2);        tBsplineDeriv(4,1) = (1 - aXi(0)) * (-1.0) * aXi(2);        tBsplineDeriv(4,2) = (1 - aXi(0)) * (1 - aXi(1)) * 1.0;
            tBsplineDeriv(5,0) = 1.0  * (1 - aXi(1)) * aXi(2);        tBsplineDeriv(5,1) = aXi(0)       * (-1.0) * aXi(2);        tBsplineDeriv(5,2) = aXi(0)       * (1 - aXi(1)) * 1.0;
            tBsplineDeriv(6,0) = -1.0 * aXi(1)       * aXi(2);        tBsplineDeriv(6,1) = (1 - aXi(0)) * 1.0    * aXi(2);        tBsplineDeriv(6,2) = (1 - aXi(0)) * aXi(1)       * 1.0;
            tBsplineDeriv(7,0) = 1.0  * aXi(1)       * aXi(2);        tBsplineDeriv(7,1) = aXi(0)       * 1.0    * aXi(2);        tBsplineDeriv(7,2) = aXi(0)       * aXi(1)       * 1.0;
        }
        else if( aPolynomialDegree == 2 )
        {
            tBsplineDeriv(0,0)  = (aXi(0)-1.0)     * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(0 ,1)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * (aXi(1)-1.0)     * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(0 ,2)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2)-1.0);
            tBsplineDeriv(1,0)  = (1.0-2.0*aXi(0)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(1 ,1)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (aXi(1)-1.0)     * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(1 ,2)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2)-1.0);
            tBsplineDeriv(2,0)  = aXi(0)           * ((1-aXi(1))*(1-aXi(1))*0.5)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(2 ,1)  = (0.5*aXi(0)*aXi(0))                 * (aXi(1)-1.0)     * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(2 ,2)  = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2)-1.0);
            tBsplineDeriv(3,0)  = (aXi(0)-1.0)     * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(3 ,1)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * (1.0-2.0*aXi(1)) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(3 ,2)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2)-1.0);
            tBsplineDeriv(4,0)  = (1.0-2.0*aXi(0)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(4 ,1)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (1.0-2.0*aXi(1)) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(4 ,2)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2)-1.0);
            tBsplineDeriv(5,0)  = aXi(0)           * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(5 ,1)  = (0.5*aXi(0)*aXi(0))                 * (1.0-2.0*aXi(1)) * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(5 ,2)  = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2)-1.0);
            tBsplineDeriv(6,0)  = (aXi(0)-1.0)     * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(6 ,1)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * aXi(1)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(6 ,2)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * (aXi(2)-1.0);
            tBsplineDeriv(7,0)  = (1.0-2.0*aXi(0)) * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(7 ,1)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * aXi(1)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(7 ,2)  = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * (aXi(2)-1.0);
            tBsplineDeriv(8,0)  = aXi(0)           * (0.5*aXi(1)*aXi(1))                   * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(8 ,1)  = (0.5*aXi(0)*aXi(0))                 * aXi(1)           * ((1-aXi(2))*(1-aXi(2))*0.5);          tBsplineDeriv(8 ,2)  = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * (aXi(2)-1.0);
            tBsplineDeriv(9,0)  = (aXi(0)-1.0)     * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(9 ,1)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * (aXi(1)-1.0)     * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(9 ,2)  = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * (1.0-2.0*aXi(2));
            tBsplineDeriv(10,0) = (1.0-2.0*aXi(0)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(10,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * (aXi(1)-1.0)     * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(10,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (1.0-2.0*aXi(2));
            tBsplineDeriv(11,0) = aXi(0)           * ((1-aXi(1))*(1-aXi(1))*0.5)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(11,1) = (0.5*aXi(0)*aXi(0))                  * (aXi(1)-1.0)     * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(11,2) = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * (1.0-2.0*aXi(2));
            tBsplineDeriv(12,0) = (aXi(0)-1.0)     * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(12,1) = ((1-aXi(0))*(1-aXi(0))*0.5)          * (1.0-2.0*aXi(1)) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(12,2) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (1.0-2.0*aXi(2));
            tBsplineDeriv(13,0) = (1.0-2.0*aXi(0)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(13,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * (1.0-2.0*aXi(1)) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(13,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (1.0-2.0*aXi(2));
            tBsplineDeriv(14,0) = aXi(0)           * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(14,1) = (0.5*aXi(0)*aXi(0))                  * (1.0-2.0*aXi(1)) * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(14,2) = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (1.0-2.0*aXi(2));
            tBsplineDeriv(15,0) = (aXi(0)-1.0)     * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(15,1) = ((1-aXi(0))*(1-aXi(0))*0.5)          * aXi(1)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(15,2) = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * (1.0-2.0*aXi(2));
            tBsplineDeriv(16,0) = (1.0-2.0*aXi(0)) * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(16,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * aXi(1)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(16,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * (1.0-2.0*aXi(2));
            tBsplineDeriv(17,0) = aXi(0)           * (0.5*aXi(1)*aXi(1))                   * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));  tBsplineDeriv(17,1) = (0.5*aXi(0)*aXi(0))                  * aXi(1)           * (-0.5*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(17,2) = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * (1.0-2.0*aXi(2));
            tBsplineDeriv(18,0) = (aXi(0)-1.0)     * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(18,1) = ((1-aXi(0))*(1-aXi(0))*0.5)          * (aXi(1)-1.0)     * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(18,2) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2));
            tBsplineDeriv(19,0) = (1.0-2.0*aXi(0)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(19,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * (aXi(1)-1.0)     * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(19,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2));
            tBsplineDeriv(20,0) = aXi(0)           * ((1-aXi(1))*(1-aXi(1))*0.5)           * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(20,1) = (0.5*aXi(0)*aXi(0))                  * (aXi(1)-1.0)     * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(20,2) = (0.5*aXi(0)*aXi(0))                 * ((1-aXi(1))*(1-aXi(1))*0.5)           * (aXi(2));
            tBsplineDeriv(21,0) = (aXi(0)-1.0)     * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(21,1) = ((1-aXi(0))*(1-aXi(0))*0.5)          * (1.0-2.0*aXi(1)) * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(21,2) = ((1-aXi(0))*(1-aXi(0))*0.5)         * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2));
            tBsplineDeriv(22,0) = (1.0-2.0*aXi(0)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(22,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * (1.0-2.0*aXi(1)) * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(22,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2));
            tBsplineDeriv(23,0) = aXi(0)           * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(23,1) = (0.5*aXi(0)*aXi(0))                  * (1.0-2.0*aXi(1)) * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(23,2) = (0.5*aXi(0)*aXi(0))                 * ((-0.5*(2*aXi(1)*aXi(1)-2*aXi(1)-1))) * (aXi(2));
            tBsplineDeriv(24,0) = (aXi(0)-1.0)     * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(24,1) = ((1-aXi(0))*(1-aXi(0))*0.5)          * aXi(1)           * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(24,2) = ((1-aXi(0))*(1-aXi(0))*0.5)         * (0.5*aXi(1)*aXi(1))                   * (aXi(2));
            tBsplineDeriv(25,0) = (1.0-2.0*aXi(0)) * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(25,1) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1))  * aXi(1)           * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(25,2) = (-0.5*(2*aXi(0)*aXi(0)-2*aXi(0)-1)) * (0.5*aXi(1)*aXi(1))                   * (aXi(2));
            tBsplineDeriv(26,0) = aXi(0)           * (0.5*aXi(1)*aXi(1))                   * (0.5*aXi(2)*aXi(2));                  tBsplineDeriv(26,1) = (0.5*aXi(0)*aXi(0))                  * aXi(1)           * (0.5*aXi(2)*aXi(2));                   tBsplineDeriv(26,2) = (0.5*aXi(0)*aXi(0))                 * (0.5*aXi(1)*aXi(1))                   * (aXi(2));
        }
        else if( aPolynomialDegree == 3 )
        {
            tBsplineDeriv(0 ,0)  = -0.5*pow(aXi(0)-1.0,2)           * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(0 ,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-0.5*pow(aXi(1)-1.0,2))         * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(0 ,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(1 ,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(1 ,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-0.5*pow(aXi(1)-1.0,2))         * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(1 ,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(2 ,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(2 ,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-0.5*pow(aXi(1)-1.0,2))         * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(2 ,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(3 ,0)  = aXi(0)*aXi(0)/2.0                * 1.0/6.0*pow(1-aXi(1),3)                                                              * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(3 ,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-0.5*pow(aXi(1)-1.0,2))         * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(3 ,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(4 ,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(4 ,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(4 ,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(5 ,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(5 ,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(5 ,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(6 ,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(6 ,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(6 ,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(7 ,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(7 ,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(7 ,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(8 ,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(8 ,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(8 ,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(9 ,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(9 ,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(9 ,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(10,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(10,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(10,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(11,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(11,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(11,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(12,0)  = -0.5*pow(aXi(0)-1.0,2)           * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(12,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)/2.0                * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(12,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(13,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(13,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)/2.0                * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(13,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(14,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(14,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)/2.0                * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(14,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(15,0)  = aXi(0)*aXi(0)/2.0                * aXi(1)*aXi(1)*aXi(1)/6                                                               * 1.0/6.0*pow(1-aXi(2),3);                                                                 tBsplineDeriv(15,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)/2.0                * 1.0/6.0*pow(1-aXi(2),3);                                                                tBsplineDeriv(15,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-0.5*pow(aXi(2)-1.0,2));
            tBsplineDeriv(16,0)  = -0.5*pow(aXi(0)-1.0,2)           * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(16,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(16,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(17,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(17,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(17,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(18,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(18,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(18,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(19,0)  = aXi(0)*aXi(0)/2.0                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(19,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(19,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(20,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(20,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(20,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(21,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(21,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(21,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(22,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(22,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(22,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(23,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(23,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(23,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(24,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(24,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(24,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(25,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(25,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(25,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(26,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(26,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(26,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(27,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(27,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(27,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(28,0)  = -0.5*pow(aXi(0)-1.0,2)           * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(28,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)/2.0                * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(28,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(29,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(29,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)/2.0                * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(29,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(30,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(30,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)/2.0                * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(30,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(31,0)  = aXi(0)*aXi(0)/2.0                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));    tBsplineDeriv(31,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)/2.0                * (1.0/6.0*pow(1-aXi(2),2)*(aXi(2)+2)-1.0/6.0*(2-aXi(2))*(2*aXi(2)*aXi(2)-2*aXi(2)-1));   tBsplineDeriv(31,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.5*aXi(2)*aXi(2)-2.0*aXi(2));
            tBsplineDeriv(32,0)  = -0.5*pow(aXi(0)-1.0,2)           * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(32,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(32,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(33,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(33,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(33,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(34,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(34,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(34,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(35,0)  = aXi(0)*aXi(0)/2.0                * 1.0/6.0*pow(1-aXi(1),3)                                                              * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(35,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-0.5*pow(aXi(1)-1.0,2))         * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(35,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(36,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(36,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(36,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(37,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(37,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(37,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(38,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(38,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(38,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(39,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(39,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(39,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(40,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(40,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(40,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(41,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(41,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(41,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(42,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(42,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(42,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(43,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(43,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(43,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(44,0)  = -0.5*pow(aXi(0)-1.0,2)           * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(44,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)/2.0                * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(44,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(45,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(45,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)/2.0                * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(45,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(46,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(46,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)/2.0                * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(46,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(47,0)  = aXi(0)*aXi(0)/2.0                * aXi(1)*aXi(1)*aXi(1)/6                                                               * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));      tBsplineDeriv(47,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)/2.0                * (1.0/6.0*(3-aXi(2))*aXi(2)*aXi(2)-1.0/6.0*(aXi(2)+1)*(2*aXi(2)*aXi(2)-2*aXi(2)-1));     tBsplineDeriv(47,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * (-1.5*aXi(2)*aXi(2)+aXi(2)+0.5);
            tBsplineDeriv(48,0)  = -0.5*pow(aXi(0)-1.0,2)           * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(48,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-0.5*pow(aXi(1)-1.0,2))         * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(48,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)/2;
            tBsplineDeriv(49,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(49,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-0.5*pow(aXi(1)-1.0,2))         * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(49,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)/2;
            tBsplineDeriv(50,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(50,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-0.5*pow(aXi(1)-1.0,2))         * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(50,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)/2;
            tBsplineDeriv(51,0)  = aXi(0)*aXi(0)/2.0                * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(51,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-0.5*pow(aXi(1)-1.0,2))         * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(51,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * 1.0/6.0*pow(1-aXi(1),3)                                                              * aXi(2)*aXi(2)/2;
            tBsplineDeriv(52,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(52,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(52,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(53,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(53,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(53,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(54,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(54,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(54,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(55,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(55,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.5*aXi(1)*aXi(1)-2.0*aXi(1))   * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(55,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*pow(1-aXi(1),2)*(aXi(1)+2)-1.0/6.0*(2-aXi(1))*(2*aXi(1)*aXi(1)-2*aXi(1)-1)) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(56,0)  = -0.5*pow(aXi(0)-1.0,2)           * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(56,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(56,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(57,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(57,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(57,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(58,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(58,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(58,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(59,0)  = aXi(0)*aXi(0)/2.0                * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(59,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (-1.5*aXi(1)*aXi(1)+aXi(1)+0.5)  * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(59,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * (1.0/6.0*(3-aXi(1))*aXi(1)*aXi(1)-1.0/6.0*(aXi(1)+1)*(2*aXi(1)*aXi(1)-2*aXi(1)-1)  ) * aXi(2)*aXi(2)/2;
            tBsplineDeriv(60,0)  = -0.5*pow(aXi(0)-1.0,2)           * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(60,1) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)/2.0                * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(60,2) = 1.0/6.0*pow(1-aXi(0),3)                                                                * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)/2;
            tBsplineDeriv(61,0)  = (1.5*aXi(0)*aXi(0)-2.0*aXi(0))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(61,1) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)/2.0                * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(61,2) = (1.0/6.0*pow(1-aXi(0),2)*(aXi(0)+2)-1.0/6.0*(2-aXi(0))*(2*aXi(0)*aXi(0)-2*aXi(0)-1))   * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)/2;
            tBsplineDeriv(62,0)  = (-1.5*aXi(0)*aXi(0)+aXi(0)+0.5)  * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(62,1) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)/2.0                * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(62,2) = (1.0/6.0*(3-aXi(0))*aXi(0)*aXi(0)-1.0/6.0*(aXi(0)+1)*(2*aXi(0)*aXi(0)-2*aXi(0)-1))     * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)/2;
            tBsplineDeriv(63,0)  = aXi(0)*aXi(0)/2.0                * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)*aXi(2)/6;                                                                  tBsplineDeriv(63,1) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)/2.0                * aXi(2)*aXi(2)*aXi(2)/6;                                                                 tBsplineDeriv(63,2) = aXi(0)*aXi(0)*aXi(0)/6                                                                 * aXi(1)*aXi(1)*aXi(1)/6                                                               * aXi(2)*aXi(2)/2;
        }
        else
        {
            Mat<real> tKnotVector(2*(aPolynomialDegree+1),1,0);
            for(uint i = 0; i < tKnotVector.length(); i++)
                tKnotVector(i) = -(real)aPolynomialDegree + (real)i;
            Mat<uint> tElement_flag(aDim,1,0);
            tBsplineDeriv = build_spline_deriv_uniform_nd(aXi,tKnotVector,aPolynomialDegree,tElement_flag,aDim);
        }
    }
    return tBsplineDeriv;
}

Mat<real>
Bspline::build_spline_deriv_uniform_nd(Mat<real> const & aXi,
        Mat<real> & aKnotVector,
        uint const & aPolynomialDegree,
        Mat<uint> const & aElement_flag,
        uint const & aDim)
{
    Mat<real> bspline_uniform_deriv(pow(aPolynomialDegree+1,aDim),aDim,0); // Initialize the size of the vector
    Mat<real> tbspline;
    Mat<real> tbspline1;
    Mat<real> tbspline2;
    Mat<real> tbspline3;
    Mat<real> ttBsplineDeriv;
    Mat<real> ttBsplineDeriv1;
    Mat<real> ttBsplineDeriv2;
    Mat<real> ttBsplineDeriv3;
    uint count = 0; // Counter is needed for a loop

    if(aDim>=2 && aDim<=3)
    {
        tbspline = Bspline::build_spline_1d(aXi(0),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_1
        tbspline1 = tbspline.cols( aElement_flag(0), aElement_flag(0)+aPolynomialDegree);
        ttBsplineDeriv = Bspline::build_spline_deriv_1d(aXi(0),aKnotVector,aPolynomialDegree);  //Derivative of the B-Splines for the direction Xi_1
        ttBsplineDeriv1 = ttBsplineDeriv.rows( aElement_flag(0), aElement_flag(0)+aPolynomialDegree);
        tbspline = Bspline::build_spline_1d(aXi(1),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_2
        tbspline2 = tbspline.cols( aElement_flag(1), aElement_flag(1)+aPolynomialDegree);
        ttBsplineDeriv = Bspline::build_spline_deriv_1d(aXi(1),aKnotVector,aPolynomialDegree);  //Derivative of the B-Splines for the direction Xi_2
        ttBsplineDeriv2 = ttBsplineDeriv.rows( aElement_flag(1), aElement_flag(1)+aPolynomialDegree);
        if(aDim==3)
        {
            tbspline = Bspline::build_spline_1d(aXi(2),aKnotVector,aPolynomialDegree);  //B-Splines for the direction Xi_3
            tbspline3 = tbspline.cols( aElement_flag(2), aElement_flag(2)+aPolynomialDegree);
            ttBsplineDeriv = Bspline::build_spline_deriv_1d(aXi(2),aKnotVector,aPolynomialDegree);  //Derivative of the B-Splines for the direction Xi_2
            ttBsplineDeriv3 = ttBsplineDeriv.rows( aElement_flag(2), aElement_flag(2)+aPolynomialDegree);
        }
    }
    else
    {
        MORIS_LOG_ERROR << "Dimension is out of range 1 < dim < 4 ";
    }

    if(aDim==2)
    {
        for (uint  j = 0; j<aPolynomialDegree+1;j++) //Loop over basis functions in Xi_2 direction
        {
            for (uint  k = 0; k<aPolynomialDegree+1;k++) //Loop over basis functions in Xi_1 direction
            {
                bspline_uniform_deriv(count,0) = ttBsplineDeriv1(k)*tbspline2(j);
                bspline_uniform_deriv(count,1) = tbspline1(k)*ttBsplineDeriv2(j);
                count++;
            }
        }
    }
    else if(aDim==3)
    {
        for (uint  i = 0; i<aPolynomialDegree+1;i++) //Loop over basis functions in Xi_3 direction
        {
            for (uint  j = 0; j<aPolynomialDegree+1;j++) //Loop over basis functions in Xi_2 direction
            {
                for (uint  k = 0; k<aPolynomialDegree+1;k++) //Loop over basis functions in Xi_1 direction
                {
                    bspline_uniform_deriv(count,0) = ttBsplineDeriv1(k)*tbspline2(j)*tbspline3(i);
                    bspline_uniform_deriv(count,1) = tbspline1(k)*ttBsplineDeriv2(j)*tbspline3(i);
                    bspline_uniform_deriv(count,2) = tbspline1(k)*tbspline2(j)*ttBsplineDeriv3(i);
                    count++;
                }
            }
        }
    }
    return bspline_uniform_deriv;
}

