/*
 * UT_AD_Scratch.cpp
 *
 *  Created on: Aug 22, 2019
 *      Author: doble
 */


#include "catch.hpp"

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>
#include "cl_SphereBox_AD.hpp"
#include "cl_Multi_Geometry_KS_AD.hpp"
#include "cl_Geometry.hpp"
#include "cl_Geometry_AD.hpp"

template <typename ScalarType>
using ScalarArray3DT = typename Kokkos::View<ScalarType***>;


struct SimplexFadTypes
{

    using StateFad   = Sacado::Fad::SFad<moris::real, 12>;
    using ControlFad = Sacado::Fad::SFad<moris::real, 4>;
    using ConfigFad  = Sacado::Fad::SFad<moris::real, 12>;
};


template<moris::moris_index SpaceDim>
class ComputeSurfaceArea
{
private:
public:
    ComputeSurfaceArea(){}



    template<typename ScalarType>
    inline void
    operator()(moris::moris_index         aCellOrdinal,
               int*                       aCellLocalNodeOrdinals,
               ScalarArray3DT<ScalarType> config,
               ScalarType& sideArea) const;
};


template<>
template<typename ScalarType>
inline void
ComputeSurfaceArea<3>::operator()(moris::moris_index aCellOrdinal,
                                  int*                aCellLocalNodeOrdinals,
                                  ScalarArray3DT<ScalarType> config,
                                  ScalarType & sideArea) const
{
    ScalarType ab0 = config(aCellOrdinal,aCellLocalNodeOrdinals[0],0) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],0);
    ScalarType ab1 = config(aCellOrdinal,aCellLocalNodeOrdinals[0],1) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],1);
    ScalarType ab2 = config(aCellOrdinal,aCellLocalNodeOrdinals[0],2) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],2);

    ScalarType bc0 = config(aCellOrdinal,aCellLocalNodeOrdinals[1],0) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],0);
    ScalarType bc1 = config(aCellOrdinal,aCellLocalNodeOrdinals[1],1) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],1);
    ScalarType bc2 = config(aCellOrdinal,aCellLocalNodeOrdinals[1],2) - config(aCellOrdinal,aCellLocalNodeOrdinals[2],2);


    ScalarType Cross0 = ab1 * bc2 - ab2 * bc1;
    ScalarType Cross1 = ab2 * bc0 - ab0 * bc2;
    ScalarType Cross2 = ab0 * bc1 - ab1 * bc0;

    std::cout<<"Cross0 = "<<Cross0<<std::endl;
    std::cout<<"Cross1 = "<<Cross1<<std::endl;
    std::cout<<"Cross2 = "<<Cross2<<std::endl;


    sideArea = sqrt(Cross0*Cross0 + Cross1*Cross1 + Cross2*Cross2)/2;
}

namespace moris
{

TEST_CASE("Automatic Differentiation","[AD_SURF]")
{
    // nodes on side
    int tNodeInds[3];
    tNodeInds[0] = 3;
    tNodeInds[1] = 1;
    tNodeInds[2] = 0;

    // Compute surface area with automatic differentiation (ordered by cell)
    ScalarArray3DT<Sacado::Fad::SFad<double, 12>> tADNodeCoordinates("no ad coords", 1, 4, 3);
    tADNodeCoordinates(0,0,0) =  0.15;     tADNodeCoordinates(0,0,1) =  0.1;     tADNodeCoordinates(0,0,2) =  0.05;
    tADNodeCoordinates(0,1,0) =  0.5;     tADNodeCoordinates(0,1,1) =  -0.3;     tADNodeCoordinates(0,1,2) =  -0.2;
    tADNodeCoordinates(0,2,0) =  -0.4;     tADNodeCoordinates(0,2,1) =  0.9;     tADNodeCoordinates(0,2,2) =  0.2;
    tADNodeCoordinates(0,3,0) =  0.3;     tADNodeCoordinates(0,3,1) =  -0.3;     tADNodeCoordinates(0,3,2) =  1.4;

    tADNodeCoordinates(0,0,0).diff(0,12);
    tADNodeCoordinates(0,0,1).diff(1,12);
    tADNodeCoordinates(0,0,2).diff(2,12);
    tADNodeCoordinates(0,1,0).diff(3,12);
    tADNodeCoordinates(0,1,1).diff(4,12);
    tADNodeCoordinates(0,1,2).diff(5,12);
    tADNodeCoordinates(0,2,0).diff(6,12);
    tADNodeCoordinates(0,2,1).diff(7,12);
    tADNodeCoordinates(0,2,2).diff(8,12);
    tADNodeCoordinates(0,3,0).diff(9,12);
    tADNodeCoordinates(0,3,1).diff(10,12);
    tADNodeCoordinates(0,3,2).diff(11,12);

    Sacado::Fad::SFad<double, 12> tSurfaceAreaAD = 0.0;
    ComputeSurfaceArea<3> tSurfacAreaFunc;
    tSurfacAreaFunc(0,tNodeInds,tADNodeCoordinates,tSurfaceAreaAD);

    auto tDx0 = tSurfaceAreaAD.dx(0);


    // perturb node 0 in the x direction
    moris::real aLength = -0.001;
    tADNodeCoordinates(0,1,0) =  aLength + tADNodeCoordinates(0,1,0).val() ;


    Sacado::Fad::SFad<double, 12> tSurfaceAreaPerturbed = 0.0;
    tSurfacAreaFunc(0,tNodeInds,tADNodeCoordinates,tSurfaceAreaPerturbed);
    std::cout<<"tSurfaceArea  AutoDifferentiation : "<<tSurfaceAreaAD<<std::endl;
    std::cout<<"tSurfaceArea  Perturbed           : "<<tSurfaceAreaPerturbed<<std::endl;
    std::cout<<"Expected                          : "<<(tSurfaceAreaPerturbed - tSurfaceAreaAD )/ (aLength)<<std::endl;

}


template<moris::moris_index SpaceDim>
class ADFunc
{
private:

public:
    ADFunc(){}

    template<typename ScalarType>
    inline void
    operator()(ScalarType aRad,
               ScalarType aXc,
               ScalarType aYc,
               ScalarType& aPhi) const;
};


template<>
template<typename ScalarType>
inline void
ADFunc<2>::operator()(ScalarType aRad,
                        ScalarType aXc,
                        ScalarType aYc,
                        ScalarType& aPhi) const
{
    aPhi = aRad*aRad - aXc*aXc - aYc*aYc*aYc;
}


TEST_CASE("Automatic Differentiation Circle","[AD_Circle]")
{
    // parameters
    Sacado::Fad::SFad<double, 3> tRadius = 1.0;
    Sacado::Fad::SFad<double, 3> tXc     = 1.0;
    Sacado::Fad::SFad<double, 3> tYc     = 1.0;

    // declare independent variables
    tRadius.diff(0,3);
    tXc.diff(1,3);
    tYc.diff(2,3);

    ADFunc<2> tCircleAd;

    Sacado::Fad::SFad<double, 3> tPhi = 0.0;
    tCircleAd(tRadius,tXc,tYc,tPhi);

    std::cout<<"tPhi = "<<tPhi<<std::endl;


}

TEST_CASE("AD Sphere Box","[AD_Spherebox]")
{
    // setup design variables
    Sacado::Fad::DFad<moris::real> Width  = 6.00;
    Sacado::Fad::DFad<moris::real> Depth  = 6.00;
    Sacado::Fad::DFad<moris::real> Height = 6.00;
    Sacado::Fad::DFad<moris::real> H_thk  = 1.10;
    Sacado::Fad::DFad<moris::real> V_thk  = 1.10;
    Sacado::Fad::DFad<moris::real> slot_w = 0.75;

    // declare as independent variables
    Width.diff(0,5);
    Depth.diff(1,5);
    Height.diff(2,5);
    H_thk.diff(3,6);
    V_thk.diff(4,6);
    slot_w.diff(5,6);

    std::cout<<"Number of derivatives = "<<Width.size()<<std::endl;

    std::cout<<"Design Variables:"<<std::endl;
    std::cout<<"\nWidth  = " << Width <<std::endl;
    std::cout<<"Depth  = " << Depth <<std::endl;
    std::cout<<"Height = " << Height<<std::endl;
    std::cout<<"H_thk  = " << H_thk <<std::endl;
    std::cout<<"V_thk  = " << V_thk <<std::endl;
    std::cout<<"slot_w = " << slot_w<<std::endl;

    // Set parameters which do not have derivatives
    // bottom  back left corner off set
    moris::real tXOff = 2.0;
    moris::real tYOff = 2.0;
    moris::real tZOff = 4.02;

    // box sphere exponent
    moris::real tNExp = 10;

    // Bottom Plate
    Sacado::Fad::DFad<moris::real>  tSx = Depth/2;
    Sacado::Fad::DFad<moris::real>  tSy = Width/2;
    Sacado::Fad::DFad<moris::real>  tSz = V_thk/2;
    Sacado::Fad::DFad<moris::real>  tCx = Depth/2 + tXOff;
    Sacado::Fad::DFad<moris::real>  tCy = Width/2 + tYOff;
    Sacado::Fad::DFad<moris::real>  tCz = V_thk/2 + tZOff;

    xtk::Sphere_Box_AD tBottomPlate(tSx,tSy,tSz,tCx,tCy,tCz,tNExp);
    tBottomPlate.print();


    // right plate (pos y face)
    Sacado::Fad::DFad<moris::real> tSx2 = Depth/2;
    Sacado::Fad::DFad<moris::real> tSy2 = H_thk/2;
    Sacado::Fad::DFad<moris::real> tSz2 = Height/2-V_thk*0.75;
    Sacado::Fad::DFad<moris::real> tCx2 = Depth/2 + tXOff;
    Sacado::Fad::DFad<moris::real> tCy2 = Width - H_thk/2 + tYOff;
    Sacado::Fad::DFad<moris::real> tCz2 = Height/2 + tZOff;
    xtk::Sphere_Box_AD tRightPlate(tSx2,tSy2,tSz2,tCx2,tCy2,tCz2,tNExp);
    tRightPlate.print();

    moris::Matrix<moris::DDRMat> tCoords {{ 5, 3, 6.2}, {5,5,4.57}};

    Sacado::Fad::DFad<moris::real>  tPhiA;
    tBottomPlate.ad_sensitivity_evaluation(0,tCoords,tPhiA);

    std::cout<<"tPhiA = "<<tPhiA<<std::endl;

    Sacado::Fad::DFad<moris::real>  tPhiB;
    tBottomPlate.ad_sensitivity_evaluation(1,tCoords,tPhiB);
    std::cout<<"tPhiB = "<<tPhiB<<std::endl;


    Sacado::Fad::DFad<moris::real> xsi = (tPhiA + tPhiB) / (tPhiA-tPhiB);
    std::cout<<"xsi = "<<xsi<<std::endl;
    Sacado::Fad::DFad<moris::real> tN1 = 0.5*(1-xsi);
    Sacado::Fad::DFad<moris::real> tN2 = 0.5*(1+xsi);

    std::cout<<"tXGamma1 = "<<tN1<<std::endl;
    std::cout<<"tXGamma2 = "<<tN2<<std::endl;

    Sacado::Fad::DFad<moris::real> tXGamma1 = tN1*tCoords(0,0) + tN2*tCoords(1,0);
    Sacado::Fad::DFad<moris::real> tXGamma2 = tN1*tCoords(0,1) + tN2*tCoords(1,1);
    Sacado::Fad::DFad<moris::real> tXGamma3 = tN1*tCoords(0,2) + tN2*tCoords(1,2);

    std::cout<<"tXGamma1 = "<<tXGamma1<<std::endl;
    std::cout<<"tXGamma2 = "<<tXGamma2<<std::endl;
    std::cout<<"tXGamma3 = "<<tXGamma3<<std::endl;


    std::cout<<"\n----------------------------------------------"<<std::endl;
    std::cout<<"KS - Geometry 2 Plates:"<<std::endl;
    // create a KS multigeometry
    moris::real tBeta = 50;
    moris::Cell<xtk::Geometry_AD*> tGeomVect = {&tBottomPlate,&tRightPlate};

    xtk::Multi_Geometry_KS_AD tMultiGeometryAD(tGeomVect,tBeta,6);

    tMultiGeometryAD.ad_sensitivity_evaluation(0,tCoords,tPhiA);

    std::cout<<"tPhiA = "<<tPhiA<<std::endl;

    tMultiGeometryAD.ad_sensitivity_evaluation(1,tCoords,tPhiB);
    std::cout<<"tPhiB = "<<tPhiB<<std::endl;

    moris::Matrix<moris::DDRMat> tDPhiADp = tMultiGeometryAD.evaluate_sensitivity_dphi_dp_with_coordinate(0,tCoords);

    print(tDPhiADp,"tDPhiADp");

    xsi = (tPhiA + tPhiB) / (tPhiA-tPhiB);
    std::cout<<"xsi = "<<xsi<<std::endl;
    tN1 = 0.5*(1-xsi);
    tN2 = 0.5*(1+xsi);

    std::cout<<"tN1 = "<<tN1<<std::endl;
    std::cout<<"tN2 = "<<tN2<<std::endl;

    tXGamma1 = tN1*tCoords(0,0) + tN2*tCoords(1,0);
    tXGamma2 = tN1*tCoords(0,1) + tN2*tCoords(1,1);
    tXGamma3 = tN1*tCoords(0,2) + tN2*tCoords(1,2);

    std::cout<<"tXGamma1 = "<<tXGamma1<<std::endl;
    std::cout<<"tXGamma2 = "<<tXGamma2<<std::endl;
    std::cout<<"tXGamma3 = "<<tXGamma3<<std::endl;


}

}

