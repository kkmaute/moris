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
}

