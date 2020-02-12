/*
 * cl_SphereBox.hpp
 *
 *  Created on: Aug 8, 2019
 *      Author: doble
 */
#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_BOX_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_BOX_HPP_


#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry.hpp"

namespace moris
{
namespace ge
{
class Sphere_Box : public GEN_Geometry
{
public:
    Sphere_Box(){}

    Sphere_Box( moris::real aSX,
                moris::real aSY,
                moris::real aSZ,
                moris::real aXCenter,
                moris::real aYCenter,
                moris::real aZCenter,
                moris::real aNexp) :
                mSX(aSX),
                mSY(aSY),
                mSZ(aSZ),
                mXCenter(aXCenter),
                mYCenter(aYCenter),
                mZCenter(aZCenter),
                mNexp(aNexp)
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                                     moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {

        // Find distance from point to center of inclusion
        moris::real dx = aCoordinates(aRowIndex,0) - mXCenter;
        moris::real dy = aCoordinates(aRowIndex,1) - mYCenter;
        moris::real dz = aCoordinates(aRowIndex,2) - mZCenter;

//        moris::real avgRad = (mSX + mSY + mSZ)/((mSX !=0.0) + (mSY !=0.0) + (mSZ !=0.0)); // only count radii that are nonzero
        moris::real avgRad = 1;


        moris::real LsValue = avgRad - std::pow(   std::pow( dx/mSX*avgRad, mNexp )
                                                 + std::pow( dy/mSY*avgRad, mNexp )
                                                 + std::pow( dz/mSZ*avgRad, mNexp ), 1.0/mNexp );

        //                    std::printf(" ... phi: %+1.6f\n", LsValue);

        // Makes diamonds
        // @link http://en.wikipedia.org/wiki/Square#Coordinates_and_equations
        // @link http://polymathprogrammer.com/2010/03/01/answered-can-you-describe-a-square-with-1-equation/
        // moris::real LsValue = radius - ( std::abs( dx ) + std::abs( dy ) + std::abs( dz ) );


        return LsValue;
    }


    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                                                 moris::Matrix< moris::DDRMat > const & aCoordinates) const
            {
        return moris::Matrix< moris::DDRMat >(0,0);
            }

    void
    print()
    {
        std::cout<<std::setw(9)<<std::left<<"X Center"<<" = "<<std::setw(10)<<mXCenter<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"Y Center"<<" = "<<std::setw(10)<<mYCenter<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"Z Center"<<" = "<<std::setw(10)<<mZCenter<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"X Semi Axis"<<" = "<<std::setw(10)<<mSX<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"Y Semi Axis"<<" = "<<std::setw(10)<<mSY<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"Z Semi Axis"<<" = "<<std::setw(10)<<mSZ<<std::endl;
        std::cout<<std::setw(9)<<std::left<<"Nexp"<<" = "<<std::setw(10)<<mNexp<<std::endl;

    }

private:
    moris::real mSX;
    moris::real mSY;
    moris::real mSZ;
    moris::real mXCenter;
    moris::real mYCenter;
    moris::real mZCenter;
    moris::real mNexp;
};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_BOX_HPP_ */
