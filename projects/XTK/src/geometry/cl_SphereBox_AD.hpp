/*
 * cl_SphereBox.hpp
 *
 *  Created on: Aug 8, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_SPHEREBOX_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_SPHEREBOX_HPP_


#include "cl_Matrix.hpp"
#include "cl_Geometry_AD.hpp"
#include <Sacado.hpp>

namespace xtk
{
class Sphere_Box_AD : public Geometry_AD
{
public:
    Sphere_Box_AD(){}

    Sphere_Box_AD( Sacado::Fad::DFad<moris::real>  aSX,
                   Sacado::Fad::DFad<moris::real>  aSY,
                   Sacado::Fad::DFad<moris::real>  aSZ,
                   Sacado::Fad::DFad<moris::real>  aXCenter,
                   Sacado::Fad::DFad<moris::real>  aYCenter,
                   Sacado::Fad::DFad<moris::real>  aZCenter,
                   moris::real  aNexp) :
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
        moris::real dx = aCoordinates(aRowIndex,0) - mXCenter.val();
        moris::real dy = aCoordinates(aRowIndex,1) - mYCenter.val();
        moris::real dz = aCoordinates(aRowIndex,2) - mZCenter.val();

                moris::real avgRad = (mSX.val() + mSY.val() + mSZ.val())/((mSX.val() !=0.0) + (mSY.val() !=0.0) + (mSZ.val() !=0.0)); // only count radii that are nonzero
//        moris::real avgRad = 1;


        moris::real  LsValue = avgRad - std::pow(   std::pow( dx/mSX.val()*avgRad, mNexp )
                                      + std::pow( dy/mSY.val()*avgRad, mNexp )
                                      + std::pow( dz/mSZ.val()*avgRad, mNexp ), 1.0/mNexp );



        return LsValue;
    }


    moris::Matrix< moris::DDRMat >
    evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                  moris::Matrix< moris::DDRMat > const & aCoordinates) const
            {

        Sacado::Fad::DFad<moris::real> tLsVal;
        this->ad_sensitivity_evaluation(aRowIndex,aCoordinates,tLsVal);
        return moris::Matrix< moris::DDRMat >(0,0);
            }


    void
    ad_sensitivity_evaluation(moris::size_t const & aRowIndex,
                              moris::Matrix< moris::DDRMat > const & aCoordinates,
                              Sacado::Fad::DFad<moris::real> & aLsVal) const
    {


        Sacado::Fad::DFad<moris::real>  dx = aCoordinates(aRowIndex,0) - mXCenter;
        Sacado::Fad::DFad<moris::real>  dy = aCoordinates(aRowIndex,1) - mYCenter;
        Sacado::Fad::DFad<moris::real>  dz = aCoordinates(aRowIndex,2) - mZCenter;

        //        moris::real avgRad = (mSX + mSY + mSZ)/((mSX !=0.0) + (mSY !=0.0) + (mSZ !=0.0)); // only count radii that are nonzero
        moris::real avgRad = 1;


        aLsVal = avgRad - std::pow(   pow( dx/mSX*avgRad, mNexp )
        + std::pow( dy/mSY*avgRad, mNexp )
        + std::pow( dz/mSZ*avgRad, mNexp ), 1.0/mNexp );

    }



    void
    print()
    {
        std::cout<<"--------------------------------------------------------"<<std::endl;
        std::cout<<"X Center   : "<<std::setw(10)<<mXCenter<<std::endl;
        std::cout<<"Y Center   : "<<std::setw(10)<<mYCenter<<std::endl;
        std::cout<<"Z Center   : "<<std::setw(10)<<mZCenter<<std::endl;
        std::cout<<"X Semi Axis: "<<std::setw(10)<<mSX<<std::endl;
        std::cout<<"Y Semi Axis: "<<std::setw(10)<<mSY<<std::endl;
        std::cout<<"Z Semi Axis: "<<std::setw(10)<<mSZ<<std::endl;
        std::cout<<"Nexp       : "<<std::setw(10)<<mNexp<<std::endl;
        std::cout<<"--------------------------------------------------------"<<std::endl;
    }

public:
    Sacado::Fad::DFad<moris::real> mSX;
    Sacado::Fad::DFad<moris::real> mSY;
    Sacado::Fad::DFad<moris::real> mSZ;
    Sacado::Fad::DFad<moris::real> mXCenter;
    Sacado::Fad::DFad<moris::real> mYCenter;
    Sacado::Fad::DFad<moris::real> mZCenter;
    moris::real mNexp;
};
}



#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_SPHEREBOX_HPP_ */
