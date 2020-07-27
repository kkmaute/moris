#ifndef MORIS_CL_GEN_SPHERE_BOX_HPP
#define MORIS_CL_GEN_SPHERE_BOX_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
namespace ge
{
class Sphere_Box : public Geometry, public Field_Analytic
{
public:
    Sphere_Box( moris::real aSX,
                moris::real aSY,
                moris::real aSZ,
                moris::real aXCenter,
                moris::real aYCenter,
                moris::real aZCenter,
                moris::real aNexp)
                : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aNexp}})),
                  Geometry(0, -1, -1, -1.0, 1.0)
    {
        mSX = aSX;
        mSY = aSY;
        mSZ = aSZ;
        mXCenter = aXCenter;
        mYCenter = aYCenter;
        mZCenter = aZCenter;
        mNexp = aNexp;
    }

    //------------------------------------ these will be deleted/modified later ----------------------------------------
    virtual real evaluate_field_value(const moris::Matrix<moris::DDRMat> &aCoordinates)
    {
        return this->evaluate_field_value_with_coordinate(0,aCoordinates);
    }

    real
    operator()(const moris::Matrix<moris::DDRMat> & aCoordinates)
    {
        return this->evaluate_field_value_with_coordinate(0,aCoordinates);
    }

    void evaluate_all_sensitivities(const moris::Matrix<moris::DDRMat> &aCoordinates, Matrix<DDRMat>& aSensitivities)
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

#endif /* MORIS_CL_GEN_SPHERE_BOX_HPP */
