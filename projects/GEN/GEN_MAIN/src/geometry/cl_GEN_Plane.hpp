/*
 * cl_GEN_Plane.hpp
 *
 *  Created on: Feb 1, 2018
 *      Author: ktdoble
 */
#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_PLANE_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_PLANE_HPP_

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace ge
{
template< moris::uint SpatialDim >
class Plane : public Geometry_Analytic
{
public:
    Plane(Matrix<moris::DDRMat> const & aCenters,
          Matrix<moris::DDRMat> const & aNormals):
            Geometry_Analytic(Matrix<DDRMat>(1, 1, 0.0))
    {
        mCenters = aCenters;
        mNormals = aNormals;
    }

    //------------------------------------ these will be deleted/modified later ----------------------------------------
    virtual real evaluate_field_value(const moris::Matrix<moris::DDRMat> &aCoordinates)
    {
        return 0;
    }

    virtual moris::Matrix<moris::DDRMat> evaluate_sensitivity(const moris::Matrix<moris::DDRMat> &aCoordinates)
    {
        return Matrix<DDRMat>(1, 1, 0.0);
    }






    bool is_analytic() const
    {
        return true;
    }

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 6;
        aNumCols = 3;
    }


    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                                     moris::Matrix< moris::DDRMat > const & aCoordinates) const;


private:
    Matrix<moris::DDRMat> mCenters;
    Matrix<moris::DDRMat> mNormals;
};

template<>
moris::real Plane<3>::evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                                 moris::Matrix< moris::DDRMat > const & aCoordinates) const
{
    moris::real tDist = mNormals(0)*(aCoordinates(aRowIndex,0)-mCenters(0)) + mNormals(1)*(aCoordinates(aRowIndex,1)-mCenters(1)) + mNormals(2)*(aCoordinates(aRowIndex,2)-mCenters(2));
    return tDist;
}
template<>
moris::real Plane<2>::evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                                 moris::Matrix< moris::DDRMat > const & aCoordinates) const
{
    moris::real tDist = mNormals(0)*(aCoordinates(aRowIndex,0)-mCenters(0)) + mNormals(1)*(aCoordinates(aRowIndex,1)-mCenters(1));
    return tDist;
}

}
}


#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_PLANE_HPP_ */
