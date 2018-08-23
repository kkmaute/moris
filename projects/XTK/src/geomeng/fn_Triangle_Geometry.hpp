/*
 * fn_Triangle_Geometry.hpp
 *
 *  Created on: Jan 22, 2018
 *      Author: doble
 */

#ifndef SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_
#define SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_


#include "linalg/cl_XTK_Matrix.hpp"
#include "tools/fn_tet_volume.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{

/*
 * Compute the surface normal of a triangle (Note this hasn't been fully optimized)
 * @param[in] aTriangleNodes - Nodes of the triangle
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
void compute_tri_surface_normal( Mat<Integer, Integer_Matrix> const & aTriangleNodes,
                                 Mat<Real, Real_Matrix> const &       aNodeCoordinates,
                                 Mat<Real,Real_Matrix> &              aSurfaceNormal,
                                 bool aPosDirection = false)
{

    Mat<Real,Real_Matrix> tVec1(1,3,0);
    Mat<Real,Real_Matrix> tVec2(1,3,0);


    Integer tN1 = aTriangleNodes(0,0);
    Integer tN2 = aTriangleNodes(0,1);
    Integer tN3 = aTriangleNodes(0,2);

    tVec1(0,0) = aNodeCoordinates(tN3,0) - aNodeCoordinates(tN1,0);
    tVec1(0,1) = aNodeCoordinates(tN3,1) - aNodeCoordinates(tN1,1);
    tVec1(0,2) = aNodeCoordinates(tN3,2) - aNodeCoordinates(tN1,2);


    tVec2(0,0) = aNodeCoordinates(tN2,0) - aNodeCoordinates(tN1,0);
    tVec2(0,1) = aNodeCoordinates(tN2,1) - aNodeCoordinates(tN1,1);
    tVec2(0,2) = aNodeCoordinates(tN2,2) - aNodeCoordinates(tN1,2);


    // Compute cross products
    aSurfaceNormal(0,0) = tVec1(0,1) * tVec2(0,2) - tVec1(0,2) * tVec2(0,1);
    aSurfaceNormal(1,0) = tVec1(0,2) * tVec2(0,0) - tVec1(0,0) * tVec2(0,2);
    aSurfaceNormal(2,0) = tVec1(0,0) * tVec2(0,1) - tVec1(0,1) * tVec2(0,0);

    // Normalized
    Real tLenSquared3 = std::pow( aSurfaceNormal(0,0) , 2) + std::pow( aSurfaceNormal(1,0) , 2) + std::pow( aSurfaceNormal(2,0) , 2);
    Real tLen3 = std::pow(tLenSquared3,0.5);

    XTK_ASSERT(tLen3>0.000005,"Dividing by near zero value");

    Real tLenInv = 1/tLen3;

    if(!aPosDirection)
    {
        aSurfaceNormal.matrix_data() = tLenInv*aSurfaceNormal.matrix_data();
    }
    else
    {
        aSurfaceNormal.matrix_data() = tLenInv*component_wise_abs(aSurfaceNormal.matrix_data());
    }


}

real
compute_volume_for_multiple_tets(Mat<real, Default_Matrix_Real> const & aAllNodeCoords,
                                 Mat<size_t, Default_Matrix_Integer> const & aElementToNodeConnectivity)
{
    size_t tNumTets = aElementToNodeConnectivity.get_num_rows();
    Mat<real, Default_Matrix_Real> tCoords(aElementToNodeConnectivity.get_num_columns(),3);
    Mat<real, Default_Matrix_Real> tCoordRow(1,3);

    real tTotalVolume = 0;
    real tSingleVolume;
    size_t k = 0;

    for(size_t i = 0; i<tNumTets; i++)
    {
        for(size_t j = 0; j<aElementToNodeConnectivity.get_num_columns();j++)
        {
            tCoordRow = aAllNodeCoords.get_row(aElementToNodeConnectivity(i,j));
            tCoords.set_row(k,tCoordRow);
            k++;
        }
        k=0;
        tSingleVolume=xtk::vol_tetrahedron(tCoords);

        tTotalVolume+=tSingleVolume;
    }

    return tTotalVolume;
}

template< typename Real, typename Real_Matrix>
Mat<Real,Real_Matrix>
compute_val_at_tet_centroid(Mat<Real,Real_Matrix> const & aNodeVals)
{
    Real_Matrix tVal = 1/4*(aNodeVals.matrix_data());
    Mat<Real,Real_Matrix> tValMat(tVal);

    return tValMat;
}

}


#endif /* SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_ */
