/*
 * fn_Triangle_Geometry.hpp
 *
 *  Created on: Jan 22, 2018
 *      Author: doble
 */

#ifndef SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_
#define SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_


#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_comp_abs.hpp"
#include "tools/fn_tet_volume.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{

/*
 * Compute the surface normal of a triangle (Note this hasn't been fully optimized)
 * @param[in] aTriangleNodes - Nodes of the triangle
 */
template<typename Real_Matrix, typename Integer_Matrix>
void compute_tri_surface_normal( moris::Matrix< Integer_Matrix > const & aTriangleNodes,
                                 moris::Matrix< Real_Matrix > const &       aNodeCoordinates,
                                 moris::Matrix< Real_Matrix > &              aSurfaceNormal,
                                 bool aPosDirection = false)
{

    moris::Matrix< Real_Matrix > tVec1(1,3,0);
    moris::Matrix< Real_Matrix > tVec2(1,3,0);


    typename moris::Matrix< Integer_Matrix >::Data_Type tN1 = aTriangleNodes(0,0);
    typename moris::Matrix< Integer_Matrix >::Data_Type tN2 = aTriangleNodes(0,1);
    typename moris::Matrix< Integer_Matrix >::Data_Type tN3 = aTriangleNodes(0,2);

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
    typename moris::Matrix< Real_Matrix >::Data_Type tLenSquared3 = std::pow( aSurfaceNormal(0,0) , 2) + std::pow( aSurfaceNormal(1,0) , 2) + std::pow( aSurfaceNormal(2,0) , 2);
    typename moris::Matrix< Real_Matrix >::Data_Type tLen3 = std::pow(tLenSquared3,0.5);

    XTK_ASSERT(tLen3>0.000005,"Dividing by near zero value");

    typename moris::Matrix< Real_Matrix >::Data_Type tLenInv = 1/tLen3;

    if(!aPosDirection)
    {
        aSurfaceNormal = tLenInv*aSurfaceNormal.matrix_data();
    }
    else
    {
        aSurfaceNormal = tLenInv*comp_abs(aSurfaceNormal);
    }


}

typename moris::Matrix< Default_Matrix_Real >::Data_Type
compute_volume_for_multiple_tets(moris::Matrix< Default_Matrix_Real > const & aAllNodeCoords,
                                 moris::Matrix< Default_Matrix_Integer > const & aElementToNodeConnectivity)
{
    size_t tNumTets = aElementToNodeConnectivity.n_rows();
    moris::Matrix< Default_Matrix_Real > tCoords(aElementToNodeConnectivity.n_cols(),3);
    moris::Matrix< Default_Matrix_Real > tCoordRow(1,3);

    moris::Matrix< Default_Matrix_Real >::Data_Type tTotalVolume = 0;
    moris::Matrix< Default_Matrix_Real >::Data_Type tSingleVolume;
    size_t k = 0;

    for(size_t i = 0; i<tNumTets; i++)
    {
        for(size_t j = 0; j<aElementToNodeConnectivity.n_cols();j++)
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

template< typename Real_Matrix>
moris::Matrix< Real_Matrix >
compute_val_at_tet_centroid(moris::Matrix< Real_Matrix > const & aNodeVals)
{
    Real_Matrix tVal = 1/4*(aNodeVals.matrix_data());
    moris::Matrix< Real_Matrix > tValMat(tVal);

    return tValMat;
}

}


#endif /* SRC_GEOMENG_FN_TRIANGLE_GEOMETRY_HPP_ */
