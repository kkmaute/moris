/*
 * cl_MTK_Interpolation_Mesh.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_

#include "cl_MTK_Mesh_Core.hpp"

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace mtk
{
class Interpolation_Mesh: public virtual Mesh_Core
{
    // Functions only valid for interpolation mIntegrationMeshes


    /*
     * Get elements interpolated into by a basis function. For a Lagrange mesh,
     * the elements in support of basis is equivalent to the elements connected
     * to a node. Therefore, a call to get_elements
     */
    virtual
    Matrix< IndexMat >
    get_elements_in_support_of_basis(moris_index aBasisIndex)
    {
        MORIS_ERROR( false, "get_elements_in_support_of_basis() not implemented for this mesh" );
        return Matrix<IndexMat>(0,0);
    }

    //FIXME: IMPLEMENT THIS FUNCTION IN STK,XTK
    /*
     * Get number of B-Spline coefficients
     */
    virtual uint
    get_num_coeffs(const uint aOrder) const
    {
        MORIS_ERROR( false, "get_num_coeffs() not implemented for this mesh" );
        return 0;
    }

    //------------------------------------------------------------------------------
    //FIXME: IMPLEMENT THIS FUNCTION IN STK
    virtual const Matrix< DDRMat > &
    get_t_matrix_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank )
            {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return mDummyMatrix;
            }

    //FIXME: IMPLEMENT THIS FUNCTION IN XTK
    virtual Matrix< IndexMat >
    get_bspline_inds_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank )
            {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return Matrix<IndexMat>(0,0);
            }

    /*
     * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
     * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
     */
    virtual
    uint
    get_num_basis_functions()
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (get_num_basis_functions function is not implemented)");
        return 0;
    }

    //FIXME: Rename or use get loc entity id from global entity id
    void
    virtual
    get_adof_map( const uint aOrder,
                  map< moris_id, moris_index > & aAdofMap ) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
    }

    /**
     * return the interpolation order of this field
     */
    virtual uint
    get_order_of_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank )
    {
        MORIS_ERROR( false ,"get_order_of_field() not implemented" );
        return 0;
    }

    Matrix< DDRMat > mDummyMatrix;
};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_ */
