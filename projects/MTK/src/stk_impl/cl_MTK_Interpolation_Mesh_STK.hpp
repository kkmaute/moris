/*
 * cl_MTK_Interpolation_Mesh_STK.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
namespace moris
{
namespace mtk
{
class Interpolation_Mesh_STK : public Mesh_Core_STK, public Interpolation_Mesh
{
    // Functions only valid for interpolation mIntegrationMeshes

public:
    Interpolation_Mesh_STK(std::shared_ptr<Mesh_Data_STK> aSTKMeshData):
        Mesh_Core_STK(aSTKMeshData)
    {

    }

    Interpolation_Mesh_STK(
            std::string    aFileName,
            MtkMeshData*   aSuppMeshData,
            const bool     aCreateFacesAndEdges = true ):
                Mesh_Core_STK(aFileName,aSuppMeshData,aCreateFacesAndEdges)

    {

    }

    Interpolation_Mesh_STK(MtkMeshData & aMeshData ):
                Mesh_Core_STK(aMeshData)

    {

    }


    std::shared_ptr<Mesh_Data_STK>
    get_stk_data_shared_pointer()
    {
        return mSTKMeshData;
    }

//    /*
//     * Get elements interpolated into by a basis function. For a Lagrange mesh,
//     * the elements in support of basis is equivalent to the elements connected
//     * to a node. Therefore, a call to get_elements
//     */
//    Matrix< IndexMat >
//    get_elements_in_support_of_basis(moris_index aBasisIndex);
//
//    //FIXME: IMPLEMENT THIS FUNCTION IN STK,XTK
//    /*
//     * Get number of B-Spline coefficients
//     */
//    uint
//    get_num_coeffs(const uint aOrder) const;
//
//    //------------------------------------------------------------------------------
//    const Matrix< DDRMat > &
//    get_t_matrix_of_node_loc_ind(
//            const moris_index aNodeIndex,
//            const EntityRank  aBSplineRank );
//
//    Matrix< IndexMat >
//    get_bspline_inds_of_node_loc_ind(
//            const moris_index aNodeIndex,
//            const EntityRank  aBSplineRank );
//
//    /*
//     * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
//     * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
//     */
//    uint
//    get_num_basis_functions();
//
//    //FIXME: Rename or use get loc entity id from global entity id
//    void
//    get_adof_map( const uint aOrder,
//                  map< moris_id, moris_index > & aAdofMap ) const;
//
//    /**
//     * return the interpolation order of this field
//     */
//    uint
//    get_order_of_field(
//            const moris_index     aFieldIndex,
//            const enum EntityRank aEntityRank );

};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_ */
