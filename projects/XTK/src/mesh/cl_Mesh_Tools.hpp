/*
 * cl_Mesh_Tools.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_TOOLS_HPP_
#define SRC_MESH_CL_MESH_TOOLS_HPP_


#include"tools/fn_tet_volume.hpp"
#include"linalg/cl_XTK_Matrix.hpp"




using namespace moris;

namespace mesh
{
class Mesh_Helper
{
public:
//    template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
//    static moris::Matrix< Integer_Matrix >
//    get_glb_entity_id_from_entity_loc_index_range(Mesh_Data<Real,Integer, Real_Matrix, Integer_Matrix> const & aMeshData,
//                                                  moris::Matrix< Integer_Matrix > const & tEntityIndices,
//                                                  enum EntityRank aEntityRank)
//    {
//        Integer tNumEntities = tEntityIndices.n_cols();
//        moris::Matrix< Integer_Matrix > tEntityIds(1,tNumEntities);
//
//        for(Integer i =0; i<tNumEntities; i++)
//        {
//            tEntityIds(0,i) = aMeshData.get_glb_entity_id_from_entity_loc_index(tEntityIndices(0,i),aEntityRank);
//        }
//        return tEntityIds;
//    }

};

}

#endif /* SRC_MESH_CL_MESH_TOOLS_HPP_ */
