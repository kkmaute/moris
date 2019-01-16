/*
 * cl_MGE_Geometry_Object_Manager.hpp
 *
 *  Created on: Nov 21, 2017
 *      Author: doble
 */

#ifndef SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_MANAGER_HPP_
#define SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_MANAGER_HPP_

#include <unordered_map>


#include"linalg/cl_XTK_Matrix_Base.hpp"

namespace xtk
{
class Geometry_Object_Manager
{
public:
    Geometry_Object_Manager()
    {

    }

    /*
     * Stores geometry objects in the geometry object manager associated with nodes
     */
    void
    store_geometry_objects(moris::Matrix< moris::IndexMat > const & aNodeIndices,
                           Cell<Geometry_Object>            const & aGeometryObjects)
    {
        moris::size_t tNumExistingGeometryObjects = mGeometryObjects.size();
        moris::size_t tNumNewGeometryObjects = aNodeIndices.n_cols();

        XTK_ASSERT(tNumNewGeometryObjects == aGeometryObjects.size(),
                   "Number of geometry objects does not match number of node indices provided.");

        // append the geometry object cell
        mGeometryObjects.append(aGeometryObjects);

        // Add to map
        for(moris::size_t i = 0; i< tNumNewGeometryObjects; i++)
        {
            mNodeToGeomObjectMap[aNodeIndices(0,i)] = i + tNumExistingGeometryObjects;
        }
    }

    /*
     * Returns the geometry object associated with the specified node index
     */
    Geometry_Object &
    get_geometry_object_from_manager(moris::moris_index const & aNodeIndex)
    {
        XTK_ASSERT(mNodeToGeomObjectMap.find(aNodeIndex)!=mNodeToGeomObjectMap.end(),
                   "Node index does not have an associated geometry object");

        moris::moris_index tGOIndex = mNodeToGeomObjectMap[aNodeIndex];

        return mGeometryObjects(tGOIndex);
    }

     void
     link_to_node_to_another_nodes_geometry_object(moris::moris_index aNodeIndexWithGeomObj,
                                                   moris::moris_index aNodeIndexToLink)
     {
         // Geometry object index
         moris::moris_index tGOIndex = mNodeToGeomObjectMap[aNodeIndexWithGeomObj];

         // Link new node by putting it's index in map with same tGOIndex as the aNodeIndexWithGeomObj has
         mNodeToGeomObjectMap[aNodeIndexToLink] = tGOIndex;

     }


     /*
      * Returns the geometry object associated with the specified node index
      * Const version of above
      */
     Geometry_Object const &
     get_geometry_object_from_manager(moris::moris_index const & aNodeIndex) const
     {
       auto  tIter = mNodeToGeomObjectMap.find(aNodeIndex);
         XTK_ASSERT(tIter!=mNodeToGeomObjectMap.end(),"Node index does not have an associated geometry object");

         moris::moris_index tGOIndex = tIter->second;

         return mGeometryObjects(tGOIndex);
     }

private:
    // Geometry objects
    Cell<Geometry_Object> mGeometryObjects;

    // Node to Geometry Object Map (key - node index, val - geometry object index)
    std::unordered_map<moris::moris_index, moris::moris_index> mNodeToGeomObjectMap;

};
}


#endif /* SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_MANAGER_HPP_ */
