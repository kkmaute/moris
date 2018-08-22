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
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
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
    store_geometry_objects(Mat<Integer,Integer_Matrix> const & aNodeIndices,
                           Cell<Geometry_Object<Real, Integer, Real_Matrix,Integer_Matrix>> const & aGeometryObjects)
    {
        Integer tNumExistingGeometryObjects = mGeometryObjects.size();
        Integer tNumNewGeometryObjects = aNodeIndices.get_num_columns();

        XTK_ASSERT(tNumNewGeometryObjects == aGeometryObjects.size(),"Number of geometry objects does not match number of node indices provided.");

        // append the geometry object cell
        mGeometryObjects.append(aGeometryObjects);

        // Add to map
        for(Integer i = 0; i< tNumNewGeometryObjects; i++)
        {
            mNodeToGeomObjectMap[aNodeIndices(0,i)] = i + tNumExistingGeometryObjects;
        }
    }

    /*
     * Returns the geometry object associated with the specified node index
     */
    Geometry_Object<Real,Integer,Real_Matrix,Integer_Matrix> &
    get_geometry_object_from_manager(Integer const & aNodeIndex)
    {
        XTK_ASSERT(mNodeToGeomObjectMap.find(aNodeIndex)!=mNodeToGeomObjectMap.end(),"Node index does not have an associated geometry object");

        Integer tGOIndex = mNodeToGeomObjectMap[aNodeIndex];

        return mGeometryObjects(tGOIndex);
    }

     /*
      * Returns the geometry object associated with the specified node index
      * Const version of above
      */
     Geometry_Object<Real,Integer,Real_Matrix,Integer_Matrix> const &
     get_geometry_object_from_manager(Integer const & aNodeIndex) const
     {
       auto  tIter = mNodeToGeomObjectMap.find(aNodeIndex);
         XTK_ASSERT(tIter!=mNodeToGeomObjectMap.end(),"Node index does not have an associated geometry object");

         Integer tGOIndex = tIter->second;

         return mGeometryObjects(tGOIndex);
     }

private:
    // Geometry objects
    Cell<Geometry_Object<Real, Integer, Real_Matrix,Integer_Matrix>> mGeometryObjects;

    // Node to Geometry Object Map (key - node index, val - geometry object index)
    std::unordered_map<Integer, Integer> mNodeToGeomObjectMap;

};
}


#endif /* SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_MANAGER_HPP_ */
