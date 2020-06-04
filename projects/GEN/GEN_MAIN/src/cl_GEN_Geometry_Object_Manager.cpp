#include "cl_GEN_Geometry_Object_Manager.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Object_Manager::Geometry_Object_Manager()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Object_Manager::store_geometry_objects(moris::Matrix<moris::IndexMat> const &aNodeIndices, moris::Cell<GEN_Geometry_Object> const &aGeometryObjects)
        {
            moris::size_t tNumExistingGeometryObjects = mGeometryObjects.size();
            moris::size_t tNumNewGeometryObjects = aNodeIndices.n_cols();

            MORIS_ASSERT(tNumNewGeometryObjects == aGeometryObjects.size(),
                    "Number of geometry objects does not match number of node indices provided.");

            // append the geometry object cell
            mGeometryObjects.append(aGeometryObjects);

            // Add to map
            for (moris::size_t i = 0; i < tNumNewGeometryObjects; i++)
            {
                mNodeToGeomObjectMap[aNodeIndices(0, i)] = i + tNumExistingGeometryObjects;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        GEN_Geometry_Object& Geometry_Object_Manager::get_geometry_object_from_manager(moris::moris_index const &aNodeIndex)
        {
            MORIS_ASSERT(mNodeToGeomObjectMap.find(aNodeIndex) != mNodeToGeomObjectMap.end(),
                    "Node index does not have an associated geometry object");

            moris::moris_index tGOIndex = mNodeToGeomObjectMap[aNodeIndex];

            return mGeometryObjects(tGOIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_Object_Manager::link_to_node_to_another_nodes_geometry_object(moris::moris_index aNodeIndexWithGeomObj, moris::moris_index aNodeIndexToLink)
        {
            // Geometry object index
            moris::moris_index tGOIndex = mNodeToGeomObjectMap[aNodeIndexWithGeomObj];

            // Link new node by putting it's index in map with same tGOIndex as the aNodeIndexWithGeomObj has
            mNodeToGeomObjectMap[aNodeIndexToLink] = tGOIndex;

        }

        //--------------------------------------------------------------------------------------------------------------

        const GEN_Geometry_Object& Geometry_Object_Manager::get_geometry_object_from_manager(moris::moris_index const &aNodeIndex) const
        {
            auto tIter = mNodeToGeomObjectMap.find(aNodeIndex);
            MORIS_ASSERT(tIter != mNodeToGeomObjectMap.end(),
                    "Node index does not have an associated geometry object");

            moris::moris_index tGOIndex = tIter->second;

            return mGeometryObjects(tGOIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Geometry_Object_Manager::get_num_geometry_objects()
        {
            return mGeometryObjects.size();
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}