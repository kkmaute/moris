#ifndef MORIS_CL_GEN_GEOMETRY_OBJECT_MANAGER_HPP_
#define MORIS_CL_GEN_GEOMETRY_OBJECT_MANAGER_HPP_

#include <unordered_map>

#include "cl_GEN_Geometry_Object.hpp"
#include "cl_Matrix.hpp"


namespace moris
{
    namespace ge
    {
        class Geometry_Object_Manager
        {
        private:
            // Geometry objects
            moris::Cell<GEN_Geometry_Object> mGeometryObjects;

            // Node to Geometry Object Map (key - node index, val - geometry object index)
            std::unordered_map<moris::moris_index, moris::moris_index> mNodeToGeomObjectMap;

        public:
            /**
             * Constructor
             */
            Geometry_Object_Manager();

            /**
             * stores geometry objects associated with nodes in the geometry object manager
             */
            void store_geometry_objects( moris::Matrix< moris::IndexMat > const & aNodeIndices,
                                    moris::Cell<GEN_Geometry_Object> const & aGeometryObjects );

            /**
             * Returns the geometry object associated with the specified node index [DEPRECATED]
             */
            GEN_Geometry_Object& get_geometry_object_from_manager(moris::moris_index const & aNodeIndex);

            /**
             * Link a node to the geometry object belonging to a different node
             *
             * @param aNodeIndexWithGeomObj The index of the node with the geometry object
             * @param aNodeIndexToLink The index of the node to link
             */
             void link_to_node_to_another_nodes_geometry_object(moris::moris_index aNodeIndexWithGeomObj,
                                                           moris::moris_index aNodeIndexToLink);

             /**
              * Returns the geometry object associated with the specified node index
              * Const version of above
              * [DEPRECATED]
              */
             const GEN_Geometry_Object& get_geometry_object_from_manager(moris::moris_index const & aNodeIndex) const;

             /**
              * [DEPRECATED]
              */
             GEN_Geometry_Object* get_geometry_object(uint aNodeIndex)
             {
                 auto tIter = mNodeToGeomObjectMap.find(aNodeIndex);
                 if (tIter == mNodeToGeomObjectMap.end())
                 {
                     return nullptr;
                 }
                 else
                 {
                     moris::moris_index tGOIndex = tIter->second;
                     return &(mGeometryObjects(tGOIndex));
                 }
             }

             /**
              * Gets the number of geometry objects stored in the manager
              *
              * @return Number of geometry objects
              */
             uint get_num_geometry_objects();

        };
    }
}

#endif /* MORIS_CL_GEN_GEOMETRY_OBJECT_MANAGER_HPP_ */
