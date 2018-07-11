/*
 * cl_Entity_Tracker.hpp
 *
 *  Created on: Apr 3, 2017
 *      Author: doble
 */

#ifndef SRC_XTK_CL_ENTITY_TRACKER_HPP_
#define SRC_XTK_CL_ENTITY_TRACKER_HPP_

#include "cl_XTK_Enums.hpp"               // For use of the xtk specific enums
#include "cl_Mat.hpp"              // For use of moris::Mat // LNA/src
#include "cl_Mesh_Enums.hpp"         // For use of entity ranks // MTK/src
#include "cl_Cell.hpp"         // For use of moris::Cell // CON/src


namespace moris
{
    namespace xtk
    {
        class EntityTracker
        {
        public:
            EntityTracker(enum EntityRank aEntityRanktoTrack,
                          enum EntityRank aChildEntityRank,
                          moris::uint     aNumEntitiestoTrack,
                          moris::uint     aNumChildrenAllowed);

            ~EntityTracker();

            void
            mark_entity_as_used(moris::uint aEntityIndex);


            void
            set_child_entity_lcl_index(moris::uint aParentEntityIndex,
                    moris::uint aSecondaryIndex,
                    moris::uint aChildEntityIndex);

            void
            set_child_entity_glb_id(moris::uint aParentEntityIndex,
                    moris::uint aSecondaryIndex,
                    moris::uint aChildEntityId);

            bool
            is_parent_entity_used(moris::uint aEntityIndex);

            /*
             * Same as above but for when more than 1 children are allowed. Sets new request if unused
             * Returns Null ptr if entity has not been used
             *         Cell 0 - NULL if not used, otherwise it has a garbage address
             *         Cell 1 - Index Pointer
             *         Cell 2 - Id    Pointer
             *
             */
            moris::Cell<moris::uint*>
            is_parent_entity_used(moris::uint aEntityIndex,
                    moris::uint aSecondaryIndex);

            moris::uint*
            get_index_pointer(moris::uint aParentEntityIndex);

            moris::uint*
            get_id_pointer(moris::uint aParentEntityIndex);

            moris::uint
            get_request_index_from_entity_index(moris::uint aEntityIndex);

            moris::uint
            get_num_children_allowed();



            void print();
        private:
            moris::Mat<moris::uint>  mUseMarker;           // Marks how many times an entity has been used
            // Id then indices
            moris::Mat<moris::uint>  mEntityTrackerInfo;   // Requests point to a location in this matrix
            moris::Mat<moris::uint>  mRequestIndexTracker; //
            moris::uint mChildrenAllowed;                  // Number of children allowed
            moris::uint mReqCounter;


            /*
             * Returns the child index for a given parent index
             */
            moris::uint
            get_child_location(moris::uint  aParentIndex,
                               moris::uint  aSecondaryIndex);

        };
    }

}

#endif /* SRC_XTK_CL_ENTITY_TRACKER_HPP_ */
