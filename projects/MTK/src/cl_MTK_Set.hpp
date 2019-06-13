/*
 * cl_MTK_Element_Set.hpp
 *
 *  Created on: Jun 24, 2019
 *      Author: schmidt
 */

#ifndef SRC_MESH_CL_MTK_SET_HPP_
#define SRC_MESH_CL_MTK_SET_HPP_

#include <string>

#include "typedefs.hpp" //MRS/COR/src
#include "fn_unique.hpp" //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

#include "cl_MTK_Cell_Cluster.hpp" //MTK/src
#include "cl_MTK_Side_Cluster.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
        class Set
        {
        private :

//------------------------------------------------------------------------------

        protected :

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Set()
            { };

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Set(){};

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
//              virtual const moris::Matrix< DDUMat > &
//              get_list_of_block_cell_clusters() const = 0;

//------------------------------------------------------------------------------

              virtual const Cell_Cluster  *
              get_cell_clusters_by_index( moris_index aCellClusterIndex ) const
              {
                  MORIS_ASSERT(false, "get_cell_clusters_by_index() virtual base class used");
                  return nullptr;
              };

              virtual const Side_Cluster  *
              get_side_clusters_by_index( moris_index aCellClusterIndex ) const
              {
                  MORIS_ASSERT(false, "get_side_clusters_by_index() virtual base class used");
                  return nullptr;
              };

//------------------------------------------------------------------------------

              virtual const uint
              get_num_vertieces_on_set() const = 0;

//------------------------------------------------------------------------------

              virtual moris::Matrix< DDSMat >
              get_vertieces_inds_on_block() const = 0;

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_SET_HPP_ */
