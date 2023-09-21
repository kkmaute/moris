/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometric_Query_Interface.hpp
 *
 */

#ifndef MORIS_CL_GEN_Geometric_Query_Interface_HPP_
#define MORIS_CL_GEN_Geometric_Query_Interface_HPP_

#include "cl_Mesh_Enums.hpp"

namespace moris
{
    namespace ge
    {
        enum class Query_Type
        {
            INTERSECTION_NO_LOCATION,
            INTERSECTION_LOCATION,
            VALUE_AT_PARAMETRIC_LOCATION,
            INVALID
        };

        class Geometric_Query_Interface
        {
          private:

          public:
            Geometric_Query_Interface()
            {
            }

            virtual ~Geometric_Query_Interface()
            {
            }

            virtual Query_Type get_query_type() const = 0;

            virtual moris_index get_geometric_index() const = 0;

            virtual EntityRank get_query_entity_rank() const = 0;

            virtual Matrix< IndexMat > const & get_query_entity_to_vertex_connectivity() const = 0;

            virtual moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > >* get_query_indexed_coordinates() const = 0;

            virtual Matrix< DDRMat > get_vertex_local_coord_wrt_parent_entity( moris_index aVertexIndex ) const = 0;

            virtual enum EntityRank get_query_parent_entity_rank() const = 0;

            virtual Matrix< IndexMat > get_query_parent_entity_connectivity() const = 0;

            virtual Matrix< DDRMat > get_query_parent_coordinates() const = 0;

            virtual moris_index get_query_parent_entity_id() const = 0;

            virtual moris_index max_query_entity_intersection() const = 0;
        };

    }    // namespace ge

}    // namespace moris

#endif
