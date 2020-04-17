/*
 * cl_GEN_Geom_Field.hpp
 *
 *  Created on: Feb 25, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_GEOM_FIELD_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_GEOM_FIELD_HPP_

#include "cl_Matrix.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Geometry.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Field : public Geometry
        {
        private :
            GEN_Field* mField = nullptr;

        //------------------------------------------------------------------------------
        public :

            ~Geometry_Field(){}
        //------------------------------------------------------------------------------
            bool is_analytic() const
            {
                return false;
            }
        //------------------------------------------------------------------------------
            virtual moris::real access_field_value_with_entity_index( moris::moris_index     aEntityIndex,
                                                                      enum moris::EntityRank aEntityRank ) const
            {
                MORIS_ASSERT( mField != nullptr, "GEN_Geom_Field::evaluate_field_value_with_coordinate() - field pointer not set" );
                return mField->get_field_val_at_vertex( aEntityIndex );    // note: aEntityIndex corresponds to the node index in the phase table
            }
        //------------------------------------------------------------------------------
            GEN_Field* get_field_pointer()
            {
                return mField;
            }
        //------------------------------------------------------------------------------
        };

    }   /* end ge namespace */
}   /* end moris namespace */



#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_GEOM_FIELD_HPP_ */
