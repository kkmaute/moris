/*
 * cl_GEN_Geom_Field.hpp
 *
 *  Created on: Feb 25, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_GEOM_FIELD_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_GEOM_FIELD_HPP_

#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Discrete.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Field : public Geometry_Discrete
        {
        private :
            GEN_Field* mField = nullptr;

        //------------------------------------------------------------------------------
        public :

            Geometry_Field(GEN_Field* aField)
            : Geometry_Discrete(Matrix<DDRMat>(1, 1, 0.0))
            {
                mField = aField;
            }

            ~Geometry_Field(){}

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aEntityIndex the index of the field value
             * @return field value at the specified index
             */
            real evaluate_field_value(moris_index aEntityIndex)
            {
                MORIS_ASSERT( mField != nullptr, "GEN_Geom_Field::evaluate_field_value_with_coordinate() - field pointer not set" );
                return mField->get_field_val_at_vertex( aEntityIndex );    // note: aEntityIndex corresponds to the node index in the phase table
            }

            /**
             * Given an index, the discrete geometry needs to return sensitivites with respect to the field value
             *
             * @param aEntityIndex the index of the field value
             * @return matrix of sensitivities
             */
            virtual Matrix<DDRMat> evaluate_sensitivity(moris_index aEntityIndex)
            {
                MORIS_ERROR(false, "Sensitivity not implemented for Geometry_Field.");
                return Matrix<DDRMat>(1, 1, 0.0);
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
