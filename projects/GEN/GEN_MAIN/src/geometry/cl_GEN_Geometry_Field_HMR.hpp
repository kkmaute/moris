/*
 * cl_GEN_Geom_Field_HMR.hpp
 *
 *  Created on: Oct 24, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_FIELD_HMR_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_FIELD_HMR_HPP_

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete.hpp"
#include "cl_Matrix.hpp"
#include "cl_HMR_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Field_HMR : public Geometry, public Field_Discrete
        {
        public:
            Geometry_Field_HMR(std::shared_ptr<moris::hmr::Field> aFieldPtr )
                    : Field(Matrix<DDRMat>(1, 1, 0.0))
            {
                mField = aFieldPtr;
            }

            //----------------------------------------------------------------------------------------------------------

            real evaluate_field_value(uint aEntityIndex)
            {
                const moris::Matrix< moris::DDRMat > & tNodeVals = mField->get_node_values();

                return tNodeVals(aEntityIndex);
            }

            //----------------------------------------------------------------------------------------------------------

            void evaluate_all_sensitivities(uint aEntityIndex, Matrix<DDRMat>& aSensitivities)
            {
                MORIS_ERROR(false, "evaluate_sensitivity function is not implemented in HMR geometry field"); //TODO: Implement this function
            }




            moris::Matrix< moris::IndexMat >
            get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
            {

                MORIS_ERROR(0,"get_node_adv_indices function is not implemented in geometry field");

                return moris::Matrix<moris::IndexMat>(0,0);
            }

            moris::real access_field_value_with_entity_index(moris::moris_index     aEntityIndex,
                                                             enum moris::EntityRank aEntityRank) const
            {
                MORIS_ASSERT(aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported");
                const moris::Matrix< moris::DDRMat > & tNodeVals = mField->get_node_values();

                return tNodeVals(aEntityIndex);
            }

            moris::Matrix< moris::DDRMat > evaluate_sensitivity_dx_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::uint aEntityIndex, enum EntityRank aEntityRank)
            {
                //TODO: Implement this function
                moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
                MORIS_ERROR(0,"evaluate_sensitivity_dx_dp function is not implemented in geometry field");
                return tSensitivityDxDp;
            }


        private:

            std::shared_ptr<moris::hmr::Field> mField;
        };
    }
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_FIELD_HMR_HPP_ */
