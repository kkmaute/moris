/*
 * cl_Geom_Field.hpp
 *
 *  Created on: Mar 1, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_


#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"
#include "cl_HMR_Field.hpp"
namespace xtk
{
class Geom_Field : public Geometry
{
public:
    Geom_Field()
    {
    }

    Geom_Field( std::shared_ptr<moris::hmr::Field> aFieldPtr ):
        mField(aFieldPtr)
    {
    }

    bool is_analytic() const
    {
        return false;
    }

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 2;
        aNumCols = 3;
    }


    moris::Matrix< moris::IndexMat >
    get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
    {

        MORIS_ERROR(0,"get_node_adv_indices function is not implemented in geometry field");

        return moris::Matrix<moris::IndexMat>(0,0);
    }

    /**
     * This assumes you are working with the active level set mesh
     */
    moris::real access_field_value_with_entity_index(moris::moris_index     aEntityIndex,
                                                     enum moris::EntityRank aEntityRank) const
    {
        MORIS_ASSERT(aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported");
        const moris::Matrix< moris::DDRMat > & tNodeVals =mField->get_node_values();

        return tNodeVals(aEntityIndex);
    }

    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dx_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::uint aEntityIndex, enum EntityRank aEntityRank)
    {
        //TODO: Implement this function
        moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
        MORIS_ERROR(0,"evaluate_sensitivity_dx_dp function is not implemented in level set mesh");
        return tSensitivityDxDp;
    }


private:

    std::shared_ptr<moris::hmr::Field> mField;
};
}



#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_ */
