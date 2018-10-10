/*
 * cl_geom_Field.hpp
 *
 *  Created on: Oct 3, 2018
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_

// XTKL: Matrix Include

#include "cl_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"
namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Geom_Field : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Geom_Field()
    {
    }

    Geom_Field( std::shared_ptr<moris::mtk::Field> aFieldPtr ):
        mField(aFieldPtr)
    {
    }

    bool is_analytic() const
    {
        return false;
    }

    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 2;
        aNumCols = 3;
    }


    moris::Matrix< Integer_Matrix >
    get_node_adv_indices(moris::Matrix< Integer_Matrix > const & aNodeIndices)
    {

        XTK_ERROR<<"get_node_adv_indices function is not implemented in geometry field";

        return moris::Matrix<Integer_Matrix>(0,0);
    }

    /**
     * This assumes you are working with the active level set mesh
     */
    Real access_field_value_with_entity_index(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        XTK_ASSERT(aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported");
        const moris::Matrix< moris::DDRMat > & tNodeVals =mField->get_node_values();

        return tNodeVals(aEntityIndex);
    }

    moris::Matrix< Real_Matrix > evaluate_sensitivity_dx_dp(moris::Matrix< Real_Matrix > const & aLocalCoordinate, Integer aEntityIndex, enum EntityRank aEntityRank)
    {
        //TODO: Implement this function
        moris::Matrix< Real_Matrix > tSensitivityDxDp(1,1,0);
        XTK_ERROR<<"evaluate_sensitivity_dx_dp function is not implemented in level set mesh";
        return tSensitivityDxDp;
    }


private:

    std::shared_ptr<moris::mtk::Field> mField;
};
}



#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_GEOM_FIELD_HPP_ */
