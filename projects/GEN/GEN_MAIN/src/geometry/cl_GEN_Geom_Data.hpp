/*
 * cl_GEN_Geom_Data.hpp
 *
 *  Created on: Nov 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_DATA_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_DATA_HPP_

#include "cl_GEN_Geometry.hpp"
#include "cl_Matrix.hpp"
#include "cl_HMR_Field.hpp"

namespace moris
{
namespace ge
{
class GEN_Geom_Data : public GEN_Geometry
{
public:
    GEN_Geom_Data()
    {
    }

    GEN_Geom_Data( Matrix< DDRMat > aGeomData ):
        mGeomData(aGeomData)
    {
    }

    bool is_analytic() const
    {
        return false;
    }

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }


    moris::Matrix< moris::IndexMat >
    get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
    {

        MORIS_ERROR(0,"get_node_adv_indices function is not implemented in geometry data");

        return moris::Matrix<moris::IndexMat>(0,0);
    }

    /**
     * This assumes you are working with the active level set mesh
     */
    moris::real access_field_value_with_entity_index(moris::moris_index     aEntityIndex,
                                                     enum moris::EntityRank aEntityRank) const
    {
        MORIS_ASSERT( aEntityRank==EntityRank::NODE,"Only nodal levelset values are supported" );
        const moris::Matrix< moris::DDRMat > & tNodeVals = mGeomData;

        return tNodeVals(aEntityIndex);
    }

private:

    Matrix< DDRMat > mGeomData;
};

}
}


#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_GEOM_DATA_HPP_ */
