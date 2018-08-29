/*
 * cl_Mesh_Entity.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_ENTITY_HPP_
#define SRC_MESH_CL_MESH_ENTITY_HPP_

#include <limits>
#include <memory>

#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace mesh
{
template<typename Real, typename Integer, typename Real_Matrix>
class Entity
{
public:
    Entity() :
            mGlbId(std::numeric_limits<Integer>::max()), mLocInd(std::numeric_limits<Integer>::max())
    {

    }

    ~Entity()
    {

    }

    void set_entity_identifiers(Integer aGlbId, Integer aLocInd, enum EntityRank aEntityRank)
    {
        mGlbId = aGlbId;
        mLocInd = aLocInd;
        mEntityRank = aEntityRank;
    }

    void set_entity_coords(moris::Mat_New<Real,Real_Matrix> const & aCoordinates)
    {
        if (mEntityRank == EntityRank::NODE)
        {
            mEntityCoordinates = aCoordinates.copy();
        }

        else
        {
            XTK_ERROR<<"Only nodes should have coordinates in this context to avoid duplicate coordinate storage";
        }
    }

    void set_field_data(moris::Mat_New<Real,Real_Matrix> const & aFieldData)
    {
        mNumFields = aFieldData.n_cols();
        mFieldData = aFieldData.copy();
    }

    Integer get_entity_loc_index() const
    {
        XTK_ASSERT(mLocInd!=std::numeric_limits<Integer>::max(),"Index has not been set");

        return mLocInd;
    }

    Integer get_entity_glb_id() const
    {
        XTK_ASSERT(mGlbId!=std::numeric_limits<Integer>::max(),"Id has not been set");
        return mGlbId;
    }

    moris::Mat_New<Real,Real_Matrix> const &
    get_entity_coords() const
    {
        return mEntityCoordinates;
    }

    Real
    get_field_data(Integer aFieldIndex) const
    {
        XTK_ASSERT(mNumFields!=0,"Fields have not been set");
        XTK_ASSERT(aFieldIndex<mNumFields,"Field index is outside of bounds. Note this function should not be used directly but via STK_Mesh_Data only");
        return mFieldData(0,aFieldIndex);
    }



private:
    Integer mGlbId;
    Integer mLocInd;
    Integer mNumFields;
    enum EntityRank mEntityRank;
    moris::Mat_New<Real,Real_Matrix> mFieldData;
    moris::Mat_New<Real,Real_Matrix> mEntityCoordinates; // If its a node
};
}

#endif /* SRC_MESH_CL_MESH_ENTITY_HPP_ */
