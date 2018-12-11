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
#include "cl_Mesh_Enums.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace mesh
{
class Entity
{
public:
    Entity() :
            mGlbId(std::numeric_limits<moris::moris_id>::max()), mLocInd(std::numeric_limits<moris::moris_id>::max())
    {

    }

    ~Entity()
    {

    }

    void set_entity_identifiers(moris::moris_id aGlbId,
                                moris::moris_id aLocInd,
                                enum moris::EntityRank aEntityRank)
    {
        mGlbId = aGlbId;
        mLocInd = aLocInd;
        mEntityRank = aEntityRank;
    }

    void set_entity_coords(moris::Matrix< moris::DDRMat > const & aCoordinates)
    {
        if (mEntityRank == moris::EntityRank::NODE)
        {
            mEntityCoordinates = aCoordinates.copy();
        }

        else
        {
            XTK_ERROR<<"Only nodes should have coordinates in this context to avoid duplicate coordinate storage";
        }
    }

    void set_field_data(moris::Matrix< moris::DDRMat > const & aFieldData)
    {
        mNumFields = aFieldData.n_cols();
        mFieldData = aFieldData.copy();
    }

    moris::moris_index get_entity_loc_index() const
    {
        XTK_ASSERT(mLocInd!=std::numeric_limits<moris::moris_id>::max(),"Index has not been set");

        return mLocInd;
    }

    moris::moris_index get_entity_glb_id() const
    {
        XTK_ASSERT(mGlbId!=std::numeric_limits<moris::moris_id>::max(),"Id has not been set");
        return mGlbId;
    }

    moris::Matrix< moris::DDRMat > const &
    get_entity_coords() const
    {
        return mEntityCoordinates;
    }

    moris::real
    get_field_data(moris::moris_index aFieldIndex) const
    {
        XTK_ASSERT(mNumFields!=0,"Fields have not been set");
        XTK_ASSERT(aFieldIndex<(moris::moris_index)mNumFields,"Field index is outside of bounds. Note this function should not be used directly but via STK_Mesh_Data only");
        return mFieldData(0,aFieldIndex);
    }



private:
    moris::moris_id              mGlbId;
    moris::moris_id              mLocInd;
    moris::size_t                mNumFields;
    enum moris::EntityRank       mEntityRank;
    moris::Matrix< moris::DDRMat > mFieldData;
    moris::Matrix< moris::DDRMat > mEntityCoordinates; // If its a node
};
}

#endif /* SRC_MESH_CL_MESH_ENTITY_HPP_ */
