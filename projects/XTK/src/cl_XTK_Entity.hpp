/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Entity.hpp
 *
 */

#ifndef SRC_XTK_CL_XTK_ENTITY_HPP_
#define SRC_XTK_CL_XTK_ENTITY_HPP_

// Standard includes
#include <limits>

// XTKL: Logging and Assertion Includes

namespace xtk
{
template<typename Real, typename Integer>
class Entity
{
public:
    Entity()
    {

    }
    Entity(mtk::EntityRank aEntityRank) :
            mEntityRank(aEntityRank), mEntityLclInd(std::numeric_limits<Integer>::max()), mEntityPhase(std::numeric_limits<Integer>::max())
    {
    }

    ~Entity()
    {

    }

    void set_entity_lcl_index(Integer aLclInd)
    {
        MORIS_ASSERT(mEntityLclInd == std::numeric_limits<Integer>::max(), "Entity Index already set");
        mEntityLclInd = aLclInd;
    }

    Integer get_entity_lcl_index() const
    {
        MORIS_ASSERT(mEntityLclInd != std::numeric_limits<Integer>::max(), "No index set");
        return mEntityLclInd;
    }

    void set_entity_phase_index(Integer aPhase)
    {
        MORIS_ASSERT(mEntityPhase == std::numeric_limits<Integer>::max(), "Entity phase already set");
        mEntityPhase = aPhase;
    }

    Integer get_entity_phase_index() const
    {
        MORIS_ASSERT(mEntityPhase != std::numeric_limits<Integer>::max(), "No phase set");
        return mEntityPhase;
    }

private:
    // Local to Simple Mesh
    mtk::EntityRank mEntityRank;
    Integer mEntityLclInd;
    Integer mEntityPhase;
};
}

#endif /* SRC_XTK_CL_XTK_ENTITY_HPP_ */

