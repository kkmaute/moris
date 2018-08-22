/*
 * cl_XTK_Entity.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_ENTITY_HPP_
#define SRC_XTK_CL_XTK_ENTITY_HPP_

// Standard includes
#include <limits>

// XTKL: Logging and Assertion Includes
#include "assert/fn_xtk_assert.hpp"

namespace xtk
{
template<typename Real, typename Integer>
class Entity
{
public:
    Entity()
    {

    }
    Entity(enum EntityRank aEntityRank) :
            mEntityRank(aEntityRank), mEntityLclInd(std::numeric_limits<Integer>::max()), mEntityPhase(std::numeric_limits<Integer>::max())
    {
    }

    ~Entity()
    {

    }

    void set_entity_lcl_index(Integer aLclInd)
    {
        XTK_ASSERT(mEntityLclInd == std::numeric_limits<Integer>::max(), "Entity Index already set");
        mEntityLclInd = aLclInd;
    }

    Integer get_entity_lcl_index() const
    {
        XTK_ASSERT(mEntityLclInd != std::numeric_limits<Integer>::max(), "No index set");
        return mEntityLclInd;
    }

    void set_entity_phase_index(Integer aPhase)
    {
        XTK_ASSERT(mEntityPhase == std::numeric_limits<Integer>::max(), "Entity phase already set");
        mEntityPhase = aPhase;
    }

    Integer get_entity_phase_index() const
    {
        XTK_ASSERT(mEntityPhase != std::numeric_limits<Integer>::max(), "No phase set");
        return mEntityPhase;
    }

private:
    // Local to Simple Mesh
    enum EntityRank mEntityRank;
    Integer mEntityLclInd;
    Integer mEntityPhase;
};
}


#endif /* SRC_XTK_CL_XTK_ENTITY_HPP_ */
