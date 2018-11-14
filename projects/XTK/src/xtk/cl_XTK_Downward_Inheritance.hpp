/*
 * cl_XTK_Downward_Inheritance.hpp
 *
 *  Created on: Jul 21, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_DOWNWARD_INHERITANCE_HPP_
#define SRC_XTK_CL_XTK_DOWNWARD_INHERITANCE_HPP_


#include <utility>
#include <limits>
#include <unordered_map>

#include "containers/cl_XTK_Cell.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{

template<typename T1,typename T2>
class Downward_Inheritance
{
public:

    Downward_Inheritance()
    {}


    Downward_Inheritance(size_t aNumEntities):
        mInheritance(aNumEntities)
    {
    }

    /*
     * Returns a bool but not the value
     */
    bool has_inheritance(T1 const & aKey) const
    {
        bool tHasInheritance = false;

        if(mInheritance.find(aKey) != mInheritance.end())
        {
            tHasInheritance = true;
        }
        return tHasInheritance;
    }


    /*
     * Returns the Value associated with a Key (i.e. XTK Element)
     */
    T2 const & get_inheritance(T1 const & aKey )
    {
        bool tTrue = has_inheritance(aKey);
        XTK_ASSERT(tTrue,"No Key Located");
        return mInheritance[aKey];
    }

    /*
     * i - processor local index of parent entity to establish downward inheritance for
     * n - number of children built on the parent entity
     * @param[in] aKey - Index Element in XTK Mesh
     * @param[in] aValue - Index of the child mesh in Cut Mesh
     */
    void register_new_inheritance_pair(T1 const & aKey, T2 const & aValue)
    {
        if(!has_inheritance(aKey))
        {
            mInheritance[aKey] = aValue;
        }
        else
        {
            XTK_ERROR<<"Already registered\n";
        }
    }

private:
    /*
     *
     */
    std::unordered_map<T1, T2> mInheritance;

private:
};
}

#endif /* SRC_XTK_CL_XTK_DOWNWARD_INHERITANCE_HPP_ */
