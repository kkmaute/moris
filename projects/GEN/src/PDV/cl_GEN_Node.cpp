/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Node.cpp
 *
 */

#include "cl_GEN_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Node::Node( uint aIndex )
            : mIndex( aIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node::get_index() const
    {
        return mIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Node::depends_on_advs() const
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Node::get_coordinate_determining_adv_ids() const
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------
}
