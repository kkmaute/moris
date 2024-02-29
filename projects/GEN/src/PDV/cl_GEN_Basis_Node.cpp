/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Basis_Node.cpp
 *
 */

#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node::Basis_Node(
            const Node& aNode,
            real        aBasis )
            : mNode( aNode )
            , mBasis( aBasis )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node::Basis_Node(
            const Parent_Node& aParentNode,
            real               aBasis )
            : mNode( aParentNode.mNode )
            , mBasis( aBasis )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Basis_Node::get_index() const
    {
        return mNode.get_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Basis_Node::get_global_coordinates() const
    {
        return mNode.get_global_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Basis_Node::get_basis() const
    {
        return mBasis;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Basis_Node::get_locator_nodes() const
    {
        return mNode.get_locator_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Basis_Node::depends_on_advs() const
    {
        return mNode.depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Basis_Node::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        mNode.append_dcoordinate_dadv( aCoordinateSensitivities, aSensitivityFactor );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Basis_Node::get_coordinate_determining_adv_ids() const
    {
        return mNode.get_coordinate_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

}
