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
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Interpolation_Function.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node::Basis_Node(
            Node* aNode,
            real  aBasis )
            : mNode( aNode )
            , mBasis( aBasis )
    {
        MORIS_ASSERT( aNode, "A GEN Basis_Node must be created with a valid node." );
    }

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node::Basis_Node(
            const Parent_Node& aParentNode,
            real               aBasis )
            : mNode( aParentNode.get_node() )
            , mBasis( aBasis )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Basis_Node::get_index() const
    {
        return mNode->get_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Basis_Node::get_global_coordinates() const
    {
        return mNode->get_global_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Basis_Node::get_basis() const
    {
        return mBasis;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Basis_Node::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor )
    {
        mNode->append_dcoordinate_dadv( aCoordinateSensitivities, aSensitivityFactor );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Basis_Node::get_coordinate_determining_adv_ids()
    {
        return mNode->get_coordinate_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

}
