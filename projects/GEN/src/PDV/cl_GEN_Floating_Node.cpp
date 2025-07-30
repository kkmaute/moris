/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Floating_Node.cpp
 *
 */

#include "cl_GEN_Floating_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Floating_Node::Floating_Node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Matrix< DDRMat >&           aParametricCoordinates,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
            : Derived_Node(
                      aNodeIndex,
                      aBackgroundNodes,
                      aParametricCoordinates,
                      aBackgroundGeometryType,
                      aBackgroundInterpolationOrder )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Floating_Node::depends_on_advs() const
    {
        return this->get_interface_geometry().depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Floating_Node::is_on_interface( const Geometry& aGeometry ) const
    {
        return &aGeometry == &this->get_interface_geometry();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Floating_Node::get_num_pdvs()
    {
        return this->get_global_coordinates().numel();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Floating_Node::set_starting_pdv_id( moris_id aPDVStartingID )
    {
        mPDVStartingID = aPDVStartingID;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Floating_Node::get_starting_pdv_id()
    {
        return mPDVStartingID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Floating_Node::set_id( moris_id aNodeID )
    {
        mNodeID = aNodeID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Floating_Node::set_owner( moris_index aNodeOwner )
    {
        mNodeOwner = aNodeOwner;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Floating_Node::get_id() const
    {
        return mNodeID;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Floating_Node::get_owner()
    {
        return mNodeOwner;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
