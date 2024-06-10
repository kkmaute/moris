/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node.cpp
 *
 */

#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node::Intersection_Node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            real                              aLocalCoordinate,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
            : Derived_Node(
                    aNodeIndex,
                    aBackgroundNodes,
                    0.5 * ( 1.0 - aLocalCoordinate ) * aFirstParentNode.get_parametric_coordinates() + 0.5 * ( 1.0 + aLocalCoordinate ) * aSecondParentNode.get_parametric_coordinates(),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder )
            , mParentNodes( { Basis_Node( aFirstParentNode, 0.5 * ( 1.0 - aLocalCoordinate ) ), Basis_Node( aSecondParentNode, 0.5 * ( 1.0 + aLocalCoordinate ) ) } )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Basis_Node& Intersection_Node::get_first_parent_node() const
    {
        return mParentNodes( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Basis_Node& Intersection_Node::get_second_parent_node() const
    {
        return mParentNodes( 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Intersection_Node::depends_on_advs() const
    {
        return this->get_interface_geometry().depends_on_advs() or mParentNodes( 0 ).depends_on_advs() or mParentNodes( 1 ).depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Intersection_Node::get_locator_nodes() const
    {
        return mParentNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Intersection_Node::is_on_interface( const Geometry& aGeometry ) const
    {
        return &aGeometry == &this->get_interface_geometry();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Intersection_Node::parent_edge_is_intersected()
    {
        return std::abs( this->get_local_coordinate() ) <= 1.0 or this->is_first_parent_on_interface() or this->is_second_parent_on_interface();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Intersection_Node::is_first_parent_on_interface()
    {
        return this->get_interface_geometry().get_geometric_region( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() ) == Geometric_Region::INTERFACE
            or std::abs( this->get_local_coordinate() + 1.0 ) < this->get_interface_geometry().get_intersection_tolerance();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Intersection_Node::is_second_parent_on_interface()
    {
        return this->get_interface_geometry().get_geometric_region( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() ) == Geometric_Region::INTERFACE
            or std::abs( this->get_local_coordinate() - 1.0 ) < this->get_interface_geometry().get_intersection_tolerance();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Intersection_Node::get_local_coordinate() const
    {
        return 1.0 - 2.0 * this->get_first_parent_node().get_basis();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Intersection_Node::get_num_pdvs()
    {
        return this->get_global_coordinates().numel();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node::set_starting_pdv_id( moris_id aPDVStartingID )
    {
        mPDVStartingID = aPDVStartingID;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Intersection_Node::get_starting_pdv_id()
    {
        return mPDVStartingID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node::set_id( moris_id aNodeID )
    {
        mNodeID = aNodeID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node::set_owner( moris_index aNodeOwner )
    {
        mNodeOwner = aNodeOwner;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Intersection_Node::get_id() const
    {
        return mNodeID;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Intersection_Node::get_owner()
    {
        return mNodeOwner;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node::join_adv_ids(
            Vector< sint >&       aCombinedIDs,
            const Vector< sint >& aIDsToAdd )
    {
        // Resize IDs
        uint tJoinedSensitivityLength = aCombinedIDs.size();
        aCombinedIDs.resize( tJoinedSensitivityLength + aIDsToAdd.size() );

        // Join IDs
        for ( uint tAddedSensitivity = 0; tAddedSensitivity < aIDsToAdd.size(); tAddedSensitivity++ )
        {
            aCombinedIDs( tJoinedSensitivityLength + tAddedSensitivity ) = aIDsToAdd( tAddedSensitivity );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
