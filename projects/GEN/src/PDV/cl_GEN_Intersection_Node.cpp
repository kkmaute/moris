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
#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node::Intersection_Node(
            uint                        aNodeIndex,
            const Cell< Node* >&        aBaseNodes,
            const Parent_Node&          aFirstParentNode,
            const Parent_Node&          aSecondParentNode,
            real                        aLocalCoordinate,
            mtk::Geometry_Type          aBaseGeometryType,
            std::shared_ptr< Geometry > aInterfaceGeometry )
            : Derived_Node(
                    aNodeIndex,
                    aBaseNodes,
                    0.5 * ( 1.0 - aLocalCoordinate ) * aFirstParentNode.get_parametric_coordinates() + 0.5 * ( 1.0 + aLocalCoordinate ) * aSecondParentNode.get_parametric_coordinates(),
                    aBaseGeometryType )
            , mParentNodes( { Basis_Node( aFirstParentNode, 0.5 * ( 1.0 - aLocalCoordinate ) ), Basis_Node( aSecondParentNode, 0.5 * ( 1.0 + aLocalCoordinate ) ) } )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node::initialize()
    {
        mParentVector  = trans( this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates() );
    }

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node& Intersection_Node::get_first_parent_node()
    {
        return mParentNodes( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Basis_Node& Intersection_Node::get_second_parent_node()
    {
        return mParentNodes( 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Intersection_Node::depends_on_advs()
    {
        return mInterfaceGeometry->depends_on_advs() or this->get_first_parent_node().depends_on_advs() or this->get_second_parent_node().depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Basis_Node >& Intersection_Node::get_locator_nodes()
    {
        return mParentNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Intersection_Node::is_on_interface( Geometry* aGeometry )
    {
        return aGeometry == mInterfaceGeometry.get();
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
        return mInterfaceGeometry->get_geometric_region( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() ) == Geometric_Region::INTERFACE;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Intersection_Node::is_second_parent_on_interface()
    {
        return mInterfaceGeometry->get_geometric_region( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() ) == Geometric_Region::INTERFACE;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node::get_local_coordinate()
    {
        return 1.0 - 2.0 * this->get_first_parent_node().get_basis();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node::get_coordinate_value( uint aCoordinateIndex )
    {
        return this->get_global_coordinates()( aCoordinateIndex );
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
        mPDVStartingID    = aPDVStartingID;
        mPDVStartingIDSet = true;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Intersection_Node::get_starting_pdv_id()
    {
        MORIS_ASSERT( mPDVStartingIDSet, "PDV Starting ID must be set for an intersection." );
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
    Intersection_Node::get_id()
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

    void
    Intersection_Node::join_adv_ids( const Matrix< DDSMat >& aIDsToAdd )
    {
        // Resize IDs
        uint tJoinedSensitivityLength = mCoordinateDeterminingADVIDs.n_cols();
        mCoordinateDeterminingADVIDs.resize( 1, tJoinedSensitivityLength + aIDsToAdd.length() );

        // Join IDs
        for ( uint tAddedSensitivity = 0; tAddedSensitivity < aIDsToAdd.length(); tAddedSensitivity++ )
        {
            mCoordinateDeterminingADVIDs( tJoinedSensitivityLength + tAddedSensitivity ) = aIDsToAdd( tAddedSensitivity );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
