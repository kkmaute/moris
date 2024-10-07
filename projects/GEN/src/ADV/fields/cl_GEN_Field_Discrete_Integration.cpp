/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Discrete_Integration.cpp
 *
 */

#include "cl_GEN_Field_Discrete_Integration.hpp"

#include <utility>
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Basis_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            const Vector< ADV >& aADVs,
            std::string          aName )
            : Field( aADVs, std::move( aName ) )
            , mMeshPair( nullptr, nullptr )    // FIXME
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            const Vector< sint >& aSharedADVIds,
            const mtk::Mesh_Pair& aMeshPair,
            std::string           aName )
            : Field( aSharedADVIds, std::move( aName ) )
            , mMeshPair( aMeshPair )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field_Discrete_Integration::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return this->get_field_value( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field_Discrete_Integration::get_field_value(
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        // Return result
        return this->get_interpolated_field_value(
                aDerivedNode.get_locator_nodes(),
                aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Field_Discrete_Integration::get_dfield_dadvs(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        return this->get_dfield_dadvs( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_dfield_dadvs(
            Matrix< DDRMat >&   aSensitivities,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        this->append_interpolated_dfield_dadvs(
                aSensitivities,
                aDerivedNode.get_locator_nodes(),
                aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Field_Discrete_Integration::get_determining_adv_ids(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        return this->get_determining_adv_ids( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_determining_adv_ids(
            Vector< sint >&   aDeterminingADVIDs,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        this->append_interpolated_determining_adv_ids(
                aDeterminingADVIDs,
                aDerivedNode.get_locator_nodes(),
                aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_dfield_dcoordinates(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ERROR( false, "Discrete dfield_dcoordinates is right now only handled by the level set geometry." );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Field_Discrete_Integration::get_determining_adv_ids( uint aNodeIndex )
    {
        return Field::get_determining_adv_ids( aNodeIndex, {{}} );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Field_Discrete_Integration::get_mtk_field()
    {
        return this->create_mtk_field( mMeshPair );
    }

    //--------------------------------------------------------------------------------------------------------------

}
