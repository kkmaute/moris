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
#include "cl_MTK_Interpolation_Function_Factory.hpp"
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
        return aNodeIndex == MORIS_INDEX_MAX ? std::numeric_limits< real >::quiet_NaN() : this->get_field_value( aNodeIndex );
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

    const Matrix< DDRMat >& Field_Discrete_Integration::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
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
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return this->get_determining_adv_ids( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_determining_adv_ids(
            Vector< sint >&     aDeterminingADVIDs,
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
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ERROR( false, "Discrete dfield_dcoordinates is right now only handled by the level set geometry." );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_dfield_dcoordinates(
            const Derived_Node& aDerivedNode,
            Matrix< DDRMat >&   aSensitivities )
    {
        // Get the background nodes (they are Basis_Node objects because this is a derived node, the basis just stores the basis value)
        Vector< const Background_Node* > tBackgroundNodes = aDerivedNode.get_background_nodes();

        // Get the parametric coordinates of the derived node within its background element
        const Matrix< DDRMat >& tParametricCoords = aDerivedNode.get_parametric_coordinates();
        uint                    tNumDims          = tParametricCoords.length();

        // Build interpolator based on the background element geometry type and interpolation order
        // Note that we use the background element interpolation order, not the locator interpolation order
        mtk::Interpolation_Function_Factory tFactory;
        mtk::Interpolation_Function_Base*   tInterpolation;
        tInterpolation = tFactory.create_interpolation_function(
                aDerivedNode.get_background_element_geometry_type(),
                mtk::Interpolation_Type::LAGRANGE,
                aDerivedNode.get_background_element_interpolation_order() );

        // Get the derivative of all the basis functions
        Matrix< DDRMat > tBasisGradXi;
        tInterpolation->eval_dNdXi( tParametricCoords, tBasisGradXi );
        uint tNumBases = tBasisGradXi.n_cols();

        // Compute the spatial jacobian
        Matrix< DDRMat > tInvJacobian = tBasisGradXi * aDerivedNode.get_background_element_nodal_coordinates();
        for ( uint iDim = 0; iDim < tNumDims; iDim++ )
        {
            tInvJacobian( iDim, iDim ) = 1.0 / tInvJacobian( iDim, iDim );
        }

        // Loop through the bases
        for ( uint iBasis = 0; iBasis < tNumBases; iBasis++ )
        {
            // Get the field value at this node
            real tPhiBasis = this->get_field_value( tBackgroundNodes( iBasis )->get_index(), tBackgroundNodes( iBasis )->get_global_coordinates() );

            aSensitivities += tInvJacobian * tBasisGradXi.get_column( iBasis ) * tPhiBasis;
        }

        // Clean up
        delete tInterpolation;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Field_Discrete_Integration::get_determining_adv_ids( uint aNodeIndex )
    {
        return Field::get_determining_adv_ids( aNodeIndex, { {} } );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Field_Discrete_Integration::get_mtk_field()
    {
        return this->create_mtk_field( mMeshPair );
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
