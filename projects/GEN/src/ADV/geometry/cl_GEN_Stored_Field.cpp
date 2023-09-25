/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Stored_Geometry.cpp
 *
 */

#include "cl_GEN_Stored_Field.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Stored_Field::Stored_Field(
            mtk::Mesh*               aMesh,
            std::shared_ptr< Field > aField )
            : Field_Discrete_Integration( {{}}, aMesh->get_num_nodes() )
            , mField( aField )
            , mMesh( aMesh )
            , mFieldValues( aMesh->get_num_nodes(), 1 )
    {
        this->evaluate_nodal_values();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Stored_Field::get_field_value(uint aNodeIndex)
    {
        return mFieldValues( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>&
    Stored_Field::get_dfield_dadvs(uint aNodeIndex)
    {
        return mField->get_dfield_dadvs( aNodeIndex, mMesh->get_node_coordinate( aNodeIndex ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat>
    Stored_Field::get_determining_adv_ids(uint aNodeIndex)
    {
        return mField->get_determining_adv_ids( aNodeIndex, mMesh->get_node_coordinate( aNodeIndex ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Stored_Field::reset_nodal_data()
    {
        // Reset child nodes
        Field_Discrete_Integration::reset_nodal_data();

        // Re-evaluate field values
        this->evaluate_nodal_values();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Stored_Field::evaluate_nodal_values()
    {
        // Assign nodal values
        uint tNumNodes = mMesh->get_num_nodes();
        for (uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++)
        {
            mFieldValues( tNodeIndex ) = mField->get_field_value( tNodeIndex, mMesh->get_node_coordinate( tNodeIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
