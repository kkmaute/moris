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
#include "cl_MTK_Field_Discrete.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            Matrix< DDRMat > aConstants,
            uint             aNumOriginalNodes )
            : Field( aConstants, "" )
            , mMeshPair( nullptr, nullptr ) // FIXME
    {
        mNumOriginalNodes = aNumOriginalNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            const Matrix< DDSMat >& aSharedADVIds,
            mtk::Mesh_Pair          aMeshPair,
            std::string             aName )
            : Field( aSharedADVIds, aName )
            , mMeshPair( aMeshPair )
    {
        mNumOriginalNodes = mMeshPair.get_interpolation_mesh()->get_num_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field_Discrete_Integration::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            return this->get_field_value( aNodeIndex );
        }
        else
        {
            return this->get_interpolated_field_value( aNodeIndex, aCoordinates );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Field_Discrete_Integration::get_dfield_dadvs(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            return this->get_dfield_dadvs( aNodeIndex );
        }
        else
        {
            return this->get_interpolated_dfield_dadvs( aNodeIndex, aCoordinates );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat> Field_Discrete_Integration::get_determining_adv_ids(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            return this->get_determining_adv_ids( aNodeIndex );
        }
        else
        {
            return this->get_interpolated_determining_adv_ids( aNodeIndex, aCoordinates );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::get_dfield_dcoordinates(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ASSERT(aNodeIndex >= mNumOriginalNodes,
                    "Discrete field dfield_dcoordinates is only valid for intersection node indices.");
        MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                     "A discrete field dfield_dcoordinates was requested from a node that this field doesn't know. "
                     "Perhaps a child node was not added to this field?");
        mChildNodes(aNodeIndex - mNumOriginalNodes)->get_dfield_dcoordinates(this, aSensitivities);
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat> Field_Discrete_Integration::get_determining_adv_ids(uint aNodeIndex)
    {
        return Field::get_determining_adv_ids(aNodeIndex, {{}});
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
    {
        MORIS_ASSERT(aNodeIndex == mNumOriginalNodes + mChildNodes.size(),
                "Child nodes must be added to a level set field in order by node index.");
        mChildNodes.push_back(aChildNode);
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field_Discrete_Integration::reset_nodal_data()
    {
        mChildNodes.resize( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Field_Discrete_Integration::get_mtk_field()
    {
        return this->create_mtk_field( mMeshPair );
    }

    //--------------------------------------------------------------------------------------------------------------

}
