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

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            Matrix< DDRMat > aConstants,
            uint             aNumOriginalNodes )
            : Field( aConstants, "" )
    {
        mNumOriginalNodes = aNumOriginalNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Field_Discrete_Integration::Field_Discrete_Integration(
            const Matrix< DDSMat >& aSharedADVIds,
            uint                    aNumOriginalNodes,
            std::string             aName )
            : Field( aSharedADVIds, aName )
    {
        mNumOriginalNodes = aNumOriginalNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field_Discrete_Integration::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( aNodeIndex < mNumOriginalNodes )
        {
            return this->get_field_value( aNodeIndex );
        }
        else
        {
            MORIS_ASSERT( ( aNodeIndex - mNumOriginalNodes ) < mChildNodes.size(),
                    "A discrete field value was requested from a node that this field doesn't know. "
                    "Perhaps a child node was not added to this field?" );
            return mChildNodes( aNodeIndex - mNumOriginalNodes )->interpolate_field_value( this );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Field_Discrete_Integration::get_dfield_dadvs(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        if (aNodeIndex < mNumOriginalNodes)
        {
            return this->get_dfield_dadvs(aNodeIndex);
        }
        else
        {
            MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                    "A discrete field dfield_dadvs was requested from a node that this field doesn't know. "
                    "Perhaps a child node was not added to this field?");
            return mChildNodes(aNodeIndex - mNumOriginalNodes)->join_field_sensitivities(this);
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat> Field_Discrete_Integration::get_determining_adv_ids(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        if (aNodeIndex < mNumOriginalNodes)
        {
            return this->get_determining_adv_ids(aNodeIndex);
        }
        else
        {
            MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                         "Determining ADV IDs were requested from a node that this discrete field doesn't know. "
                         "Perhaps a child node was not added to this field?");
            return mChildNodes(aNodeIndex - mNumOriginalNodes)->join_determining_adv_ids(this);
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
        mChildNodes.resize(0);
    }

    //--------------------------------------------------------------------------------------------------------------

}
