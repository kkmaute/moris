/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Scalar_Field_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SCALAR_FIELD_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SCALAR_FIELD_INFO_HPP_

#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris
{
namespace mtk
{
template<typename Matrix_Type>
struct Scalar_Field_Info
{
    typedef typename Matrix<Matrix_Type>::Data_Type Field_Data_Type;

    /*
     * Default constructor
     */
    Scalar_Field_Info():
        mFieldEntityRank(EntityRank::INVALID),
        mFieldData(nullptr)
    {

    }

    /*
     * Set the field name for this field
     */
    void
    set_field_name(std::string aFieldName)
    {
        MORIS_ASSERT(!field_has_name(),"Trying to overwrite field name");
        mFieldName = aFieldName;
    }

    /*
     * Retrieve the field name
     */
    std::string
    get_field_name() const
    {
        MORIS_ASSERT(field_has_name(),"Field does not have a name associated with it");
        return mFieldName;
    }

    /*
     * Ask if the field has a name
     */
    bool
    field_has_name() const
    {
        return !mFieldName.empty();
    }

    // -----------------------------------------------------------------------

    /*
     * Set the part name (i.e. block set this field is on)
     */
    void
    set_part_name(std::string aFieldPart)
    {
        MORIS_ASSERT(!field_has_part_name(),"Field already has a part");
        mFieldPart = aFieldPart;
    }

    /*
     * Get the part name
     */
    std::string
    get_part_name() const
    {
        MORIS_ASSERT(field_has_part_name(),"Field already has part name associated with it");
        return mFieldPart;
    }

    /*
     * Ask if this mesh has a part name
     */
    bool
    field_has_part_name() const
    {
        return !mFieldPart.empty();
    }

    // -----------------------------------------------------------------------

    /*
     * Set the entity rank this field is associated with
     */
    void
    set_field_entity_rank(enum EntityRank aEntityRank)
    {
        MORIS_ASSERT(mFieldEntityRank == EntityRank::INVALID,"Field already has entity rank");
        mFieldEntityRank = aEntityRank;
    }

    enum EntityRank
    get_field_entity_rank()
    {
        MORIS_ASSERT(mFieldEntityRank!=EntityRank::INVALID,"Field has not been given a valid entity rank");
        return mFieldEntityRank;
    }

    /*
     * Add field data with field entity ids and pointer to the data
     */
    void
    add_field_data( Matrix<IdMat>*     aEntityIds,
                    Matrix<Matrix_Type>* aFieldDataPtr)
    {
        MORIS_ASSERT(!field_has_data(),"Field already has data on it");
        MORIS_ASSERT(aEntityIds->numel() == aFieldDataPtr->numel(),"Mismatch dimension between provided entity ids and field data");
        mEntityIds = aEntityIds;
        mFieldData = aFieldDataPtr;
    }

    Matrix<IdMat> const &
    get_field_entity_ids() const
    {
        return *mEntityIds;
    }

    Matrix<Matrix_Type> const &
    get_field_data() const
    {
        return *mFieldData;
    }

    /*
     * Ask if field has data
     */
    bool
    field_has_data() const
    {
        if(mFieldData==NULL || mFieldData->numel() == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    uint
    get_num_components()
    {
        return 1;
    }
    // -----------------------------------------------------------------------

private:
    // Name of the field
    std::string     mFieldName;

    // Rank of the entity this field lives on, i.e.
    // a node or element.
    enum EntityRank mFieldEntityRank;

    // Part of the mesh the field is on. i.e. a blockset
    // with the name "block_1"
    std::string     mFieldPart;

    // Field Data
    Matrix<Matrix_Type>* mFieldData;

    // Field Entity Ids associated with mFieldData (ordering here matters)
    Matrix<IdMat>* mEntityIds;
};
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_SCALAR_FIELD_INFO_HPP_ */

