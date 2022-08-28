/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Matrix_Field_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MATRIX_FIELD_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MATRIX_FIELD_INFO_HPP_

#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"

namespace moris
{
namespace mtk
{
template<typename Matrix_Type>
struct Matrix_Field_Info
{

    typedef typename Matrix<Matrix_Type>::Data_Type Field_Data_Type;

    /*
     * Default constructor
     */
    Matrix_Field_Info():
                           mNumRows(0),
                           mNumCols(0),
                           mFieldEntityRank(EntityRank::INVALID),
                           mFieldData(0)
    {

    }

    /*
     * Require all matrix fields to have the same amount of rows and columns
     * for that reason this needs to be specified.
     */
    Matrix_Field_Info(uint aNumRows,
                       uint aNumCols):
                           mNumRows(aNumRows),
                           mNumCols(aNumCols),
                           mFieldEntityRank(EntityRank::INVALID),
                           mFieldData(0)
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

    /*
     * Add field data with field entity ids and pointer to the data
     */
    void
    add_field_data( Matrix<IdMat>*                    aEntityIds,
                    moris::Cell<Matrix<Matrix_Type>*> aFieldDataPtr)
    {
        MORIS_ASSERT(!field_has_data(),"Field already has data on it");
        MORIS_ASSERT(aEntityIds->numel() == aFieldDataPtr.size(),"Mismatch dimension between provided entity ids and field data");
        mEntityIds = aEntityIds;
        mFieldData = aFieldDataPtr;
        MORIS_ASSERT(verify_matrix_size_matches_field_size(),"Field Matrix Data is not compatible with expected rows and columns in field");
    }

    /*
     * Set the entity rank this field is associated with
     */
    enum EntityRank
    get_field_entity_rank()
    {
        MORIS_ASSERT(mFieldEntityRank != EntityRank::INVALID,"Field does not have an entity rank");
        return mFieldEntityRank;
    }

    Matrix<IdMat> const &
    get_field_entity_ids() const
    {
        return *mEntityIds;
    }

    moris::Cell<Matrix<Matrix_Type>*> const &
    get_field_data() const
    {
        return mFieldData;
    }

    /*
     * Ask if field has data
     */
    bool
    field_has_data() const
    {
        if(mFieldData.size() == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    // -----------------------------------------------------------------------

    uint
    get_num_rows() const
    {
        return mNumRows;
    }

    uint
    get_num_cols() const
    {
        return mNumCols;
    }

    uint
    get_num_components()
    {
        return mNumCols*mNumRows;
    }

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
    const uint mNumCols;
    const uint mNumRows;
    moris::Cell<Matrix<Matrix_Type>*> mFieldData;

    // Field Entity Ids associated with mFieldData (ordering here matters)
    Matrix<IdMat>* mEntityIds;

    bool
    verify_matrix_size_matches_field_size()
    {
        bool tValid = true;
        for(uint i = 0; i <mFieldData.size(); i++)
        {
            if(mFieldData(i).n_rows() != mNumRows || mFieldData(i).n_cols() != mNumCols)
            {
                return false;
            }
        }

        return tValid;
    }

};
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MATRIX_FIELD_INFO_HPP_ */

