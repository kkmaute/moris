/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Fields_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Matrix_Field_Info.hpp"

namespace moris
{
namespace mtk
{

struct MtkFieldsInfo
{
    MtkFieldsInfo():
        mRealScalarFields(0),
        mRealMatrixFields(0),
        mSintScalarFields(0),
        mSintMatrixFields(0)
    {

    }

    uint
    get_num_real_scalar_fields() const
    {
        return mRealScalarFields.size();
    }

    uint
    get_num_real_matrix_fields() const
    {
        return mRealMatrixFields.size();
    }

    uint
    get_num_sint_scalar_fields() const
    {
        return mSintScalarFields.size();
    }

    uint
    get_num_sint_matrix_fields() const
    {
        return mSintMatrixFields.size();
    }

    uint
    get_num_fields() const
    {
        uint tNumFields = get_num_real_scalar_fields()+get_num_real_matrix_fields()+get_num_sint_scalar_fields()+get_num_sint_matrix_fields();
        return tNumFields;
    }

    void
    clear_fields()
    {
        mRealScalarFields.clear();
        mRealMatrixFields.clear();
        mSintScalarFields.clear();
        mSintMatrixFields.clear();
    }

    void
    combine_fields_info( MtkFieldsInfo & aOtherFields)
    {
        mRealScalarFields.append(aOtherFields.mRealScalarFields);
        mRealMatrixFields.append(aOtherFields.mRealMatrixFields);
        mSintScalarFields.append(aOtherFields.mSintScalarFields);
        mSintMatrixFields.append(aOtherFields.mSintMatrixFields);
    }

    moris::Cell<Scalar_Field_Info<DDRMat>*> mRealScalarFields;
    moris::Cell<Matrix_Field_Info<DDRMat>*> mRealMatrixFields;
    moris::Cell<Scalar_Field_Info<DDSMat>*> mSintScalarFields;
    moris::Cell<Matrix_Field_Info<DDSMat>*> mSintMatrixFields;

    void
    print()
    {
        std::cout<<"Number of Real Scalar Fields: "<<std::setw(6)<<mRealScalarFields.size()<<std::endl;
        std::cout<<"Number of Real Matrix Fields: "<<std::setw(6)<<mRealMatrixFields.size()<<std::endl;
        std::cout<<"Number of Int  Scaler Fields: "<<std::setw(6)<<mSintScalarFields.size()<<std::endl;
        std::cout<<"Number of Int  Matrix Fields: "<<std::setw(6)<<mSintMatrixFields.size()<<std::endl;

        std::cout<<" Real Scalar Field Information:"<<std::endl;
        for(moris::uint  i = 0 ; i < mRealScalarFields.size(); i++)
        {
            std::cout<<"    Name: "<<mRealScalarFields(i)->get_field_name()<<" | Entity Rank: "<<(uint)mRealScalarFields(i)->get_field_entity_rank()<<std::endl;
        }
    }

};

template<typename Field_Ptr>
inline
void
add_field_for_mesh_input(Field_Ptr aFieldToAdd,
                         MtkFieldsInfo & aFieldsInfo)
{
    MORIS_ERROR(false, "Field type specified is not supported");
}

inline
void
add_field_for_mesh_input(Scalar_Field_Info<DDRMat>* aField,
                         MtkFieldsInfo &            aFieldsInfo)
{
    aFieldsInfo.mRealScalarFields.push_back(aField);
}

}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_ */

