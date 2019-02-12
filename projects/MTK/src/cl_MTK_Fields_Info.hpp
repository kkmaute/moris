/*
 * cl_MTK_Fields_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
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


    moris::Cell<Scalar_Field_Info<DDRMat>*> mRealScalarFields;
    moris::Cell<Matrix_Field_Info<DDRMat>*> mRealMatrixFields;
    moris::Cell<Scalar_Field_Info<DDSMat>*> mSintScalarFields;
    moris::Cell<Matrix_Field_Info<DDSMat>*> mSintMatrixFields;
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
