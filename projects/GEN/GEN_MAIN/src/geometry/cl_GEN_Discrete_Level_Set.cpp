//
// Created by christopherson on 4/16/20.
//

#include "cl_GEN_Discrete_Level_Set.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Level_Set::Discrete_Level_Set(moris::mtk::Interpolation_Mesh* aMesh,
                moris::Cell<std::string> const & aFieldNames, EntityRank aEntityRank)
                : Geometry_Discrete(Matrix<DDRMat>(1, 1, 0.0)),
        {
            mFieldNames = aFieldNames;
            mMesh = aMesh;
            mEntityRank = aEntityRank;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Level_Set::evaluate_field_value(moris_index aEntityIndex);
        {
            std::string tActiveFieldName = mLevelSetFieldNames(mActiveFieldIndex);
            return mLevelSetMesh->get_entity_field_value_real_scalar({{aEntityIndex}}, tActiveFieldName, mEntityRank;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Discrete_Level_Set::evaluate_sensitivity(moris_index aEntityIndex);
        {
            MORIS_ERROR(false, "evaluate_sensitivity_dx_dp function is not implemented in level set mesh"); //TODO: Implement this function

            moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
            return tSensitivityDxDp;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Discrete_Level_Set::set_active_field(size_t aActiveFieldIndex)
        {
            mActiveFieldIndex = aActiveFieldIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}