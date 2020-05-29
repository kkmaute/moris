#include "cl_GEN_Discrete_Level_Set.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Level_Set::Discrete_Level_Set(moris::mtk::Interpolation_Mesh* aMesh,
                moris::Cell<std::string> const & aFieldNames, EntityRank aEntityRank)
                : Field(Matrix<DDRMat>(1, 1, 0.0))
        {
            mFieldNames = aFieldNames;
            mMesh = aMesh;
            mEntityRank = aEntityRank;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Level_Set::evaluate_field_value(uint aEntityIndex)
        {
            std::string tActiveFieldName = mFieldNames(mActiveFieldIndex);
            return mMesh->get_entity_field_value_real_scalar({{moris_index(aEntityIndex)}}, tActiveFieldName, mEntityRank)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Discrete_Level_Set::evaluate_all_sensitivities(uint aEntityIndex, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "evaluate_sensitivity function is not implemented in level set mesh"); //TODO: Implement this function
        }

        //--------------------------------------------------------------------------------------------------------------

        void Discrete_Level_Set::set_active_field(size_t aActiveFieldIndex)
        {
            mActiveFieldIndex = aActiveFieldIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Discrete_Level_Set::sensitivities_available()
        {
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}