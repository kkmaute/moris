#include "cl_GEN_Discrete_Level_Set.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Level_Set::Discrete_Level_Set(moris::mtk::Interpolation_Mesh* aMesh,
                                               moris::Cell<std::string> const & aFieldNames,
                                               EntityRank aEntityRank,
                                               sint aNumRefinements,
                                               sint aRefinementFunctionIndex)
                : Field(Matrix<DDRMat>(1, 1, 0.0)),
                  Geometry(aNumRefinements, aRefinementFunctionIndex),
                  mFieldNames(aFieldNames),
                  mMesh(aMesh),
                  mEntityRank(aEntityRank),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Level_Set::evaluate_field_value(uint aEntityIndex)
        {
            if (aEntityIndex < mNumOriginalNodes)
            {
                std::string tActiveFieldName = mFieldNames(mActiveFieldIndex);
                return mMesh->get_entity_field_value_real_scalar({{moris_index(aEntityIndex)}}, tActiveFieldName, mEntityRank)(0); // TODO implement getting just real value from all meshes
            }
            else
            {
                MORIS_ASSERT((aEntityIndex - mNumOriginalNodes) < mChildNodes.size(),
                        "A discrete level set field value was requested from a node that this field doesn't know. "
                        "Perhaps a child node was not added to this field?");
                return mChildNodes(aEntityIndex - mNumOriginalNodes).interpolate_geometry_field_value(this);
            }
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

        void Discrete_Level_Set::add_child_node(uint aNodeIndex, Child_Node aChildNode)
        {
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}