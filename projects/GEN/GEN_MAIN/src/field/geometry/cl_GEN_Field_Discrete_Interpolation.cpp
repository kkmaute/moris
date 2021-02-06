#include "cl_GEN_Field_Discrete_Interpolation.hpp"

#include "cl_MTK_Mesh.hpp"
namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Field_Discrete_Interpolation::get_field_value(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            if ((moris_index)aNodeIndex < mNumOriginalNodes)
            {
                return mMesh->get_entity_field_value_real_scalar({{(moris_index)mMesh->get_base_node_index(aNodeIndex)}},mFieldName,EntityRank::NODE)(0);
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                        "A discrete field value was requested from a node that this field doesn't know. "
                        "Perhaps a child node was not added to this field?");
                return mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_field_value(this);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Field_Discrete_Interpolation::get_field_sensitivities(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            if ((moris_index)aNodeIndex < mNumOriginalNodes)
            {
                return this->get_field_sensitivities(aNodeIndex);
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                        "A discrete field sensitivity was requested from a node that this field doesn't know. "
                        "Perhaps a child node was not added to this field?");
                return mChildNodes(aNodeIndex - mNumOriginalNodes)->join_field_sensitivities(this);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_determining_adv_ids(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            if ((moris_index)aNodeIndex < mNumOriginalNodes)
            {
                return this->get_determining_adv_ids(aNodeIndex);
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                             "A discrete field sensitivity was requested from a node that this field doesn't know. "
                             "Perhaps a child node was not added to this field?");
                return mChildNodes(aNodeIndex - mNumOriginalNodes)->join_determining_adv_ids(this);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_determining_adv_ids(uint aNodeIndex)
        {
            return Field::get_determining_adv_ids(aNodeIndex, {{}});
        }

        //--------------------------------------------------------------------------------------------------------------
        const Matrix<DDRMat>& 
        Field_Discrete_Interpolation::get_field_sensitivities(
            uint                  aNodeIndex)
        {
            return mFieldSensitivity;
        }
        //--------------------------------------------------------------------------------------------------------------
        
        void Field_Discrete_Interpolation::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
            MORIS_ASSERT(aNodeIndex == mNumOriginalNodes + mChildNodes.size(),
                    "Child nodes must be added to a level set field in order by node index.");
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete_Interpolation::reset_nodal_information()
        {
            mChildNodes.resize(0);
        }
    }
}
