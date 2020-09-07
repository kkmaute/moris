#include "cl_GEN_Field_Discrete.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Discrete::Field_Discrete(
                uint aNumOriginalNodes)
        : mNumOriginalNodes(aNumOriginalNodes)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Discrete::evaluate_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            if (aNodeIndex < mNumOriginalNodes)
            {
                return this->evaluate_field_value(aNodeIndex);
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

        void Field_Discrete::evaluate_all_sensitivities(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            if (aNodeIndex < mNumOriginalNodes)
            {
                this->evaluate_all_sensitivities(aNodeIndex, aSensitivities);
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                        "A discrete field sensitivity was requested from a node that this field doesn't know. "
                        "Perhaps a child node was not added to this field?");
                mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_field_sensitivity(this, aSensitivities);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete::add_child_node(
                uint                        aNodeIndex,
                std::shared_ptr<Child_Node> aChildNode)
        {
            MORIS_ASSERT(aNodeIndex == mNumOriginalNodes + mChildNodes.size(),
                    "Child nodes must be added to a level set field in order by node index.");
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete::reset_child_nodes()
        {
            mChildNodes.resize(0);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
