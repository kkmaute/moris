#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Analytic::Field_Analytic()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Analytic::get_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            if( !mInterpolateChildNodes or !mMesh )
            {
                return this->get_field_value_geometry( aNodeIndex, aCoordinates );
            }
            else
            {
                if (aNodeIndex < mNumOriginalNodes)
                {
                    return this->get_field_value_geometry( aNodeIndex, aCoordinates );
                }
                else
                {
                    MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                            "A analytical field value was requested from a node that this field doesn't know. "
                            "Perhaps a child node was not added to this field?");
                    return mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_field_value(this);
                }
            }


        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Field_Analytic::get_field_sensitivities(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_sensitivities(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Analytic::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
            if( mMesh )
            {
                MORIS_ASSERT(aNodeIndex == mNumOriginalNodes + mChildNodes.size(),
                        "Child nodes must be added to a level set field in order by node index.");
                mChildNodes.push_back(aChildNode);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Analytic::reset_child_nodes()
        {
            mChildNodes.resize(0);
        }

    }
}
