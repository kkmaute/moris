#include "cl_GEN_Mesh_Field_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Field_Geometry::Mesh_Field_Geometry(moris::mtk::Interpolation_Mesh* aMesh,
                                                 std::string aFieldName,
                                                 EntityRank aEntityRank,
                                                 sint aNumRefinements,
                                                 sint aRefinementFunctionIndex)
                : Field(Matrix<DDRMat>(1, 1, 0.0)),
                  Geometry(aNumRefinements, aRefinementFunctionIndex),
                  mFieldName(aFieldName),
                  mMesh(aMesh),
                  mEntityRank(aEntityRank),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Mesh_Field_Geometry::evaluate_field_value(uint aNodeIndex)
        {
            if (aNodeIndex < mNumOriginalNodes)
            {
                return mMesh->get_entity_field_value_real_scalar({{moris_index(aNodeIndex)}}, mFieldName, mEntityRank)(0);
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                        "A discrete level set field value was requested from a node that this field doesn't know. "
                        "Perhaps a child node was not added to this field?");
                return mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_geometry_field_value(this);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mesh_Field_Geometry::evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "evaluate_sensitivity function is not implemented in level set mesh");
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Mesh_Field_Geometry::sensitivities_available()
        {
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mesh_Field_Geometry::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}