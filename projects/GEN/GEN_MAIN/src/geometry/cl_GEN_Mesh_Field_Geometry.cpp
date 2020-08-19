#include "cl_GEN_Mesh_Field_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Field_Geometry::Mesh_Field_Geometry(mtk::Mesh*  aMesh,
                                                 std::string aFieldName,
                                                 EntityRank  aEntityRank,
                                                 sint        aNumRefinements,
                                                 sint        aRefinementFunctionIndex)
                : Field(Matrix<DDRMat>(1, 1, 0.0), aNumRefinements, aRefinementFunctionIndex, -1, -1.0, 1.0),
                  Field_Discrete(aMesh->get_num_nodes()),
                  mMesh(aMesh),
                  mFieldName(aFieldName),
                  mEntityRank(aEntityRank)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Mesh_Field_Geometry::evaluate_field_value(uint aNodeIndex)
        {
            return mMesh->get_entity_field_value_real_scalar({{moris_index(aNodeIndex)}}, mFieldName, mEntityRank)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mesh_Field_Geometry::evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "evaluate_sensitivity function is not implemented for a mesh field geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}