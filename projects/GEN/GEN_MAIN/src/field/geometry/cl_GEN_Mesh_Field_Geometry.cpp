#include "cl_GEN_Mesh_Field_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Field_Geometry::Mesh_Field_Geometry(mtk::Mesh*  aMesh,
                                                 std::string aFieldName,
                                                 EntityRank  aEntityRank,
                                                 Geometry_Field_Parameters aParameters)
                : Field(Matrix<DDRMat>(1, 1, 0.0), aParameters)
                , Geometry(aParameters)
                , Field_Discrete_Integration(aMesh->get_num_nodes())
                , mMesh(aMesh)
                , mFieldName(aFieldName)
                , mEntityRank(aEntityRank)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Mesh_Field_Geometry::get_field_value(uint aNodeIndex)
        {
            return mMesh->get_entity_field_value_real_scalar({{moris_index(aNodeIndex)}}, mFieldName, mEntityRank)(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Mesh_Field_Geometry::get_dfield_dadvs(uint aNodeIndex)
        {
            MORIS_ERROR(false, "get_dfield_dadvs function is not implemented for a mesh field geometry.");
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mesh_Field_Geometry::get_dfield_dcoordinates(
                uint            aNodeIndex,
                Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for mesh field geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
