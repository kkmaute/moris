#ifndef MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP
#define MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace ge
    {
        class Mesh_Field_Geometry: public Geometry, public Field_Discrete_Integration
        {

        private:
            mtk::Mesh* mMesh;
            std::string mFieldName;
            EntityRank mEntityRank;

        public:
            /**
             * Constructor
             *
             * @param aMesh Mesh with the level set fields
             * @param aFieldNames Names of the fields
             */
            Mesh_Field_Geometry(mtk::Mesh*  aMesh,
                                std::string aFieldName,
                                EntityRank  aEntityRank = EntityRank::NODE,
                                Geometry_Field_Parameters aParameters = {});

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Distance to this geometry
             */
            real get_field_value(uint aNodeIndex);

        private:

            /**
             * Given a node index, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables. This is currently not implemented for a mesh field geometry.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_dfield_dadvs(uint aNodeIndex);

        };
    }
}

#endif /* MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP */
