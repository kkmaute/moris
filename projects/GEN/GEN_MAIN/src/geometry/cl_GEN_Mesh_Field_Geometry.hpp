#ifndef MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP
#define MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace ge
    {
        class Mesh_Field_Geometry: public Geometry, public Field_Discrete
        {

        private:
            std::string mFieldName;
            mtk::Mesh* mMesh;
            EntityRank mEntityRank;
            uint mNumOriginalNodes;
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:
            /**
             * Constructor
             *
             * @param aMesh Mesh with the level set fields
             * @param aFieldNames Names of the fields
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             */
            Mesh_Field_Geometry(mtk::Mesh*  aMesh,
                                std::string aFieldName,
                                EntityRank  aEntityRank = EntityRank::NODE,
                                sint        aNumRefinements = 0,
                                sint        aRefinementFunctionIndex = -1);

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Distance to this geometry
             */
            real evaluate_field_value(uint aNodeIndex);

            /**
             * Given a node index, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables. This is currently not implemented for a mesh field geometry.
             *
             * @param aNodeIndex Node index
             * @param aSensitivities Vector of sensitivities
             */
            void evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities);

            /**
             * Add a new child node for evaluation
             *
             * @param aNodeIndex Index of the child node
             * @param aChildNode Contains information about how the child node was created
             */
            void add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode);

        };
    }
}

#endif /* MORIS_CL_GEN_MESH_FIELD_GEOMETRY_HPP */
