#ifndef MORIS_CL_GEN_LEVEL_SET_HPP
#define MORIS_CL_GEN_LEVEL_SET_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Level_Set : public Geometry, public Field_Discrete
        {

        private:
            mtk::Mesh* mMesh;
            uint mNumOriginalNodes;
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:
            /**
             * Constructor
             *
             * @param aMeshWithLevelSetFields Mesh with the level set fields
             * @param aFieldNames Names of the fields
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             */
            Level_Set(Matrix<DDRMat>& aADVs,
                      Matrix<DDUMat>  aGeometryVariableIndices,
                      Matrix<DDUMat>  aADVIndices,
                      Matrix<DDRMat>  aConstantParameters,
                      mtk::Mesh*      aMesh,
                      sint            aNumRefinements = 0,
                      sint            aRefinementFunctionIndex = -1);

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aNodeIndex Node index for field evaluation
             * @return field value at the specified index
             */
            real evaluate_field_value(uint aNodeIndex);

            /**
             * Given an index, returns sensitivites of the geometry with respect to input parameters
             *
             * @param aNodeIndex Node index for field evaluation
             * @param aSensitivities Matrix of sensitivities to be returned
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

#endif //MORIS_CL_GEN_LEVEL_SET_HPP
