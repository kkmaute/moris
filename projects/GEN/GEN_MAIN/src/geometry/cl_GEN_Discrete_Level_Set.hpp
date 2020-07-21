#ifndef MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP
#define MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace ge
    {
        class Discrete_Level_Set: public Geometry, public Field_Discrete
        {
        private:
            size_t mActiveFieldIndex = 0;
            Cell<std::string> mFieldNames;
            moris::mtk::Interpolation_Mesh* mMesh;
            EntityRank mEntityRank;
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
            Discrete_Level_Set(moris::mtk::Interpolation_Mesh* aMeshWithLevelSetFields,
                               moris::Cell<std::string> const & aFieldNames,
                               EntityRank aEntityRank = EntityRank::NODE,
                               sint aNumRefinements = 0,
                               sint aRefinementFunctionIndex = -1);

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aEntityIndex the index of the field value
             * @return field value at the specified index
             */
            real evaluate_field_value(uint aEntityIndex);

            /**
             * Given an index, returns sensitivities of the geometry with respect to input parameters
             *
             * @param aEntityIndex the index of the field value
             * @param aSensitivities Matrix of sensitivities to be returned
             */
            void evaluate_all_sensitivities(uint aEntityIndex, Matrix<DDRMat>& aSensitivities);

            /**
             * Sets the active field to the given value
             *
             * @param aActiveFieldIndex The index of the desired field to become active
             */
            void set_active_field(size_t aActiveFieldIndex);

            /**
             * Lets the geometry engine know if sensitivities are available, otherwise it will perform finite
             * differencing instead for intersection locations
             *
             * @return If sensitivities are implemented or not (false for discrete level set)
             */
            bool sensitivities_available();

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

#endif /* MORIS_CL_GEN_DISCRETE_LEVEL_SET_HPP */
