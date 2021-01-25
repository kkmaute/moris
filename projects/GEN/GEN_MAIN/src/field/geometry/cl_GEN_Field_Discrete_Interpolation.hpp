#ifndef MORIS_CL_GEN_FIELD_DISCRETE_INTERPOLATION_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_INTERPOLATION_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete_Interpolation : virtual public Geometry
        {
        private:
            mtk::Mesh*                        mMesh;
            std::string                       mFieldName;
            moris::moris_index                mFieldIndex;
            Matrix<DDRMat>                    mFieldSensitivity;
            moris_index                       mNumOriginalNodes;
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            template< typename Vector_Type > 
            Field_Discrete_Interpolation(
                  mtk::Mesh*       aMesh,
                  Vector_Type&     aADVs,
                  Matrix<DDUMat>   aGeometryVariableIndices,
                  Matrix<DDUMat>   aADVIndices,
                  Matrix<DDRMat>   aConstants,
                  Field_Parameters aParameters = {})
                 : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters),
                  mMesh(aMesh),
                  mFieldName(aParameters.mName),
                  mFieldIndex(mMesh->get_field_ind(mFieldName,EntityRank::NODE)),
                  mNumOriginalNodes(mMesh->get_num_entities(EntityRank::NODE))
            {
            };

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aNodeIndex Node index
             * @return Field value
             */
            real get_field_value(uint aNodeIndex)
            {
                MORIS_ERROR(0,"ISSUE");
                return 0;
            };

            /**
             * Given a node index or coordinate, returns a matrix all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& 
            get_field_sensitivities(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, returns a matrix of all sensitivities.
             *
             * @param aNodeIndex Node index
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(uint aNodeIndex);

            

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations, including child nodes.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations for non-child nodes.
             *
             * @param aNodeIndex Node index
             * @return Determining ADV IDs at this node
             */
            virtual Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

            /**
             * Add a new child node for evaluation.
             *
             * @param aNodeIndex Index of the child node
             * @param aChildNode Contains information about how the child node was created
             */
            void add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode);

            /**
             * Resets all child nodes, called when a new XTK mesh is being created.
             */
            void reset_nodal_information();

        };
    }
}


#endif //MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
