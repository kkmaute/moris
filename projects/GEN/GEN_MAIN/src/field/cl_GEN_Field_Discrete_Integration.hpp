#ifndef MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete_Integration : virtual public Field
        {

        private:
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:

            /**
             * Constructor
             *
             * @param aNumOriginalNodes Number of original nodes on the base integration mesh
             */
            Field_Discrete_Integration(uint aNumOriginalNodes);

            Field_Discrete_Integration()
            {};

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, returns the field value
             *
             * @param aNodeIndex Node index
             * @return Field value
             */
            virtual real get_field_value(uint aNodeIndex) = 0;

            /**
             * Given a node index or coordinate, returns a matrix all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_dfield_dadvs(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, returns a vector of the field derivatives with respect to its ADVs.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            virtual const Matrix<DDRMat>& get_dfield_dadvs(uint aNodeIndex) = 0;

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations, including child nodes.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations for non-child nodes.
             *
             * @param aNodeIndex Node index
             * @return Determining ADV IDs at this node
             */
            virtual Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

            /**
             * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            void get_dfield_dcoordinates(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates,
                    Matrix<DDRMat>&       aSensitivities);

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
            void reset_nodal_data();

        };
    }
}


#endif //MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
