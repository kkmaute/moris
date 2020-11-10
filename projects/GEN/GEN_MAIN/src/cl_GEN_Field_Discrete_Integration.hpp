#ifndef MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete_Integration : virtual public Field
        {
        protected:
            uint mNumOriginalNodes;

        private:
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Discrete_Integration(uint aNumOriginalNodes);

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
            virtual real get_field_value(uint aNodeIndex) = 0;

            /**
             * Given a node index or coordinate, returns a matrix all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, returns a matrix of all sensitivities.
             *
             * @param aNodeIndex Node index
             * @return Matrix of sensitivities
             */
            virtual const Matrix<DDRMat>& get_field_sensitivities(uint aNodeIndex) = 0;

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
