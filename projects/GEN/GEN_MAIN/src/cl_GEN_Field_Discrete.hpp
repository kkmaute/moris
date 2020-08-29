#ifndef MORIS_CL_GEN_FIELD_DISCRETE_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete : virtual public Field
        {

        protected:
            uint mNumOriginalNodes;

        private:
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Discrete(uint aNumOriginalNodes);

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real evaluate_field_value(uint                  aNodeIndex,
                                      const Matrix<DDRMat>& aCoordinates);

            /**
             * Add a new child node for evaluation
             *
             * @param aNodeIndex Index of the child node
             * @param aChildNode Contains information about how the child node was created
             */
            void add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode);

            /**
             * Resets all child nodes, called when a new XTK mesh is being created.
             */
            void reset_child_nodes();

        private:

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aCoordinates vector of coordinate values
             * @return Field value
             */
            virtual real evaluate_field_value(uint aNodeIndex) = 0;

            /**
             * Given a node coordinate, returns a matrix of all sensitivities
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            void evaluate_all_sensitivities(uint                  aNodeIndex,
                                            const Matrix<DDRMat>& aCoordinates,
                                            Matrix<DDRMat>&       aSensitivities);

            /**
             * Given a node index, returns a matrix of all sensitivities
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            virtual void evaluate_all_sensitivities(uint            aNodeIndex,
                                                    Matrix<DDRMat>& aSensitivities) = 0;

        };
    }
}


#endif //MORIS_CL_GEN_FIELD_DISCRETE_HPP
