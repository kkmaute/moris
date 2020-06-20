#ifndef MORIS_CL_GEN_Field_HPP
#define MORIS_CL_GEN_Field_HPP

#include "cl_Matrix.hpp"
#include "cl_GEN_Child_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Field
        {
        protected:
            Cell<real*> mFieldVariables;

        private:
            Matrix<DDUMat> mADVIndices;
            Matrix<DDRMat> mConstantParameters;
            Cell<bool> mActiveVariables;
            uint mNumActiveVariables;

        protected:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aFieldVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Field(Matrix<DDRMat>& aADVs,
                  Matrix<DDUMat> aFieldVariableIndices,
                  Matrix<DDUMat> aADVIndices,
                  Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this field
             */
            Field(Matrix<DDRMat> aConstantParameters);

        public:
	
	Field(){};
            /**
             * Destructor
             */
            ~Field();

            /**
             * Get the indices of the ADVs which this field depends on
             *
             * @return Vector of ADV indices
             */
            Matrix<DDUMat> get_adv_indices();

            /**
             * Given a node index or coordinate, returns the field value
             *
             * @param aIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            virtual real evaluate_field_value(      uint            aIndex,
                                              const Matrix<DDRMat>& aCoordinates) = 0;

            /**
             * Given a node index or coordinate, returns a matrix of relevant sensitivities
             *
             * @param aIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            void evaluate_sensitivity(      uint            aIndex,
                                      const Matrix<DDRMat>& aCoordinates,
                                            Matrix<DDRMat>& aSensitivities);

            /**
             * Add a new child node for evaluation, implemented for discrete fields
             *
             * @param aNodeIndex Index of the child node
             * @param aChildNode Contains information about how the child node was created
             */
            virtual void add_child_node(uint aNodeIndex, Child_Node aChildNode);

        private:

            /**
             * Given a node coordinate @param aCoordinates, the function returns a matrix of all sensitivities
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            virtual void evaluate_all_sensitivities(      uint            aIndex,
                                                    const Matrix<DDRMat>& aCoordinates,
                                                          Matrix<DDRMat>& aSensitivities) = 0;

        };
    }
}


#endif //MORIS_CL_GEN_Field_HPP
