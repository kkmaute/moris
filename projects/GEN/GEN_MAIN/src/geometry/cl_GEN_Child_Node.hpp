#ifndef MORIS_CL_GEN_CHILD_NODE_HPP
#define MORIS_CL_GEN_CHILD_NODE_HPP

#include "cl_Matrix.hpp"
#include "cl_XTK_Basis_Function.hpp"

namespace moris
{
    namespace ge
    {
        // Forward declaration of Field class
        class Field;

        class Child_Node
        {

        protected:
            Matrix<DDUMat>       mParentNodeIndices;
            Cell<Matrix<DDRMat>> mParentNodeCoordinates;
            Matrix<DDRMat>       mBasisValues;

        private:
            Matrix<DDRMat>       mLocalCoordinates;

        public:
            Child_Node(){};
            /**
             * Constructor
             *
             * @param aParentNodeIndices Node indices of the parent of this child node
             * @param aParentNodeCoordinates Coordinates of the parent of this child node
             * @param aBasisFunction Basis function of the parent topology
             * @param aLocalCoordinates Local coordinate of this child inside of the parent element
             */
            Child_Node(Matrix<DDUMat>             aParentNodeIndices,
                       Cell<Matrix<DDRMat>>       aParentNodeCoordinates,
                       const xtk::Basis_Function& aBasisFunction,
                       Matrix<DDRMat>             aLocalCoordinates);

            /**
             * Gets the local coordinates of this child node.
             *
             * @return Local coordinates
             */
            Matrix<DDRMat> get_local_coordinates();

            /**
             * Get the field value on the child node based on values from its parents.
             *
             * @param aField Field pointer, referenced during call from field class
             * @return Field value
             */
            virtual real interpolate_field_value(Field* aField);

            /**
             * Get the field sensitivities on the child node based on its parents.
             *
             * @param aField Field pointer, referenced during call from field class
             * @return Field sensitivities
             */
            virtual Matrix<DDRMat> interpolate_field_sensitivity(Field* aField);
        };
    }
}

#endif //MORIS_CL_GEN_CHILD_NODE_HPP
