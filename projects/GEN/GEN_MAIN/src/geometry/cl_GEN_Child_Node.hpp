#ifndef MORIS_CL_GEN_CHILD_NODE_HPP
#define MORIS_CL_GEN_CHILD_NODE_HPP

#include "cl_Matrix.hpp"
#include "cl_XTK_Basis_Function.hpp"

namespace moris
{
    namespace ge
    {
        // Forward declaration of Geometry class
        class Geometry;

        class Child_Node
        {
        private:
            Matrix<DDUMat>       mParentNodeIndices;
            Cell<Matrix<DDRMat>> mParentNodeCoordinates;
            Matrix<DDRMat>       mBasisValues;

        public:

            /**
             * Constructor
             *
             * @param aParentNodeIndices Node indices of the parent of this child node
             * @param aParentNodeCoordinates Coordinates of the parent of this child node
             * @param aBasisFunction Basis function of the parent topology
             * @param aIntersectionGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of the child node
             */
            Child_Node(Matrix<DDUMat>             aParentNodeIndices,
                       Cell<Matrix<DDRMat>>       aParentNodeCoordinates,
                       const xtk::Basis_Function& aBasisFunction,
                       Matrix<DDRMat>             aLocalCoordinates);

            /**
             * Get a geometry field value on the child node based on values from the parents
             *
             * @param aGeometry Geometry pointer, called from inside a geometry
             */
            virtual real interpolate_geometry_field_value(Geometry* aGeometry);
        };
    }
}

#endif //MORIS_CL_GEN_CHILD_NODE_HPP
