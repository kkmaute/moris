#ifndef SRC_MTK_BSPLINE_FIELD_HPP_
#define SRC_MTK_BSPLINE_FIELD_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "st_MTK_Mesh_Pair.hpp"

namespace moris
{
    namespace mtk
    {
        class BSpline_Field : public Field
        {
        private:
            Mesh_Pair      mMeshPair;
            Matrix<DDRMat> mCoefficients;
            Matrix<DDRMat> mNodalValues;

        public:

            BSpline_Field(
                    Mesh_Pair aMeshPair,
                    uint      aDiscretizationMeshIndex = 0);

            ~BSpline_Field();

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
             * Sets the B-spline coefficients which interpolate this field.
             *
             * @param aCoefficients B-spline coefficients
             */
            void set_coefficients(const Matrix<DDRMat>& aCoefficients);

            /**
             * Gets the nodal values of this B-spline field.
             *
             * @return Nodal values
             */
            const Matrix<DDRMat>& get_nodal_values();

            /**
             * Gets the mesh pair associated with this field.
             *
             * @return Mesh pair
             */
            Mesh_Pair get_mesh_pair();

        };
    }
}
#endif /* SRC_MTK_BSPLINE_FIELD_HPP_ */
