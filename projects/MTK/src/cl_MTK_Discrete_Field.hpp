#ifndef SRC_MTK_DISCRETE_FIELD_HPP_
#define SRC_MTK_DISCRETE_FIELD_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "st_MTK_Mesh_Pair.hpp"

namespace moris
{
    namespace mtk
    {
        class Discrete_Field : public Field
        {
        private:
            Mesh_Pair      mMeshPair;
            Matrix<DDRMat> mNodalValues;

        public:

            Discrete_Field(
                    Mesh_Pair aMeshPair,
                    uint      aDiscretizationMeshIndex = 0);

            ~Discrete_Field();

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
             * Sets the nodal values of this discrete field.
             *
             * @param aNodalValues Nodal values
             */
            void set_nodal_values(const Matrix<DDRMat>& aNodalValues);

            /**
             * Gets the nodal values of this discrete field.
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
#endif /* SRC_MTK_DISCRETE_FIELD_HPP_ */
