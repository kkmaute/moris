#ifndef SRC_MTK_BSPLINE_FIELD_HPP_
#define SRC_MTK_BSPLINE_FIELD_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        class BSpline_Field : public Field
        {
        private:
            std::shared_ptr<mtk::Mesh_Manager> mMeshManager;
            uint mMeshIndex;
            Matrix< DDRMat > mNodalValues;
            Matrix< DDRMat > mCoefficients;

        public:

            BSpline_Field(
                    std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                    uint                               aMeshIndex,
                    uint                               aDiscretizationMeshIndex = 0 );

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

            // FIXME
            void evaluate_nodal_values();

            // FIXME
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > get_mesh_pair();

            // FIXME
            Matrix<DDRMat>& get_nodal_values()
            {
                return mNodalValues;
            }

            //------------------------------------------------------------------------------

            void set_coefficients(Matrix<DDRMat> aCoefficients);

            //------------------------------------------------------------------------------

            void transfer_coefficients(const BSpline_Field& aField);

            // ----------------------------------------------------------------------------------------------

        };
    }
}
#endif /* SRC_MTK_BSPLINE_FIELD_HPP_ */
