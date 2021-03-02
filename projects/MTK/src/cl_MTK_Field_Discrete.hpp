/*
 * cl_MTK_Field_Discrete.hpp
 */
#ifndef SRC_MTK_FIELD_DISCRETE_HPP_
#define SRC_MTK_FIELD_DISCRETE_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        class Field_Discrete : public Field
        {
            protected:

                //! Discretization Index
                moris_index mDiscretizationMeshIndex = -1;

                //! Flag whether coefficient vector is initialized
                bool mCoefficientsAreInitialized = false;

            public :

                // ----------------------------------------------------------------------------------------------

                Field_Discrete()
                {};

                // ----------------------------------------------------------------------------------------------

                Field_Discrete(
                        mtk::Mesh_Pair * aMeshPairs,
                        uint     const & aDiscretizationMeshIndex = 0,
                        uint     const & mNumberOfFields = 1);

                // ----------------------------------------------------------------------------------------------

                ~Field_Discrete()
                {};

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief child class implementation:returns order of discretization
                 */
                uint get_discretization_order() const;

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief child class implementation: returns discretization mesh index
                 */
                moris_index get_discretization_mesh_index() const;

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief cchild class implementation: fills coefficient vector
                 *
                 * @param[in]  vector of coefficients
                 */
                void set_coefficient_vector(const Matrix< DDRMat > & aCoefficients);


                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief child class implementation: computes and stores nodal values
                 */
                virtual void compute_nodal_values();

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief child class implementation: computes derivatives of nodal values
                 */
                virtual void compute_derivatives_of_field_value(
                        Matrix< DDRMat >       & aDerivatives,
                        Matrix< IndexMat >     & aCoefIndices,
                        uint             const & aNodeIndex,
                        uint             const & aFieldIndex);

                // ----------------------------------------------------------------------------------------------

        };
    }
}
#endif /* SRC_MTK_FIELD_ANALYTIC_HPP_ */
