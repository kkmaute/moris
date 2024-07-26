/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field_Discrete.hpp
 *
 */

#ifndef SRC_MTK_FIELD_DISCRETE_HPP_
#define SRC_MTK_FIELD_DISCRETE_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace mtk
    {
        /* child class implementation of mtk::Field. The discrete field class handles cases where the Lagrange interpolation mesh
         * has an underlying discretization such as B-splines. This class provides a convenient interface to access the coefficient
         * of this discretization and for computing nodal values.
         *
         * IMPORTANT: Not every node has access to the underlying discretization, such as T-matrices of aura nodes. Therefore, when
         * computing nodal values a communication happens across processors. This is implemented via a distriuted vector. Furthermore,
         * the derivatives of a nodal value with respect of coefficients returns 0 if the node does not have access to the underlying
         * discretization.
         */

        class Field_Discrete : public Field
        {
          protected:
            //! Map from mesh indices to field coefficient indices
            Matrix< IndexMat > mMeshToFieldCoefficientIndexMap;

            //! Discretization Index
            moris_index mDiscretizationMeshIndex = -1;

            //! Maximum number of coefficients used by mesh
            uint mMaxNumberOfCoefficients = 0;

            //! Flag whether coefficient vector is initialized
            bool mCoefficientsAreInitialized = false;

            //! Distributed vector of owned nodal values
            sol::Dist_Vector* mOwnedNodalValues = nullptr;

            //! distributed vector of shared nodal values
            sol::Dist_Vector* mSharedNodalValues = nullptr;

          public:
            // ----------------------------------------------------------------------------------------------

            Field_Discrete(
                    mtk::Mesh_Pair aMeshPairs,
                    uint const &   aDiscretizationMeshIndex = 0,
                    uint const &   mNumberOfFields          = 1 );

            // ----------------------------------------------------------------------------------------------

            ~Field_Discrete();

            // ----------------------------------------------------------------------------------------------

            /**
             * @brief child class implementation: updates internal data associated with coefficients
             */

            void update_coefficient_data();

            // ----------------------------------------------------------------------------------------------

            /**
             * @brief child class implementation: returns order of discretization
             */
            uint get_discretization_order() const;

            // ----------------------------------------------------------------------------------------------

            /**
             * @brief child class implementation: returns discretization mesh index
             */
            moris_index get_discretization_mesh_index() const;

            // ----------------------------------------------------------------------------------------------

            /**
             * @brief child class implementation: fills coefficient vector
             *
             * @param[in]  vector of coefficients
             */
            void set_coefficient_vector( const Matrix< DDRMat >& aCoefficients );

            //------------------------------------------------------------------------------

            /**
             * @brief child class implementation: get IDs of used coefficients of underlying discretization
             */
            const Matrix< IdMat >& get_coefficient_id_and_owner_vector();

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
                    Matrix< DDRMat >&   aDerivatives,
                    Matrix< IndexMat >& aCoefIndices,
                    uint const &        aNodeIndex,
                    uint const &        aFieldIndex );

            // ----------------------------------------------------------------------------------------------

          private:
            void communicate_missing_owned_coefficients(
                    Matrix< IdMat >& aAllCoefIds,
                    Matrix< IdMat >& aAllCoefOwners );
        };
    }    // namespace mtk
}    // namespace moris
#endif /* SRC_MTK_FIELD_ANALYTIC_HPP_ */
