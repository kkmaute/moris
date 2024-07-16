/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field_Analytic.hpp
 *
 */

#ifndef SRC_MTK_FIELD_ANALYTIC_HPP_
#define SRC_MTK_FIELD_ANALYTIC_HPP_

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    // User-defined field functions
    typedef real ( *Analytic_Field_Function ) (
            const moris::Matrix< DDRMat > & aCoordinates,
            const moris::Matrix< DDRMat > & aParameters );   // evaluates function given coordinates and coefficients
    typedef void ( *Analytic_Derivative_Function ) (
            const moris::Matrix< DDRMat > & aCoordinates,
            const moris::Matrix< DDRMat > & aParameters,
            moris::Matrix< DDRMat >       & aReturnValue);   // evaluates derivative of function given coordinates and coefficients

    namespace mtk
    {
        /**
         * @ brief Implementation of mesh based nodally discretized scalar and vector field, based on mtk::Field class.
         *
         * For each field, a function needs to be provided that evaluates the field at a node given the nodal
         * coordinates and a coefficient vector. In addition, a function needs to be provided that evaluates the
         * derivatives of nodal field value with respect to all parameters. Finally, the number of coefficients the value
         * and derivative functions are operating on needs to be provided.
         *
         * Note: the value and derivative function operate on the same coefficient vector !
         * */

        class Field_Analytic : public Field
        {
            protected:
                Vector<Analytic_Field_Function>      mAnalyticFieldValueFunction;
                Vector<Analytic_Derivative_Function> mAnalyticDerivativeFunction;

            public :

                Field_Analytic(
                        Vector<Analytic_Field_Function>       aFunction,
                        Vector<Analytic_Derivative_Function>  aDerivativeFunction,
                        moris::Matrix<DDRMat>              const & aCoefficients,
                        Mesh_Pair                                  aMeshPairs,
                        uint                               const & aNumberOfFields = 1);

                // ----------------------------------------------------------------------------------------------

                Field_Analytic(
                        Analytic_Field_Function        aFunction,
                        Analytic_Derivative_Function   aDerivativeFunction,
                        moris::Matrix<DDRMat>  const & aCoefficients,
                        Mesh_Pair                      aMeshPairs,
                        uint                   const & aNumberOfFields = 1);

                // ----------------------------------------------------------------------------------------------

                ~Field_Analytic();

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

