/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bspline.hpp
 *
 */

#ifndef SRC_MODEL_CL_BSPLINE_HPP_
#define SRC_MODEL_CL_BSPLINE_HPP_

// MORIS library header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Mat.hpp"

namespace moris
{
    class Bspline
    {
    protected:

    public:

        /**
         * Bspline constructor
         */
        Bspline()
        {
        }
        /**
         * Bspline destructor.
         */
        ~Bspline() = default;

        /**
         * Creates the basis functions for a polynomial degree for 1D.
         *
         * @param[in] aXi                 Natural coordinate.
         * @param[in] aKnotVector         Knot vector of the b-spline.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         *
         * @param[out] Bspline            Basis functions of the b-spline.
         *
         */
        static Mat<real>
        build_spline_1d(real const & aXi,
                Mat<real> & aKnotVector,
                uint const & aPolynomialDegree);

        /**
         * Creates the derivative of the basis functions for a polynomial degree for 1D.
         *
         * @param[in] aXi                 Natural coordinate.
         * @param[in] aKnotVector         Knot vector of the b-spline.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         *
         * @param[out] Bspline_deriv      Derivative of the b-spline basis functions.
         *
         */
        static Mat<real>
        build_spline_deriv_1d(real const & aXi,
                Mat<real> & aKnotVector,
                uint const & aPolynomialDegree);

        /**
         * Creates two or three-dimensional basis functions for a polynomial degree.
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aKnotVector1        Knot vector in xi_1 of the b-spline.
         * @param[in] aKnotVector2        Knot vector in xi_2 of the b-spline.
         * @param[in] aKnotVector3        Knot vector in xi_3 of the b-spline.
         * @param[in] aPolynomialDegree   Matrix with Polynomial degrees of the b-splines for each direction [p_1, p_2, p_3].
         * @param[in] aElement_flag       Which basis functions are active in the respective element for each direction [start_1 start_2 start_3].
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline            Basis functions of the b-spline. For example in 2D and a quadratic approach , the vector creates\n
         *                                .-------------------\n
         *                                .| N7  N8  N9 |\n
         *                                .| N4  N5  N6 |\n
         *                                .| N1  N2  N3 |\n
         *                                .-------------------\n
         *                                In 3D, it creates in each plain of the domain e_1 - e_2 the same order. For example in 3D and a linear approach\n
         *                                Plain 1\n
         *                                .--------------\n
         *                                .| N3  N4 |\n
         *                                .| N1  N2 |\n
         *                                .--------------\n
         *                                Plain 2\n
         *                                .-------------\n
         *                                .| N7  N8 |\n
         *                                .| N5  N6 |\n
         *                                .-------------
         *
         */
        static Mat<real>
        build_spline_nd(Mat<real> const & aXi,
                Mat<real> & aKnotVector1,
                Mat<real> & aKnotVector2,
                Mat<real> & aKnotVector3,
                Mat<uint> const & aPolynomialDegree,
                Mat<uint> const & aElement_flag,
                uint const & aDim);

        /**
         * Creates two or three-dimensional derivative of the basis functions for a polynomial degree.
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aKnotVector1        Knot vector in xi_1 of the b-spline.
         * @param[in] aKnotVector2        Knot vector in xi_2 of the b-spline.
         * @param[in] aKnotVector3        Knot vector in xi_3 of the b-spline.
         * @param[in] aPolynomialDegree   Matrix with Polynomial degrees of the b-splines for each direction [p_1, p_2, p_3].
         * @param[in] aElement_flag       Which basis functions are active in the respective element for each direction [start_1 start_2 start_3].
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline_deriv      Derivative of the basis functions. It creates a Matrix with the derivatives in the rows and the directions in the columns.
         *                                The order is similiar to: build_spline_nd
         *
         */
        static Mat<real>
        build_spline_deriv_nd(Mat<real> const & aXi,
                Mat<real> & aKnotVector1,
                Mat<real> & aKnotVector2,
                Mat<real> & aKnotVector3,
                Mat<uint> const & aPolynomialDegree,
                Mat<uint> const & aElement_flag,
                uint const & aDim);

        /**
         * Creates two or three-dimensional uniform basis functions for a polynomial degree.
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aKnotVector         Knot vector for each direction of the b-splines.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-splines for each direction p_1 = p_2 = p_3.
         * @param[in] aElement_flag       Which basis functions are active in the respective element for each direction [start_1 start_2 start_3].
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline            B-spline basis functions. It creates a Matrix with the derivatives in the columns and the directions in the rows.
         *                                The order is similiar to: build_spline_nd
         */
        static Mat<real>
        build_spline_uniform_nd(Mat<real> const & aXi,
                Mat<real> & aKnotVector,
                uint const & aPolynomialDegree,
                Mat<uint> const & aElement_flag,
                uint const & aDim);

        /**
         * Creates the basis functions for a polynomial degree for 1D and a standard element with the element domain 0 <= xi <= 1
         *
         * @param[in] aXi                 Natural coordinate.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         *
         * @param[out] Bspline            Basis functions of the b-spline.
         *
         */
        static Mat<real>
        build_spline_uniform_1d(real const & aXi,
                uint const & aPolynomialDegree);

        /**
         * Creates the basis functions for a polynomial degree for 1D and a standard element with the element domain 0 <= xi <= 1
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aPolynomialDegree   Polynomial degree of the b-splines for each direction p_1 = p_2 = p_3.
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline            B-spline basis functions. It creates a Matrix with the derivatives in the columns and the directions in the rows.
         *                                The order is similiar to: build_spline_nd
         */

        static Mat<real>
        build_spline_uniform_nd(
                Mat<real> const& aXi,
                uint const & aPolynomialDegree,
                uint const & aDim);

        /**
         * Creates the derivative of the basis functions for a polynomial degree for 1D and a standard element with the element domain 0 <= xi <= 1
         *
         * @param[in] aXi                 Natural coordinate.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         *
         * @param[out] Bspline            Basis functions of the b-spline.
         *
         */
        static Mat<real>
        build_spline_deriv_uniform_1d(
                real const & aXi,
                uint const & aPolynomialDegree);

        /**
         * Creates two or three-dimensional derivative of the basis functions for a polynomial degree.
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aPolynomialDegree   Matrix with Polynomial degrees of the b-splines for each direction [p_1, p_2, p_3].
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline_deriv      Derivative of the basis functions. It creates a Matrix with the derivatives in the rows and the directions in the columns.
         *                                The order is similiar to: build_spline_nd
         *
         */
        static Mat<real>
        build_spline_deriv_uniform_nd(
                Mat<real> const & aXi,
                uint const & aPolynomialDegree,
                uint const & aDim);

        /**
         * Creates two or three-dimensional derivative of uniform basis functions for a polynomial degree.
         *
         * @param[in] aXi                 Matrix with natural coordinates in each direction [xi_1, xi_2, xi_3] .
         * @param[in] aKnotVector         Knot vector for each direction of the b-splines.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-splines for each direction p_1 = p_2 = p_3.
         * @param[in] aElement_flag       Which basis functions are active in the respective element for each direction [start_1 start_2 start_3].
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] Bspline_deriv      Derivative of the basis functions. It creates a Matrix with the derivatives in the rows and the directions in the columns.
         *                                The order is similiar to: build_spline_nd
         *
         */
        static Mat<real>
        build_spline_deriv_uniform_nd(Mat<real> const & aXi,
                Mat<real> & aKnotVector,
                uint const & aPolynomialDegree,
                Mat<uint> const & aElement_flag,
                uint const & aDim);
    };

}   // namespace moris

#endif /* SRC_MODEL_CL_BSPLINE_HPP_ */

