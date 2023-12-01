/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator_Test_Polynomial.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_POLYNOMIAL_FOR_TEST_HPP_
#define SRC_MTK_CL_MTK_POLYNOMIAL_FOR_TEST_HPP_

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        /**
         * \brief A special class to test the functionality of the integrator
         */
        class Integrator_Test_Polynomial
        {
            //! coefficients
            Matrix< DDRMat > mCoeffs;

            //! exponents
            Matrix< DDRMat > mExponents;

            //! number of coefficients
            uint mNumberOfCoeffs;

            //! number of dimensions
            uint mNumberOfDimensions;

            //! the solution calculated by MATLAB
            real mIntegral;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * loads a MATLAB-generates binary file and generates the
             * testing polybomial
             *
             * @param[ in ] aPath   path to binary file
             */
            Integrator_Test_Polynomial( const std::string & aPath );

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Integrator_Test_Polynomial(){}

//------------------------------------------------------------------------------

            /**
             * evaluates the polynomial at specified position
             */
            real
            eval( const Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------

            /**
             * returns the MATLAB precalculated integral
             */
            real
            get_integral();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
     } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_POLYNOMIAL_FOR_TEST_HPP_ */
