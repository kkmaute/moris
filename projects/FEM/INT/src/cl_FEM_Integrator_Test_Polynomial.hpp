/*
 * cl_FEM_Integrator_Test_Polynomial.hpp
 *
 *  Created on: Jul 20, 2018
 *      Author: messe
 */
#include <string>
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src


#ifndef SRC_FEM_CL_FEM_POLYNOMIAL_FOR_TEST_HPP_
#define SRC_FEM_CL_FEM_POLYNOMIAL_FOR_TEST_HPP_

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * \brief A special class to test the functionality of the integrator
         */
        class Integrator_Test_Polynomial
        {
            //! coefficients
            Mat< real > mCoeffs;

            //! exponents
            Mat< real > mExponents;

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
            eval( const Mat<real> & aXi );

//------------------------------------------------------------------------------

            /**
             * returns the MATLAB precalculated integral
             */
            real
            get_integral();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
     } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_POLYNOMIAL_FOR_TEST_HPP_ */
