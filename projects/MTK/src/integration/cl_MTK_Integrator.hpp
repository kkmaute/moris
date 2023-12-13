/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATOR_HPP_
#define SRC_MTK_CL_MTK_INTEGRATOR_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_MTK_Enums.hpp"                 //MTK/src
#include "cl_MTK_Integration_Rule.hpp"      //MTK/src
#include "cl_MTK_Integration_Coeffs.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        class Integrator
        {
            // pointer to space rule, if specified
            Integration_Coeffs_Base *mSpaceCoeffs = nullptr;

            // pointer to time rule, if specified
            Integration_Coeffs_Base *mTimeCoeffs = nullptr;

            // number of points in space
            uint mNumOfSpacePoints;

            // number of points in time
            uint mNumOfTimePoints;

            // matrix with space points
            Matrix< DDRMat > mSpacePoints;

            // matrix with time points
            Matrix< DDRMat > mTimePoints;

            // matrix with space weights
            Matrix< DDRMat > mSpaceWeights;

            // matrix with time weights
            Matrix< DDRMat > mTimeWeights;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * constructs an integrator from an integration rule
             **/
            Integrator( const Integration_Rule &aIntegrationRule );

            //------------------------------------------------------------------------------

            /**
             * destructor
             **/
            ~Integrator();

            //------------------------------------------------------------------------------

            /**
             * get the number of integration points
             **/
            uint get_number_of_points();

            //------------------------------------------------------------------------------

            /**
             * get the integration points as a (d x n) matrix, where d is the dimension of the space
             * and n is the number of integration points (each column contains one point)
             **/
            void             get_points( Matrix< DDRMat > &aIntegrationPoints ) const;
            Matrix< DDRMat > get_points() const;

            //------------------------------------------------------------------------------

            /**
             * get the integration point weights
             **/
            void             get_weights( Matrix< DDRMat > &aIntegrationWeights ) const;
            Matrix< DDRMat > get_weights() const;

            //------------------------------------------------------------------------------

        };    // class Integrator

        //------------------------------------------------------------------------------

    }    // namespace mtk

}    // namespace moris

#endif /* SRC_MTK_CL_MTK_INTEGRATOR_HPP_ */
