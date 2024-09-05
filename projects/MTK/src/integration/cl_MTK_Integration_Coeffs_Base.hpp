/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Base.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BASE_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BASE_HPP_

#include "moris_typedefs.hpp"     //MRS/COR/src
#include "cl_Matrix.hpp"          //LNA/src
#include "linalg_typedefs.hpp"    //LNA/src
#include "cl_MTK_Enums.hpp"       //MTK/src

namespace moris::mtk
{

    //------------------------------------------------------------------------------

    class Integration_Coeffs_Base
    {

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /* trivial constructor */
        Integration_Coeffs_Base(){};

        //------------------------------------------------------------------------------

        /* trivial destructor */
        virtual ~Integration_Coeffs_Base(){};

        //------------------------------------------------------------------------------

        /**
         * returns the number of dimensions
         */
        virtual uint
        get_number_of_dimensions() = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the number of points
         */
        virtual uint
        get_number_of_points() = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the integration weights
         *
         * @param[ in ] aIntegrationWeights
         */
        virtual void get_weights( Matrix< DDRMat > &aIntegrationWeights ) = 0;

        //------------------------------------------------------------------------------

        /**
         * writes the integration points into given Mat
         *
         * @param[ in ] aIntegrationPoints
         */
        virtual void get_points( Matrix< DDRMat > &aIntegrationPoints ) = 0;

        //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BASE_HPP_ */
