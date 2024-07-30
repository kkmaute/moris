/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HPP_

#include "assert.hpp"

#include "moris_typedefs.hpp"                    //MRS/COR/src
#include "cl_Matrix.hpp"                         //LNA/src
#include "linalg_typedefs.hpp"                   //LNA/src
#include "cl_MTK_Enums.hpp"                      //MTK/src
#include "cl_MTK_Integration_Coeffs_Base.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        template< Integration_Type T,
                Integration_Order  P >
        class Integration_Coeffs : public Integration_Coeffs_Base
        {
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Integration_Coeffs(){};

            //------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Integration_Coeffs(){};
            //------------------------------------------------------------------------------

            /**
             * tells how many dimensions this rule has
             */
            uint get_number_of_dimensions();

            //------------------------------------------------------------------------------

            /**
             * tells how many integration points this rule used
             */
            uint get_number_of_points();

            //------------------------------------------------------------------------------
            /**
             * returns the integration weights
             *
             * @param[ in ] aIntegrationWeights
             */
            void get_weights( Matrix< DDRMat > &aIntegrationWeights );

            //------------------------------------------------------------------------------

            /**
             * writes the integration points into given Mat
             *
             * @param[ in ] aIntegrationPoints
             */
            void get_points( Matrix< DDRMat > &aIntegrationPoints );
        };

        //------------------------------------------------------------------------------

        template< Integration_Type T, Integration_Order P >
        uint
        Integration_Coeffs< T, P >::get_number_of_dimensions()
        {
            MORIS_ERROR( false,
                    "get_number_of_dimensions() not implemented for this integration rule." );
            return 0;
        }

        //------------------------------------------------------------------------------

        template< Integration_Type T, Integration_Order P >
        uint
        Integration_Coeffs< T, P >::get_number_of_points()
        {
            MORIS_ERROR( false,
                    "get_number_of_points() not implemented for this integration rule." );
            return 0;
        }

        //------------------------------------------------------------------------------

        template< Integration_Type T, Integration_Order P >
        void
        Integration_Coeffs< T, P >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            MORIS_ERROR( false,
                    "get_weights() not implemented for this rule." );
        }

        //------------------------------------------------------------------------------

        template< Integration_Type T, Integration_Order P >
        void
        Integration_Coeffs< T, P >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            MORIS_ERROR( false,
                    "get_points() not implemented for this rule." );
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HPP_ */
