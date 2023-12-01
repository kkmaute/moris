/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Rule.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRAITION_RULE_HPP_
#define SRC_MTK_CL_MTK_INTEGRAITION_RULE_HPP_

#include "cl_MTK_Enums.hpp"                   //MTK/src

#include "cl_MTK_Enums.hpp"                   //MTK/src
#include "cl_MTK_Integration_Coeffs_Base.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
    /**
     *  /brief a container that defines an interpolation function
     */
    class Integration_Rule
    {
        const Geometry_Type        mGeometryType;
        const Integration_Type          mSpaceIntegrationType;
        const Integration_Order         mSpaceIntegrationOrder;

        const Geometry_Type        mTimeGeometryType;
        const Integration_Type          mTimeIntegrationType;
        const Integration_Order         mTimeIntegrationOrder;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructs an integration rule
         *
         * @param[ in ] aGeometryType             eg. QUAD, HEX ...
         * @param[ in ] aSpaceIntegrationType     eg. GAUSS, CONSTANT ...
         * @param[ in ] aSpaceIntegrationOrder    eg. 2x2, 3x3 ...
         * @param[ in ] aTimeIntegrationType      eg. GAUSS, CONSTANT ...
         * @param[ in ] aTimeIntegrationOrder     eg. 2x2, 3x3 ...
         *
         */
        Integration_Rule( const Geometry_Type  & aGeometryType,
                          const Integration_Type    & aSpaceIntegrationType,
                          const Integration_Order   & aSpaceIntegrationOrder,
                          const Integration_Type    & aTimeIntegrationType,
                          const Integration_Order   & aTimeIntegrationOrder );

        Integration_Rule( const Geometry_Type  & aGeometryType,
                          const Integration_Type    & aSpaceIntegrationType,
                          const Integration_Order   & aSpaceIntegrationOrder,
                          const Geometry_Type  & aTimeGeometryType,
                          const Integration_Type    & aTimeIntegrationType,
                          const Integration_Order   & aTimeIntegrationOrder );

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Integration_Rule(){};

//------------------------------------------------------------------------------
        /**
         * returns the integration primitive
         */
        Geometry_Type get_geometry_type() const
        {
            return mGeometryType;
        }

//------------------------------------------------------------------------------
        /**
         * returns the integration primitive
         */
        Geometry_Type get_time_geometry_type() const
        {
            return mTimeGeometryType;
        }
//------------------------------------------------------------------------------
        /**
         * returns the space integration order
         */
        Integration_Order get_space_integration_order() const
        {
            return mSpaceIntegrationOrder;
        }

//------------------------------------------------------------------------------
        /**
         * returns the time integration order
         */
        Integration_Order get_time_integration_order() const
        {
            return mTimeIntegrationOrder;
        }

//------------------------------------------------------------------------------
        /**
         * returns the space integration type
         */
        Integration_Type get_space_integration_type() const
        {
            return mSpaceIntegrationType;
        }

//------------------------------------------------------------------------------
        /**
         * returns the time integration type
         */
        Integration_Type get_time_integration_type() const
        {
            return mTimeIntegrationType;
        }

//------------------------------------------------------------------------------
        /**
         * create the space integration coefficients
         */
        Integration_Coeffs_Base * create_space_coeffs() const;

//------------------------------------------------------------------------------
        /**
         * create the time integration coefficients
         */
        Integration_Coeffs_Base * create_time_coeffs() const;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
        /**
         * the factory function
         */
        Integration_Coeffs_Base * create_coeffs( const Geometry_Type & aGeometryType,
                                                 const Integration_Type   & aIntegrationType,
                                                 const Integration_Order  & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_point( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_bar( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_quad( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_hex( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_tri( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base * create_coeffs_gauss_tet( const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTEGRATION_RULE_CPP_ */
