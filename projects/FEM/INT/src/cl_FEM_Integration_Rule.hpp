/*
 * cl_FEM_Integration_Rule.cpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRAITION_RULE_HPP_
#define SRC_FEM_CL_FEM_INTEGRAITION_RULE_HPP_

#include "cl_MTK_Enums.hpp"                   //MTK/src

#include "cl_FEM_Enums.hpp"                   //FEM/INT/src
#include "cl_FEM_Integration_Coeffs_Base.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
    /**
     *  /brief a container that defines an interpolation function
     */
    class Integration_Rule
    {
        const mtk::Geometry_Type        mGeometryType;

        const Integration_Type          mSpaceIntegrationType;
        const Integration_Order         mSpaceIntegrationOrder;

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
        Integration_Rule( const mtk::Geometry_Type  & aGeometryType,
                          const Integration_Type    & aSpaceIntegrationType,
                          const Integration_Order   & aSpaceIntegrationOrder,
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
        mtk::Geometry_Type get_geometry_type() const
        {
            return mGeometryType;
        }

//------------------------------------------------------------------------------
        /**
         * returns the space integration order
         */
        fem::Integration_Order get_space_integration_order() const
        {
            return mSpaceIntegrationOrder;
        }

//------------------------------------------------------------------------------
        /**
         * returns the time integration order
         */
        fem::Integration_Order get_time_integration_order() const
        {
            return mTimeIntegrationOrder;
        }

//------------------------------------------------------------------------------
        /**
         * returns the space integration type
         */
        fem::Integration_Type get_space_integration_type() const
        {
            return mSpaceIntegrationType;
        }

//------------------------------------------------------------------------------
        /**
         * returns the time integration type
         */
        fem::Integration_Type get_time_integration_type() const
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
        Integration_Coeffs_Base * create_coeffs( const mtk::Geometry_Type & aGeometryType,
                                                 const Integration_Type   & aIntegrationType,
                                                 const Integration_Order  & aIntegrationOrder ) const;

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
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRAITION_RULE_CPP_ */
