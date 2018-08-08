/*
 * cl_FEM_Integration_Rule.cpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRAITION_RULE_HPP_
#define SRC_FEM_CL_FEM_INTEGRAITION_RULE_HPP_

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp" //FEM/INT/src
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
        const mtk::Geometry_Type             mGeometryType;

        const Integration_Type          mSpaceIntegrationType;
        const Integration_Order         mSpaceIntegrationOrder;

        const Integration_Type          mTimeIntegrationType;
        const Integration_Order         mTimeIntegrationOrder;

        const Integration_Type          mSpaceTimeIntegrationType;
        const Integration_Order         mSpaceTimeIntegrationOrder;

        //! flag telling if integration rule is a combination of two
        const bool                      mHasTwoRulesFlag;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * constructs an integration rule
         *
         * @param[ in ] aIntegrationType          eg. GAUSS, CONSTANT ...
         * @param[ in ] aInterpolationPrimitive   eg. QUAD, HEX ...
         * @param[ in ] aIntegrationOrder         eg. 2x2, 3x3 ...
         *
         */
        Integration_Rule(
                const mtk::Geometry_Type       & aGeometryType,
                const Integration_Type    & aSpaceTimeIntegrationType,
                const Integration_Order   & aSpaceTimeIntegrationOrder );


        Integration_Rule(
                const mtk::Geometry_Type       & aGeometryType,
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
         * returns the integration type
         */
        /* auto
        get_type() const -> decltype ( mIntegrationType )
        {
            return mIntegrationType;
        } */

//------------------------------------------------------------------------------

        /**
         * tells if this rule is combined by space and time
         */
        auto
        has_two_rules() const -> decltype ( mHasTwoRulesFlag )
        {
            return mHasTwoRulesFlag;
        }

//------------------------------------------------------------------------------

        /**
         * returns the interpolation primitive
         */
        auto
        get_geometry_type() const -> decltype ( mGeometryType )
        {
            return mGeometryType;
        }

//------------------------------------------------------------------------------

        /**
         * returns the interpolation order
         */
        auto
        get_order() const -> decltype ( mSpaceTimeIntegrationOrder )
        {
            return mSpaceTimeIntegrationOrder;
        }

//------------------------------------------------------------------------------

        Integration_Coeffs_Base *
        create_space_coeffs() const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base *
        create_time_coeffs() const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base *
        create_space_time_coeffs() const;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
        /**
         * the factory function
         */
        Integration_Coeffs_Base *
        create_coeffs(
                const mtk::Geometry_Type       & aGeometryType,
                const Integration_Type    & aIntegrationType,
                const Integration_Order   & aIntegrationOrder ) const;

//------------------------------------------------------------------------------


       Integration_Coeffs_Base *
       create_coeffs_gauss_bar(
               const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base *
        create_coeffs_gauss_quad(
                const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

        Integration_Coeffs_Base *
        create_coeffs_gauss_hex(
                const Integration_Order & aIntegrationOrder ) const;

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRAITION_RULE_CPP_ */
