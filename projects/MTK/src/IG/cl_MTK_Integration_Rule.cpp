//MTK/src
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integration_Coeffs.hpp"
#include "cl_MTK_Integration_Coeffs_Point.hpp"
//bar
#include "cl_MTK_Integration_Coeffs_Bar_1.hpp"
#include "cl_MTK_Integration_Coeffs_Bar_2.hpp"
#include "cl_MTK_Integration_Coeffs_Bar_3.hpp"
#include "cl_MTK_Integration_Coeffs_Bar_4.hpp"
#include "cl_MTK_Integration_Coeffs_Bar_5.hpp"
#include "cl_MTK_Integration_Coeffs_Bar_6.hpp"
//quad
#include "cl_MTK_Integration_Coeffs_Quad_2x2.hpp"
#include "cl_MTK_Integration_Coeffs_Quad_3x3.hpp"
#include "cl_MTK_Integration_Coeffs_Quad_4x4.hpp"
#include "cl_MTK_Integration_Coeffs_Quad_5x5.hpp"
//tri
#include "cl_MTK_Integration_Coeffs_Tri_1.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_3.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_4.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_6.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_7.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_12.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_13.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_16.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_19.hpp"
#include "cl_MTK_Integration_Coeffs_Tri_25.hpp"
//tet
#include "cl_MTK_Integration_Coeffs_Tet_1.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_4.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_5.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_11.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_15.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_20.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_35.hpp"
#include "cl_MTK_Integration_Coeffs_Tet_56.hpp"
//hex
#include "cl_MTK_Integration_Coeffs_Hex_2x2x2.hpp"
#include "cl_MTK_Integration_Coeffs_Hex_3x3x3.hpp"
#include "cl_MTK_Integration_Coeffs_Hex_4x4x4.hpp"
#include "cl_MTK_Integration_Coeffs_Hex_5x5x5.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Integration_Rule::Integration_Rule(
                const Geometry_Type & aGeometryType,
                const Integration_Type   & aSpaceIntegrationType,
                const Integration_Order  & aSpaceIntegrationOrder,
                const Integration_Type   & aTimeIntegrationType,
                const Integration_Order  & aTimeIntegrationOrder )
        : mGeometryType( aGeometryType ),
          mSpaceIntegrationType( aSpaceIntegrationType ),
          mSpaceIntegrationOrder( aSpaceIntegrationOrder ),
          mTimeGeometryType( Geometry_Type::LINE ),
          mTimeIntegrationType( aTimeIntegrationType ),
          mTimeIntegrationOrder( aTimeIntegrationOrder )
        {}

        Integration_Rule::Integration_Rule(
                const Geometry_Type & aGeometryType,
                const Integration_Type   & aSpaceIntegrationType,
                const Integration_Order  & aSpaceIntegrationOrder,
                const Geometry_Type & aTimeGeometryType,
                const Integration_Type   & aTimeIntegrationType,
                const Integration_Order  & aTimeIntegrationOrder )
        : mGeometryType( aGeometryType ),
          mSpaceIntegrationType( aSpaceIntegrationType ),
          mSpaceIntegrationOrder( aSpaceIntegrationOrder ),
          mTimeGeometryType( aTimeGeometryType ),
          mTimeIntegrationType( aTimeIntegrationType ),
          mTimeIntegrationOrder( aTimeIntegrationOrder )
        {}

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_space_coeffs() const
        {
            return this->create_coeffs(
                    mGeometryType,
                    mSpaceIntegrationType,
                    mSpaceIntegrationOrder );
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_time_coeffs() const
        {
            return this->create_coeffs(
                    mTimeGeometryType,
                    mTimeIntegrationType,
                    mTimeIntegrationOrder );
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs(
                const Geometry_Type & aGeometryType,
                const Integration_Type   & aIntegrationType,
                const Integration_Order  & aIntegrationOrder ) const
        {
            switch( aIntegrationType )
            {
                case Integration_Type::GAUSS :
                {
                    switch( aGeometryType )
                    {
                        case Geometry_Type::POINT :
                            return this->create_coeffs_gauss_point( aIntegrationOrder );

                        case Geometry_Type::LINE :
                            return this->create_coeffs_gauss_bar( aIntegrationOrder );

                        case Geometry_Type::QUAD :
                            return this->create_coeffs_gauss_quad( aIntegrationOrder );

                        case Geometry_Type::HEX :
                            return this->create_coeffs_gauss_hex( aIntegrationOrder );

                        case Geometry_Type::TRI :
                            return this->create_coeffs_gauss_tri( aIntegrationOrder );

                        case Geometry_Type::TET :
                            return this->create_coeffs_gauss_tet( aIntegrationOrder );

                        default :
                            MORIS_ERROR( false, " Integration_Rule::create_coeffs - unknown geometry type. ");
                            return nullptr;

                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Integration_Rule::create_coeffs - unknown integration type. ");
                    return nullptr;
                }
            }
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_point(
                const Integration_Order & aIntegrationOrder ) const
        {
            return new Integration_Coeffs<
                    Integration_Type::GAUSS,
                    Integration_Order::POINT >();
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_bar(
                const Integration_Order & aIntegrationOrder ) const
        {
            switch( aIntegrationOrder )
            {
                case Integration_Order::BAR_1 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_1 >();

                case Integration_Order::BAR_2 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_2 >();

                case Integration_Order::BAR_3 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_3 >();

                case Integration_Order::BAR_4 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_4 >();

                case Integration_Order::BAR_5 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_5 >();

                case Integration_Order::BAR_6 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::BAR_6 >();

                default :
                    MORIS_ERROR( false, "Integration_Rule::create_coeffs_gauss_bar - integration order not implemented/allowed for bar.");
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_quad(
                const Integration_Order & aIntegrationOrder ) const
        {
            switch( aIntegrationOrder )
            {
                case Integration_Order::QUAD_2x2 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::QUAD_2x2 >();

                case Integration_Order::QUAD_3x3 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::QUAD_3x3 >();

                case Integration_Order::QUAD_4x4 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::QUAD_4x4 >();

                case Integration_Order::QUAD_5x5 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::QUAD_5x5 >();

                default :
                    MORIS_ERROR( false, "Integration_Rule::create_coeffs_gauss_quad - integration order not implemented/allowed for quad.");
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_hex
        ( const Integration_Order & aIntegrationOrder ) const
        {
            switch( aIntegrationOrder )
            {
                case Integration_Order::HEX_2x2x2 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::HEX_2x2x2 >();

                case Integration_Order::HEX_3x3x3 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::HEX_3x3x3 >();

                case Integration_Order::HEX_4x4x4 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::HEX_4x4x4 >();

                case Integration_Order::HEX_5x5x5 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::HEX_5x5x5 >();

                default :
                    MORIS_ERROR( false, " Integration_Rule::create_coeffs_gauss_hex - integration order not implemented/allowed for hex.");
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_tri(
                const Integration_Order & aIntegrationOrder ) const
        {
            switch( aIntegrationOrder )
            {
                // for polynomial order 1
                case Integration_Order::TRI_1 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_1 >();

                    // for polynomial order 2
                case Integration_Order::TRI_3 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_3 >();

                    // for polynomial order 3
                case Integration_Order::TRI_4 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_4 >();

                    // for polynomial order 4
                case Integration_Order::TRI_6 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_6 >();

                    // for polynomial order 5
                case Integration_Order::TRI_7 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_7 >();

                    // for polynomial order 6
                case Integration_Order::TRI_12 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_12 >();

                    // for polynomial order 7
                case Integration_Order::TRI_13 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_13 >();

                    // for polynomial order 8
                case Integration_Order::TRI_16 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_16 >();

                    // for polynomial order 9
                case Integration_Order::TRI_19 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_19 >();

                    // for polynomial order 10
                case Integration_Order::TRI_25 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TRI_25 >();

                default :
                    MORIS_ERROR( false, " Integration_Rule::create_coeffs_gauss_tri - integration order not implemented/allowed for tri.");
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Integration_Coeffs_Base * Integration_Rule::create_coeffs_gauss_tet(
                const Integration_Order & aIntegrationOrder ) const
        {
            switch( aIntegrationOrder )
            {
                // for polynomial order 1
                case Integration_Order::TET_1 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_1 >();

                    // for polynomial order 2
                case Integration_Order::TET_4 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_4 >();

                    // for polynomial order 3
                case Integration_Order::TET_5 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_5 >();

                    // for polynomial order 4
                case Integration_Order::TET_11 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_11 >();

                    // for polynomial order 5
                case Integration_Order::TET_15 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_15 >();

                    // for polynomial order 6
                case Integration_Order::TET_20 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_20 >();

                    // for polynomial order 7
                case Integration_Order::TET_35 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_35 >();

                    // for polynomial order 9
                case Integration_Order::TET_56 :
                    return new Integration_Coeffs<
                            Integration_Type::GAUSS,
                            Integration_Order::TET_56 >();

                default :
                    MORIS_ERROR( false, " Integration_Rule::create_coeffs_gauss_tet - integration order not implemented/allowed for tet.");
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */