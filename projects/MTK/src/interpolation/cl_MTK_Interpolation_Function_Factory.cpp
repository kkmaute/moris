/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Factory.cpp
 *
 */

#include "assert.hpp"
#include "cl_MTK_Enums.hpp"                             //MTK/src
#include "cl_MTK_Interpolation_Rule.hpp"                //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"            //MTK/src

#include "cl_MTK_Interpolation_Function_Lagrange_Bar1.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Bar2.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Bar3.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Bar4.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Quad4.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Quad8.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Quad9.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Quad16.hpp"    //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Hex8.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Hex20.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Hex27.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Hex64.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tri3.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tri6.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tri10.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tet4.hpp"      //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tet10.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Lagrange_Tet20.hpp"     //MTK/src

#include "cl_MTK_Interpolation_Function_Constant_Bar2.hpp"     //MTK/src
#include "cl_MTK_Interpolation_Function_Constant_Point.hpp"    //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_interpolation_function(
            const Geometry_Type&       aGeometryType,
            const Interpolation_Type&  aInterpolationType,
            const Interpolation_Order& aInterpolationOrder )
    {
        // select type
        switch ( aInterpolationType )
        {
            case ( Interpolation_Type::LAGRANGE ):
            {
                switch ( aGeometryType )
                {
                    case ( Geometry_Type::LINE ):
                    {
                        return this->create_lagrange_bar( aInterpolationOrder );
                        break;
                    }
                    case ( Geometry_Type::QUAD ):
                    {
                        return this->create_lagrange_quad( aInterpolationOrder );
                        break;
                    }
                    case ( Geometry_Type::HEX ):
                    {
                        return this->create_lagrange_hex( aInterpolationOrder );
                        break;
                    }
                    case ( Geometry_Type::TRI ):
                    {
                        return this->create_lagrange_tri( aInterpolationOrder );
                        break;
                    }
                    case ( Geometry_Type::TET ):
                    {
                        return this->create_lagrange_tet( aInterpolationOrder );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Interpolation_Function_Factory::create_interpolation_function - unknown element geometry type" );
                        return nullptr;
                        break;
                    }
                }
                break;
            }

            case ( Interpolation_Type::CONSTANT ):
            {
                switch ( aGeometryType )
                {
                    case ( Geometry_Type::POINT ):
                    {
                        return this->create_constant_point( aInterpolationOrder );
                        break;
                    }
                    case ( Geometry_Type::LINE ):
                    {
                        return this->create_constant_bar( aInterpolationOrder );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Interpolation_Function_Factory::create_interpolation_function - unknown element geometry type" );
                        return nullptr;
                        break;
                    }
                }
                break;
            }

            default:
            {
                MORIS_ERROR( false, "Interpolation_Function_Factory::create_interpolation_function - unknown interpolation type" );
                return nullptr;
                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_lagrange_quad( const Interpolation_Order& aInterpolationOrder )
    {
        switch ( aInterpolationOrder )
        {
            case ( Interpolation_Order::LINEAR ):
            {
                return new Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >();
                break;
            }
            case ( Interpolation_Order::SERENDIPITY ):
            {
                return new Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >();
                break;
            }
            case ( Interpolation_Order::QUADRATIC ):
            {
                return new Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >();
                break;
            }
            case ( Interpolation_Order::CUBIC ):
            {
                return new Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_quad - unknown interpolation order for QUAD" );
                return nullptr;
                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_lagrange_hex( const Interpolation_Order& aInterpolationOrder )
    {

        switch ( aInterpolationOrder )
        {
            case ( Interpolation_Order::LINEAR ):
            {
                return new Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >();
                break;
            }
            case ( Interpolation_Order::SERENDIPITY ):
            {
                return new Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >();
                break;
            }
            case ( Interpolation_Order::QUADRATIC ):
            {
                return new Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27 >();
                break;
            }
            case ( Interpolation_Order::CUBIC ):
            {
                return new Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 64 >();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_hex - unknown interpolation order for HEX" );
                return nullptr;
                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_lagrange_bar( const Interpolation_Order& aInterpolationOrder )
    {
        switch ( aInterpolationOrder )
        {
            case ( Interpolation_Order::LINEAR ):
            {
                // bar2
                return new Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >();
                break;
            }
            case ( Interpolation_Order::SERENDIPITY ):
            {
                // bar3
                return new Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 3 >();
                break;
            }
            case ( Interpolation_Order::QUADRATIC ):
            {
                // bar3
                return new Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 3 >();
                break;
            }
            case ( Interpolation_Order::CUBIC ):
            {
                // bar4
                return new Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >();
                break;
            }
            default:
            {
                MORIS_ERROR( false, " Interpolation_Function_Factory::create_lagrange_bar - unknown number of nodes" );
                return nullptr;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_lagrange_tri( const Interpolation_Order& aInterpolationOrder )
    {
        switch ( aInterpolationOrder )
        {
            case ( Interpolation_Order::LINEAR ):
            {
                return new Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >();
                break;
            }
            case ( Interpolation_Order::QUADRATIC ):
            {
                return new Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >();
                break;
            }
            case ( Interpolation_Order::CUBIC ):
            {
                return new Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_tri - unknown interpolation order for TRI" );
                return nullptr;
                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_lagrange_tet( const Interpolation_Order& aInterpolationOrder )
    {
        switch ( aInterpolationOrder )
        {
            case ( Interpolation_Order::LINEAR ):
            {
                return new Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >();
                break;
            }
            case ( Interpolation_Order::QUADRATIC ):
            {
                return new Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >();
                break;
            }
            case ( Interpolation_Order::CUBIC ):
            {
                return new Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_tet - unknown interpolation order for TET" );
                return nullptr;
                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_constant_bar( const Interpolation_Order& aInterpolationOrder )
    {
        // bar 2 with one constant shape function
        return new Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::CONSTANT, 1, 1 >();
    }

    //------------------------------------------------------------------------------

    Interpolation_Function_Base*
    Interpolation_Function_Factory::create_constant_point( const Interpolation_Order& aInterpolationOrder )
    {
        // point with one constant shape function
        return new Interpolation_Function< Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1 >();
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk
