#include "assert.hpp"
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad4.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad8.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad16.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex20.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex64.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Interpolation_Function_Base *
		Interpolation_Function_Factory::create_interpolation_function( const mtk::Geometry_Type       & aGeometryType,
                                                                       const Interpolation_Type       & aInterpolationType,
                                                                       const mtk::Interpolation_Order & aInterpolationOrder )
        {
            // select type
            switch ( aInterpolationType )
            {
                case( Interpolation_Type::LAGRANGE ) :
                {
                    switch ( aGeometryType )
                    {
                        case( mtk::Geometry_Type::LINE ) :
                        {
                            return this->create_lagrange_bar( aInterpolationOrder );
                            break;
                        }
                        case( mtk::Geometry_Type::QUAD ) :
                        {
                           return this->create_lagrange_quad( aInterpolationOrder );
                           break;
                        }
                        case( mtk::Geometry_Type::HEX ) :
                        {
                        	return this->create_lagrange_hex( aInterpolationOrder );
                        	break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, "unknown element geometry type" );
                            return nullptr;
                            break;
                        }
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown interpolation type" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        Interpolation_Function_Base * Interpolation_Function_Factory::create_lagrange_quad( const mtk::Interpolation_Order  & aInterpolationOrder )
        {
            switch ( aInterpolationOrder )
            {
                case ( mtk::Interpolation_Order::LINEAR ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4 >();
                    break;
                }
                case ( mtk::Interpolation_Order::SERENDIPITY ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8 >();
                    break;
                }
                case ( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 9 >();
                    break;
                }
                case ( mtk::Interpolation_Order::CUBIC ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown interpolation order for QUAD" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        Interpolation_Function_Base *
        Interpolation_Function_Factory::create_lagrange_hex( const mtk::Interpolation_Order & aInterpolationOrder )
        {

            switch ( aInterpolationOrder )
            {
                case ( mtk::Interpolation_Order::LINEAR ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8 >();
                    break;
                }
                case ( mtk::Interpolation_Order::SERENDIPITY ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >();
                    break;
                }
                case ( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27 >();
                    break;
                }
                case ( mtk::Interpolation_Order::CUBIC ) :
                {
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown interpolation order for HEX" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        Interpolation_Function_Base *
        Interpolation_Function_Factory::create_lagrange_bar( const mtk::Interpolation_Order & aInterpolationOrder)
        {
            switch ( aInterpolationOrder )
            {
                case(mtk::Interpolation_Order::CONSTANT ) :
                {
                	// bar1
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1 >();
                    break;
                }
                case( mtk::Interpolation_Order::LINEAR ) :
                {
                	// bar2
                	return new Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >();
                	break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ) :
                {
                	// bar3
                    return new Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "create_lagrange_bar: unknown number of nodes" );
                    return nullptr;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
