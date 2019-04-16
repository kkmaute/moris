#include "assert.hpp"
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Bar4.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad4.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad8.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Quad16.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex20.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Hex64.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Constant_Bar2.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tri3.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tri6.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tri10.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tet4.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tet10.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Lagrange_Tet20.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    Interpolation_Function_Base * Interpolation_Function_Factory::create_interpolation_function
        ( const mtk::Geometry_Type       & aGeometryType,
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
                        case( mtk::Geometry_Type::TRI ) :
                        {
                            return this->create_lagrange_tri( aInterpolationOrder );
                            break;
                        }
                        case( mtk::Geometry_Type::TET ) :
                        {
                            return this->create_lagrange_tet( aInterpolationOrder );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, " Interpolation_Function_Factory::create_interpolation_function - unknown element geometry type" );
                            return nullptr;
                            break;
                        }
                    }
                    break;
                }

                case( Interpolation_Type::CONSTANT ) :
                {
                    switch ( aGeometryType )
                    {
                        case( mtk::Geometry_Type::LINE ) :
                        {
                            return this->create_constant_bar( aInterpolationOrder );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, " Interpolation_Function_Factory::create_interpolation_function - unknown element geometry type" );
                            return nullptr;
                            break;
                        }
                    }
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, "Interpolation_Function_Factory::create_interpolation_function - unknown interpolation type" );
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
                    return new Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >();
                    break;
                }
                case ( mtk::Interpolation_Order::SERENDIPITY ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >();
                    break;
                }
                case ( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >();
                    break;
                }
                case ( mtk::Interpolation_Order::CUBIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_quad - unknown interpolation order for QUAD" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

    Interpolation_Function_Base * Interpolation_Function_Factory::create_lagrange_hex( const mtk::Interpolation_Order  & aInterpolationOrder )
        {

            switch ( aInterpolationOrder )
            {
                case ( mtk::Interpolation_Order::LINEAR ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >();
                    break;
                }
                case ( mtk::Interpolation_Order::SERENDIPITY ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >();
                    break;
                }
                case ( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27 >();
                    break;
                }
                case ( mtk::Interpolation_Order::CUBIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 64 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_hex - unknown interpolation order for HEX" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

    Interpolation_Function_Base * Interpolation_Function_Factory::create_lagrange_bar( const mtk::Interpolation_Order  & aInterpolationOrder)
        {
            switch ( aInterpolationOrder )
            {
                case(mtk::Interpolation_Order::CONSTANT ) :
                {
                    // bar1
                    return new Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >();
                    break;
                }
                case( mtk::Interpolation_Order::LINEAR ) :
                {
                    // bar2
                    return new Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >();
                    break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    // bar3
                    return new Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 3 >();
                    break;
                }
                case( mtk::Interpolation_Order::CUBIC ) :
                {
                    // bar4
                    return new Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Interpolation_Function_Factory::create_lagrange_bar - unknown number of nodes" );
                    return nullptr;
                }
            }
        }

//------------------------------------------------------------------------------

    Interpolation_Function_Base * Interpolation_Function_Factory::create_lagrange_tri
        ( const mtk::Interpolation_Order  & aInterpolationOrder )
        {
            switch ( aInterpolationOrder )
            {
                case ( mtk::Interpolation_Order::LINEAR ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >();
                    break;
                }
                case ( mtk::Interpolation_Order::QUADRATIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >();
                    break;
                }
                case ( mtk::Interpolation_Order::CUBIC ) :
                {
                    return new Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_tri - unknown interpolation order for TRI" );
                    return nullptr;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        Interpolation_Function_Base * Interpolation_Function_Factory::create_lagrange_tet
            ( const mtk::Interpolation_Order  & aInterpolationOrder )
            {
                switch ( aInterpolationOrder )
                {
                    case ( mtk::Interpolation_Order::LINEAR ) :
                    {
                        return new Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >();
                        break;
                    }
                    case ( mtk::Interpolation_Order::QUADRATIC ) :
                    {
                        return new Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >();
                        break;
                    }
                    case ( mtk::Interpolation_Order::CUBIC ) :
                    {
                        return new Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >();
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Interpolation_Function_Factory::create_lagrange_tet - unknown interpolation order for TET" );
                        return nullptr;
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------

    Interpolation_Function_Base *
    Interpolation_Function_Factory::create_constant_bar( const mtk::Interpolation_Order  & aInterpolationOrder)
    {
        // bar 2 with one constant shape fucntion
        return new Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::CONSTANT, 1, 1 >();
    }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
