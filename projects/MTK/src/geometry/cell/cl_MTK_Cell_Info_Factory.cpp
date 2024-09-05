/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Factory.cpp
 *
 */

#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "cl_MTK_Cell_Info_Tri6.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell_Info_Quad8.hpp"
#include "cl_MTK_Cell_Info_Quad9.hpp"
#include "cl_MTK_Cell_Info_Quad16.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Cell_Info_Tet10.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Hex20.hpp"
#include "cl_MTK_Cell_Info_Hex27.hpp"
#include "cl_MTK_Cell_Info_Hex64.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::mtk
{
    moris::mtk::Cell_Info*
    Cell_Info_Factory::create_cell_info( enum CellTopology aCellTopo )
    {
        moris::mtk::Cell_Info* tConn = nullptr;
        switch ( aCellTopo )
        {
            case ( CellTopology::TRI3 ):
            {
                tConn = new Cell_Info_Tri3();
                break;
            }
            case ( CellTopology::TRI6 ):
            {
                tConn = new Cell_Info_Tri6();
                break;
            }
                // case( CellTopology::TRI10  ):{ tConn = new Cell_Info_Tri10();  break; }

            case ( CellTopology::QUAD4 ):
            {
                tConn = new Cell_Info_Quad4();
                break;
            }
            case ( CellTopology::QUAD8 ):
            {
                tConn = new Cell_Info_Quad8();
                break;
            }
            case ( CellTopology::QUAD9 ):
            {
                tConn = new Cell_Info_Quad9();
                break;
            }
            case ( CellTopology::QUAD16 ):
            {
                tConn = new Cell_Info_Quad16();
                break;
            }

            case ( CellTopology::TET4 ):
            {
                tConn = new Cell_Info_Tet4();
                break;
            }
            case ( CellTopology::TET10 ):
            {
                tConn = new Cell_Info_Tet10();
                break;
            }
                // case( CellTopology::TET20 ):{ tConn = new Cell_Info_Tet20(); break; }

            case ( CellTopology::HEX8 ):
            {
                tConn = new Cell_Info_Hex8();
                break;
            }
            case ( CellTopology::HEX20 ):
            {
                tConn = new Cell_Info_Hex20();
                break;
            }
            case ( CellTopology::HEX27 ):
            {
                tConn = new Cell_Info_Hex27();
                break;
            }
            case ( CellTopology::HEX64 ):
            {
                tConn = new Cell_Info_Hex64();
                break;
            }

            default:
            {
                MORIS_ERROR( 0, "Invalid cell topology specified for Cell_Info factory (this could be because the connecitivty class is not implemented or because aCellTopo is invalid" );
                break;
            }
        }

        return tConn;
    }

    std::shared_ptr< moris::mtk::Cell_Info >
    Cell_Info_Factory::create_cell_info_sp( enum CellTopology aCellTopo )
    {
        std::shared_ptr< moris::mtk::Cell_Info > tConn = nullptr;
        switch ( aCellTopo )
        {
            case ( CellTopology::TRI3 ):
            {
                tConn = std::make_shared< Cell_Info_Tri3 >();
                break;
            }
            case ( CellTopology::TRI6 ):
            {
                tConn = std::make_shared< Cell_Info_Tri6 >();
                break;
            }
                // case( CellTopology::TRI10 ):{ tConn = std::make_shared< Cell_Info_Tri10 >(); break; }

            case ( CellTopology::QUAD4 ):
            {
                tConn = std::make_shared< Cell_Info_Quad4 >();
                break;
            }
            case ( CellTopology::QUAD8 ):
            {
                tConn = std::make_shared< Cell_Info_Quad8 >();
                break;
            }
            case ( CellTopology::QUAD9 ):
            {
                tConn = std::make_shared< Cell_Info_Quad9 >();
                break;
            }
            case ( CellTopology::QUAD16 ):
            {
                tConn = std::make_shared< Cell_Info_Quad16 >();
                break;
            }

            case ( CellTopology::TET4 ):
            {
                tConn = std::make_shared< Cell_Info_Tet4 >();
                break;
            }
            case ( CellTopology::TET10 ):
            {
                tConn = std::make_shared< Cell_Info_Tet10 >();
                break;
            }
                // case( CellTopology::TET20 ):{ tConn = std::make_shared< Cell_Info_Tet20>(); break; }

            case ( CellTopology::HEX8 ):
            {
                tConn = std::make_shared< Cell_Info_Hex8 >();
                break;
            }
            case ( CellTopology::HEX20 ):
            {
                tConn = std::make_shared< Cell_Info_Hex20 >();
                break;
            }
            case ( CellTopology::HEX27 ):
            {
                tConn = std::make_shared< Cell_Info_Hex27 >();
                break;
            }
            case ( CellTopology::HEX64 ):
            {
                tConn = std::make_shared< Cell_Info_Hex64 >();
                break;
            }

            default:
            {
                MORIS_ERROR( 0, "Invalid cell topology specified for Cell_Info factory (this could be because the connecitivty class is not implemented or because aCellTopo is invalid" );
                break;
            }
        }

        return tConn;
    }

    moris::mtk::Cell_Info*
    Cell_Info_Factory::create_cell_info( enum Geometry_Type aCellGeom,
            enum Interpolation_Order                        aInterpOrder )
    {
        moris::mtk::Cell_Info* tConn = nullptr;
        switch ( aCellGeom )
        {
            case ( Geometry_Type::HEX ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = new Cell_Info_Hex8();
                        break;
                    }
                    case ( Interpolation_Order::SERENDIPITY ):
                    {
                        tConn = new Cell_Info_Hex20();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = new Cell_Info_Hex27();
                        break;
                    }
                    case ( Interpolation_Order::CUBIC ):
                    {
                        tConn = new Cell_Info_Hex64();
                        break;
                    }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid hex interpolation order" );
                        break;
                    }
                }

                break;
            }
            case ( Geometry_Type::QUAD ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = new Cell_Info_Quad4();
                        break;
                    }
                    case ( Interpolation_Order::SERENDIPITY ):
                    {
                        tConn = new Cell_Info_Quad8();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = new Cell_Info_Quad9();
                        break;
                    }
                    case ( Interpolation_Order::CUBIC ):
                    {
                        tConn = new Cell_Info_Quad16();
                        break;
                    }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid quad interpolation order" );
                        break;
                    }
                }
                break;
            }
            case ( Geometry_Type::TET ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = new Cell_Info_Tet4();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = new Cell_Info_Tet10();
                        break;
                    }
                        // case(Interpolation_Order::CUBIC):{  tConn = new Cell_Info_Tet20(); break; }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid tet interpolation order" );
                        break;
                    }
                }
                break;
            }
            case ( Geometry_Type::TRI ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = new Cell_Info_Tri3();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = new Cell_Info_Tri6();
                        break;
                    }
                        // case(Interpolation_Order::CUBIC):{  tConn = new Cell_Info_Tri10(); break; }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid tri interpolation order" );
                        break;
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Invalid geometry type" );
                break;
            }
        }

        return tConn;
    }

    std::shared_ptr< moris::mtk::Cell_Info >
    Cell_Info_Factory::create_cell_info_sp( enum Geometry_Type aCellGeom,
            enum Interpolation_Order                           aInterpOrder )
    {
        std::shared_ptr< moris::mtk::Cell_Info > tConn = nullptr;
        switch ( aCellGeom )
        {
            case ( Geometry_Type::HEX ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = std::make_shared< Cell_Info_Hex8 >();
                        break;
                    }
                    case ( Interpolation_Order::SERENDIPITY ):
                    {
                        tConn = std::make_shared< Cell_Info_Hex20 >();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Hex27 >();
                        break;
                    }
                    case ( Interpolation_Order::CUBIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Hex64 >();
                        break;
                    }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid hex interpolation order" );
                        break;
                    }
                }

                break;
            }
            case ( Geometry_Type::QUAD ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = std::make_shared< Cell_Info_Quad4 >();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Quad9 >();
                        break;
                    }
                    case ( Interpolation_Order::SERENDIPITY ):
                    {
                        tConn = std::make_shared< Cell_Info_Quad8 >();
                        break;
                    }
                    case ( Interpolation_Order::CUBIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Quad16 >();
                        break;
                    }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid quad interpolation order" );
                        break;
                    }
                }
                break;
            }
            case ( Geometry_Type::TET ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = std::make_shared< Cell_Info_Tet4 >();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Tet10 >();
                        break;
                    }
                        // case(Interpolation_Order::CUBIC):{  tConn = std::make_shared< Cell_Info_Tet20 >(); break; }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid tet interpolation order" );
                        break;
                    }
                }
                break;
            }
            case ( Geometry_Type::TRI ):
            {
                switch ( aInterpOrder )
                {
                    case ( Interpolation_Order::LINEAR ):
                    {
                        tConn = std::make_shared< Cell_Info_Tri3 >();
                        break;
                    }
                    case ( Interpolation_Order::QUADRATIC ):
                    {
                        tConn = std::make_shared< Cell_Info_Tri6 >();
                        break;
                    }
                        // case(Interpolation_Order::CUBIC):{  tConn = std::make_shared< Cell_Info_Tri10 >(); break; }

                    default:
                    {
                        MORIS_ERROR( 0, "Invalid tri interpolation order" );
                        break;
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Invalid geometry type" );
                break;
            }
        }

        return tConn;
    }

}    // namespace moris::mtk
