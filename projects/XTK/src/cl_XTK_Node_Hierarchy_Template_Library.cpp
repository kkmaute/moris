#include "cl_XTK_Node_Hierarchy_Interface.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_Pairing.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include <algorithm>// std::sort, std::stable_sort
#include <numeric>
#include "cl_Tracer.hpp"
#include <chrono>
namespace xtk
{

// ----------------------------------------------------------------------------------

void
Node_Hierachy_Template_Library::load_template(
    moris_index              aSpatialDim,
    moris_index              aTemplateId,
    Node_Hierarchy_Template* aNodeHierTemplate )
{
    if ( aSpatialDim == 2 )
    {
        this->load_2d_template( aTemplateId, aNodeHierTemplate );
    }
    else if ( aSpatialDim == 3 )
    {
        this->load_3d_template( aTemplateId, aNodeHierTemplate );
    }
    else
    {
        MORIS_ERROR( false, "Node_Hierachy_Template_Library::load_template() - Number of spatial dimensions not 2 or 3." );
    }
}

// ----------------------------------------------------------------------------------

void
Node_Hierachy_Template_Library::load_2d_template(
    moris_index              aTemplateId,
    Node_Hierarchy_Template* aNodeHierTemplate )
{
    aNodeHierTemplate->mCellTopology = CellTopology::TRI3;
    switch ( aTemplateId )
    {
        case 1:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 2 }, { 3, 1, 4 },{ 3, 4, 2 } } );
            aNodeHierTemplate->mNumCells          = 3;
            break;
        }

        case 2:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 4 }, { 3, 2, 4 }, { 3, 1, 2 } } );
            aNodeHierTemplate->mNumCells          = 3;
            break;
        }

        case 3:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 3 }, { 0, 3, 4 }, { 3, 2, 4 } } );
            aNodeHierTemplate->mNumCells          = 3;
            break;
        }

        case 10:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 2 }, { 3, 1, 2 } } );
            aNodeHierTemplate->mNumCells          = 2;
            break;
        }

        case 11:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 3 }, { 0, 3, 2 } } );
            aNodeHierTemplate->mNumCells          = 2;
            break;
        }

        case 12:
        {
            aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 3 }, { 1, 2, 3 } } );
            aNodeHierTemplate->mNumCells          = 2;
            break;
        }

        default:
        {
            std::cout << "Node_Hierachy_Template_Library::load_2d_template(): \n" << std::flush;
            std::cout << "WARNING INVALID TEMPLATE ID: " << aTemplateId << std::endl;
            break;
        }
    }
}

// ----------------------------------------------------------------------------------

void
Node_Hierachy_Template_Library::load_3d_template(
    moris_index              aTemplateId,
    Node_Hierarchy_Template* aNodeHierTemplate )
{
    aNodeHierTemplate->mCellTopology = CellTopology::TET4;
    switch ( aTemplateId )
    {
    case ( 320 ):
    {
        // Permutation 320
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 32 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }


    case ( 203 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 251 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 512 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 125 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 140 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 401 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 14 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 453 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 534 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 345 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 230 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 302 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 23 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 521 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 152 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 215 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 410 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 41 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 104 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 543 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 354 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 435 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 5420 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5240 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 425 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3501 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1503 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3051 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1053 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4312 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2314 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4132 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2134 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 245 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4502 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4052 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1243 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2504 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2054 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5310 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5130 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 135 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 315 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3421 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3241 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1423 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4250 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2450 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4205 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2405 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5031 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5013 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 531 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 513 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3124 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1342 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1324 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3142 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5024 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 524 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5042 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 542 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3150 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1350 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3105 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1305 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4231 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2431 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4213 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2413 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4520 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2540 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4025 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2045 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5301 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5103 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 351 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 153 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3412 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3214 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1432 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1234 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5402 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 452 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5204 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 254 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3510 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1530 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3015 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1035 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4321 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2341 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4123 ):
    {

        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2143 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 10000 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 3 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10001 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 0, 4, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10002 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10003 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 4 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10004 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 3 }, { 0, 1, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10005 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 0, 1, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;


    case ( 10250 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 0, 4, 2, 5 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10520 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 0, 4, 5, 3 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10260 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 0, 1, 5, 3 }, { 0, 4, 2, 5 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10620 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 0, 5, 4, 3 }, { 0, 5, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10560 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 0, 1, 2, 5 }, { 0, 4, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10650 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 0, 5, 2, 4 }, { 0, 5, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10640 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 4, 3 }, { 1, 2, 5, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10460 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 0, 1, 5, 4 }, { 4, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10230 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 4, 2, 3 }, { 1, 4, 5, 3 }, { 0, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10320 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 2, 4, 3 }, { 0, 5, 4, 3 }, { 0, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10130 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 3 }, { 5, 4, 1, 3 }, { 0, 4, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10310 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 3 }, { 4, 5, 2, 3 }, { 0, 5, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10120 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 3 }, { 0, 4, 5, 3 }, { 4, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10210 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 3 }, { 5, 4, 2, 3 }, { 5, 1, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10140 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 5, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10410 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 4 }, { 4, 5, 2, 3 }, { 5, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10150 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10510 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 4 }, { 5, 4, 2, 3 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10450 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 4 }, { 0, 1, 2, 5 }, { 4, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10540 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 2, 4 }, { 5, 4, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10360 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 0, 1, 4, 5 }, { 1, 2, 4, 5 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10630 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 1, 4, 5, 3 }, { 1, 2, 5, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10340 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 4, 1, 2, 5 }, { 1, 2, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10430 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 4, 1, 5, 3 }, { 1, 2, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    default:
    {
        std::cout << "Node_Hierachy_Template_Library::load_2d_template(): \n" << std::flush;
        std::cout << "WARNING INVALID TEMPLATE ID: " << aTemplateId << std::endl;
        break;
    }
    }
}

// ----------------------------------------------------------------------------------

}// namespace xtk
