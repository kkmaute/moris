/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Factory.cpp
 *
 */

#include "cl_HMR_Factory.hpp"               //HMR/src

#include "cl_HMR_Background_Mesh.hpp"       //HMR/src
#include "cl_HMR_Background_Mesh_2D.hpp"    //HMR/src
// #include "cl_HMR_Background_Mesh_3D.hpp" //HMR/src
#include "cl_HMR_BSpline_Element.hpp"           //HMR/src
#include "cl_HMR_BSpline_Element_Hex27.hpp"     //HMR/src
#include "cl_HMR_BSpline_Element_Hex64.hpp"     //HMR/src
#include "cl_HMR_BSpline_Element_Hex8.hpp"      //HMR/src
#include "cl_HMR_BSpline_Element_Quad16.hpp"    //HMR/src
#include "cl_HMR_BSpline_Element_Quad4.hpp"     //HMR/src
#include "cl_HMR_BSpline_Element_Quad9.hpp"     //HMR/src
#include "cl_HMR_BSpline_Mesh.hpp"              //HMR/src
#include "cl_HMR_Lagrange_Edge2.hpp"            //HMR/src
#include "cl_HMR_Lagrange_Edge3.hpp"            //HMR/src
#include "cl_HMR_Lagrange_Edge4.hpp"            //HMR/src
#include "cl_HMR_Lagrange_Element.hpp"          //HMR/src
#include "cl_HMR_Lagrange_Facet_Line2.hpp"      //HMR/src
#include "cl_HMR_Lagrange_Facet_Line3.hpp"      //HMR/src
#include "cl_HMR_Lagrange_Facet_Line4.hpp"      //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad16.hpp"     //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad4.hpp"      //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad9.hpp"      //HMR/src
#include "cl_HMR_Lagrange_Mesh.hpp"             //HMR/src
#include "cl_HMR_T_Matrix.hpp"

namespace moris::hmr
{
    //-------------------------------------------------------------------------------

    Factory::Factory( const Parameters* aParameters )
            : mParameters( aParameters )
    {
    }

    //-------------------------------------------------------------------------------

    Background_Mesh_Base*
    Factory::create_background_mesh()
    {
        // create background mesh object
        Background_Mesh_Base* aMesh;

        // get number of dimensions from settings
        uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

        switch ( tNumberOfDimensions )
        {
            case ( 1 ):
            {
                // create mesh object
                aMesh = new Background_Mesh< 1 >( mParameters );
                break;
            }
            case ( 2 ):
            {
                aMesh = new Background_Mesh< 2 >( mParameters );
                break;
            }
            case ( 3 ):
            {
                aMesh = new Background_Mesh< 3 >( mParameters );
                break;
            }
            default:
            {
                std::fprintf(
                        stdout,
                        "hmr::Factory::create_background_mesh(): unknown number of dimensions %u",
                        (unsigned int)tNumberOfDimensions );
                exit( -1 );
            }
        }

        // reset main patterns of this mesh
        // fixme: this should be its own function
        aMesh->reset_pattern( mParameters->get_bspline_input_pattern() );
        aMesh->reset_pattern( mParameters->get_lagrange_input_pattern() );
        aMesh->reset_pattern( mParameters->get_bspline_output_pattern() );
        aMesh->reset_pattern( mParameters->get_lagrange_output_pattern() );

        return aMesh;
    }

    //-------------------------------------------------------------------------------

    Lagrange_Mesh_Base*
    Factory::create_lagrange_mesh(
            Background_Mesh_Base*      aBackgroundMesh,
            Vector< BSpline_Mesh_Base* > aBSplineMeshes,
            uint                       aActivationPattern,
            uint                       aPolynomialDegree,
            uint                       aMeshIndex )
    {
        // get number of dimensions from settings
        uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

        switch ( tNumberOfDimensions )
        {
            case ( 2 ):
            {
                switch ( aPolynomialDegree )
                {
                    case ( 1 ):
                    {
                        return new Lagrange_Mesh< 2, 1 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 2 ):
                    {
                        return new Lagrange_Mesh< 2, 2 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 3 ):
                    {
                        return new Lagrange_Mesh< 2, 3 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u",
                                (unsigned int)aPolynomialDegree,
                                (unsigned int)tNumberOfDimensions );
                        return nullptr;
                    }
                }
            }
            case ( 3 ):
            {
                switch ( aPolynomialDegree )
                {
                    case ( 1 ):
                    {
                        return new Lagrange_Mesh< 3, 1 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 2 ):
                    {
                        return new Lagrange_Mesh< 3, 2 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 3 ):
                    {
                        return new Lagrange_Mesh< 3, 3 >(
                                mParameters,
                                aBackgroundMesh,
                                aBSplineMeshes,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u",
                                (unsigned int)aPolynomialDegree,
                                (unsigned int)tNumberOfDimensions );
                        return nullptr;
                    }
                }
            }
            default:
            {
                MORIS_ERROR(
                        false,
                        "hmr::Factory::create_lagrange_mesh(): unknown number of dimensions %u",
                        (unsigned int)tNumberOfDimensions );
                return nullptr;
            }
        }
    }

    //-------------------------------------------------------------------------------

#define SWITCH_ORDER_Z( x, y, z )                                                   \
    switch ( z )                                                                    \
    {                                                                               \
        case 0:                                                                     \
            return new BSpline_Mesh< x, y, 0 >(                                     \
                    mParameters,                                                    \
                    aBackgroundMesh,                                                \
                    aPattern,                                                       \
                    aMeshIndex );                                                   \
        case 1:                                                                     \
            return new BSpline_Mesh< x, y, 1 >(                                     \
                    mParameters,                                                    \
                    aBackgroundMesh,                                                \
                    aPattern,                                                       \
                    aMeshIndex );                                                   \
        case 2:                                                                     \
            return new BSpline_Mesh< x, y, 2 >(                                     \
                    mParameters,                                                    \
                    aBackgroundMesh,                                                \
                    aPattern,                                                       \
                    aMeshIndex );                                                   \
        case 3:                                                                     \
            return new BSpline_Mesh< x, y, 3 >(                                     \
                    mParameters,                                                    \
                    aBackgroundMesh,                                                \
                    aPattern,                                                       \
                    aMeshIndex );                                                   \
        default:                                                                    \
            MORIS_ERROR( false, "Cannot create B-spline mesh with z order %d", z ); \
            return nullptr;                                                         \
    }

#define SWITCH_ORDER_Y( x, y, z )                                                   \
    switch ( y )                                                                    \
    {                                                                               \
        case 0:                                                                     \
            SWITCH_ORDER_Z( x, 0, z )                                               \
        case 1:                                                                     \
            SWITCH_ORDER_Z( x, 1, z )                                               \
        case 2:                                                                     \
            SWITCH_ORDER_Z( x, 2, z )                                               \
        case 3:                                                                     \
            SWITCH_ORDER_Z( x, 3, z )                                               \
        default:                                                                    \
            MORIS_ERROR( false, "Cannot create B-spline mesh with y order %d", y ); \
            return nullptr;                                                         \
    }

    //-------------------------------------------------------------------------------

    BSpline_Mesh_Base*
    Factory::create_bspline_mesh(
            Background_Mesh_Base* aBackgroundMesh,
            uint                  aPattern,
            uint                  aPolynomialDegree,
            uint                  aMeshIndex )
    {
        // get number of dimensions from settings
        uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

        // Set B-spline order vector
        Vector< uint > tBSplineOrders;
        if ( tNumberOfDimensions == 2 )
        {
            tBSplineOrders = { aPolynomialDegree, aPolynomialDegree, 0 };
        }
        else
        {
            tBSplineOrders = { aPolynomialDegree, aPolynomialDegree, aPolynomialDegree };
        }

        switch ( tBSplineOrders( 0 ) )
        {
            case ( 1 ):
                SWITCH_ORDER_Y( 1, tBSplineOrders( 1 ), tBSplineOrders( 2 ) )
            case ( 2 ):
                SWITCH_ORDER_Y( 2, tBSplineOrders( 1 ), tBSplineOrders( 2 ) )
            case ( 3 ):
                SWITCH_ORDER_Y( 3, tBSplineOrders( 1 ), tBSplineOrders( 2 ) )
            default:
            {
                MORIS_ERROR( false, "Cannot create B-spline mesh with x order %d", tBSplineOrders( 0 ) );
                return nullptr;
            }
        }
    }

    //-------------------------------------------------------------------------------

    BSpline_Mesh_Base*
    Factory::create_dummy_bspline_mesh(
            Background_Mesh_Base* aBackgroundMesh,
            uint                  aOrder,
            uint                  aPattern )
    {
        // get number of dimensions from settings
        uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

        switch ( tNumberOfDimensions )
        {
            case ( 2 ):
            {
                switch ( aOrder )
                {
                    case ( 1 ):
                    {
                        return new BSpline_Mesh< 1, 1, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 2 ):
                    {
                        return new BSpline_Mesh< 2, 2, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 3 ):
                    {
                        return new BSpline_Mesh< 3, 3, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 4 ):
                    {
                        return new BSpline_Mesh< 4, 4, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 5 ):
                    {
                        return new BSpline_Mesh< 5, 5, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_bspline_mesh(): unsupported polynomial degree %u for dimension %u",
                                (unsigned int)aOrder,
                                (unsigned int)tNumberOfDimensions );
                        return nullptr;
                    }
                }
            }
            case ( 3 ):
            {
                switch ( aOrder )
                {
                    case ( 1 ):
                    {
                        return new BSpline_Mesh< 1, 1, 1 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 2 ):
                    {
                        return new BSpline_Mesh< 2, 2, 2 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 3 ):
                    {
                        return new BSpline_Mesh< 3, 3, 3 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 4 ):
                    {
                        return new BSpline_Mesh< 4, 4, 4 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    case ( 5 ):
                    {
                        return new BSpline_Mesh< 5, 5, 5 >(
                                mParameters,
                                aBackgroundMesh,
                                aPattern,
                                MORIS_UINT_MAX );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_bspline_mesh(): unsupported polynomial degree %u for dimension %u",
                                (unsigned int)aOrder,
                                (unsigned int)tNumberOfDimensions );
                        return nullptr;
                    }
                }
            }
            default:
            {
                MORIS_ERROR(
                        false,
                        "hmr::Factory::create_bspline_mesh(): unknown number of dimensions %u",
                        (unsigned int)tNumberOfDimensions );
                return nullptr;
            }
        }
    }

    //-------------------------------------------------------------------------------

    T_Matrix_Base*
    Factory::create_t_matrix(
            Lagrange_Mesh_Base* aLagrangeMesh,
            BSpline_Mesh_Base*  aBSplineMesh,
            Lagrange_Mesh_Base* aLagrangeMeshFine )
    {
        switch ( mParameters->get_number_of_dimensions() )
        {
            case 2:
                return create_t_matrix< 2 >( aLagrangeMesh, aBSplineMesh, aLagrangeMeshFine );
            case 3:
                return create_t_matrix< 3 >( aLagrangeMesh, aBSplineMesh, aLagrangeMeshFine );
            default:
            {
                MORIS_ERROR(
                        false,
                        "hmr::create_t_matrix()::create_t_matrix(): unknown number of dimensions %u",
                        (unsigned int)mParameters->get_number_of_dimensions() );
                return nullptr;
            }
        }
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
