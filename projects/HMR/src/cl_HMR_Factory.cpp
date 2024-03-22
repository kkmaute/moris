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
                break;
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
            luint                      aPolynomialDegree,
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

    BSpline_Mesh_Base*
    Factory::create_bspline_mesh(
            Background_Mesh_Base* aBackgroundMesh,
            uint                  aActivationPattern,
            luint                 aPolynomialDegree,
            uint                  aMeshIndex )
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
                        return new BSpline_Mesh< 1, 1, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 2 ):
                    {
                        return new BSpline_Mesh< 2, 2, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 3 ):
                    {
                        return new BSpline_Mesh< 3, 3, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 4 ):
                    {
                        return new BSpline_Mesh< 4, 4, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 5 ):
                    {
                        return new BSpline_Mesh< 5, 5, 0 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_bspline_mesh(): unsupported polynomial degree %u for dimension %u",
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
                        return new BSpline_Mesh< 1, 1, 1 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 2 ):
                    {
                        return new BSpline_Mesh< 2, 2, 2 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 3 ):
                    {
                        return new BSpline_Mesh< 3, 3, 3 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 4 ):
                    {
                        return new BSpline_Mesh< 4, 4, 4 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    case ( 5 ):
                    {
                        return new BSpline_Mesh< 5, 5, 5 >(
                                mParameters,
                                aBackgroundMesh,
                                aActivationPattern,
                                aMeshIndex );
                    }
                    default:
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Factory::create_bspline_mesh(): unsupported polynomial degree %u for dimension %u",
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
                        "hmr::Factory::create_bspline_mesh(): unknown number of dimensions %u",
                        (unsigned int)tNumberOfDimensions );
                return nullptr;
            }
        }
    }

    //-------------------------------------------------------------------------------

#define CREATE_BSPLINE_MESH( x, y, z ) \
    switch ( z )                                                                         \
    {                                                                                    \
        case 0:                                                                          \
            return new BSpline_Mesh< x, y, 0 >(                                          \
                    mParameters,                                                         \
                    aBackgroundMesh,                                                     \
                    tPattern,                                                            \
                    aMeshIndex );                                                        \
        case 1:                                                                          \
            return new BSpline_Mesh< x, y, 1 >(                                          \
                    mParameters,                                                         \
                    aBackgroundMesh,                                                     \
                    tPattern,                                                            \
                    aMeshIndex );                                                        \
        case 2:                                                                          \
            return new BSpline_Mesh< x, y, 2 >(                                          \
                    mParameters,                                                         \
                    aBackgroundMesh,                                                     \
                    tPattern,                                                            \
                    aMeshIndex );                                                        \
        case 3:                                                                          \
            return new BSpline_Mesh< x, y, 3 >(                                          \
                    mParameters,                                                         \
                    aBackgroundMesh,                                                     \
                    tPattern,                                                            \
                    aMeshIndex );                                                        \
        default:                                                                         \
            MORIS_ERROR( false, "creating B-spline mesh with incorrect z order %d", z ); \
            return nullptr;                                                              \
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
