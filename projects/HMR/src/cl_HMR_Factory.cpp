/*
 * cl_HMR_Factory.cpp
 *
 *  Created on: May 7, 2018
 *      Author: messe
 */


#include "cl_HMR_Factory.hpp" //HMR/src

#include "cl_HMR_Background_Mesh.hpp" //HMR/src
//#include "cl_HMR_Background_Mesh_1D.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_2D.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_3D.hpp" //HMR/src

#include "cl_HMR_Lagrange_Element.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Quad4.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Quad9.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Quad16.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Hex8.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Hex27.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element_Hex64.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh.hpp" //HMR/src

#include "cl_HMR_BSpline_Element.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad4.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad9.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad16.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex8.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex27.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex64.hpp" //HMR/src

#include "cl_HMR_BSpline_Mesh.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {

//-------------------------------------------------------------------------------

    /**
     * creates a backgroud mesh depending on the number of dimensions set
     */
    Background_Mesh_Base*
    Factory::create_background_mesh( const Parameters * aParameters )
    {

        // create background mesh object
        Background_Mesh_Base* aMesh;

        // get number of dimensions from settings
        auto tNumberOfDimensions =  aParameters->get_number_of_dimensions();

        switch( tNumberOfDimensions )
        {
        case( 1 ) :
        {
            // create mesh object
            aMesh = new Background_Mesh< 1 >( aParameters );
            break;
        }
        case( 2 ) :
        {
            aMesh = new Background_Mesh< 2 >( aParameters );
            break;
        }
        case( 3 ) :
        {
            aMesh = new Background_Mesh< 3 >( aParameters );
            break;
        }
        default :
        {
            std::fprintf( stdout,
                    "create_background_mesh(): unknown number of dimensions %u\n",
                    ( unsigned int ) tNumberOfDimensions );
            exit(-1);
            break;
        }
        }

        return aMesh;
    }
//-------------------------------------------------------------------------------

        Lagrange_Mesh_Base*
        Factory::create_lagrange_mesh(
                const Parameters      * aParameters,
                Background_Mesh_Base  * aBackgroundMesh,
                BSpline_Mesh_Base     * aBSplineMesh,
                const  uint           & aActivePattern,
                const luint           & aPolynomialDegree )
        {
            Lagrange_Mesh_Base* aMesh;

            // get number of dimensions from settings
            auto tNumberOfDimensions =  aParameters->get_number_of_dimensions();

            switch( tNumberOfDimensions )
            {
            case( 2 ) :
            {
                switch ( aPolynomialDegree )
                {
                case( 1 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 1 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 2 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 3 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 4 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 5 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                default:
                {
                    std::fprintf(
                            stdout,
                            "create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u\n",
                             ( unsigned int )  aPolynomialDegree,
                             ( unsigned int )  tNumberOfDimensions );
                    exit(-1);
                    break;
                }
                }
                break;
            }
            case( 3 ) :
            {
                switch ( aPolynomialDegree )
                {
                case( 1 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 1 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 2 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 3 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 4 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 5 >( aParameters, aBackgroundMesh, aBSplineMesh );
                    break;
                }
                default:
                {
                    std::fprintf(
                            stdout,
                            "create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
                    exit(-1);
                    break;
                }
                }
                break;
            }
            default :
            {
                std::fprintf( stdout,
                        "create_lagrange_mesh(): unknown number of dimensions %u\n",
                        ( unsigned int )  tNumberOfDimensions );
                exit(-1);
                break;
            }
            }

            aMesh->set_activation_pattern( aActivePattern );

            return aMesh;
        }
//-------------------------------------------------------------------------------

        BSpline_Mesh_Base*
        Factory::create_bspline_mesh(
                const Parameters      * aParameters,
                Background_Mesh_Base  * aBackgroundMesh,
                const  uint           & aActivePattern,
                const luint           & aPolynomialDegree )
        {
            BSpline_Mesh_Base* aMesh;

            // get number of dimensions from settings
            auto tNumberOfDimensions =  aParameters->get_number_of_dimensions();

            switch( tNumberOfDimensions )
            {
            case( 2 ) :
            {
                switch ( aPolynomialDegree )
                {
                case( 1 ):
                {
                    aMesh = new BSpline_Mesh< 2, 1 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new BSpline_Mesh< 2, 2 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new BSpline_Mesh< 2, 3 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new BSpline_Mesh< 2, 4 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new BSpline_Mesh< 2, 5 >( aParameters, aBackgroundMesh );
                    break;
                }
                default:
                {
                    std::fprintf(
                            stdout,
                            "create_bspline_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
                    exit(-1);
                    break;
                }
            }
            }
            break;
            case( 3 ) :
            {
                switch ( aPolynomialDegree )
                {
                case( 1 ):
                {
                    aMesh = new BSpline_Mesh< 3, 1 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new BSpline_Mesh< 3, 2 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new BSpline_Mesh< 3, 3 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new BSpline_Mesh< 3, 4 >( aParameters, aBackgroundMesh );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new BSpline_Mesh< 3, 5 >( aParameters, aBackgroundMesh );
                    break;
                }
                default:
                {
                    std::fprintf(
                            stdout,
                            "create_bspline_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
                    exit(-1);
                    break;
                }
           }
           break;
           }
           default :
           {
               std::fprintf( stdout,
                        "create_bspline_mesh(): unknown number of dimensions %u\n",
                        ( unsigned int )  tNumberOfDimensions );
                exit(-1);
                break;
            }
            }

            aMesh->set_activation_pattern( aActivePattern );

            return aMesh;
        }

//-------------------------------------------------------------------------------

    }  /* namespace hmr */
} /* namespace moris */
