/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Factory.cpp
 *
 */

#include "cl_HMR_Factory.hpp" //HMR/src

#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_2D.hpp" //HMR/src
// #include "cl_HMR_Background_Mesh_3D.hpp" //HMR/src
#include "cl_HMR_BSpline_Element.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex27.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex64.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Hex8.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad16.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad4.hpp" //HMR/src
#include "cl_HMR_BSpline_Element_Quad9.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh.hpp" //HMR/src
#include "cl_HMR_Lagrange_Edge2.hpp" //HMR/src
#include "cl_HMR_Lagrange_Edge3.hpp" //HMR/src
#include "cl_HMR_Lagrange_Edge4.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Line2.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Line3.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Line4.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad16.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad4.hpp" //HMR/src
#include "cl_HMR_Lagrange_Facet_Quad9.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {

//-------------------------------------------------------------------------------

    /**
     * creates a backgroud mesh depending on the number of dimensions set
     */
    Background_Mesh_Base * Factory::create_background_mesh( const Parameters * aParameters )
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

        // reset main patterns of this mesh
        // fixme: this should be its own funcion
        aMesh->reset_pattern( aParameters->get_bspline_input_pattern() );
        aMesh->reset_pattern( aParameters->get_lagrange_input_pattern() );
        aMesh->reset_pattern( aParameters->get_bspline_output_pattern() );
        aMesh->reset_pattern( aParameters->get_lagrange_output_pattern() );

        return aMesh;
    }
//-------------------------------------------------------------------------------

        Lagrange_Mesh_Base * Factory::create_lagrange_mesh( const Parameters                  * aParameters,
                                                                  Background_Mesh_Base        * aBackgroundMesh,
                                                                  Cell< BSpline_Mesh_Base * >   aBSplineMeshes,
                                                            const uint                        & aActivationPattern,
                                                            const luint                       & aPolynomialDegree )
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
                    aMesh = new Lagrange_Mesh< 2, 1 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 2 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new Lagrange_Mesh< 2, 3 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
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
                    aMesh = new Lagrange_Mesh< 3, 1 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 2 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new Lagrange_Mesh< 3, 3 >( aParameters,
                                                       aBackgroundMesh,
                                                       aBSplineMeshes,
                                                       aActivationPattern );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "create_lagrange_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
                    break;
                }
                }
                break;
            }
            default :
            {
                MORIS_ERROR( false, "create_lagrange_mesh(): unknown number of dimensions %u\n",
                        ( unsigned int )  tNumberOfDimensions );
                break;
            }
            }

            return aMesh;
        }
//-------------------------------------------------------------------------------

        BSpline_Mesh_Base * Factory::create_bspline_mesh( const Parameters           * aParameters,
                                                                Background_Mesh_Base * aBackgroundMesh,
                                                          const uint                 & aActivationPattern,
                                                          const luint                & aPolynomialDegree )
        {
            BSpline_Mesh_Base * aMesh;

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
                    aMesh = new BSpline_Mesh< 2, 1 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new BSpline_Mesh< 2, 2 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new BSpline_Mesh< 2, 3 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new BSpline_Mesh< 2, 4 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new BSpline_Mesh< 2, 5 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "create_bspline_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
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
                    aMesh = new BSpline_Mesh< 3, 1 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 2 ):
                {
                    aMesh = new BSpline_Mesh< 3, 2 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 3 ):
                {
                    aMesh = new BSpline_Mesh< 3, 3 >( aParameters,
                                                     aBackgroundMesh,
                                                     aActivationPattern );
                    break;
                }
                case( 4 ):
                {
                    aMesh = new BSpline_Mesh< 3, 4 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                case( 5 ):
                {
                    aMesh = new BSpline_Mesh< 3, 5 >( aParameters,
                                                      aBackgroundMesh,
                                                      aActivationPattern );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "create_bspline_mesh(): unsupported polynomial degree %u for dimension %u\n",
                            ( unsigned int )  aPolynomialDegree,
                            ( unsigned int )  tNumberOfDimensions );
                    break;
                }
           }
           break;
           }
           default :
           {
               MORIS_ERROR( false, "create_bspline_mesh(): unknown number of dimensions %u\n",
                       ( unsigned int )  tNumberOfDimensions );
               break;
           }
           }

           return aMesh;
        }

//-------------------------------------------------------------------------------

    }  /* namespace hmr */
} /* namespace moris */

