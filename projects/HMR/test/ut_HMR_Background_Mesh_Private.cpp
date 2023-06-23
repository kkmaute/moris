/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Mesh_Private.cpp
 *
 */

#include <catch.hpp>

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Map.hpp" //CNT/srcf
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp"
#define protected public
#define private   public
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Background_Element.hpp" //HMR/src
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#undef protected
#undef private

TEST_CASE("HMR_Background_Mesh_Private", "[moris],[mesh],[hmr][Background_Mesh_private]")
{

//-------------------------------------------------------------------------------
    if(  moris::par_size() == 1  ||  moris::par_size() == 2  || moris::par_size() == 4 )
    {
//-------------------------------------------------------------------------------

        SECTION( "Background mesh 2D: test calc_child_index() function")
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {6}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer sizes to zero
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            tBackgroundMesh->collect_active_elements();

            // create alternating pedigree path
            moris::Matrix< moris::DDUMat > tInputTree( moris::hmr::gMaxNumberOfLevels-1, 1 );
            moris::uint tChildIndex = 0;
            for( moris::uint l=0; l<moris::hmr::gMaxNumberOfLevels-1; ++l )
            {
                tInputTree( l ) = tChildIndex++;
                if( tChildIndex == 4 )
                {
                    tChildIndex = 0;
                }
            }

            // get first active element on proc
            moris::hmr::Background_Element_Base * tElement = tBackgroundMesh->get_element( 0 );

            // refine element up to max possible level
            for( moris::uint l=0; l<moris::hmr::gMaxNumberOfLevels-1; ++l )
            {
                // call refinement procedure
                tBackgroundMesh->refine_element( tElement, false );

                // get child of element
                tElement = tElement->get_child( tInputTree( l ) );

                // get ij position of element
                const moris::luint* tIJ = tElement->get_ijk();

                // calculate child index
                tChildIndex = tBackgroundMesh->calc_child_index( tIJ[ 0 ], tIJ[ 1 ] );

                // make sure that index is correct
                REQUIRE( tChildIndex ==  tInputTree( l ) );
            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }

//-------------------------------------------------------------------------------

     SECTION( "Background mesh 3D: test calc_child_index() function")
     {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {6}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh
            = tFactory.create_background_mesh();

            // update element table
            tBackgroundMesh->collect_active_elements();

            // create alternating pedigree path
            moris::Matrix< moris::DDUMat > tInputTree( moris::hmr::gMaxNumberOfLevels-1, 1 );
            moris::uint tChildIndex = 0;
            for( moris::uint l=0; l<moris::hmr::gMaxNumberOfLevels-1; ++l )
            {
                tInputTree( l ) = tChildIndex++;
                if( tChildIndex == 8 )
                {
                    tChildIndex = 0;
                }
            }

            // get first active element on proc
            moris::hmr::Background_Element_Base * tElement = tBackgroundMesh->get_element( 0 );

            // refine element up to max possible level
            for( moris::uint l=0; l<moris::hmr::gMaxNumberOfLevels-1; ++l )
            {
                // call refinement procedure
                tBackgroundMesh->refine_element( tElement, false );

                // get child of element
                tElement = tElement->get_child( tInputTree( l ) );

                // get ijk position of element
                const moris::luint* tIJK = tElement->get_ijk();

                // calculate child index
                tChildIndex = tBackgroundMesh->calc_child_index( tIJK[ 0 ], tIJK[ 1 ], tIJK[ 2 ] );

                // make sure that index is correct
                REQUIRE( tChildIndex ==  tInputTree( l ) );
            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
//-------------------------------------------------------------------------------

       SECTION( "Background mesh 2D: test neighborhood calculation")
        {
            //
            // In this test, a 6x4 mesh is generated and uniformly refined
            // two times. We pick elements near the border of the proc domain,
            // manually calculate the subdomain ID and test if the neighbors
            // were found correctly.
            //

            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {6}, {4} };

            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            tBackgroundMesh->collect_active_elements();

            // max level to which we refine to
            moris::uint tLevel = 2;

            // refine the whole mesh two times
            for( moris::uint l=0; l<tLevel; ++l  )
            {
                auto tNumberOfElements = tBackgroundMesh->get_number_of_active_elements_on_proc();

                // refine all 16 elements
                for( moris::luint k=0; k<tNumberOfElements; ++k )
                {
                    // get element
                    moris::hmr::Background_Element_Base* tElement = tBackgroundMesh->get_element( k );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                // refine mesh
                tBackgroundMesh->perform_refinement( 0 );
            }

            // collect neighbors
            tBackgroundMesh->collect_neighbors();

            // get number of active elements
            auto tNumberOfActiveElements = tBackgroundMesh->get_number_of_active_elements_on_proc();

            tBackgroundMesh->save_to_vtk( "BackgorundMesh_6x4.vtk");

            // create element map so that we can access elements over global ID
            moris::map<moris::luint, moris::luint> tElementMap;

            // loop over all elements on proc
            for( moris::luint k=0; k<tNumberOfActiveElements; ++k )
            {
                // pick element
                moris::hmr::Background_Element_Base* tElement = tBackgroundMesh->get_element( k );

                // write ID into map
                tElementMap[ tElement->get_hmr_id() ] = k;
            }

            // now we pick some elements and store their local indices
            // in this matrix
            moris::Matrix< moris::DDLUMat > tElements;

            if( moris::par_size() == 1 )
            {
                // pick an element in the middle
                tElements.set_size( 1, 1, tElementMap.find( 1138 ) );
            }
            else if( moris::par_size() == 2 )
            {
                // pick two elements near a point where both procs meet
                tElements.set_size( 2, 1 );
                if( moris::par_rank() == 0 )
                {
                    tElements( 0 ) = tElementMap.find( 1138 );
                    tElements( 1 ) = tElementMap.find( 1139 );
                }
                else
                {
                    tElements( 0 ) = tElementMap.find( 1140 );
                    tElements( 1 ) = tElementMap.find( 1142 );
                }
            }
            else if( moris::par_size() == 4 )
            {
                // pick four elements at the corner where all four procs meet
                tElements.set_size( 4, 1 );
                if( moris::par_rank() == 0 )
                {
                    tElements( 0 ) = tElementMap.find(  978 );
                    tElements( 1 ) = tElementMap.find(  979 );
                    tElements( 2 ) = tElementMap.find( 1018 );
                    tElements( 3 ) = tElementMap.find( 1019 );
                }
                else if( moris::par_rank() == 1 )
                {
                    tElements( 0 ) = tElementMap.find( 1058 );
                    tElements( 1 ) = tElementMap.find( 1059 );
                    tElements( 2 ) = tElementMap.find( 1098 );
                    tElements( 3 ) = tElementMap.find( 1099 );
                }
                else if( moris::par_rank() == 2 )
                {
                    tElements( 0 ) = tElementMap.find(  980 );
                    tElements( 1 ) = tElementMap.find(  981 );
                    tElements( 2 ) = tElementMap.find( 1020 );
                    tElements( 3 ) = tElementMap.find( 1021 );
                }
                else
                {
                    tElements( 0 ) = tElementMap.find( 1060 );
                    tElements( 1 ) = tElementMap.find( 1061 );
                    tElements( 2 ) = tElementMap.find( 1100 );
                    tElements( 3 ) = tElementMap.find( 1101 );
                }
            }

            // loop over all selected elements
            for( moris::luint k=0; k<tElements.length(); ++k )
            {
                // pick an element
                moris::hmr::Background_Element_Base* tElement
                    = tBackgroundMesh->get_element( tElements( k ) );

                // get ijk position of element
                const moris::luint* tElIJ  = tElement->get_ijk();

                // calculate ij positions for neighbor elements
                moris::Matrix< moris::DDLUMat > tIJ( 2, 8 );

                // neighbor 0
                tIJ( 0, 0 ) = tElIJ[ 0 ];
                tIJ( 1, 0 ) = tElIJ[ 1 ] - 1;

                // neighbor 1
                tIJ( 0, 1 ) = tElIJ[ 0 ] + 1;
                tIJ( 1, 1 ) = tElIJ[ 1 ];

                // neighbor 2
                tIJ( 0, 2 ) = tElIJ[ 0 ];
                tIJ( 1, 2 ) = tElIJ[ 1 ] + 1;

                // neighbor 3
                tIJ( 0, 3 ) = tElIJ[ 0 ] - 1;
                tIJ( 1, 3 ) = tElIJ[ 1 ];

                // neighbor 4
                tIJ( 0, 4 ) = tElIJ[ 0 ] - 1;
                tIJ( 1, 4 ) = tElIJ[ 1 ] - 1;

                // neighbor 5
                tIJ( 0, 5 ) = tElIJ[ 0 ] + 1;
                tIJ( 1, 5 ) = tElIJ[ 1 ] - 1;

                // neighbor 6
                tIJ( 0, 6 ) = tElIJ[ 0 ] + 1;
                tIJ( 1, 6 ) = tElIJ[ 1 ] + 1;

                // neighbor 7
                tIJ( 0, 7 ) = tElIJ[ 0 ] - 1;
                tIJ( 1, 7 ) = tElIJ[ 1 ] + 1;

                // loop over all neighbors
                for( moris::uint e=0; e<8; ++e )
                {
                    // get neighbor of element
                    moris::hmr::Background_Element_Base * tNeighbor = tElement->get_neighbor( e );

                    // calculate expected id of element
                    moris::luint tDomainID = tBackgroundMesh->calc_domain_id_of_element( tLevel,
                                                                                         tIJ( 0 , e ),
                                                                                         tIJ( 1 , e ) );

                    // test if domain ID is correct
                    REQUIRE( tDomainID == tNeighbor->get_hmr_id() );
                }
            }

            // delete mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }

//-------------------------------------------------------------------------------

       SECTION( "Background mesh 3D: test neighborhood calculation")
       {
                // create settings object
                auto tParameters = new moris::hmr::Parameters;

                // set number of elements
                moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {2}, {3} };

                tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

                // set buffer size to zero
                tParameters->set_refinement_buffer( 0 );
                tParameters->set_staircase_buffer( 0 );

                // deactivate truncation
                tParameters->set_bspline_truncation( false );

                // create factory
                moris::hmr::Factory tFactory( tParameters );

                // create background mesh object
                moris::hmr::Background_Mesh_Base * tBackgroundMesh = tFactory.create_background_mesh();

                // update element table
                tBackgroundMesh->collect_active_elements();

                // max level to which we refine to
                moris::uint tLevel = 3;

                // refine the whole mesh three times
                for( moris::uint l=0; l<tLevel; ++l  )
                {
                    // ask background mesh for number of elements
                    auto tNumberOfElements =  tBackgroundMesh->get_number_of_active_elements_on_proc();

                    // refine every element
                    for( moris::luint k=0; k<tNumberOfElements; ++k )
                    {
                        // get element
                        moris::hmr::Background_Element_Base* tElement = tBackgroundMesh->get_element( k );

                        //flag element for refinement
                        tElement->put_on_refinement_queue();
                    }

                    // refine mesh
                    tBackgroundMesh->perform_refinement( 0 );
                }

                // ask background mesh for number of elements on proc
                auto tNumberOfElementsOnProc = tBackgroundMesh->get_number_of_active_elements_on_proc_including_aura();

                // pointer to element that is to be tested
                moris::hmr::Background_Element_Base * tElement = tBackgroundMesh->get_element_from_proc_domain_including_aura( 0 );

                // find element 40271 ( this is a bit brute force, but for the test, it is fast enough )
                // in parallel, this element has neighbors in the aura.

                for ( moris::luint k=0; k<tNumberOfElementsOnProc; ++k )
                {
                    tElement = tBackgroundMesh->get_element_from_proc_domain_including_aura( k );

                    if ( tElement->get_hmr_id() == 40271 )
                    {
                        break;
                    }
                }

                // calculate neighbor stencil
                const moris::luint* tElIJK = tElement->get_ijk();

                // extract i j and k
                moris::luint tI = tElIJK[ 0 ];
                moris::luint tJ = tElIJK[ 1 ];
                moris::luint tK = tElIJK[ 2 ];

                // calculate positions of neighbor elements
                moris::Matrix< moris::DDLUMat > tIJK( 3, 26 );

                // neighbor 0
                tIJK( 0,  0 ) = tI;
                tIJK( 1,  0 ) = tJ - 1;
                tIJK( 2,  0 ) = tK;

                // neighbor 1
                tIJK( 0,  1 ) = tI + 1;
                tIJK( 1,  1 ) = tJ;
                tIJK( 2,  1 ) = tK;

                // neighbor 2
                tIJK( 0,  2 ) = tI;
                tIJK( 1,  2 ) = tJ + 1;
                tIJK( 2,  2 ) = tK;

                // neighbor 3
                tIJK( 0,  3 ) = tI - 1;
                tIJK( 1,  3 ) = tJ;
                tIJK( 2,  3 ) = tK;

                // neighbor 4
                tIJK( 0,  4 ) = tI;
                tIJK( 1,  4 ) = tJ;
                tIJK( 2,  4 ) = tK - 1;

                // neighbor 5
                tIJK( 0,  5 ) = tI;
                tIJK( 1,  5 ) = tJ;
                tIJK( 2,  5 ) = tK + 1;

                // neighbor 6
                tIJK( 0,  6 ) = tI;
                tIJK( 1,  6 ) = tJ - 1;
                tIJK( 2,  6 ) = tK - 1;

                // neighbor 7
                tIJK( 0,  7 ) = tI + 1;
                tIJK( 1,  7 ) = tJ;
                tIJK( 2,  7 ) = tK - 1;

                // neighbor 8
                tIJK( 0,  8 ) = tI;
                tIJK( 1,  8 ) = tJ + 1;
                tIJK( 2,  8 ) = tK - 1;

                // neighbor 9
                tIJK( 0,  9 ) = tI - 1;
                tIJK( 1,  9 ) = tJ;
                tIJK( 2,  9 ) = tK - 1;

                // neighbor 10
                tIJK( 0, 10 ) = tI - 1;
                tIJK( 1, 10 ) = tJ - 1;
                tIJK( 2, 10 ) = tK;

                // neighbor 11
                tIJK( 0, 11 ) = tI + 1;
                tIJK( 1, 11 ) = tJ - 1;
                tIJK( 2, 11 ) = tK;

                // neighbor 12
                tIJK( 0, 12 ) = tI + 1;
                tIJK( 1, 12 ) = tJ + 1;
                tIJK( 2, 12 ) = tK;

                // neighbor 13
                tIJK( 0, 13 ) = tI - 1;
                tIJK( 1, 13 ) = tJ + 1;
                tIJK( 2, 13 ) = tK;

                // neighbor 14
                tIJK( 0, 14 ) = tI ;
                tIJK( 1, 14 ) = tJ - 1;
                tIJK( 2, 14 ) = tK + 1;

                // neighbor 15
                tIJK( 0, 15 ) = tI + 1;
                tIJK( 1, 15 ) = tJ;
                tIJK( 2, 15 ) = tK + 1;

                // neighbor 16
                tIJK( 0, 16 ) = tI;
                tIJK( 1, 16 ) = tJ + 1;
                tIJK( 2, 16 ) = tK + 1;

                // neighbor 17
                tIJK( 0, 17 ) = tI - 1;
                tIJK( 1, 17 ) = tJ;
                tIJK( 2, 17 ) = tK + 1;

                // neighbor 18
                tIJK( 0, 18 ) = tI - 1;
                tIJK( 1, 18 ) = tJ - 1;
                tIJK( 2, 18 ) = tK - 1;

                // neighbor 19
                tIJK( 0, 19 ) = tI + 1;
                tIJK( 1, 19 ) = tJ - 1;
                tIJK( 2, 19 ) = tK - 1;

                // neighbor 20
                tIJK( 0, 20 ) = tI + 1;
                tIJK( 1, 20 ) = tJ + 1;
                tIJK( 2, 20 ) = tK - 1;

                // neighbor 21
                tIJK( 0, 21 ) = tI - 1;
                tIJK( 1, 21 ) = tJ + 1;
                tIJK( 2, 21 ) = tK - 1;

                // neighbor 22
                tIJK( 0, 22 ) = tI - 1;
                tIJK( 1, 22 ) = tJ - 1;
                tIJK( 2, 22 ) = tK + 1;

                // neighbor 23
                tIJK( 0, 23 ) = tI + 1;
                tIJK( 1, 23 ) = tJ - 1;
                tIJK( 2, 23 ) = tK + 1;

                // neighbor 24
                tIJK( 0, 24 ) = tI + 1;
                tIJK( 1, 24 ) = tJ + 1;
                tIJK( 2, 24 ) = tK + 1;

                // neighbor 25
                tIJK( 0, 25 ) = tI - 1;
                tIJK( 1, 25 ) = tJ + 1;
                tIJK( 2, 25 ) = tK + 1;

                // loop over all neighbors
                for ( moris::uint e = 0; e<26; ++e )
                {
                   // get neighbor of element
                   moris::hmr::Background_Element_Base * tNeighbor = tElement->get_neighbor( e );

                    // calculate expected id of element
                    moris::luint tDomainID = tBackgroundMesh->calc_domain_id_of_element( tLevel,
                                                                                         tIJK( 0, e ),
                                                                                         tIJK( 1, e ),
                                                                                         tIJK( 2, e ) );

                    // test if domain ID is correct
                    REQUIRE( tDomainID == tNeighbor->get_hmr_id() );
                }

                // delete mesh
                delete tBackgroundMesh;

                // delete settings object
                delete tParameters;
        }
//-------------------------------------------------------------------------------
    }  // end if par_rank() == 1, 2, or 4

    if(  moris::par_size() == 1 )
    {
//-------------------------------------------------------------------------------

        SECTION( "Background mesh 2D: test get_neighbors_from_same_level() function")
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {6}, {6} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to one
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh
            = tFactory.create_background_mesh();

            // get number of elements
            moris::luint tNumberOfActiveElements
                = tBackgroundMesh->get_number_of_active_elements_on_proc();

            // refine all elements on first level
            for( moris::luint k=0; k<tNumberOfActiveElements; ++k )
            {
                tBackgroundMesh->get_element( k )->put_on_refinement_queue();
            }

            // perform refinement
            tBackgroundMesh->perform_refinement( 0 );

            // pick element in the middle
            moris::hmr::Background_Element_Base * tElement = tBackgroundMesh->get_element( 81 );

            moris::Cell< moris::hmr::Background_Element_Base * > tNeighbors;

            for( uint tOrder=1; tOrder<=2; ++tOrder )
            {
                // number of expected neighbors
                moris::luint tNumberOfNeighbors = std::pow( 2*tOrder+1, tParameters->get_number_of_dimensions() ) - 1;

                // get IJK of this element
                const moris::luint* tIJK = tElement->get_ijk();

                // initialize counter
                moris::luint tCount = 0;

                moris::Matrix< moris::DDLUMat > tExpectedIDs( tNumberOfNeighbors, 1 );

                // loop over all expected neighbors and calculate domain id
                for( uint j=tIJK[ 1 ]-tOrder; j<=tIJK[ 1 ]+tOrder; ++j )
                {
                    for( uint i=tIJK[ 0 ]-tOrder; i<=tIJK[ 0 ]+tOrder; ++i )
                    {
                        if( !( i == tIJK[ 0 ] && j == tIJK[ 1 ] ) )
                        {
                            tExpectedIDs( tCount++ ) = tBackgroundMesh->calc_domain_id_of_element( 1, i, j );
                        }
                    }
                }

                // sort
                moris::Matrix< moris::DDLUMat >tExpectedUniqueIDs;
                moris::unique( tExpectedIDs, tExpectedUniqueIDs );

                moris::Matrix< moris::DDLUMat > tCalculatedIDs( tNumberOfNeighbors, 1 );

                tElement->get_neighbors_from_same_level( tOrder, tNeighbors );

                // reset counter
                tCount = 0;

                // copy domain IDs from neighbors
                for( auto tNeighbor: tNeighbors )
                {
                    tCalculatedIDs( tCount++ ) = tNeighbor->get_hmr_id();
                }

                // test for correct size
                REQUIRE( tCount == tNumberOfNeighbors );

                // make result unique
                moris::Matrix< moris::DDLUMat >tCalculatedUniqueIDs;
                moris::unique( tCalculatedIDs, tCalculatedUniqueIDs );

                // make sure that IDs are correct
                bool tOK = true;
                for( moris::luint k=0; k<tNumberOfNeighbors; ++k )
                {
                    tOK = tOK &&  tCalculatedUniqueIDs( k ) == tExpectedUniqueIDs( k );
                }

                // request correct IDs
                REQUIRE( tOK );
            }
            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
//-------------------------------------------------------------------------------

       SECTION( "Background mesh 3D: test get_neighbors_from_same_level() function")
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {6}, {6}, {6} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to one
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh
            = tFactory.create_background_mesh();

            // get number of elements
            moris::luint tNumberOfActiveElements = tBackgroundMesh->get_number_of_active_elements_on_proc();

            // refine all elements on first level
            for( moris::luint k=0; k<tNumberOfActiveElements; ++k )
            {
                tBackgroundMesh->get_element( k )->put_on_refinement_queue();
            }

            // perform refinement
            tBackgroundMesh->perform_refinement( 0 );

            // pick an element in the middle
            moris::hmr::Background_Element_Base * tElement = tBackgroundMesh->get_element( 1025 );

            moris::Cell< moris::hmr::Background_Element_Base *  > tNeighbors;

            for( uint tOrder=1; tOrder<=2; ++tOrder )
            {
                // number of expected neighbors
                moris::luint tNumberOfNeighbors = std::pow( 2*tOrder+1, tParameters->get_number_of_dimensions() ) - 1;

                // get IJK of this element
                const moris::luint* tIJK = tElement->get_ijk();

                // initialize counter
                moris::luint tCount = 0;

                moris::Matrix< moris::DDLUMat > tExpectedIDs( tNumberOfNeighbors, 1 );

                // loop over all expected neighbors and calculate domain id
                for( uint k=tIJK[ 2 ]-tOrder; k<=tIJK[ 2 ]+tOrder; ++k )
                {
                    for( uint j=tIJK[ 1 ]-tOrder; j<=tIJK[ 1 ]+tOrder; ++j )
                    {
                        for( uint i=tIJK[ 0 ]-tOrder; i<=tIJK[ 0 ]+tOrder; ++i )
                        {
                            if( !( i == tIJK[ 0 ] && j == tIJK[ 1 ] && k == tIJK[ 2 ] ) )
                            {
                                tExpectedIDs( tCount++ ) = tBackgroundMesh->calc_domain_id_of_element( 1, i, j, k );
                            }
                        }
                    }
                }

                // sort
                moris::Matrix< moris::DDLUMat > tUniqueIDs;
                moris::unique( tExpectedIDs, tUniqueIDs );

                moris::Matrix< moris::DDLUMat > tCalculatedIDs( tNumberOfNeighbors, 1 );

                tElement->get_neighbors_from_same_level( tOrder, tNeighbors );

                // reset counter
                tCount = 0;

                // copy domain IDs from neighbors
                for( auto tNeighbor: tNeighbors )
                {
                    tCalculatedIDs( tCount++ ) = tNeighbor->get_hmr_id();
                }

                // test for correct size
                REQUIRE(  tCount == tNumberOfNeighbors );

                // make result unique
                moris::Matrix< moris::DDLUMat > tCalculatedUniqueIDs;
                moris::unique( tCalculatedIDs,  tCalculatedUniqueIDs );

                // make sure that IDs are correct
                bool tOK = true;
                for( moris::luint k=0; k<tNumberOfNeighbors; ++k )
                {
                    tOK = tOK &&  tCalculatedUniqueIDs( k ) == tUniqueIDs( k );
                }

                REQUIRE( tOK );
            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
    } // end serial case
//-------------------------------------------------------------------------------
} // end test

