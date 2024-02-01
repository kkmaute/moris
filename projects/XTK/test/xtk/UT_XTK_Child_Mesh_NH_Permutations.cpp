/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Child_Mesh_NH_Permutations.cpp
 *
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_plus.hpp"

#include "cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "fn_verify_tet_topology.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "fn_equal_to.hpp"

#include <chrono>
#include <thread>
#include <iostream>
#include "cl_XTK_Child_Mesh.hpp"

#include "fn_GEN_Triangle_Geometry.hpp"

namespace xtk
{
    class Permutations
    {
      public:
        Permutations( std::string tType = "" )

        {
            if ( tType.compare( "" ) == 0 )
            {
                this->load_all_permutations();
            }

            else if ( tType.compare( "pa" ) == 0 )
            {
                this->load_pos_a_permutations();
            }

            else if ( tType.compare( "na" ) == 0 )
            {
                this->load_neg_a_permutations();
            }

            else if ( tType.compare( "pb" ) == 0 )
            {
                this->load_pos_b_permutations();
            }
            else if ( tType.compare( "nb" ) == 0 )
            {
                this->load_neg_b_permutations();
            }

            else if ( tType.compare( "pc" ) == 0 )
            {
                this->load_pos_c_permutations();
            }

            else if ( tType.compare( "nc" ) == 0 )
            {
                this->load_neg_c_permutations();
            }
            else if ( tType.compare( "3p" ) == 0 )
            {
                this->load_pos_3_permutations();
            }
            else if ( tType.compare( "3n" ) == 0 )
            {
                this->load_neg_3_permutations();
            }
            else if ( tType.compare( "2_node" ) == 0 )
            {
                this->load_2_node_permutations();
            }
            else
            {
                std::cout << "String Not Recognized" << std::endl;
            }
            mNumPermutations = mPermutations.size();
        }

        size_t
        get_num_permutations()
        {
            return mNumPermutations;
        }

        Cell< size_t >
        get_permutation( size_t aIndex )
        {
            return mPermutations( aIndex );
        }

      private:
        size_t                 mNumPermutations;
        Cell< Cell< size_t > > mPermutations;

      private:
        void
        load_all_permutations()
        {
            mPermutations = {
                { 0, 2, 4, 5 },
                { 0, 2, 5, 4 },
                { 0, 4, 2, 5 },
                { 0, 4, 5, 2 },
                { 0, 5, 2, 4 },
                { 0, 5, 4, 2 },

                { 2, 0, 4, 5 },
                { 2, 0, 5, 4 },
                { 2, 4, 0, 5 },
                { 2, 4, 5, 0 },
                { 2, 5, 0, 4 },
                { 2, 5, 4, 0 },

                { 4, 0, 2, 5 },
                { 4, 0, 5, 2 },
                { 4, 2, 0, 5 },
                { 4, 2, 5, 0 },
                { 4, 5, 0, 2 },
                { 4, 5, 2, 0 },

                { 5, 0, 2, 4 },
                { 5, 0, 4, 2 },
                { 5, 2, 0, 4 },
                { 5, 2, 4, 0 },
                { 5, 4, 0, 2 },
                { 5, 4, 2, 0 },

                { 0, 1, 3, 5 },
                { 0, 1, 5, 3 },
                { 0, 3, 1, 5 },
                { 0, 3, 5, 1 },
                { 0, 5, 1, 3 },
                { 0, 5, 3, 1 },

                { 1, 0, 3, 5 },
                { 1, 0, 5, 3 },
                { 1, 3, 0, 5 },
                { 1, 3, 5, 0 },
                { 1, 5, 0, 3 },
                { 1, 5, 3, 0 },

                { 3, 0, 1, 5 },
                { 3, 0, 5, 1 },
                { 3, 1, 0, 5 },
                { 3, 1, 5, 0 },
                { 3, 5, 0, 1 },
                { 3, 5, 1, 0 },

                { 5, 0, 1, 3 },
                { 5, 0, 3, 1 },
                { 5, 1, 0, 3 },
                { 5, 1, 3, 0 },
                { 5, 3, 0, 1 },
                { 5, 3, 1, 0 },

                { 2, 1, 3, 4 },
                { 2, 1, 4, 3 },
                { 2, 3, 1, 4 },
                { 2, 3, 4, 1 },
                { 2, 4, 1, 3 },
                { 2, 4, 3, 1 },

                { 1, 2, 3, 4 },
                { 1, 2, 4, 3 },
                { 1, 3, 2, 4 },
                { 1, 3, 4, 2 },
                { 1, 4, 2, 3 },
                { 1, 4, 3, 2 },

                { 3, 2, 1, 4 },
                { 3, 2, 4, 1 },
                { 3, 1, 2, 4 },
                { 3, 1, 4, 2 },
                { 3, 4, 2, 1 },
                { 3, 4, 1, 2 },

                { 4, 2, 1, 3 },
                { 4, 2, 3, 1 },
                { 4, 1, 2, 3 },
                { 4, 1, 3, 2 },
                { 4, 3, 2, 1 },
                { 4, 3, 1, 2 }
            };
        }

        void
        load_pos_a_permutations()
        {

            mPermutations = { { 0, 2, 4, 5 },    // 5420
                { 0, 4, 2, 5 },                  // 5240
                { 5, 4, 2, 0 },                  // 245
                { 5, 2, 4, 0 },                  // 425
                { 1, 0, 5, 3 },                  // 3501
                { 3, 0, 5, 1 },                  // 1503
                { 1, 5, 0, 3 },                  // 3051
                { 3, 5, 0, 1 },                  // 1053
                { 2, 1, 3, 4 },                  // 4312
                { 4, 1, 3, 2 },                  // 2314
                { 2, 3, 1, 4 },                  // 4132
                { 4, 3, 1, 2 } };                // 2134
        }

        void
        load_neg_a_permutations()
        {
            mPermutations = { { 2, 0, 5, 4 },
                { 2, 5, 0, 4 },
                { 4, 0, 5, 2 },
                { 4, 5, 0, 2 },
                { 0, 1, 3, 5 },
                { 0, 3, 1, 5 },
                { 5, 3, 1, 0 },
                { 5, 1, 3, 0 },
                { 1, 2, 4, 3 },
                { 1, 4, 2, 3 },
                { 3, 2, 4, 1 },
                { 3, 4, 2, 1 } };
        }

        void
        load_pos_b_permutations()
        {
            mPermutations = { { 0, 5, 2, 4 },
                { 0, 5, 4, 2 },
                { 5, 0, 2, 4 },
                { 5, 0, 4, 2 },
                { 1, 3, 0, 5 },
                { 3, 1, 0, 5 },
                { 1, 3, 5, 0 },
                { 3, 1, 5, 0 },
                { 4, 2, 1, 3 },
                { 2, 4, 3, 1 },
                { 4, 2, 3, 1 },
                { 2, 4, 1, 3 } };
        }

        void
        load_neg_b_permutations()
        {
            mPermutations = { { 4, 2, 0, 5 },
                { 4, 2, 5, 0 },
                { 2, 4, 0, 5 },
                { 2, 4, 5, 0 },
                { 0, 5, 1, 3 },
                { 0, 5, 3, 1 },
                { 5, 0, 1, 3 },
                { 5, 0, 3, 1 },
                { 1, 3, 2, 4 },
                { 1, 3, 4, 2 },
                { 3, 1, 2, 4 },
                { 3, 1, 4, 2 } };
        }

        void
        load_pos_c_permutations()
        {

            mPermutations = { { 0, 2, 5, 4 },
                { 0, 4, 5, 2 },
                { 5, 2, 0, 4 },
                { 5, 4, 0, 2 },
                { 1, 0, 3, 5 },
                { 3, 0, 1, 5 },
                { 1, 5, 3, 0 },
                { 3, 5, 1, 0 },
                { 2, 1, 4, 3 },
                { 4, 1, 2, 3 },
                { 2, 3, 4, 1 },
                { 4, 3, 2, 1 } };
        }

        void
        load_neg_c_permutations()
        {
            mPermutations = { { 2, 0, 4, 5 },
                { 2, 5, 4, 0 },
                { 4, 0, 2, 5 },
                { 4, 5, 2, 0 },
                { 0, 1, 5, 3 },
                { 0, 3, 5, 1 },
                { 5, 1, 0, 3 },
                { 5, 3, 0, 1 },
                { 1, 2, 3, 4 },
                { 1, 4, 3, 2 },
                { 3, 2, 1, 4 },
                { 3, 4, 1, 2 } };
        }

        void
        load_pos_3_permutations()
        {
            mPermutations = { { 0, 2, 3, 0 },
                { 2, 3, 0, 0 },
                { 3, 0, 2, 0 },
                { 1, 5, 2, 0 },
                { 2, 1, 5, 0 },
                { 5, 2, 1, 0 },
                { 0, 4, 1, 0 },
                { 1, 0, 4, 0 },
                { 4, 1, 0, 0 },
                { 3, 5, 4, 0 },
                { 4, 3, 5, 0 },
                { 5, 4, 3, 0 } };
        }

        void
        load_neg_3_permutations()
        {

            mPermutations = { { 0, 3, 2, 0 },
                { 2, 0, 3, 0 },
                { 3, 2, 0, 0 },
                { 1, 2, 5, 0 },
                { 2, 5, 1, 0 },
                { 5, 1, 2, 0 },
                { 0, 1, 4, 0 },
                { 1, 4, 0, 0 },
                { 4, 0, 1, 0 },
                { 3, 4, 5, 0 },
                { 4, 5, 3, 0 },
                { 5, 3, 4, 0 } };
        }

        void
        load_2_node_permutations()
        {

            mPermutations = { { 1, 4 },
                { 4, 1 },
                { 1, 5 },
                { 5, 1 },
                { 4, 5 },
                { 5, 4 },
                { 5, 3 },
                { 3, 5 },
                { 1, 2 },
                { 2, 1 },
                { 0, 2 },
                { 2, 0 },
                { 0, 1 },
                { 1, 0 },
                { 0, 3 },
                { 3, 0 },
                { 0, 4 },
                { 4, 0 },
                { 3, 4 },
                { 4, 3 },
                { 2, 5 },
                { 5, 2 },
                { 2, 3 },
                { 3, 2 } };
        }
    };

    void
    get_midside_coordinate( moris::moris_index const &aEdgeIndex,
            moris::Matrix< moris::DDRMat >           &aMidEdgeCoordinate )
    {
        if ( aEdgeIndex == 0 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.5;
            aMidEdgeCoordinate( 0, 1 ) = 0.0;
            aMidEdgeCoordinate( 0, 2 ) = 0.0;
        }

        else if ( aEdgeIndex == 1 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.5;
            aMidEdgeCoordinate( 0, 1 ) = 0.5;
            aMidEdgeCoordinate( 0, 2 ) = 0.0;
        }

        else if ( aEdgeIndex == 2 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.0;
            aMidEdgeCoordinate( 0, 1 ) = 0.5;
            aMidEdgeCoordinate( 0, 2 ) = 0.0;
        }

        else if ( aEdgeIndex == 3 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.0;
            aMidEdgeCoordinate( 0, 1 ) = 0.0;
            aMidEdgeCoordinate( 0, 2 ) = 0.5;
        }

        else if ( aEdgeIndex == 4 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.5;
            aMidEdgeCoordinate( 0, 1 ) = 0.0;
            aMidEdgeCoordinate( 0, 2 ) = 0.5;
        }

        else if ( aEdgeIndex == 5 )
        {
            aMidEdgeCoordinate( 0, 0 ) = 0.0;
            aMidEdgeCoordinate( 0, 1 ) = 0.5;
            aMidEdgeCoordinate( 0, 2 ) = 0.5;
        }
        else
        {
            std::cout << "UNDEFINED EDGE" << std::endl;
        }
    }

    void
    setup_node_coordinates_3_node( moris::moris_index const &tEdgeL,
            moris::moris_index const                        &tEdgeM,
            moris::moris_index const                        &tEdgeH,
            moris::Matrix< moris::DDRMat >                  &aNodeCoordinates )
    {
        moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates( 1, 3, 75 );

        aNodeCoordinates         = moris::Matrix< moris::DDRMat >( 7, 3 );
        aNodeCoordinates( 0, 0 ) = 0.0;
        aNodeCoordinates( 0, 1 ) = 0.0;
        aNodeCoordinates( 0, 2 ) = 0.0;
        aNodeCoordinates( 1, 0 ) = 1.0;
        aNodeCoordinates( 1, 1 ) = 0.0;
        aNodeCoordinates( 1, 2 ) = 0.0;
        aNodeCoordinates( 2, 0 ) = 0.0;
        aNodeCoordinates( 2, 1 ) = 1.0;
        aNodeCoordinates( 2, 2 ) = 0.0;
        aNodeCoordinates( 3, 0 ) = 0.0;
        aNodeCoordinates( 3, 1 ) = 0.0;
        aNodeCoordinates( 3, 2 ) = 1.0;

        get_midside_coordinate( tEdgeL, tEdgeNodeCoordinates );
        aNodeCoordinates( 4, 0 ) = tEdgeNodeCoordinates( 0, 0 );
        aNodeCoordinates( 4, 1 ) = tEdgeNodeCoordinates( 0, 1 );
        aNodeCoordinates( 4, 2 ) = tEdgeNodeCoordinates( 0, 2 );

        get_midside_coordinate( tEdgeM, tEdgeNodeCoordinates );
        aNodeCoordinates( 5, 0 ) = tEdgeNodeCoordinates( 0, 0 );
        aNodeCoordinates( 5, 1 ) = tEdgeNodeCoordinates( 0, 1 );
        aNodeCoordinates( 5, 2 ) = tEdgeNodeCoordinates( 0, 2 );

        get_midside_coordinate( tEdgeH, tEdgeNodeCoordinates );
        aNodeCoordinates( 6, 0 ) = tEdgeNodeCoordinates( 0, 0 );
        aNodeCoordinates( 6, 1 ) = tEdgeNodeCoordinates( 0, 1 );
        aNodeCoordinates( 6, 2 ) = tEdgeNodeCoordinates( 0, 2 );
    }

    void
    setup_node_coordinates_4_node( size_t const &tEdgeL,
            size_t const                        &tEdgeML,
            size_t const                        &tEdgeMH,
            size_t const                        &tEdgeH,
            moris::Matrix< moris::DDRMat >      &aNodeCoordinates )
    {
        moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates( 1, 3, 75 );

        aNodeCoordinates         = moris::Matrix< moris::DDRMat >( 8, 3 );
        aNodeCoordinates( 0, 0 ) = 0.0;
        aNodeCoordinates( 0, 1 ) = 0.0;
        aNodeCoordinates( 0, 2 ) = 0.0;
        aNodeCoordinates( 1, 0 ) = 1.0;
        aNodeCoordinates( 1, 1 ) = 0.0;
        aNodeCoordinates( 1, 2 ) = 0.0;
        aNodeCoordinates( 2, 0 ) = 0.0;
        aNodeCoordinates( 2, 1 ) = 1.0;
        aNodeCoordinates( 2, 2 ) = 0.0;
        aNodeCoordinates( 3, 0 ) = 0.0;
        aNodeCoordinates( 3, 1 ) = 0.0;
        aNodeCoordinates( 3, 2 ) = 1.0;

        get_midside_coordinate( tEdgeL, tEdgeNodeCoordinates );
        aNodeCoordinates( 4, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 4, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 4, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );

        get_midside_coordinate( tEdgeML, tEdgeNodeCoordinates );
        aNodeCoordinates( 5, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 5, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 5, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );

        get_midside_coordinate( tEdgeMH, tEdgeNodeCoordinates );
        aNodeCoordinates( 6, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 6, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 6, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );

        get_midside_coordinate( tEdgeH, tEdgeNodeCoordinates );
        aNodeCoordinates( 7, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 7, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 7, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );
    }

    void
    setup_node_coordinates_2_node( size_t const &tEdgeL,
            size_t const                        &tEdgeH,
            moris::Matrix< moris::DDRMat >      &aNodeCoordinates )
    {
        moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates( 1, 3, 75 );

        aNodeCoordinates         = moris::Matrix< moris::DDRMat >( 6, 3 );
        aNodeCoordinates( 0, 0 ) = 0.0;
        aNodeCoordinates( 0, 1 ) = 0.0;
        aNodeCoordinates( 0, 2 ) = 0.0;
        aNodeCoordinates( 1, 0 ) = 1.0;
        aNodeCoordinates( 1, 1 ) = 0.0;
        aNodeCoordinates( 1, 2 ) = 0.0;
        aNodeCoordinates( 2, 0 ) = 0.0;
        aNodeCoordinates( 2, 1 ) = 1.0;
        aNodeCoordinates( 2, 2 ) = 0.0;
        aNodeCoordinates( 3, 0 ) = 0.0;
        aNodeCoordinates( 3, 1 ) = 0.0;
        aNodeCoordinates( 3, 2 ) = 1.0;

        get_midside_coordinate( tEdgeL, tEdgeNodeCoordinates );
        aNodeCoordinates( 4, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 4, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 4, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );

        get_midside_coordinate( tEdgeH, tEdgeNodeCoordinates );
        aNodeCoordinates( 5, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 5, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 5, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );
    }

    void
    setup_node_coordinates_bisected( size_t const &tEdgeOrd,
            moris::Matrix< moris::DDRMat >        &aNodeCoordinates )
    {
        moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates( 1, 3, 75 );

        aNodeCoordinates         = moris::Matrix< moris::DDRMat >( 8, 3 );
        aNodeCoordinates( 0, 0 ) = 0.0;
        aNodeCoordinates( 0, 1 ) = 0.0;
        aNodeCoordinates( 0, 2 ) = 0.0;
        aNodeCoordinates( 1, 0 ) = 1.0;
        aNodeCoordinates( 1, 1 ) = 0.0;
        aNodeCoordinates( 1, 2 ) = 0.0;
        aNodeCoordinates( 2, 0 ) = 0.0;
        aNodeCoordinates( 2, 1 ) = 1.0;
        aNodeCoordinates( 2, 2 ) = 0.0;
        aNodeCoordinates( 3, 0 ) = 0.0;
        aNodeCoordinates( 3, 1 ) = 0.0;
        aNodeCoordinates( 3, 2 ) = 1.0;

        get_midside_coordinate( tEdgeOrd, tEdgeNodeCoordinates );
        aNodeCoordinates( 4, 0 ) = (tEdgeNodeCoordinates)( 0, 0 );
        aNodeCoordinates( 4, 1 ) = (tEdgeNodeCoordinates)( 0, 1 );
        aNodeCoordinates( 4, 2 ) = (tEdgeNodeCoordinates)( 0, 2 );
    }

    bool
    verify_edge_is_on_parent_edge( moris::Matrix< moris::IndexMat > const &aParentEdgeNodes,
            moris::Matrix< moris::DDRMat >                                &aNodeCoordinates,
            moris::Matrix< moris::IndexMat >                              &aEdgeNodes )
    {
        bool tEdgeIsOnParent = false;

        // Collect information for two point line form in 3d
        moris::real tX1 = aNodeCoordinates( aParentEdgeNodes( 0 ), 0 );
        moris::real tY1 = aNodeCoordinates( aParentEdgeNodes( 0 ), 1 );
        moris::real tZ1 = aNodeCoordinates( aParentEdgeNodes( 0 ), 2 );
        moris::real tX2 = aNodeCoordinates( aParentEdgeNodes( 1 ), 0 );
        moris::real tY2 = aNodeCoordinates( aParentEdgeNodes( 1 ), 1 );
        moris::real tZ2 = aNodeCoordinates( aParentEdgeNodes( 1 ), 2 );

        // Avg value
        moris::real tXAvg = ( tX1 + tX2 ) / 2;
        moris::real tYAvg = ( tY1 + tY2 ) / 2;
        moris::real tZAvg = ( tZ1 + tZ2 ) / 2;

        //    std::cout<<"X1 = " <<tX1<<","<<tY1<<","<<tZ1<<std::endl;
        //    std::cout<<"X2 = " <<tX2<<","<<tY2<<","<<tZ2<<std::endl;

        moris::real tVal = 0.0;
        for ( moris::uint i = 0; i < aEdgeNodes.numel(); i++ )
        {
            moris::real tX = aNodeCoordinates( aEdgeNodes( i ), 0 );
            moris::real tY = aNodeCoordinates( aEdgeNodes( i ), 1 );
            moris::real tZ = aNodeCoordinates( aEdgeNodes( i ), 2 );
            //        std::cout<<"X = " <<tX<<","<<tY<<","<<tZ<<std::endl;

            if ( moris::equal_to( tX1, tX ) && moris::equal_to( tY1, tY ) && moris::equal_to( tZ1, tZ ) )
            {
                continue;
            }

            else if ( moris::equal_to( tX2, tX ) && moris::equal_to( tY2, tY ) && moris::equal_to( tZ2, tZ ) )
            {
                continue;
            }

            else if ( moris::equal_to( tXAvg, tX ) && moris::equal_to( tYAvg, tY ) && moris::equal_to( tZAvg, tZ ) )
            {
                continue;
            }
            else
            {
                tVal = tVal + 10.0;
            }
        }

        if ( moris::equal_to( tVal, 0 ) )
        {
            tEdgeIsOnParent = true;
        }

        return tEdgeIsOnParent;
    }

    TEST_CASE( "Node Hierarchy Template 3 Node Case Permutations", "[3_NODE]" )
    {
        // Tests:
        // Floodfill which checks whether the element to element connectivity is traversable
        // Topology which checks whether the new tets have the correct topology
        // Face ancestry check using surface normals of tri
        // Number of child faces created on parent face

        moris::Matrix< moris::IndexMat > tNodeIndex( { { 0, 1, 2, 3 } } );

        moris::Matrix< moris::IdMat > tNodeIds( { { 1, 2, 3, 4, 5, 6, 7 } } );

        moris::Matrix< moris::IndexMat > tElementsAncestry( { { 0 } } );    // Not used
        moris::Matrix< moris::DDSTMat >  tElementNodeParentRanks( 1, 4, 0 );
        moris::Matrix< moris::IndexMat > tParentEdgeInds( { { 0, 1, 2, 3, 4, 5 } } );
        moris::Matrix< moris::DDSTMat >  tParentEdgeRanks( 1, 6, 1 );
        moris::Matrix< moris::IndexMat > tParentFaceInds( { { 0, 1, 2, 3 } } );
        moris::Matrix< moris::DDSTMat >  tParentFaceRanks( 1, 4, 2 );

        Cell< std::string > tCaseStrings = { "3p", "3n" };
        for ( size_t iCase = 0; iCase < tCaseStrings.size(); iCase++ )
        {
            /*
             * Initialize Permutations
             */
            std::string &tCases = tCaseStrings( iCase );
            Permutations tPermutations( tCases );
            size_t       tNumPermutations = tPermutations.get_num_permutations();
            for ( size_t iPerm = 0; iPerm < tNumPermutations; iPerm++ )
            {
                // Initialize Template
                Mesh_Modification_Template tMeshTemplate( tElementsAncestry( 0, 0 ),
                        0,
                        tNodeIndex,
                        tNodeIndex,
                        tElementNodeParentRanks,
                        tParentEdgeInds,
                        tParentEdgeRanks,
                        tParentFaceInds,
                        tParentFaceRanks,
                        TemplateType::TET_4 );

                // Initialize child mesh with template (in this case a tet4)
                Child_Mesh tChildMesh( tMeshTemplate );

                tChildMesh.add_new_geometry_interface( 0 );

                // add new node indices
                tChildMesh.add_node_indices( { { 4, 5, 6 } } );
                tChildMesh.add_node_ids( tNodeIds );

                // select template
                Cell< size_t >     tCurrentPermutation = tPermutations.get_permutation( iPerm );
                moris::moris_index tEdgeL              = tCurrentPermutation( 0 );
                moris::moris_index tEdgeM              = tCurrentPermutation( 1 );
                moris::moris_index tEdgeH              = tCurrentPermutation( 2 );

                // Set up node coordinates
                moris::Matrix< moris::DDRMat > tNodeCoords;
                setup_node_coordinates_3_node( tEdgeL, tEdgeM, tEdgeH, tNodeCoords );

                // Compute base tet volume
                real tTetVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );

                // Compute base element surface normals (parent faces)
                moris::Matrix< moris::IndexMat > const &tParentFaceToNode = tChildMesh.get_face_to_node();
                size_t                                  tNumParentFaces   = tParentFaceToNode.n_rows();
                moris::Matrix< moris::DDRMat >          tParentFaceNormals( 3, tNumParentFaces );

                // Iterate through and compute all face normals
                for ( size_t iF = 0; iF < tNumParentFaces; iF++ )
                {
                    // Get face iF nodes
                    moris::Matrix< moris::IndexMat > tFaceNodes = tParentFaceToNode.get_row( iF );

                    moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                    gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                    tParentFaceNormals.set_column( iF, tFaceNormal );
                }

                tChildMesh.add_entity_to_intersect_connectivity( 4, tEdgeL, true );
                tChildMesh.add_entity_to_intersect_connectivity( 5, tEdgeM, true );
                tChildMesh.add_entity_to_intersect_connectivity( 6, tEdgeH, true );

                tChildMesh.modify_child_mesh( TemplateType::HIERARCHY_TET4 );

                // Verify that if we set each element to the same bulk phase we can traverse the element to element connectivity
                moris::moris_index               tMax       = std::numeric_limits< moris::moris_index >::max();
                moris::size_t                    tNumPhases = 2;
                moris::Matrix< moris::IndexMat > tActiveElements( { { 0, 1, 2, 3 } } );
                moris::Matrix< moris::IndexMat > tIncludedElementMarker( 1, 4, 1 );
                moris::Matrix< moris::IndexMat > tElementPhase( 1, 4, 0 );
                moris::moris_index               tMaxFloodFill = 0;
                moris::Matrix< moris::IndexMat > tElementSubphase =
                        flood_fill( tChildMesh.get_element_to_element(),
                                tElementPhase,
                                tActiveElements,
                                tIncludedElementMarker,
                                tNumPhases,
                                tMax,
                                tMaxFloodFill,
                                true );

                moris::Matrix< moris::IndexMat > tExpElementSubphase( 4, 1, 0 );
                CHECK( equal_to( tExpElementSubphase, tElementSubphase ) );

                // Verify that the tets created have correct topology
                bool tValidTopo = verify_tet4_topology( tChildMesh.get_element_to_node(),
                        tChildMesh.get_element_to_edge(),
                        tChildMesh.get_element_to_face(),
                        tChildMesh.get_edge_to_node(),
                        tChildMesh.get_face_to_node() );

                CHECK( tValidTopo );

                // verify volume is conserved
                real tTotalChildVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );
                ;

                CHECK( approximate( tTetVol, tTotalChildVol ) );

                // Check ancestry of faces
                moris::Matrix< moris::IndexMat > const &tFaceToNode      = tChildMesh.get_face_to_node();
                moris::Matrix< moris::IndexMat > const &tFaceParentInds  = tChildMesh.get_face_parent_inds();
                moris::Matrix< moris::DDSTMat > const  &tFaceParentRanks = tChildMesh.get_face_parent_ranks();
                size_t                                  tNumFaces        = tFaceToNode.n_rows();
                moris::Matrix< moris::DDRMat >          tFaceNormals( 3, tNumFaces );

                // Iterate through and compute all face normals
                for ( size_t iF = 0; iF < tNumFaces; iF++ )
                {
                    // Get face iF nodes
                    moris::Matrix< moris::IndexMat > tFaceNodes = tFaceToNode.get_row( iF );

                    moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                    gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                    tFaceNormals.set_column( iF, tFaceNormal );
                }

                // Check to see if the surface normal matches parent face normals
                // Count for number of children faces born on a parent face
                size_t tChildFacewithParentFaceRank = 0;
                for ( size_t iF = 0; iF < tNumFaces; iF++ )
                {
                    if ( tFaceParentRanks( 0, iF ) == 2 )
                    {
                        moris::Matrix< moris::DDRMat > tChildFaceNormal  = tFaceNormals.get_column( iF );
                        moris::Matrix< moris::DDRMat > tParentFaceNormal = tParentFaceNormals.get_column( tFaceParentInds( 0, iF ) );
                        tChildFacewithParentFaceRank++;
                        CHECK( equal_to( tChildFaceNormal, tParentFaceNormal ) );
                    }
                }
                CHECK( tChildFacewithParentFaceRank == 10 );

                // Check Edge Ancestry
                moris::mtk::Cell_Info_Tet4              tConn;
                moris::Matrix< moris::IndexMat >        tParentEdgeToNodeMap = tConn.get_node_to_edge_map();
                moris::Matrix< moris::IndexMat > const &tEdgeToNode          = tChildMesh.get_edge_to_node();
                moris::Matrix< moris::IndexMat > const &tEdgeParentInds      = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat > const  &tEdgeParentRanks     = tChildMesh.get_edge_parent_ranks();
                for ( size_t iEdge = 0; iEdge < tEdgeToNode.n_rows(); iEdge++ )
                {
                    // Verify all edges with a parent rank of 1 have nodes which belong on that edge
                    if ( tEdgeParentRanks( iEdge ) == 1 )
                    {
                        moris::Matrix< moris::IndexMat > tParentEdgeNodes = tParentEdgeToNodeMap.get_row( tEdgeParentInds( iEdge ) );
                        moris::Matrix< moris::IndexMat > tChildEdgeNodes  = tEdgeToNode.get_row( iEdge );

                        CHECK( verify_edge_is_on_parent_edge( tParentEdgeNodes, tNodeCoords, tChildEdgeNodes ) );
                    }
                }
            }
        }
    }

    TEST_CASE( "Node Hierarchy Template 4 Node Case Permutations", "[4_NODE]" )
    {
        // Tests:
        // Floodfill which checks whether the element to element connectivity is traversable
        // Topology which checks whether the new tets have the correct topology

        moris::Matrix< moris::IndexMat > tNodeIndex( { { 0, 1, 2, 3 } } );

        moris::Matrix< moris::IdMat > tNodeIds( { { 1, 2, 3, 4, 5, 6, 7, 8 } } );

        moris::Matrix< moris::IndexMat > tElementsAncestry( { { 0 } } );    // Not used
        moris::Matrix< moris::DDSTMat >  tElementNodeParentRanks( 1, 4, 0 );
        moris::Matrix< moris::IndexMat > tParentEdgeInds( { { 0, 1, 2, 3, 4, 5 } } );
        moris::Matrix< moris::DDSTMat >  tParentEdgeRanks( 1, 6, 1 );
        moris::Matrix< moris::IndexMat > tParentFaceInds( { { 0, 1, 2, 3 } } );
        moris::Matrix< moris::DDSTMat >  tParentFaceRanks( 1, 4, 2 );

        Cell< std::string > tCaseNames = { "pa", "na", "pb", "nb", "pc", "nc" };
        for ( size_t iCase = 0; iCase < tCaseNames.size(); iCase++ )
        {
            /*
             * Initialize Permutations
             */
            std::string &tCases = tCaseNames( iCase );
            Permutations tPermutations( tCases );
            size_t       tNumPermutations = tPermutations.get_num_permutations();

            for ( size_t iPerm = 0; iPerm < tNumPermutations; iPerm++ )
            {

                // Initialize Template
                Mesh_Modification_Template tMeshTemplate( tElementsAncestry( 0, 0 ),
                        0,
                        tNodeIndex,
                        tNodeIndex,
                        tElementNodeParentRanks,
                        tParentEdgeInds,
                        tParentEdgeRanks,
                        tParentFaceInds,
                        tParentFaceRanks,
                        TemplateType::TET_4 );

                // Initialize child mesh with template (in this case a tet4)
                Child_Mesh tChildMesh( tMeshTemplate );

                tChildMesh.add_new_geometry_interface( 0 );

                // add new node indices
                tChildMesh.add_node_indices( { { 4, 5, 6, 7 } } );
                tChildMesh.add_node_ids( tNodeIds );

                // select template
                Cell< size_t >     tCurrentPermutation = tPermutations.get_permutation( iPerm );
                moris::moris_index tEdgeL              = tCurrentPermutation( 0 );
                moris::moris_index tEdgeML             = tCurrentPermutation( 1 );
                moris::moris_index tEdgeMH             = tCurrentPermutation( 2 );
                moris::moris_index tEdgeH              = tCurrentPermutation( 3 );

                // Set up node coordinates
                moris::Matrix< moris::DDRMat > tNodeCoords;
                setup_node_coordinates_4_node( tEdgeL, tEdgeML, tEdgeMH, tEdgeH, tNodeCoords );
                // Compute base tet volume
                real tTetVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );

                // Compute base element surface normals (parent faces)
                moris::Matrix< moris::IndexMat > const &tParentFaceToNode = tChildMesh.get_face_to_node();
                size_t                                  tNumParentFaces   = tParentFaceToNode.n_rows();
                moris::Matrix< moris::DDRMat >          tParentFaceNormals( 3, tNumParentFaces );

                // Iterate through and compute all face normals
                for ( size_t iF = 0; iF < tNumParentFaces; iF++ )
                {
                    // Get face iF nodes
                    moris::Matrix< moris::IndexMat > tFaceNodes = tParentFaceToNode.get_row( iF );

                    moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                    gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                    tParentFaceNormals.set_column( iF, tFaceNormal );
                }

                tChildMesh.add_entity_to_intersect_connectivity( 4, tEdgeL, true );
                tChildMesh.add_entity_to_intersect_connectivity( 5, tEdgeML, true );
                tChildMesh.add_entity_to_intersect_connectivity( 6, tEdgeMH, true );
                tChildMesh.add_entity_to_intersect_connectivity( 7, tEdgeH, true );

                tChildMesh.modify_child_mesh( TemplateType::HIERARCHY_TET4 );

                // Verify that if we set each element to the same bulk phase we can traverse the element to element connectivity
                moris::moris_index               tMax       = std::numeric_limits< moris::moris_index >::max();
                size_t                           tNumPhases = 2;
                moris::Matrix< moris::IndexMat > tActiveElements( { { 0, 1, 2, 3, 4, 5 } } );
                moris::Matrix< moris::IndexMat > tIncludedElementMarker( 1, 6, 1 );
                moris::Matrix< moris::IndexMat > tElementPhase( 1, 6, 0 );
                moris::moris_index               tMaxFloodFill = 0;
                moris::Matrix< moris::IndexMat > tElementSubphase =
                        flood_fill( tChildMesh.get_element_to_element(),
                                tElementPhase,
                                tActiveElements,
                                tIncludedElementMarker,
                                tNumPhases,
                                tMax,
                                tMaxFloodFill,
                                true );

                moris::Matrix< moris::IndexMat > tExpElementSubphase( 6, 1, 0 );
                CHECK( equal_to( tExpElementSubphase, tElementSubphase ) );

                // Verify that the tets created have correct topology
                bool tValidTopo = verify_tet4_topology( tChildMesh.get_element_to_node(),
                        tChildMesh.get_element_to_edge(),
                        tChildMesh.get_element_to_face(),
                        tChildMesh.get_edge_to_node(),
                        tChildMesh.get_face_to_node() );

                CHECK( tValidTopo );

                // verify volume is conserved
                real tTotalChildVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );
                ;

                CHECK( approximate( tTetVol, tTotalChildVol ) );

                // Check ancestry of faces
                moris::Matrix< moris::IndexMat > const &tFaceToNode      = tChildMesh.get_face_to_node();
                moris::Matrix< moris::IndexMat > const &tFaceParentInds  = tChildMesh.get_face_parent_inds();
                moris::Matrix< moris::DDSTMat > const  &tFaceParentRanks = tChildMesh.get_face_parent_ranks();
                size_t                                  tNumFaces        = tFaceToNode.n_rows();
                moris::Matrix< moris::DDRMat >          tFaceNormals( 3, tNumFaces );

                // Iterate through and compute all face normals
                for ( size_t iF = 0; iF < tNumFaces; iF++ )
                {
                    // Get face iF nodes
                    moris::Matrix< moris::IndexMat > tFaceNodes = tFaceToNode.get_row( iF );

                    moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                    gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                    tFaceNormals.set_column( iF, tFaceNormal );
                }

                // Check to see if the surface normal matches parent face normals
                // Count for number of children faces born on a parent face
                size_t tChildFacewithParentFaceRank = 0;
                for ( size_t iF = 0; iF < tNumFaces; iF++ )
                {
                    if ( tFaceParentRanks( 0, iF ) == 2 )
                    {
                        moris::Matrix< moris::DDRMat > tChildFaceNormal  = tFaceNormals.get_column( iF );
                        moris::Matrix< moris::DDRMat > tParentFaceNormal = tParentFaceNormals.get_column( tFaceParentInds( 0, iF ) );
                        tChildFacewithParentFaceRank++;
                        CHECK( equal_to( tChildFaceNormal, tParentFaceNormal ) );
                    }
                }
                CHECK( tChildFacewithParentFaceRank == 12 );

                // Check Edge Ancestry
                moris::mtk::Cell_Info_Tet4              tConn;
                moris::Matrix< moris::IndexMat >        tParentEdgeToNodeMap = tConn.get_node_to_edge_map();
                moris::Matrix< moris::IndexMat > const &tEdgeToNode          = tChildMesh.get_edge_to_node();
                moris::Matrix< moris::IndexMat > const &tEdgeParentInds      = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat > const  &tEdgeParentRanks     = tChildMesh.get_edge_parent_ranks();
                for ( size_t iEdge = 0; iEdge < tEdgeToNode.n_rows(); iEdge++ )
                {
                    // Verify all edges with a parent rank of 1 have nodes which belong on that edge
                    if ( tEdgeParentRanks( iEdge ) == 1 )
                    {
                        moris::Matrix< moris::IndexMat > tParentEdgeNodes = tParentEdgeToNodeMap.get_row( tEdgeParentInds( iEdge ) );
                        moris::Matrix< moris::IndexMat > tChildEdgeNodes  = tEdgeToNode.get_row( iEdge );

                        CHECK( verify_edge_is_on_parent_edge( tParentEdgeNodes, tNodeCoords, tChildEdgeNodes ) );
                    }
                }
            }
        }
    }

    TEST_CASE( "Bisected Tetrahedral Template", "[BISECT_TEMPLATE]" )
    {
        moris::Matrix< moris::IndexMat > tNodeIndex( { { 0, 1, 2, 3 } } );

        moris::Matrix< moris::IdMat > tNodeIds( { { 1, 2, 3, 4, 5 } } );

        moris::Matrix< moris::IndexMat > tElementsAncestry( { { 0 } } );    // Not used
        moris::Matrix< moris::DDSTMat >  tElementNodeParentRanks( 1, 4, 0 );
        moris::Matrix< moris::IndexMat > tParentEdgeInds( { { 0, 1, 2, 3, 4, 5 } } );
        moris::Matrix< moris::DDSTMat >  tParentEdgeRanks( 1, 6, 1 );
        moris::Matrix< moris::IndexMat > tParentFaceInds( { { 0, 1, 2, 3 } } );
        moris::Matrix< moris::DDSTMat >  tParentFaceRanks( 1, 4, 2 );

        moris::moris_index tNumPermutations = 6;

        for ( moris::moris_index iEdge = 0; iEdge < tNumPermutations; iEdge++ )
        {
            // Initialize Template
            Mesh_Modification_Template tMeshTemplate( tElementsAncestry( 0, 0 ),
                    0,
                    tNodeIndex,
                    tNodeIndex,
                    tElementNodeParentRanks,
                    tParentEdgeInds,
                    tParentEdgeRanks,
                    tParentFaceInds,
                    tParentFaceRanks,
                    TemplateType::TET_4 );

            // Initialize child mesh with template (in this case a tet4)
            Child_Mesh tChildMesh( tMeshTemplate );

            tChildMesh.add_new_geometry_interface( 0 );

            // add new node indices
            tChildMesh.add_node_indices( { { 4 } } );
            tChildMesh.add_node_ids( tNodeIds );

            // Set up node coordinates
            moris::Matrix< moris::DDRMat > tNodeCoords;
            setup_node_coordinates_bisected( iEdge, tNodeCoords );

            // Compute base tet volume
            real tTetVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );

            // Compute base element surface normals (parent faces)
            moris::Matrix< moris::IndexMat > const &tParentFaceToNode = tChildMesh.get_face_to_node();
            size_t                                  tNumParentFaces   = tParentFaceToNode.n_rows();
            moris::Matrix< moris::DDRMat >          tParentFaceNormals( 3, tNumParentFaces );

            // Iterate through and compute all face normals
            for ( size_t iF = 0; iF < tNumParentFaces; iF++ )
            {
                // Get face iF nodes
                moris::Matrix< moris::IndexMat > tFaceNodes = tParentFaceToNode.get_row( iF );

                moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                tParentFaceNormals.set_column( iF, tFaceNormal );
            }

            // Initialize/set  intersection connectivity in child mehs
            tChildMesh.add_entity_to_intersect_connectivity( 4, iEdge, true );

            tChildMesh.modify_child_mesh( TemplateType::HIERARCHY_TET4 );

            // Verify that if we set each element to the same bulk phase we can traverse the element to element connectivity
            moris::moris_index               tMax       = std::numeric_limits< moris::moris_index >::max();
            size_t                           tNumPhases = 2;
            moris::Matrix< moris::IndexMat > tActiveElements( { { 0, 1 } } );
            moris::Matrix< moris::IndexMat > tIncludedElementMarker( 1, 2, 1 );
            moris::Matrix< moris::IndexMat > tElementPhase( 1, 2, 0 );

            moris::moris_index               tMaxFloodFill = 0;
            moris::Matrix< moris::IndexMat > tElementSubphase =
                    flood_fill( tChildMesh.get_element_to_element(),
                            tElementPhase,
                            tActiveElements,
                            tIncludedElementMarker,
                            tNumPhases,
                            tMax,
                            tMaxFloodFill,
                            true );

            moris::Matrix< moris::IndexMat > tExpElementSubphase( 2, 1, 0 );
            CHECK( equal_to( tExpElementSubphase, tElementSubphase ) );

            // Verify that the tets created have correct topology
            bool tValidTopo = verify_tet4_topology( tChildMesh.get_element_to_node(),
                    tChildMesh.get_element_to_edge(),
                    tChildMesh.get_element_to_face(),
                    tChildMesh.get_edge_to_node(),
                    tChildMesh.get_face_to_node() );

            CHECK( tValidTopo );

            // verify volume is conserved
            real tTotalChildVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );
            ;

            CHECK( approximate( tTetVol, tTotalChildVol ) );

            // Check ancestry of faces
            moris::Matrix< moris::IndexMat > const &tFaceToNode      = tChildMesh.get_face_to_node();
            moris::Matrix< moris::IndexMat > const &tFaceParentInds  = tChildMesh.get_face_parent_inds();
            moris::Matrix< moris::DDSTMat > const  &tFaceParentRanks = tChildMesh.get_face_parent_ranks();
            size_t                                  tNumFaces        = tFaceToNode.n_rows();
            moris::Matrix< moris::DDRMat >          tFaceNormals( 3, tNumFaces );

            // Iterate through and compute all face normals
            for ( size_t iF = 0; iF < tNumFaces; iF++ )
            {
                // Get face iF nodes
                moris::Matrix< moris::IndexMat > tFaceNodes = tFaceToNode.get_row( iF );

                moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                tFaceNormals.set_column( iF, tFaceNormal );
            }

            // Check to see if the surface normal matches parent face normals
            // Count for number of children faces born on a parent face
            size_t tChildFacewithParentFaceRank = 0;
            for ( size_t iF = 0; iF < tNumFaces; iF++ )
            {
                if ( tFaceParentRanks( 0, iF ) == 2 )
                {
                    moris::Matrix< moris::DDRMat > tChildFaceNormal  = tFaceNormals.get_column( iF );
                    moris::Matrix< moris::DDRMat > tParentFaceNormal = tParentFaceNormals.get_column( tFaceParentInds( 0, iF ) );
                    tChildFacewithParentFaceRank++;
                    CHECK( equal_to( tChildFaceNormal, tParentFaceNormal ) );
                }
            }
            CHECK( tChildFacewithParentFaceRank == 6 );

            // Check Edge Ancestry
            moris::mtk::Cell_Info_Tet4              tConn;
            moris::Matrix< moris::IndexMat >        tParentEdgeToNodeMap = tConn.get_node_to_edge_map();
            moris::Matrix< moris::IndexMat > const &tEdgeToNode          = tChildMesh.get_edge_to_node();
            moris::Matrix< moris::IndexMat > const &tEdgeParentInds      = tChildMesh.get_edge_parent_inds();
            moris::Matrix< moris::DDSTMat > const  &tEdgeParentRanks     = tChildMesh.get_edge_parent_ranks();
            for ( size_t iEdge = 0; iEdge < tEdgeToNode.n_rows(); iEdge++ )
            {
                // Verify all edges with a parent rank of 1 have nodes which belong on that edge
                if ( tEdgeParentRanks( iEdge ) == 1 )
                {
                    moris::Matrix< moris::IndexMat > tParentEdgeNodes = tParentEdgeToNodeMap.get_row( tEdgeParentInds( iEdge ) );
                    moris::Matrix< moris::IndexMat > tChildEdgeNodes  = tEdgeToNode.get_row( iEdge );

                    CHECK( verify_edge_is_on_parent_edge( tParentEdgeNodes, tNodeCoords, tChildEdgeNodes ) );
                }
            }
        }
    }

    TEST_CASE( "2 Edge intersected Tetrahedral Template", "[2_NODE]" )
    {
        moris::Matrix< moris::IndexMat > tNodeIndex( { { 0, 1, 2, 3 } } );

        moris::Matrix< moris::IdMat > tNodeIds( { { 1, 2, 3, 4, 5, 6 } } );

        moris::Matrix< moris::IndexMat > tElementsAncestry( { { 0 } } );    // Not used
        moris::Matrix< moris::DDSTMat >  tElementNodeParentRanks( 1, 4, 0 );
        moris::Matrix< moris::IndexMat > tParentEdgeInds( { { 0, 1, 2, 3, 4, 5 } } );
        moris::Matrix< moris::DDSTMat >  tParentEdgeRanks( 1, 6, 1 );
        moris::Matrix< moris::IndexMat > tParentFaceInds( { { 0, 1, 2, 3 } } );
        moris::Matrix< moris::DDSTMat >  tParentFaceRanks( 1, 4, 2 );

        Permutations tPermutations( "2_node" );
        size_t       tNumPermutations = tPermutations.get_num_permutations();

        for ( size_t iPerm = 0; iPerm < tNumPermutations; iPerm++ )
        {
            // select template
            Cell< size_t > tCurrentPermutation = tPermutations.get_permutation( iPerm );

            moris::moris_index tEdgeL = tCurrentPermutation( 0 );
            moris::moris_index tEdgeH = tCurrentPermutation( 1 );

            // Initialize Template
            Mesh_Modification_Template tMeshTemplate( tElementsAncestry( 0, 0 ),
                    0,
                    tNodeIndex,
                    tNodeIndex,
                    tElementNodeParentRanks,
                    tParentEdgeInds,
                    tParentEdgeRanks,
                    tParentFaceInds,
                    tParentFaceRanks,
                    TemplateType::TET_4 );

            // Initialize child mesh with template (in this case a tet4)
            Child_Mesh tChildMesh( tMeshTemplate );

            tChildMesh.add_new_geometry_interface( 0 );

            // add new node indices
            tChildMesh.add_node_indices( { { 4, 5 } } );
            tChildMesh.add_node_ids( tNodeIds );

            // Set up node coordinates
            moris::Matrix< moris::DDRMat > tNodeCoords;
            setup_node_coordinates_2_node( tEdgeL, tEdgeH, tNodeCoords );

            // Compute base tet volume
            real tTetVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );
            // Compute base element surface normals (parent faces)
            moris::Matrix< moris::IndexMat > const &tParentFaceToNode = tChildMesh.get_face_to_node();
            size_t                                  tNumParentFaces   = tParentFaceToNode.n_rows();
            moris::Matrix< moris::DDRMat >          tParentFaceNormals( 3, tNumParentFaces );

            // Iterate through and compute all face normals
            for ( size_t iF = 0; iF < tNumParentFaces; iF++ )
            {
                // Get face iF nodes
                moris::Matrix< moris::IndexMat > tFaceNodes = tParentFaceToNode.get_row( iF );

                moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                tParentFaceNormals.set_column( iF, tFaceNormal );
            }

            // Initialize/set  intersection connectivity in child mehs
            tChildMesh.add_entity_to_intersect_connectivity( 4, tEdgeL, true );
            tChildMesh.add_entity_to_intersect_connectivity( 5, tEdgeH, true );

            tChildMesh.modify_child_mesh( TemplateType::HIERARCHY_TET4 );

            // Verify that if we set each element to the same bulk phase we can traverse the element to element connectivity
            moris::moris_index               tMax       = std::numeric_limits< moris::moris_index >::max();
            size_t                           tNumPhases = 2;
            moris::Matrix< moris::IndexMat > tActiveElements( { { 0, 1, 2 } } );
            moris::Matrix< moris::IndexMat > tIncludedElementMarker( 1, 3, 1 );
            moris::Matrix< moris::IndexMat > tElementPhase( 1, 3, 0 );
            moris::moris_index               tMaxFloodFill = 0;
            moris::Matrix< moris::IndexMat > tElementSubphase =
                    flood_fill( tChildMesh.get_element_to_element(),
                            tElementPhase,
                            tActiveElements,
                            tIncludedElementMarker,
                            tNumPhases,
                            tMax,
                            tMaxFloodFill,
                            true );

            moris::Matrix< moris::IndexMat > tExpElementSubphase( 3, 1, 0 );
            CHECK( equal_to( tExpElementSubphase, tElementSubphase ) );

            // Verify that the tets created have correct topology
            bool tValidTopo = verify_tet4_topology( tChildMesh.get_element_to_node(),
                    tChildMesh.get_element_to_edge(),
                    tChildMesh.get_element_to_face(),
                    tChildMesh.get_edge_to_node(),
                    tChildMesh.get_face_to_node() );

            CHECK( tValidTopo );

            // verify volume is conserved
            real tTotalChildVol = gen::compute_volume_for_multiple_tets( tNodeCoords, tChildMesh.get_element_to_node() );

            CHECK( approximate( tTetVol, tTotalChildVol ) );

            // Check ancestry of faces
            moris::Matrix< moris::IndexMat > const &tFaceToNode      = tChildMesh.get_face_to_node();
            moris::Matrix< moris::IndexMat > const &tFaceParentInds  = tChildMesh.get_face_parent_inds();
            moris::Matrix< moris::DDSTMat > const  &tFaceParentRanks = tChildMesh.get_face_parent_ranks();
            size_t                                  tNumFaces        = tFaceToNode.n_rows();
            moris::Matrix< moris::DDRMat >          tFaceNormals( 3, tNumFaces );

            // Iterate through and compute all face normals
            for ( size_t iF = 0; iF < tNumFaces; iF++ )
            {
                // Get face iF nodes
                moris::Matrix< moris::IndexMat > tFaceNodes = tFaceToNode.get_row( iF );

                moris::Matrix< moris::DDRMat > tFaceNormal( 3, 1, 9.0 );
                gen::compute_tri_surface_normal( tFaceNodes, tNodeCoords, tFaceNormal, true );

                tFaceNormals.set_column( iF, tFaceNormal );
            }

            // Check to see if the surface normal matches parent face normals
            // Count for number of children faces born on a parent face
            size_t tChildFacewithParentFaceRank = 0;
            for ( size_t iF = 0; iF < tNumFaces; iF++ )
            {
                if ( tFaceParentRanks( 0, iF ) == 2 )
                {
                    moris::Matrix< moris::DDRMat > tChildFaceNormal  = tFaceNormals.get_column( iF );
                    moris::Matrix< moris::DDRMat > tParentFaceNormal = tParentFaceNormals.get_column( tFaceParentInds( 0, iF ) );
                    tChildFacewithParentFaceRank++;
                    CHECK( equal_to( tChildFaceNormal, tParentFaceNormal ) );
                }
            }
            CHECK( tChildFacewithParentFaceRank == 8 );

            // Check Edge Ancestry
            moris::mtk::Cell_Info_Tet4              tConn;
            moris::Matrix< moris::IndexMat >        tParentEdgeToNodeMap = tConn.get_node_to_edge_map();
            moris::Matrix< moris::IndexMat > const &tEdgeToNode          = tChildMesh.get_edge_to_node();
            moris::Matrix< moris::IndexMat > const &tEdgeParentInds      = tChildMesh.get_edge_parent_inds();
            moris::Matrix< moris::DDSTMat > const  &tEdgeParentRanks     = tChildMesh.get_edge_parent_ranks();
            for ( size_t iEdge = 0; iEdge < tEdgeToNode.n_rows(); iEdge++ )
            {
                // Verify all edges with a parent rank of 1 have nodes which belong on that edge
                if ( tEdgeParentRanks( iEdge ) == 1 )
                {
                    moris::Matrix< moris::IndexMat > tParentEdgeNodes = tParentEdgeToNodeMap.get_row( tEdgeParentInds( iEdge ) );
                    moris::Matrix< moris::IndexMat > tChildEdgeNodes  = tEdgeToNode.get_row( iEdge );

                    CHECK( verify_edge_is_on_parent_edge( tParentEdgeNodes, tNodeCoords, tChildEdgeNodes ) );
                }
            }
        }
    }

    // TEST_CASE("Recursive modification of child mesh","[NH_RECURSION_4_NODE]")
    //{
    //     // Tests:
    //     // Floodfill which checks whether the element to element connectivity is traversable
    //     // Topology which checks whether the new tets have the correct topology
    //
    //     moris::Matrix< moris::IndexMat > tNodeIndex({{0,1,2,3}});
    //
    //     moris::Matrix< moris::IdMat > tNodeIds({{1,2,3,4,5,6,7,8}});
    //
    //
    //     moris::Matrix< moris::IndexMat > tElementsAncestry({{0}}); // Not used
    //     moris::Matrix< moris::IndexMat > tParentEdgeInds({{0,1,2,3,4,5}});
    //     moris::Matrix< moris::DDSTMat > tParentEdgeRanks(1,6,1);
    //     moris::Matrix< moris::IndexMat > tParentFaceInds({{0,1,2,3}});
    //     moris::Matrix< moris::DDSTMat > tParentFaceRanks(1,4,2);
    //
    //     Cell<std::string> tCaseNames = {"pa","na","pb","nb","pc","nc"};
    ////    Cell<std::string> tCaseNames = {"pa"};
    //
    //    for(size_t iCase = 0; iCase<tCaseNames.size(); iCase++)
    //    {
    //
    //        std::string & tCases = tCaseNames(iCase);
    //        Permutations tPermutations(tCases);
    //        size_t tNumPermutations = tPermutations.get_num_permutations();
    ////        size_t tNumPermutations = 1;
    //
    //        // Second template insertion permutation
    //        for(size_t iPerm2 = 0; iPerm2<tNumPermutations; iPerm2++)
    //        {
    //            // select template
    //            Cell<size_t> tCurrentPermutation = tPermutations.get_permutation(iPerm2);
    //
    //            // Edge ordinals
    //            moris::moris_index tEdgeOrdL2  = tCurrentPermutation(0);
    //            moris::moris_index tEdgeOrdML2 = tCurrentPermutation(1);
    //            moris::moris_index tEdgeOrdMH2 = tCurrentPermutation(2);
    //            moris::moris_index tEdgeOrdH2  = tCurrentPermutation(3);
    //
    //            size_t tPermutationId2 = 1000*tEdgeOrdH2 + 100*tEdgeOrdMH2 + 10 * tEdgeOrdML2 + tEdgeOrdL2;
    //
    //            // First template insertion permutation
    //            for(size_t iPerm1 = 0; iPerm1<tNumPermutations; iPerm1++)
    //            {
    //
    //               // Initialize Template
    //                Mesh_Modification_Template tMeshTemplate(tElementsAncestry(0,0),
    //                                                                                                                 0,
    //                                                                                                                 tNodeIndex,
    //                                                                                                                 tParentEdgeInds,
    //                                                                                                                 tParentEdgeRanks,
    //                                                                                                                 tParentFaceInds,
    //                                                                                                                 tParentFaceRanks,
    //                                                                                                                 TemplateType::TET_4);
    //
    //                // Initialize child mesh with template (in this case a tet4)
    //                Child_Mesh_Test tChildMesh(tMeshTemplate);
    //
    //
    //                // add new node indices
    //                tChildMesh.add_node_indices({{4,5,6,7}});
    //                tChildMesh.add_node_ids(tNodeIds);
    //
    //                // select template
    //                Cell<size_t> tCurrentPermutation = tPermutations.get_permutation(iPerm1);
    //                moris::moris_index tEdgeL1  = tCurrentPermutation(0);
    //                moris::moris_index tEdgeML1 = tCurrentPermutation(1);
    //                moris::moris_index tEdgeMH1 = tCurrentPermutation(2);
    //                moris::moris_index tEdgeH1  = tCurrentPermutation(3);
    //
    //                moris::moris_index tPermutationId1 = 1000*tEdgeH1 + 100*tEdgeMH1 + 10 * tEdgeML1 + tEdgeL1;
    //
    //                // Set up node coordinates
    //                moris::Matrix< moris::DDRMat > tNodeCoords;
    //                setup_node_coordinates_4_node(tEdgeL1,tEdgeML1,tEdgeMH1,tEdgeH1,tNodeCoords);
    //
    //                // ------------------------------------------------
    //                // Compute Base Tet Information
    //                // ------------------------------------------------
    //
    //                // Compute base tet volume
    //                real tTetVol = gen::compute_volume_for_multiple_tets(tNodeCoords,tChildMesh.get_element_to_node());
    //
    //                // Compute base element surface normals (parent faces)
    //                moris::Matrix< moris::IndexMat > const & tParentFaceToNode = tChildMesh.get_face_to_node();
    //                size_t tNumParentFaces = tParentFaceToNode.n_rows();
    //                moris::Matrix< moris::DDRMat > tParentFaceNormals(3,tNumParentFaces);
    //
    //                // Iterate through and compute all face normals
    //                for( size_t iF = 0; iF<tNumParentFaces; iF++)
    //                {
    //                    // Get face iF nodes
    //                    moris::Matrix< moris::IndexMat > tFaceNodes = tParentFaceToNode.get_row(iF);
    //
    //                    moris::Matrix< moris::DDRMat > tFaceNormal(3,1,9.0);
    //                    compute_tri_surface_normal( tFaceNodes,tNodeCoords, tFaceNormal, true);
    //
    //                    tParentFaceNormals.set_column(iF,tFaceNormal);
    //                }
    //
    //                // ------------------------------------------------
    //                // Insert First Template
    //                // ------------------------------------------------
    //
    //                // Initialize/set  intersection connectivity in child mesh for first level modification
    //                moris::Matrix< moris::IndexMat > tIntersectConn({{4,4,5,6,7,INT_MAX,INT_MAX,tEdgeL1,tEdgeML1,tEdgeMH1,tEdgeH1,INT_MAX,INT_MAX}});
    //                tChildMesh.set_intersect_connectivity(tIntersectConn);
    //                tChildMesh.modify_child_mesh(TemplateType::HIERARCHY_TET4);
    //
    //                // ------------------------------------------------
    //                // Second Template Insertion
    //                // ------------------------------------------------
    //                // Replacing element 0
    //                // Node Ids for this template
    //                moris::Matrix< moris::IndexMat > tNodeInds2({{11,10,9,8}});
    //                moris::Matrix< moris::IdMat > tNodeIds2({{12,11,10,9}});
    //                tChildMesh.add_node_indices(tNodeInds2);
    //                tChildMesh.add_node_ids(tNodeIds);
    //
    //                // Initialize intersection connectivity
    //                tChildMesh.init_intersect_connectivity();
    //
    //
    //                moris::Matrix< moris::IndexMat > const & tElemToEdge = tChildMesh.get_element_to_edge();
    //                moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();
    //                moris::moris_index  tEdgeL2 = tElemToEdge(0,tEdgeOrdL2);
    //                moris::moris_index  tEdgeML2 = tElemToEdge(0,tEdgeOrdML2);
    //                moris::moris_index  tEdgeMH2 = tElemToEdge(0,tEdgeOrdMH2);
    //                moris::moris_index  tEdgeH2 = tElemToEdge(0,tEdgeOrdH2);
    //
    //                tChildMesh.add_entity_to_intersect_connectivity(8, tEdgeL2 , 1);
    //                tChildMesh.add_entity_to_intersect_connectivity(9, tEdgeML2, 1);
    //                tChildMesh.add_entity_to_intersect_connectivity(10, tEdgeMH2, 1);
    //                tChildMesh.add_entity_to_intersect_connectivity(11, tEdgeH2 , 1);
    //
    //                // Set up node coordinates for second template insertion
    //                tNodeCoords.resize(tNodeCoords.n_rows() + 4,tNodeCoords.n_cols());
    //
    //                moris::Matrix< moris::DDRMat > tMidSideCoord(1,3);
    //
    //                moris::Matrix< moris::DDRMat > tN1Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeL2,0));
    //                moris::Matrix< moris::DDRMat > tN2Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeL2,1));
    //
    //                tMidSideCoord =  0.5 * ( tN1Coord + tN2Coord ) ;
    //                tNodeCoords.set_row(11,tMidSideCoord);
    //
    //                tN1Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeML2,0));
    //                tN2Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeML2,1));
    //                tMidSideCoord =  0.5 * ( tN1Coord + tN2Coord ) ;
    //                tNodeCoords.set_row(10,tMidSideCoord);
    //
    //                tN1Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeMH2,0));
    //                tN2Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeMH2,1));
    //                tMidSideCoord =  0.5 * ( tN1Coord + tN2Coord ) ;
    //                tNodeCoords.set_row(9,tMidSideCoord);
    //
    //
    //                tN1Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeH2,0));
    //                tN2Coord = tNodeCoords.get_row(tEdgeToNode(tEdgeH2,1));
    //                tMidSideCoord =  0.5 * ( tN1Coord + tN2Coord ) ;
    //                tNodeCoords.set_row(8,tMidSideCoord);
    //
    ////                print_to_matlab(tNodeCoords.matrix_base(),"tNodeCoords");
    //
    //                tChildMesh.modify_child_mesh(TemplateType::HIERARCHY_TET4);
    ////                print_to_matlab(tChildMesh.get_element_to_node().matrix_base(),"tElemNode");
    //
    //                // ------------------------------------------------
    //                // Test the Generated Mesh
    //                // ------------------------------------------------
    //
    //                // Verify that the tets created have correct topology
    //                bool tValidTopo = verify_tet4_topology(tChildMesh.get_element_to_node(),
    //                                                       tChildMesh.get_element_to_edge(),
    //                                                       tChildMesh.get_element_to_face(),
    //                                                       tChildMesh.get_edge_to_node(),
    //                                                       tChildMesh.get_face_to_node());
    //
    //                CHECK(tValidTopo);
    //
    //                // verify volume is conserved
    //                real tTotalChildVol = gen::compute_volume_for_multiple_tets(tNodeCoords,tChildMesh.get_element_to_node());
    //                CHECK(approximate(tTetVol,tTotalChildVol));
    //
    //
    //                // Check ancestry of faces
    //                moris::Matrix< moris::IndexMat > const & tFaceToNode      = tChildMesh.get_face_to_node();
    //                moris::Matrix< moris::IndexMat > const & tFaceParentInds  = tChildMesh.get_face_parent_inds();
    //                moris::Matrix< moris::DDSTMat > const & tFaceParentRanks = tChildMesh.get_face_parent_ranks();
    //                size_t tNumFaces = tFaceToNode.n_rows();
    //                moris::Matrix< moris::DDRMat > tFaceNormals(3,tNumFaces);
    //
    //                // Iterate through and compute all face normals
    //                for( size_t iF = 0; iF<tNumFaces; iF++)
    //                {
    //                    // Get face iF nodes
    //                    moris::Matrix< moris::IndexMat > tFaceNodes = tFaceToNode.get_row(iF);
    //
    //                    moris::Matrix< moris::DDRMat > tFaceNormal(3,1,9.0);
    //                    compute_tri_surface_normal( tFaceNodes,tNodeCoords, tFaceNormal, true);
    //
    //                    tFaceNormals.set_column(iF,tFaceNormal);
    //                }
    //
    //                // Check to see if the surface normal matches parent face normals
    //                // Count for number of children faces born on a parent face
    //                size_t tChildFacewithParentFaceRank = 0;
    //                for( size_t iF = 0; iF<tNumFaces; iF++)
    //                {
    //                    if( tFaceParentRanks(0,iF) == 2)
    //                    {
    //                        moris::Matrix< moris::DDRMat > tChildFaceNormal  = tFaceNormals.get_column(iF);
    //                        moris::Matrix< moris::DDRMat > tParentFaceNormal = tParentFaceNormals.get_column(tFaceParentInds(0,iF));
    //                        tChildFacewithParentFaceRank++;
    //
    //                        CHECK(equal_to(tChildFaceNormal, tParentFaceNormal));
    //                    }
    //                }
    //
    //            }
    //        }
    //    }
    //}

}    // namespace xtk
