/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * fn_HMR_refinement_transition_locations.hpp  
 * 
 */
#ifndef SRC_fn_HMR_refinement_transition_locations
#define SRC_fn_HMR_refinement_transition_locations

#include "cl_Cell.hpp"

namespace moris::hmr
{

    //-----------------------------------------------------------------------------
    // Ordering of child element ordinals in the function below
    //
    //            6 ------ 7
    //           /|       /|
    //          / |      / |
    //         4 ------ 5  |
    //         |  |     |  |
    //         |  2 ----|- 3
    //         | /      | /
    //  z      |/       |/
    //  |  y   0 ------ 1
    //  | /
    //  |/____ x
    //
    //-----------------------------------------------------------------------------

    /**
     * @brief Get the refinement transition location from a big to a small element knowing only the small element's child ordinal
     *
     * @param aBigElementSideOrdinal Side ordinal originating from the big element
     * @param aSmallElementChildOrdinal child ordinal of the small element
     * @param aNumSpatialDims number of spatial dimenstions
     * @return moris_index transition location from big to small element (=hypothetical child ordinal of the big element in the place of the transition)
     */
    moris_index
    get_refinement_transition_location_for_neighbor_child_ordinal(
        moris_index aBigElementSideOrdinal,
        moris_index aSmallElementChildOrdinal,
        uint        aNumSpatialDims )
    {

        switch( aNumSpatialDims )
        {
            // =====================================
            // 2D
            case 2 :
            {
                // Side Ordinal
                switch( aBigElementSideOrdinal )
                {
                    // ------------------------------
                    // aCenterSideOrdinal = 0
                    case 0 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 2 :
                            {
                                return 0;
                                break;
                            }
                            case 3 :
                            {
                                return 1;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 0" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 1
                    case 1 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 0 :
                            {
                                return 1;
                                break;
                            }
                            case 2 :
                            {
                                return 3;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 1" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 2
                    case 2 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 0 :
                            {
                                return 3;
                                break;
                            }
                            case 1 :
                            {
                                return 2;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 2" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 3
                    case 3 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 1 :
                            {
                                return 2;
                                break;
                            }
                            case 3 :
                            {
                                return 0;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 3" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    default:
                    {
                        MORIS_ERROR( false,
                            "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                            "side ordinal out of bounds for 2D" );
                    }
                    // ------------------------------
                }

                break;
            }
            // =====================================

            // =====================================
            // 3D
            case 3 :
            {
                // Side Ordinal
                switch( aBigElementSideOrdinal )
                {
                    // ------------------------------
                    // aCenterSideOrdinal = 0
                    case 0 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 2 :
                            {
                                return 0;
                                break;
                            }
                            case 3 :
                            {
                                return 1;
                                break;
                            }
                            case 6 :
                            {
                                return 5;
                                break;
                            }
                            case 7 :
                            {
                                return 4;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 0" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 1
                    case 1 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 0 :
                            {
                                return 1;
                                break;
                            }
                            case 2 :
                            {
                                return 3;
                                break;
                            }
                            case 4 :
                            {
                                return 7;
                                break;
                            }
                            case 6 :
                            {
                                return 5;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 1" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 2
                    case 2 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 0 :
                            {
                                return 3;
                                break;
                            }
                            case 1 :
                            {
                                return 2;
                                break;
                            }
                            case 4 :
                            {
                                return 6;
                                break;
                            }
                            case 5 :
                            {
                                return 7;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 2" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 3
                    case 3 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 1 :
                            {
                                return 2;
                                break;
                            }
                            case 3 :
                            {
                                return 0;
                                break;
                            }
                            case 5 :
                            {
                                return 4;
                                break;
                            }
                            case 7 :
                            {
                                return 6;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 3" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 4
                    case 4 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 4 :
                            {
                                return 2;
                                break;
                            }
                            case 5 :
                            {
                                return 3;
                                break;
                            }
                            case 6 :
                            {
                                return 1;
                                break;
                            }
                            case 7 :
                            {
                                return 0;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 4" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    // aCenterSideOrdinal = 5
                    case 5 :
                    {
                        switch( aSmallElementChildOrdinal )
                        {
                            case 0 :
                            {
                                return 4;
                                break;
                            }
                            case 1 :
                            {
                                return 5;
                                break;
                            }
                            case 2 :
                            {
                                return 7;
                                break;
                            }
                            case 3 :
                            {
                                return 6;
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                                    "there is no center child ordinal connected to this neighbor child ordinal - 5" );
                            }
                        }

                        break;
                    }
                    // ------------------------------

                    // ------------------------------
                    default:
                    {
                        MORIS_ERROR( false,
                            "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                            "side ordinal out of bounds for 3D" );
                    }
                    // ------------------------------
                }

                break;
            }
            // =====================================

            // =====================================
            // invalid number of spatial dimensions
            default:
            {
                MORIS_ERROR( false,
                    "hmr::get_refinement_transition_location_for_neighbor_child_ordinal() - "
                    "Number of spatial dimensions must be 2 or 3" );
            }
            // =====================================

        }

        // default / failed
        return -1;

    }  // end: function

    //-----------------------------------------------------------------------------

    /**
     * @brief Get the refinement transition location from a big to a small element knowing only the small element's child ordinal
     *
     * @param aSmallElementSideOrdinal Side ordinal originating from the small element
     * @param aSmallElementChildOrdinal child ordinal of the small element
     * @param aNumSpatialDims number of spatial dimenstions
     * @return moris_index transition location from big to small element (=hypothetical child ordinal of the big element in the place of the transition)
     */
    moris_index
    get_refinement_transition_location_for_neighbor_side_and_child_ordinal(
        moris_index aSmallElementSideOrdinal,
        moris_index aSmallElementChildOrdinal,
        uint        aNumSpatialDims )
    {
        // convert to the small element's child ordinal to that of the big element
        moris_index tBigElementSideOrdinal = -1;

        switch( aNumSpatialDims )
        {
            case 2 :
            {
                moris::Cell< moris_index > tConvert = { 2, 3, 0, 1 };
                tBigElementSideOrdinal = tConvert( aSmallElementSideOrdinal );
                break;
            }
            case 3 :
            {
                moris::Cell< moris_index > tConvert = { 2, 3, 0, 1, 5, 4 };
                tBigElementSideOrdinal = tConvert( aSmallElementSideOrdinal );
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                    "hmr::get_refinement_transition_location_for_neighbor_side_and_child_ordinal() - "
                    "Number of spatial dimensions must be 2 or 3" );
            }
        }

        // retrieve the value using the above function
        return get_refinement_transition_location_for_neighbor_child_ordinal( tBigElementSideOrdinal, aSmallElementChildOrdinal, aNumSpatialDims );

    } // end: function

    //-----------------------------------------------------------------------------

} // namespace moris

#endif /* SRC_fn_HMR_refinement_transition_locations */